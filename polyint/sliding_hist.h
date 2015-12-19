/*
 * sliding_hist.h
 *
 *  Created on: 2010-10-26
 *      Author: zhangyuting
 */


#ifndef SLIDING_HIST_H_
#define SLIDING_HIST_H_

#include "config.h"

#include "zutils.h"
#include "polyint.h"
#include "hashspecial.h"
#include "polygroup.h"
#include <climits>
#include <vector>
#include <memory>


/* --------------------- NumericTraits ---------------- */

template <typename T>
class piNumericTraitsBase {
public:
	enum {
		Unknown = 0,
		IntegerType = 1,
		FloatType = 2
	};
	int GetType() { return Unknown; }
	typedef piZeroScheme::Basic<T> ZeroScheme;
};

template <typename T>
class piNumericTraits : public piNumericTraitsBase<T> {
};

#define PI_NUMERIC_TRAITS_SPECIFY2INTEGER( TYPE ) \
	template <>	\
	class piNumericTraits<TYPE> : public piNumericTraitsBase<TYPE> {	\
	public:	\
		int GetType() { return IntegerType; }	\
	};


#define PI_NUMERIC_TRAITS_SPECIFY2FLOAT( TYPE ) \
	template <>	\
	class piNumericTraits<TYPE> : public piNumericTraitsBase<TYPE> {	\
	public:	\
		int GetType() { return FloatType; }	\
		typedef piZeroScheme::Epsilon<TYPE> ZeroScheme;	\
	};


PI_NUMERIC_TRAITS_SPECIFY2INTEGER(int);
//PI_NUMERIC_TRAITS_SPECIFY2INTEGER(piFloat);
PI_NUMERIC_TRAITS_SPECIFY2FLOAT(piFloat);



/* ------------- Histogram types ------------ */

#define piArrHist	array1d

template <typename T, class ZERO_SCHEME = typename piNumericTraits<T>::ZeroScheme >
class piHashHist :	public piHashSpecial<T,ZERO_SCHEME> {
public:
	typedef piHashHist<T,ZERO_SCHEME> 		MyType;
	typedef piHashSpecial<T,ZERO_SCHEME>	HashType;

	typedef typename HashType::ListType			ListType;
	typedef typename HashType::iterator			iterator;
	typedef	typename HashType::const_iterator	const_iterator;
public:
	explicit piHashHist( size_t table_size, size_t list_size = 0, const ZERO_SCHEME zero_scheme = ZERO_SCHEME() )
		: HashType( table_size, list_size, zero_scheme ) {}
	piHashHist( const MyType& _orig ) : HashType( _orig ) {}

	const MyType& operator += ( const MyType& r ) {
		for ( typename MyType::const_iterator iter = r.begin(); iter!= r.end(); ++iter )
			this->inc_s(iter->idx, iter->val );
		return *this;
	}
	const MyType& operator += ( const ListType& r ) {
		for ( typename ListType::const_iterator iter = r.begin(); iter!= r.end(); iter.vecAdvance() )
			this->inc(iter->idx, iter->val );
		return *this;
	}
	const MyType& operator -= ( const MyType& r ) {
		for ( typename MyType::const_iterator iter = r.begin(); iter!= r.end(); ++iter )
			this->inc_s(iter->idx, - iter->val );
		return *this;
	}
	const MyType& operator -= ( const ListType& r ) {
		for ( typename ListType::const_iterator iter = r.begin(); iter!= r.end(); iter.vecAdvance() )
			this->dec(iter->idx, iter->val );
		return *this;
	}
	template <typename ARR_T>
	void add2array( ARR_T& arr ) const {
		for ( const_iterator iter = this->begin(); iter!= this->end(); ++iter )
			arr[iter->idx] += iter->val;
	}
};

/* ------------- Histogram Functions ------------ */

//Please do NOT use SIZE_T as T_IN, it will results in overflow
//Try to find a way to forbid the use of unsigned type during compiling
template <typename T_OUT, typename T_IN, typename BIN_TYPE>
class piHistFunConcept {
public:
	typedef T_OUT output_type;
	typedef T_IN  input_type;
	typedef BIN_TYPE bin_type;	//Generally bin_type is size_t
	/*
	output_type operator () (bin_type, input_type) {}
	output_type zero() {}
	*/
};

//Function for histogram output
class piFun4Hist : public piHistFunConcept<piHashHist<piFloat>, piFloat, size_t> {
	size_t _length;
public:
	piFun4Hist( size_t length ) : _length(length) {}
	output_type operator () (bin_type b, input_type n) {
		output_type f(_length, 1);
		f.WriteTable(b, n);
		return f;
	}
	output_type zero() { return output_type(_length,_length); }
};

class piFun4TestSum : public piHistFunConcept<piFloat, piFloat, size_t> {
public:
	piFun4TestSum() {}
	output_type operator () (bin_type b, input_type n) {
		return (piFloat)(n-1)*(n-2);
	}
	output_type zero() { return 0.0; }
};




/* ------------------ Middle-level data  -------------- */

//LineSeg class (slope, length)
//LineSeg position
//LineSeg group
//Group   The last LineSeg (The reverse index)
//Y span of UnitSeg
//Y span of LongSeg


//Detected duplicated LongSeg analysis
//Detected duplicated UnitSeg analysis

//UnitSeg Map
//LongSeg Map
//Vertical Incremental Scheme:
////Precise:	reconstruction, r_d incremental
////Rasterized:	reconstruction, r_d incremental, vertical, diagonal
//Function value buffer array (one for horizontal, one for vertical)


//--------------------------------------------------------

class piSlidingWindows : public boost::noncopyable {
public:
	typedef array1d<iPoint>		PixelArr;
	typedef std::vector<iPoint>	PixelVec;
	typedef size_t				IdxType;
	typedef std::vector<IdxType>	SegIdxVec;
	typedef ptrvec<SegIdxVec>		SegIdxVecPVec;
	typedef array1d<SegIdxVecPVec>	SegIdxVecPVecArr;
	typedef std::vector<int>		SpanVec;
	typedef std::vector<int>		FracVec;
	typedef std::vector<int>		CoordinateVec;
	typedef std::vector<const piRasterizedLineSeg*>	LineSegPtrVec;
	typedef ptrvec<PixelArr>		PixelArrPVec;
	typedef piRasterizedLineSeg::PortionType	PortionType;
	typedef array1d<PortionType>				PortionArr;
	typedef std::vector<PortionType>			PortionVec;
	typedef ptrvec<PortionArr>					PortionArrPVec;
	typedef std::vector<iSpan>					iSpanVec;
private:
	const piPolyGroupArr& _shapes;
	const iPoint _stride;
	size_t _winNum;
	array1d<iISize>		_winSize;
	array1d<size_t>		_polyNum;

	ptrvec<PixelArr>	_interiorPixels;
	ptrvec<piPixelIntersection> _bodyIntersect;
	ptrvec<piPixelIntersection> _bodyIntersectC;
	ptrvec<PixelVec>		_boundaryPixels;
	ptrvec<PortionVec>		_boundaryPortions;

	SegIdxVecPVecArr	_segIdx;
	SegIdxVec			_seg2unit;
	FracVec				_segFrac;

	SegIdxVec			_unit2seg;
	PixelArrPVec		_unitSegPixels;
	PortionArrPVec		_unitSegPixelPortions;

	CoordinateVec	_segYlowest;
	SpanVec			_segYspan;
	SpanVec			_unitYspan;
	iSpanVec		_segXspan;
	iSpanVec		_unitXspan;
	LineSegPtrVec	_idx2seg;


public:
	const piPolyGroupArr& 		Shapes;
	const size_t 				WinNum;
	const array1d<iISize>&		WinSize;
	const array1d<size_t>&		PolyNum;

	const ptrvec<PixelArr>&		InteriorPixels;
	const ptrvec<piPixelIntersection>& BodyIntersect;
	const ptrvec<piPixelIntersection>& BodyIntersectC;
	const ptrvec<PixelVec>&		BoundaryPixels;
	const ptrvec<PortionVec>&	BoundaryPortions;

	const SegIdxVecPVecArr& 	SegIdx;
	const SegIdxVec&			Seg2Unit;
	const FracVec&				SegFrac;

	const SegIdxVec&			Unit2Seg;
	const PixelArrPVec&			UnitSegPixels;
	const PortionArrPVec&		UnitSegPixelPortions;

	const CoordinateVec&		SegYLowest;
	const SpanVec&				SegYSpan;
	const SpanVec&				UnitSegYSpan;
	const iSpanVec&				SegXSpan;
	const iSpanVec&				UnitSegXSpan;
	const LineSegPtrVec&		Idx2Seg;
	const iPoint&				Stride;
public:
	explicit piSlidingWindows( const piPolyGroupArr& shapes, const iPoint stride = iPoint(1,1) );
private:
	void _initialize( const piPolyGroupArr& shapes, const iPoint stride );
};

//--------------------------------------------------------------------

struct swUnitPixel {
	typedef piRasterizedLineSeg::PortionType PortionType;
	iPoint		pos;
	PortionType	por;
};


struct swSegStruct;

struct swUnitStruct {
	typedef array1d<swUnitPixel>	UnitPixelArr;
	typedef std::vector<swSegStruct*>	SegPtrVec;
	SegPtrVec		seg;
	UnitPixelArr	pixels;
	iPoint			vec;
	int				Yspan;
	iSpan			Xispan;
	const piUnitRasterizedLineSeg*	rdata;
	mutable iSpan		slideSpan;
	mutable void*		buffer;
	mutable void*		tag_buffer;
	swUnitStruct() : pixels(0, 0, UnitPixelArr::WithoutConstruction() ), buffer(NULL) {}
};

struct swSegStruct {
	enum {
		SS_Vertical				=	0,
		SS_UnitSeg				=	1,
		SS_NoBuffer				=	2,
		SS_Buffer				=	3,
	};
	swUnitStruct*	unit;
	int			lowY;
	int			Yspan;
	iSpan		Xispan;
	int			frac;
	int			fracY4sparse_buffer;
	iPoint		vec;
	iPoint		inc_vec;
	int			refreshType;
	mutable iSpan		slideSpan;
	mutable int			Xcycle;
	mutable void*		buffer;
	mutable void*		tag_buffer;
	swSegStruct() : buffer(NULL) {}
};

struct swEdge {
	swSegStruct*	seg;
	lineseg_t		line;
	bool			is_right;
};

class swBetterSlidingWindows;

struct swPolygonGroup {
public:
	typedef array1d<swEdge>		EdgeArr;
	typedef array1d<iPoint>		PixelArr;
	typedef array1d<swUnitPixel>	UPixelArr;
	struct ContinuousLine {
		int	  num;
		iSpan span;
	};
	typedef std::vector<ContinuousLine> CLVec;
	iISize			size;
	EdgeArr			edges;
	PixelArr		interior;
	CLVec			interiorRL;		//Run length encoded interior
private:
	PixelArr		interior_v_r1;
	PixelArr		interior_v_r2;
	PixelArr		interior_h_r1;
	PixelArr		interior_h_r2;
	bool			_is_naiveInc;	//0 - edge inc, 1 - naive inc
public:
	typedef	pair_t<PixelArr*>	PArrPPair;
	PArrPPair		bodyResidual_V;
	PArrPPair		bodyResidual_H;
	UPixelArr		boundary;
	bool			is_naiveInc;	//0 - edge inc, 1 - naive inc
	swPolygonGroup() : edges(0, 0, EdgeArr::WithoutConstruction() ), 
		interior(0, 0, PixelArr::WithoutConstruction()),
		interior_v_r1(0, 0, PixelArr::WithoutConstruction()),
		interior_v_r2(0, 0, PixelArr::WithoutConstruction()),
		interior_h_r1(0, 0, PixelArr::WithoutConstruction()),
		interior_h_r2(0, 0, PixelArr::WithoutConstruction()),
		bodyResidual_V(&interior_v_r1, &interior_v_r2), 
		bodyResidual_H(&interior_h_r1, &interior_h_r2), 
		boundary(0, 0, UPixelArr::WithoutConstruction()) {}
	friend class swBetterSlidingWindows;
};

struct swConf4BetterSlidingWindows {
	enum IncMethod {
		Inc_Auto	=	0,
		Inc_Edge	=	1,
		Inc_Naive	=	-1,
	};
	size_t	buffer_pixel_threshold;
	int		inc_method;
	iPoint	stride;
	swConf4BetterSlidingWindows( const size_t _buffer_pixel_threshold = DEFAULT_BUFFER_PIXEL_THRESHOLD,
		const int _inc_method = Inc_Auto, const iPoint _stride = iPoint(1,1) )
		: buffer_pixel_threshold(_buffer_pixel_threshold), inc_method(_inc_method), stride(_stride) {}
	swConf4BetterSlidingWindows( const swConf4BetterSlidingWindows& _orig ) 
		: buffer_pixel_threshold(_orig.buffer_pixel_threshold), inc_method(_orig.inc_method), stride(_orig.stride) {}
};

class swBetterSlidingWindows {
public:
	typedef swConf4BetterSlidingWindows	Conf;
private:
	piPolyGroupArr	_shapes;
	Conf	_conf;
	piSlidingWindows*	_psw;
public:
	typedef array1d<swUnitStruct>	UnitStructArr;
	UnitStructArr	unit;
	typedef array1d<swSegStruct>	SegStructArr;
	SegStructArr	seg;
	typedef array1d<swPolygonGroup>	WindowArr;
	WindowArr		win;
	int			unitDupFrac4sparse_buffer;
private:
	void _initialize();
public:
	swBetterSlidingWindows( const piPolyGroupArr& shapes, 
		const Conf& conf = Conf() );
	swBetterSlidingWindows( const swBetterSlidingWindows& _orig );
	~swBetterSlidingWindows();
	void ChangeIncMethod( int incMethod = Conf::Inc_Auto );
	iPoint GetStride() const { return _conf.stride; }
	int GetXStride() const { return _conf.stride.x; }
	int GetYStride() const { return _conf.stride.y; }
	size_t GetBufferPixelThreshold() const { return _conf.buffer_pixel_threshold; };
	const piSlidingWindows& OldSlidingWins() const { return *_psw; }
};


/* -------------- BIN MAP Related --------------- */

class piBinMapTraits {
public:
	//array[y][x]
	class DoubleBrackets {
	public:
		template <class BIN_MAP, class bin_t>
		static bin_t read( const BIN_MAP& _bmap, int x, int y ) {
			return _bmap[y][x];
		}
		template <class BIN_MAP>
		static typename BIN_MAP::value_type read( const BIN_MAP& _bmap, int x, int y ) {
			return _bmap[y][x];
		}
	};
	//array(x,y)
	class FunctionClass  {
	public:
		template <class BIN_MAP, class bin_t>
		static bin_t read( const BIN_MAP& _bmap, int x, int y ) {
			return _bmap(x,y);
		}
		template <class BIN_MAP>
		static typename BIN_MAP::value_type read( const BIN_MAP& _bmap, int x, int y ) {
			return _bmap(x,y);
		}
	};
};


/* ------------------ Bin Map -------------- */

class piUnitWeightArr {
public:
	typedef piUnitWeightArr	self_type;
	typedef short int	value_type;
	typedef int	idx_type;
	
	value_type operator () ( idx_type, idx_type ) const {
		return value_type(1);
	}

	class FakeRowType {
	public:
		piUnitWeightArr::value_type operator [] ( idx_type ) const {
			return piUnitWeightArr::value_type(1);
		}
	};
	const FakeRowType operator [] ( idx_type ) const {
		return FakeRowType();
	}
};

template <class BIN_MAP>
class piWeightRowType4BinMap {
	const BIN_MAP&	_bmap;
	int			_row;
public:
	piWeightRowType4BinMap( const BIN_MAP& bmap, int row ) 
		: _bmap(bmap), _row(row) {
	}
	typename BIN_MAP::weight_type operator [] ( int col ) const {
		return _bmap[ iPoint(col,_row) ];
	}
};

//a shell for arr_read
template <class ARR_T, class W_ARR_T = piUnitWeightArr, 
	class T = typename ARR_T::value_type, class W_T = typename W_ARR_T::value_type, 
	class BIN_MAP_TRAITS = typename piBinMapTraits::DoubleBrackets>
class piBinMap : public boost::noncopyable {
public:
	typedef piBinMap<ARR_T, W_ARR_T, T, W_T, BIN_MAP_TRAITS>	self_type;
	typedef self_type	MyType;
	typedef T	value_type;		//generally size_t or unsigned char
	typedef W_T	weight_type;
public:
	const ARR_T&	_arr;
	const W_ARR_T&	_warr;
	uISize _size;
	size_t _binnum;
private:
	piUnitWeightArr _unit_weight;
public:
	explicit piBinMap( const ARR_T& arr, uISize size, size_t bin_num ) 
		: _arr(arr), _warr(_unit_weight), _size(size), _binnum(bin_num) {}
	explicit piBinMap( const ARR_T& arr, const W_ARR_T& warr, const uISize size, size_t bin_num ) 
		: _arr(arr), _warr(warr), _size(size), _binnum(bin_num) {}
	//Parentheses for reading label
	T operator () ( int x, int y ) const {
		return BIN_MAP_TRAITS::read(_arr, x, y);
	}
	T operator () ( iPoint p ) const {
		return (*this)( p.x, p.y );
	}
	//Brackets for reading weight
	typedef piWeightRowType4BinMap<MyType>	WeightRowType;
	WeightRowType operator [] ( int y ) const {
		return WeightRowType( *this, y );
	}
	W_T	operator [] ( iPoint p ) const {
		return BIN_MAP_TRAITS::read( _warr, p.x, p.y );
	}
	const uISize& size() const { return _size; }
	size_t  binnum() const { return _binnum; }
};

//------------------------------------------------------------------------

template <class HIST_FUN>
class piBuffer4HistSlider {
public:
	typedef piBuffer4HistSlider<HIST_FUN>	MyType;
	typedef typename HIST_FUN::bin_type		bin_t;		//generally size_t or unsigned char
	typedef typename HIST_FUN::input_type	count_t;	//int or piFloat
	typedef typename HIST_FUN::output_type	result_t;	//ANY
	typedef piNumericTraits< typename HIST_FUN::input_type >	NumericTraits;
	typedef typename NumericTraits::ZeroScheme	ZeroScheme;
public:
	static const int stepX = 1;		//Can be changed for sparse grid

	typedef piHashHist<count_t, ZeroScheme >		ColumnHist;
	typedef array2d< ColumnHist > 	ColumnHistArr;
	typedef ptrvec<ColumnHistArr> 	ColumnHistArrPVec;

	typedef typename ColumnHist::ListType	ColumnList;
	typedef array2d< ColumnList > 	ColumnListArr;
	typedef ptrvec<ColumnListArr> 	ColumnListArrPVec;
private:
	typedef std::auto_ptr<ColumnHistArrPVec> ColumnHistArrPVec_Ptr;
	typedef std::auto_ptr<ColumnListArrPVec> ColumnListArrPVec_Ptr;
	ColumnHistArrPVec_Ptr	_seg_arr;
	ColumnListArrPVec_Ptr	_unit_arr;
	size_t _map_width;
	size_t _labelPerPixel;
	size_t _binNum;
	const swBetterSlidingWindows& _windows;
private:
	void _initialize(size_t map_width, size_t binNum, size_t labelPerPixel);
public:
	piBuffer4HistSlider( const swBetterSlidingWindows& windows, const size_t binNum = 0, const size_t labelPerPixel = 1 ) :
		_map_width(0), _labelPerPixel(binNum), _binNum(binNum), _windows(windows) {}
	void initialize( const size_t map_width, const size_t binNum = 0, const size_t labelPerPixel = 1 ) {
		_initialize( map_width, binNum, labelPerPixel );
	}
	void release() {
		_seg_arr.reset();	_unit_arr.reset();
		_map_width = 0;		_binNum    = 0;
	}
	ColumnHistArrPVec& SegArr()  { return *_seg_arr; }
	ColumnListArrPVec& UnitArr() { return *_unit_arr; }
	const swBetterSlidingWindows& Windows() { return _windows; }
	MyType* GetMe( const size_t map_width, const size_t binNum = 0, const size_t labelPerPixel = 1 ) {
		_initialize( map_width, binNum, labelPerPixel );
		return this;
	}
};


template <class HIST_FUN,  class BIN_MAP>
class piHistSlider {
public:
	typedef piHistSlider<HIST_FUN,BIN_MAP>	MyType;
	typedef BIN_MAP		BinMap;
	typedef typename HIST_FUN::bin_type		bin_t;		//generally size_t or unsigned char
	typedef typename BIN_MAP::weight_type	weight_t;
	typedef typename HIST_FUN::input_type	count_t;	//piFloat
	typedef typename HIST_FUN::output_type	result_t;	//ANY
	typedef piNumericTraits< typename HIST_FUN::input_type >	NumericTraits;
	typedef typename NumericTraits::ZeroScheme	ZeroScheme;
	typedef piBuffer4HistSlider<HIST_FUN>	MyBuffer;
public:
	typedef array1d<result_t>		ResultArr;
	typedef array2d<result_t>		AllResultArr;
	typedef array1d<iPoint>			PixelArr;
	typedef swPolygonGroup::UPixelArr	UPixelArr;
	typedef swPolygonGroup::PArrPPair	PArrPPair;

	typedef piHashHist<count_t, ZeroScheme >		ColumnHist;
	typedef array2d< ColumnHist > 	ColumnHistArr;
	typedef ptrvec<ColumnHistArr> 	ColumnHistArrPVec;

	typedef typename ColumnHist::ListType	ColumnList;
	typedef array2d< ColumnList > 	ColumnListArr;
	typedef ptrvec<ColumnListArr> 	ColumnListArrPVec;

	typedef array1d<count_t>		InternalHist;
	typedef array1d<InternalHist>	InternalHistArr;

	typedef piHashSpecialBase<count_t, ZeroScheme >		ColumnHistBase;
	typedef typename ColumnHistBase::TableType			ColumnHistTable;

	typedef bool					TagType;
	typedef array1d<bool>			TagArr;
	typedef ptrvec<TagArr>			TagArrPVec;

	typedef array1d<BIN_MAP>		BinMapArr;

public:
	static const int stepX = piBuffer4HistSlider<HIST_FUN>::stepX;		//Can be changed for sparse grid
private:
	const BIN_MAP&  		_binmap;
	const BIN_MAP*  		_binmaps;
	const size_t			_binmap_num;
	HIST_FUN*  _hist_fun;
	const swBetterSlidingWindows& _windows;

	ColumnHist			_buffer_hist;	//temporarily of no use

	InternalHistArr		_tmp_hist;
	InternalHistArr		_tmp_hist_r;
	InternalHistArr		_fixed_hist;
	ColumnHistArr		_v_boundary;
	ColumnHistArr		_h_boundary;
	int					_rbIdx;
	int					_rb_h_Idx;
	ColumnHist*			_fixed_v_boundary;

	ColumnHistTable		_tmp_hashtable;

	array1d<iISize>		_slide_range;
	array1d<iISize>		_result_size;
	iISize				_max_slide_range;
	ResultArr			_fixed_results;
	ResultArr			_fun_results;
	AllResultArr		_dim_results;
	ResultArr			_fun_results_r;
	AllResultArr		_dim_results_r;
	iPoint 				_current;
	iPoint				_currentRP;	//Result position

	typedef std::auto_ptr<MyBuffer> MyBuffer_Ptr;
	MyBuffer_Ptr		_arr_buffer_holder;
	MyBuffer*			_arr_buffer;
	ColumnHistArrPVec&	_seg_arr;
	ColumnListArrPVec&	_unit_arr;

	TagArrPVec			_seg_arr_tag;
	TagArrPVec			_unit_arr_tag;

	TagArr	_result_mask;
	size_t	_avail_result_num;
public:
	const iISize&				MaxSlidingRange;
	const ResultArr&			Results;
	const TagArr&				ResultMask;
	const size_t&				AvailResultNum;
	const InternalHistArr&		HistFixed;
	const ResultArr&			ResultsFixed;
	const array1d<iISize>&		SlidingRange;
	const array1d<iISize>&		ResultSize;
private:
	void _initialize();
	const ColumnList& _get_unitseg( swUnitStruct* u, iPoint p );
	const ColumnHist& _get_seg( swSegStruct* s, iPoint p );
	void _inc_fun_refresh( result_t& r, result_t* dim_r, InternalHist& hist, const ColumnHist& inc  );
	void _full_cal_fun( InternalHistArr& H, ResultArr& ra, AllResultArr* dim_ra_p=NULL );
	void _add_bodyInc2hist( const PArrPPair& bodyInc, iPoint p, ColumnHist& h );
	void _add_boundary2hist( const UPixelArr& boundary, iPoint p, ColumnHist& h  );
	void _fixed_position_hist( iPoint p, InternalHistArr& H, ColumnHist* BoundaryHist );
	void _initial_slide();
	void _shadow_first_column_results();
	void _row_slide();

	void _column_slide_naiveInc( iPoint ps, const swPolygonGroup& w, ColumnHist& oBoundary, ColumnHist& nBoundary );
	void _column_slide_edgeInc( iPoint ps, const swPolygonGroup& w );
	void _column_slide();
private:
	typedef void (MyType::*SLIDE_INTERNAL) ();
	SLIDE_INTERNAL _column_slide_c;
public:
	explicit piHistSlider(const BIN_MAP& binmap, HIST_FUN* hist_fun, const swBetterSlidingWindows& windows, MyBuffer* buffer = NULL );
	explicit piHistSlider(const BinMapArr& binmaps, HIST_FUN* hist_fun, const swBetterSlidingWindows& windows, MyBuffer* buffer = NULL );
	~piHistSlider() {}
	void histAt( iPoint p ) {
		_fixed_position_hist( p, _fixed_hist, _fixed_v_boundary );
	}
	void histAt( int x, int y ) { histAt(iPoint(x,y));}
	void funAt( iPoint p ) {
		histAt(p);
		_full_cal_fun( _fixed_hist, _fixed_results );
	}
	void funAt( int x,int y ) { funAt(iPoint(x,y)); }
	const iPoint& current() const { return _current; }
	const iPoint& currentRP() const { return _currentRP; }
	int stepSlide();
	//-------------------------------------------------------
public:
	typedef array2d<result_t>		ResultOutput;
	typedef ptrvec<ResultOutput>	ResultOutputPVec;
private:
	void	_generateResultAllocation ( ResultOutputPVec& rp );
	void	_generateResultAllocation ( ResultOutputPVec& rp, int );
	void	_generateFastResults(  ResultOutputPVec& r);
	void	_generateNaiveResults( ResultOutputPVec& r );
public:
	void	generateFastResults( ResultOutputPVec& rp );
	void	generateFastResults( ResultOutputPVec& rp, int );		//Without pre-construction
	void	generateNaiveResults( ResultOutputPVec& rp );
	void	generateNaiveResults( ResultOutputPVec& rp, int );		//Without pre-construction

/*--------------------------------------------------------
 * FOR Integral Image!!!!
 *--------------------------------------------------------
*/

private:		//Change to private
	typedef array2d< count_t >		 IntegralType1;		//for slope and bin
	typedef array2d< IntegralType1 > IntegralType2;		//for x and y
	typedef std::auto_ptr<IntegralType2> IMapPtr;
	IMapPtr		_imap;
	int			_intAtY;
	int			_intZeroSlope;
	int			_intInfSlope;
	iISize		_intSize;
	int			_intD;
private:		//Change to private
	void		_makeIntegral( int Y );
	void		_intHistInternal( iPoint p );
public:
	void InitializeIntegralImage();
	int IntegralSlide();
	void	generateIntegralResults( ResultOutputPVec& rp );
	void	generateIntegralResults( ResultOutputPVec& rp, int );		//Without pre-construction
private:
	void	_generateIntegralResults( ResultOutputPVec& rp );
};

//----------------------------------------------------

#include "sliding_hist_t.hpp"

#endif
