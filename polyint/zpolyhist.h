
#ifndef Z_POLYHIST_H_
#define Z_POLYHIST_H_

#include <cctype>
#include <vector>

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((visibility("default")))
    #else
      #define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #else
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((visibility("default")))
    #else
      #define DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
	#endif
  #endif
  #define DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define DLL_PUBLIC __attribute__ ((visibility("default")))
    #define DLL_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define DLL_PUBLIC
    #define DLL_LOCAL
  #endif
#endif


namespace zRawArray2D {

template < class ELEMENT_T >
struct Array2D;

template < class ELEMENT_T >
class _Array2DReleaser {
public:
	typedef Array2D<ELEMENT_T>	Array2DType;
	static void CascadeRelease( Array2DType& _p  );
};

template < class SUB_ELEMENT_T >
class _Array2DReleaser< Array2D<SUB_ELEMENT_T> > {
public:
	typedef Array2D<SUB_ELEMENT_T>	SubArray2DType;
	typedef Array2D<SubArray2DType>	Array2DType;
	static void CascadeRelease( Array2DType& _p  );
};


template < class ELEMENT_T >
struct Array2D {
	typedef Array2D<ELEMENT_T>	self_type;
	typedef ELEMENT_T			value_type;

	value_type*	data;
	size_t		width;
	size_t		height;

	Array2D() { }
	Array2D( value_type*	_data, size_t _width, size_t _height )
		: data(_data), width(_width), height(_height) { }
	Array2D( const self_type& r )
		: data(r.data), width(r.width), height(r.height) { }
	virtual ~Array2D() {}

	const self_type& operator = ( const self_type& r ) {
		data   = r.data;
		width  = r.width;
		height = r.height;
		return *this;
	}

	size_t size() const {
		return width * height;
	}

	value_type* operator [] ( size_t r_idx ) {
		return data+r_idx*width;
	}

	value_type* operator [] ( size_t r_idx ) const {
		return data+r_idx*width;
	}

	value_type& operator () ( size_t c_idx, size_t r_idx ) {
		return (*this)[r_idx][c_idx];
	}

	value_type& operator () ( size_t c_idx, size_t r_idx ) const {
		return (*this)[r_idx][c_idx];
	}

	value_type& operator () ( size_t idx ) {
		return data[idx];
	}

	value_type& operator () ( size_t idx ) const {
		return data[idx];
	}

	bool operator == ( const self_type& r ) const {
		return ( data==r.data && width == r.width && height == r.height );
	}

	static self_type Create( size_t _width, size_t _height ) {
		return self_type( new value_type[_width*_height], _width, _height );
	}

	static self_type Create( size_t _width ) {
		return Create( _width, 1 );
	}

	static void Release( self_type& _p ) {	//Can't release const object
		if ( _p.data )
			delete[] _p.data;
	}

	void Release() {
		Release( *this );
	}

	static void CascadeRelease( self_type& _p ) {
		_Array2DReleaser<value_type>::CascadeRelease( _p );
	};

	void CascadeRelease() {
		CascadeRelease( *this );
	}

	self_type clone_content() const {
		self_type l = Create( this->width, this->height );
		value_type* _eptr = this->data + this->width*this->height;
		value_type* p = this->data,
			      * q = l.data;
		while ( p!=_eptr )
			*q++ = *p++;
		return l;
	}

};

template < class ELEMENT_T >
void _Array2DReleaser<ELEMENT_T>::CascadeRelease( Array2D<ELEMENT_T>& _p ) {
	Array2DType::Release(_p);
}

template < class SUB_ELEMENT_T >
void _Array2DReleaser< Array2D<SUB_ELEMENT_T> >::CascadeRelease ( Array2D< Array2D<SUB_ELEMENT_T> >& _p ) {
	for ( size_t y=0; y<_p.height; ++y ) {
		for ( size_t x=0; x<_p.width; ++x ) {
			SubArray2DType::Release( _p[y][x] );
		}
	}
	Array2DType::Release(_p);
}

template < class ELEMENT_T >
void ReleaseArray2D( Array2D<ELEMENT_T>& _p ) {
	Array2D<ELEMENT_T>::Release(_p);
}

template < class ELEMENT_T >
void CascadeReleaseArray2D( Array2D<ELEMENT_T>& _p ) {
	Array2D<ELEMENT_T>::CascadeRelease(_p);
}


}


namespace zPolyHist {

	typedef double			FloatType;
	typedef unsigned int	UIntType;

	struct SparseHistElement {
		UIntType	idx;
		FloatType	count;
	};

	using namespace zRawArray2D;

	typedef Array2D<UIntType>		UIntArray2D;
	typedef Array2D<FloatType>		FloatArray2D;
	typedef Array2D<FloatArray2D>	OutputArr;
	typedef std::vector<SparseHistElement>	SparseHistType;
	typedef Array2D<SparseHistType>	SparseHistArr;
	typedef Array2D<SparseHistArr>	SparseHistOutputArr;
	typedef Array2D<UIntArray2D>	ArrayOfUIntArray2D;
	typedef Array2D<FloatArray2D>	ArrayOfFloatArray2D;

	// A temporary definition for point
	// Might be replaced by the utility set

	template< class TYPE >
	struct PointTemplate {
		TYPE x;
		TYPE y;

		PointTemplate() {}
		PointTemplate(TYPE _x, TYPE _y) : x(_x), y(_y) {}
		PointTemplate(const PointTemplate& _p ) : x(_p.x), y(_p.y) {}
		
		const PointTemplate& operator = ( const PointTemplate& r ) {
			x = r.x; y = r.y;
			return *this;
		}

		bool operator == ( const PointTemplate& r ) {
			return (x==r.x) && (y==r.y);
		}
	};

	typedef PointTemplate<int> IntPoint;

	//--------------------------------------------

	class OutputFactory;

	struct ShapesContent;
	class DLL_PUBLIC Shapes {
		ShapesContent* _content;
	public:
		class Exception_OperateOnNoncompletedShapes {};
		class Exception_EditCompletedShapes {};
		Shapes();
		Shapes( const Shapes& r );
		const Shapes& operator = ( const Shapes& r );
		~Shapes();
		void AddVertex2Current( IntPoint p );
		void NewPolygon();
		void NewShape();
		void Complete();
		bool IsCompleted() const;
		size_t size() const;	//How many shape
		size_t size( size_t idx ) const;	//How many polygons in the idx shape
		IntPoint lefttop( size_t idx ) const;
		IntPoint rightbottom( size_t idx ) const;

		template<class InputIterator>
		void AddMultipleVertice2Current( InputIterator _start, InputIterator _end ) {
			for ( InputIterator iter = _start; iter!=_end; ++iter ) {
				this->AddVertex2Current( *iter );
			}
		}

		template<class InputIterator>
		Shapes( InputIterator _start, InputIterator _end ) {
			AddMultipleVertice2Current( _start, _end );
			Complete();
		}

		void SetMargin4CurrentShape( int top, int right, int bottom, int left );

		friend class OutputFactory;

	private:

		bool operator == ( const Shapes& r );
	};

	//Histogram function ===========================
	class HistFunBase {
	public:
		typedef	UIntType  bin_type;
		typedef FloatType input_type;
		typedef FloatType output_type;
		HistFunBase() {}
		//Do not support custom ZERO
		virtual output_type operator () (bin_type b, input_type n) = 0;
		virtual ~HistFunBase() {};
	};

	class NullHistFun : public HistFunBase {
	public:
		typedef	UIntType  bin_type;
		typedef FloatType input_type;
		typedef FloatType output_type;
		NullHistFun() {}
		//Do not support custom ZERO
		virtual output_type operator () (bin_type b, input_type n) { return 0.0; };
		virtual ~NullHistFun() {};
	};


	//Output factory ===============================

	struct OutputFactoryContent;
	class DLL_PUBLIC OutputFactory {
		OutputFactoryContent* _content;
	public:
		typedef std::vector<size_t> SizeTypeVec;
		enum Methods {
			METHOD_FAST_AUTO  = 0,
			METHOD_FAST_EDGE  = 1,
			METHOD_FAST_NAIVE = 2,
			METHOD_NORMAL     = -2,
			METHOD_INTEGRAL   = -5
		};
		OutputFactory( const Shapes& shapes, HistFunBase& fun, size_t bin_num, IntPoint stride = IntPoint(1,1) );
		~OutputFactory();
		const OutputArr ComputeResponses( const UIntArray2D label_map, int method = METHOD_FAST_AUTO ) const;
		const OutputArr ComputeResponses( const ArrayOfUIntArray2D label_maps, int method = METHOD_FAST_AUTO ) const;
		const OutputArr ComputeResponses( const UIntArray2D label_map, const FloatArray2D weight_map, int method = METHOD_FAST_AUTO ) const;
		const OutputArr ComputeResponses( const ArrayOfUIntArray2D label_maps, const ArrayOfFloatArray2D weight_maps, int method = METHOD_FAST_AUTO ) const;

		const SparseHistOutputArr ComputeHistogram( const UIntArray2D label_map, int method = METHOD_FAST_AUTO, bool need_sorting = false  ) const;
		const SparseHistOutputArr ComputeHistogram( const ArrayOfUIntArray2D label_maps, int method = METHOD_FAST_AUTO, bool need_sorting = false  ) const;
		const SparseHistOutputArr ComputeHistogram( const UIntArray2D label_map, const FloatArray2D weight_map, int method = METHOD_FAST_AUTO, bool need_sorting = false  ) const;
		const SparseHistOutputArr ComputeHistogram( ArrayOfUIntArray2D label_maps, const ArrayOfFloatArray2D weight_maps, int method = METHOD_FAST_AUTO, bool need_sorting = false  ) const;
		const SizeTypeVec& HistgramNonzeroCount;

		void SetHistFun( HistFunBase& fun );
		void SetBinNumber( size_t bin_num );
	};
}

#endif
