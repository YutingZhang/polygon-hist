/*
 * sliding_hist_t.hpp
 *
 *  Created on: 2010-10-28
 *      Author: zhangyuting
 */

#include <algorithm>
#include <cmath>

template <class HIST_FUN>
void piBuffer4HistSlider<HIST_FUN>::_initialize(size_t map_width, size_t binNum, size_t labelPerPixel) {
	if (map_width<=_map_width && binNum<=_binNum && labelPerPixel<=_labelPerPixel )
		return;

	using namespace std;

	_labelPerPixel = labelPerPixel;
	_binNum = binNum;
	_map_width = map_width;

	//LineSeg Buffers
	_seg_arr = ColumnHistArrPVec_Ptr(new ColumnHistArrPVec);
	for ( size_t j=0; j<_windows.seg.size(); ++j ) {
		const swSegStruct& s = _windows.seg[j];
		size_t nPixel = min(s.unit->pixels.size()*s.frac*_labelPerPixel, _binNum);
		ColumnHist	zero_hist(_binNum, nPixel );
		int Xspan = ((int)_map_width-s.Xispan.e2)-(s.Xispan.e1+_windows.GetXStride())+1;
		int Xsize = (Xspan-1) / _windows.GetXStride()+1;
		int Ysize;
		switch ( s.refreshType ) {
			case swSegStruct::SS_Buffer:
			case swSegStruct::SS_NoBuffer:
				Ysize = s.Yspan;
				break;
			default:
				Ysize = 1;
				break;
		}
		ColumnHistArr* ha = new ColumnHistArr( uISize( Xsize, Ysize ),
			zero_hist, typename ColumnHistArr::WithFixedConstructor() );
		s.buffer = ha;
		_seg_arr->push_back( ha );
	}


	//UnitLineSeg Buffers
	_unit_arr = ColumnListArrPVec_Ptr(new ColumnListArrPVec);
	for ( size_t j=0; j<_windows.unit.size(); ++j ) {
		const swUnitStruct& u = _windows.unit[j];
		bool  need_buffer = false;
		for ( size_t k=0; k<u.seg.size(); ++k ) {
			const swSegStruct& s = *(u.seg[k]);
			if ( s.refreshType == swSegStruct::SS_Buffer ) {
				need_buffer = true;
				break;
			}
		}
		if ( !need_buffer )
			continue;
		size_t nPixel = min(u.pixels.size()*_windows.unitDupFrac4sparse_buffer*_labelPerPixel,_binNum);

		int Xspan = (map_width - u.Xispan.e2) - (u.Xispan.e1 + _windows.GetXStride()) + 1;
		int Xsize = (Xspan-1)/_windows.GetXStride() + 1;

		ColumnListArr* ha = new ColumnListArr( uISize(Xsize,u.Yspan),
			nPixel, typename ColumnListArr::WithFixedConstructor() );
		_unit_arr->push_back( ha );
		u.buffer = ha;
	}
}

//-------------------------------------------------------------------------------

#define PI_HISTSLIDER_CONSTRUCT_FUNC \
		  _hist_fun(hist_fun), _windows(windows),	\
		  _buffer_hist( _binmap.binnum() ),	\
		  _tmp_hist		( _windows.win.size(), _binmap.binnum(), typename InternalHistArr::WithFixedConstructor() ),	\
		  _tmp_hist_r	( _windows.win.size(), _binmap.binnum(), typename InternalHistArr::WithFixedConstructor() ),	\
		  _fixed_hist	( _windows.win.size(), _binmap.binnum(), typename InternalHistArr::WithFixedConstructor() ),	\
		  _v_boundary	( uISize(_windows.win.size(),3), _binmap.binnum(), typename ColumnHistArr::WithFixedConstructor() ),	\
		  _h_boundary	( uISize(_windows.win.size(),2), _binmap.binnum(), typename ColumnHistArr::WithFixedConstructor() ),	\
		  _rbIdx(0), _rb_h_Idx(0),	\
		  _fixed_v_boundary( _v_boundary[2] ),	\
			\
		  _tmp_hashtable( _binmap.binnum() ),	\
			\
		  _slide_range(_windows.win.size()),	\
		  _result_size(_windows.win.size()),	\
			\
		  _fixed_results( _windows.win.size(), hist_fun->zero(), typename ResultArr::WithFixedConstructor() ),	\
		  _fun_results( _windows.win.size(), hist_fun->zero(), typename ResultArr::WithFixedConstructor() ),	\
		  _dim_results( uISize(_binmap.binnum(),_windows.win.size()), hist_fun->zero(), typename AllResultArr::WithFixedConstructor() ),	\
		  _fun_results_r( _windows.win.size(), hist_fun->zero(), typename ResultArr::WithFixedConstructor() ),	\
		  _dim_results_r( uISize(_binmap.binnum(),_windows.win.size()), hist_fun->zero(), typename AllResultArr::WithFixedConstructor() ),	\
			\
		  _arr_buffer_holder( (buffer)?(NULL):(new MyBuffer(_windows)) ),		\
		  _arr_buffer( (buffer)?( buffer->GetMe(_binmap.size().width, _binmap.binnum()) ):		\
				  (_arr_buffer_holder->GetMe(_binmap.size().width, _binmap.binnum()) ) ),		\
			\
		  _seg_arr( _arr_buffer->SegArr() ),	\
		  _unit_arr( _arr_buffer->UnitArr() ),	\
			\
		  _result_mask( _windows.win.size() ),	\
			\
		  MaxSlidingRange( _max_slide_range ),	\
		  Results(_fun_results), 	\
		  ResultMask(_result_mask), 	\
		  AvailResultNum(_avail_result_num),	\
		  HistFixed(_fixed_hist),	\
		  ResultsFixed(_fixed_results),	\
		  SlidingRange(_slide_range),	\
		  ResultSize(_result_size),	\
			\
		  _column_slide_c(&MyType::_initial_slide){	\
			_tmp_hashtable.set( typename ColumnHistBase::iterator(NULL) );	\
			_initialize();	\
	}

template <typename HIST_FUN,  typename BIN_MAP>
piHistSlider<HIST_FUN,BIN_MAP>::piHistSlider(const BIN_MAP& binmap, HIST_FUN* hist_fun, const swBetterSlidingWindows& windows, MyBuffer* buffer)
	: _binmap(binmap), _binmaps(&binmap), _binmap_num(1),
	  PI_HISTSLIDER_CONSTRUCT_FUNC;

template <typename HIST_FUN,  typename BIN_MAP>
piHistSlider<HIST_FUN,BIN_MAP>::piHistSlider(const BinMapArr& binmaps, HIST_FUN* hist_fun, const swBetterSlidingWindows& windows, MyBuffer* buffer)
	: _binmap(binmaps[0]), _binmaps(binmaps), _binmap_num( binmaps.size() ),
	  PI_HISTSLIDER_CONSTRUCT_FUNC;


template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_initialize() {

	
	//Initialize the Sliding range
	_max_slide_range = iISize(0,0);
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		//_slide_range[i] = iISize(_binmap.size() - _windows.win[i].size - (_windows.GetStride()-1) );
		_slide_range[i] = iISize(_binmap.size() - _windows.win[i].size );
		PI_CONDITIONAL_REFRESH( _max_slide_range.width , _slide_range[i].width , < );
		PI_CONDITIONAL_REFRESH( _max_slide_range.height, _slide_range[i].height, < );
		if ( _slide_range[i].width<0 || _slide_range[i].height<0 ) {
			_result_size[i] = iISize(0,0);
		} else {
			_result_size[i] = iISize(0,0)+(iPoint(_slide_range[i].width,_slide_range[i].height)/_windows.GetStride()+1);
		}

	}

	//LineSeg Buffers
	for ( size_t j=0; j<_windows.seg.size(); ++j ) {
		const swSegStruct& s = _windows.seg[j];
		s.slideSpan.e1 = s.Xispan.e1 + _windows.GetXStride();
		s.slideSpan.e2 = _binmap.size().width - s.Xispan.e2;
		int Xspan = s.slideSpan.e2-s.slideSpan.e1+1;
		s.Xcycle = (Xspan-1) / _windows.GetXStride()+1;
		// No need to Plus 1, because I use the left-top end as the index, and y length is at least 1
		TagArr* ta = new TagArr(_binmap.size().height, false, TagArr::WithFixedConstructor());
		s.tag_buffer = ta;
		_seg_arr_tag.push_back( ta );
	}

	//UnitLineSeg Buffers
	for ( size_t j=0; j<_windows.unit.size(); ++j ) {
		const swUnitStruct& u = _windows.unit[j];
		u.slideSpan.e1 = u.Xispan.e1 + _windows.GetXStride();
		u.slideSpan.e2 = _binmap.size().width - u.Xispan.e2;
		TagArr* ta = new TagArr( _binmap.size().height, false, TagArr::WithFixedConstructor() );
		_unit_arr_tag.push_back( ta );
		u.tag_buffer = ta;
	}

}

template<typename T1, typename T2>
inline bool IsInMap( const zUtilities_PolyHistInternal::zPointT<T1>& p, const zUtilities_PolyHistInternal::zISizeT<T2>& size ) {
	/* x>=XStep (=1), y doesn't need to >=SomeStep, but >=0 */
	return ( p.x>0 && p.x<=(T1)size.width && p.y>=0 && p.y<=(T1)size.height );
}

template <typename HIST_FUN,  typename BIN_MAP>
const typename piHashHist<typename HIST_FUN::input_type, typename piNumericTraits< typename HIST_FUN::input_type >::ZeroScheme >::ListType&
piHistSlider<HIST_FUN,BIN_MAP>::_get_unitseg( swUnitStruct* u, iPoint p ) {
	ColumnListArr& larr = *(ColumnListArr*)(u->buffer);
	TagType& tag = (*(TagArr*)(u->tag_buffer))[p.y];
	size_t cycle_height = larr.size().height;
	size_t cy = (p.y/_windows.GetYStride()) % cycle_height;

	if ( !tag ) {
		tag = true;
		//Refresh row

		ColumnList* lrow = larr[cy];
		iPoint vec = u->vec;
		int xi = 0;
		for ( int cx = u->slideSpan.e1; cx<=u->slideSpan.e2; cx+=_windows.GetXStride() ) {
			ColumnList& l = lrow[xi++];
			l.clear();

			iPoint p1 = iPoint(cx,p.y);
			{
				ColumnHistBase hh( _tmp_hashtable, l );	//The table & list consistency is gauranteed
				for ( int i=0; i<_windows.unitDupFrac4sparse_buffer; ++i ) {
					for (swUnitStruct::UnitPixelArr::iterator iter = u->pixels.begin();
						iter != u->pixels.end(); ++iter ) {
							iPoint pp = iter->pos+p1;
							for ( size_t v=0; v<_binmap_num ; ++v ) {
								bin_t b=_binmaps[v](pp);
								weight_t	wi = _binmaps[v][pp];
								hh.inc( b, iter->por*wi );
							}
					}
					p1 = p1 + vec;
				}
				hh.clearTab();
			}
		}
	}
	return larr[cy][(p.x-u->slideSpan.e1)/_windows.GetXStride()];
}

inline int pMod( int a, int b) {
	int c = a % b;
	if (c<0) c+=b;
	return c;
}

template <typename HIST_FUN,  typename BIN_MAP>
const piHashHist<typename HIST_FUN::input_type, typename piNumericTraits< typename HIST_FUN::input_type >::ZeroScheme >&
piHistSlider<HIST_FUN,BIN_MAP>::_get_seg( swSegStruct* s, iPoint p ) {
	ColumnHistArr& harr = *(ColumnHistArr*)(s->buffer);
	TagType& tag = (*(TagArr*)(s->tag_buffer))[p.y];
	size_t cycle_height = harr.size().height;

	size_t stored_cycle_width = s->Xcycle;	//harr.size().width;
	int cy = (p.y/_windows.GetYStride()) % cycle_height;	//cycle_height == uvec.y

	iPoint uvec = s->unit->vec;
	//Compute array x from image x
	//const int x_start = stepX;
	int deltaX = (p.y/s->inc_vec.y) * (s->inc_vec.x/_windows.GetXStride());
	int startSpan = _windows.GetXStride()+s->Xispan.e1;
	int cx1 = pMod( (p.x-startSpan)/_windows.GetXStride()-deltaX, stored_cycle_width );

	if ( !tag ) {
		tag = true;

		//Refresh row
		swUnitStruct* u = s->unit;
//--------------------------------------------------------------------------------------------
#define NO_BUFFER_INITIAL \
	hh.clear();		\
	iPoint pt(p1);	\
	for ( int i=0; i<s->frac; ++i ) {	\
		for (swUnitStruct::UnitPixelArr::iterator iter = u->pixels.begin();	\
			iter != u->pixels.end(); ++iter ) {	\
				iPoint pp = iter->pos+pt;	\
				for ( size_t v=0; v<_binmap_num ; ++v ) {		\
					bin_t  b  = _binmaps[v](pp);	\
					weight_t	wi = _binmaps[v][pp];	\
					hh.inc( b, iter->por*wi );	\
				}	\
		}	\
		pt = pt + uvec;	\
	}
//-----------------------------------------------------------
		
		ColumnHist* hrow = harr[cy];

		int px = pMod(deltaX, stored_cycle_width)*_windows.GetXStride()+startSpan;
		int pxMax = (stored_cycle_width-1)*_windows.GetXStride()+startSpan;

		for ( size_t cx = 0; cx<stored_cycle_width; ++cx ) {
			//--------------------------
			ColumnHist& hh = hrow[cx];
			//int px = pMod(cx+deltaX, stored_cycle_width)*_windows.GetXStride()+startSpan;
			iPoint p1 = iPoint(px,p.y);
			iPoint p2 = p1 + s->vec;
			iPoint p0( p1-uvec * _windows.unitDupFrac4sparse_buffer );
			iPoint p15(p2-uvec * _windows.unitDupFrac4sparse_buffer );
			switch (s->refreshType) {
				case swSegStruct::SS_Vertical:
					{
						int x0 = p1.x-_windows.GetXStride();
						if ( p0.y>=s->lowY ) {
							int y15= p2.y - _windows.GetYStride();
							//Incremental refresh
							for ( int y=p1.y - _windows.GetYStride(); y<p1.y; ++y ) {
								//dec
								bin_t b;
								weight_t	wi;
								for ( int x = x0; x<p1.x; ++x ) { 
									for ( size_t v=0; v<_binmap_num ; ++v ) {
										b  = _binmaps[v](x,y);
										wi = _binmaps[v][y][x];
										hh.dec(b,wi);
									}
								}
								//inc
								for ( int x = x0; x<p1.x; ++x ) { 
									for (size_t v = 0; v < _binmap_num; ++v) {
										b = _binmaps[v](x, y15);
										wi = _binmaps[v][y15][x];
										hh.inc(b, wi);
									}
								}
								++y15;
							}
						} else {
							//Initialize;
							hh.clear();
							int y = p1.y;
							int Y = y+s->frac;
							while ( y<Y ) {
								for ( int x = x0; x<p1.x; ++x ) {
									for (size_t v = 0; v < _binmap_num; ++v) {
										bin_t b = _binmaps[v](x,y);
										weight_t wi = _binmaps[v][y][x];
										hh.inc(b,wi);
									}
								}
								++y;
							}
						}
					}
					break;
				case swSegStruct::SS_NoBuffer:
					{
						swUnitStruct* u = s->unit;
						
						if ( (p0.x>=s->slideSpan.e1 && p0.x<=s->slideSpan.e2) && p0.y>=s->lowY ) {
							iPoint pt0(p0), pt15(p15);
							for ( int i=0; i<_windows.unitDupFrac4sparse_buffer; ++i ) {
								//dec
								for (swUnitStruct::UnitPixelArr::iterator iter = u->pixels.begin();
									iter != u->pixels.end(); ++iter ) {
									for (size_t v = 0; v < _binmap_num; ++v) {
										iPoint pp = iter->pos+pt0;
										bin_t  b  = _binmaps[v](pp);
										weight_t	wi = _binmaps[v][pp];
										hh.dec( b, iter->por*wi );
									}
								}
								//inc
								for (swUnitStruct::UnitPixelArr::iterator iter = u->pixels.begin();
									iter != u->pixels.end(); ++iter ) {
									for (size_t v = 0; v < _binmap_num; ++v) {
										iPoint pp = iter->pos+pt15;
										bin_t  b  = _binmaps[v](pp);
										weight_t	wi = _binmaps[v][pp];
										hh.inc( b, iter->por*wi );
									}
								}
								pt0 = pt0 + uvec;
								pt15= pt15+ uvec;
							}
						} else {
							//initialize
							NO_BUFFER_INITIAL
						}
					}
					break;
				case swSegStruct::SS_UnitSeg:
					{
						NO_BUFFER_INITIAL
					}
					break;
				//case swSegStruct::SS_Buffer:
				default:
					if ( (p0.x>=s->slideSpan.e1 && p0.x<=s->slideSpan.e2) && p0.y>=s->lowY ) {
						//Incremental refresh
						const ColumnList& _dec = _get_unitseg( s->unit, p0 );
						const ColumnList& _inc = _get_unitseg( s->unit, p15 );
						//The order of _dec and _inc couldn't be swapped, 
						//because the preallocated memory is not enough for _inc first
						hh-=_dec;
						hh+=_inc;
					} else {
						//Initialize;
						hh.clear();
						iPoint pp(p1);
						for ( int i=0; i<s->fracY4sparse_buffer; ++i ) {
							const ColumnList& _inc = _get_unitseg( s->unit, pp );
							hh += _inc;
							pp = pp + uvec*_windows.unitDupFrac4sparse_buffer;
						}
					}
			}
			px += _windows.GetXStride();
			if (px > pxMax) {
				px = startSpan;
			}
		}

	}
	return harr[cy][cx1];
}


inline bool IsPointLessThanRange( iPoint p, iISize size ) {
	return (p.x<=size.width && p.y <=size.height);
}

//---------------------------------DEBUG
#if 0
	template <typename T>
	inline void PrintNonZero( T, size_t b, size_t i ) {
		;
	}

	template <>
	inline void PrintNonZero<piFloat>( piFloat a, size_t b, size_t i  ) {
		if (a!=0.0)
			std::cerr << "!" << i << "\t" << b << "\t" << a << std::endl;
	}
#endif
//-------------------------------------


template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_full_cal_fun( InternalHistArr& H, ResultArr& ra, AllResultArr* dim_ra_p ) {
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		if (!_result_mask[i]) 
			continue;
		const InternalHist& h = H[i];

		//calculate fun
		result_t& f = ra[i];
		f = _hist_fun->zero();
		if (dim_ra_p) {
			result_t* dr_r  =(*dim_ra_p)[i];
			for ( bin_t b=0; b<_binmap.binnum() ; ++b ) {
				dr_r[b] = _hist_fun->operator()( b, h[b] );
				f += dr_r[b];
			}
		} else {
			for ( bin_t b=0; b<_binmap.binnum() ; ++b ) {
				f += _hist_fun->operator()( b, h[b] );
			}
		}
	}
}

#define TEST_AVAILABILITY( POINT, IDX, ACC ) 	\
	if ( !IsPointLessThanRange(POINT, _slide_range[IDX]) ) { \
		_fun_results[IDX] = _hist_fun->zero(); \
		_result_mask[IDX] = false;	\
		continue;	\
	} else {	\
		_result_mask[IDX] = true;	\
		++ACC;	\
	}

template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_add_boundary2hist( const UPixelArr& boundary, iPoint p, ColumnHist& h ) {
	for ( UPixelArr::const_iterator iter = boundary.begin(); 
		iter!=boundary.end(); ++iter ) {
			iPoint	pp = iter->pos+p;
			for (size_t v = 0; v < _binmap_num; ++v) {
				bin_t b = _binmaps[v]( pp );
				weight_t wi = _binmaps[v][pp];
				h.inc(b, iter->por*wi);
			}
	}
}

template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_fixed_position_hist( iPoint p, InternalHistArr& H, ColumnHist* BoundaryHist ) {
	_avail_result_num = 0;
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		TEST_AVAILABILITY(p,i,_avail_result_num);

		InternalHist& h = H[i];
		h.set(count_t(0));
		//Interior
		const swPolygonGroup& w = _windows.win[i];

		for ( swPolygonGroup::CLVec::const_iterator iter = w.interiorRL.begin();
			iter != w.interiorRL.end(); ++iter ) {
				int x0 = p.x+iter->span.e1,
					x1 = p.x+iter->span.e2;
				int y = iter->num + p.y;
				for ( int x=x0; x<=x1; ++x ) {
					for (size_t v = 0; v < _binmap_num; ++v) {
						bin_t b = _binmaps[v](x,y);
						weight_t	wi = _binmaps[v][y][x];
						h[b] += wi;
					}
				}
		}

		//Boundary
		ColumnHist& theBoundary = BoundaryHist[i];
		theBoundary.clear();
		_add_boundary2hist( w.boundary, p, theBoundary );
		theBoundary.add2array( h );
	}
}

template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_shadow_first_column_results() {
	ColumnHist* v_boundary = _v_boundary[_rbIdx],
		* h_boudary = _h_boundary[_rb_h_Idx];
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		result_t* dr  =_dim_results[i];
		result_t* dr_r=_dim_results_r[i];
		for ( bin_t b=0; b<_binmap.binnum() ; ++b ) {
			dr[b] = dr_r[b];
		}
		_fun_results[i] = _fun_results_r[i];
		_tmp_hist[i]    = _tmp_hist_r[i];
		h_boudary[i]	= v_boundary[i];
	}
}


template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_initial_slide() {
	_column_slide_c = &MyType::_column_slide;

	_current = iPoint(0,0);
	_currentRP = iPoint(0,0);

	_rbIdx = 0;
	_fixed_position_hist( _current, _tmp_hist_r, _v_boundary[0] );

	_avail_result_num = 0;
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		TEST_AVAILABILITY(_current,i,_avail_result_num);
	}

	_full_cal_fun(_tmp_hist_r, _fun_results_r, &_dim_results_r);

	//-----------------DEBUG
#if 0
	for (size_t i=0; i<_windows.win.size(); ++i) {
		cerr << "No.\t" << i << endl;
		cerr << "dim_r:" << endl;
		for ( bin_t b=0; b<_binmap.binnum() ; ++b ) {
			PrintNonZero(_dim_results_r[i][b],b,i);
		}

		cerr << "th:" << endl;
		int k=0;
		for ( InternalHist::iterator iter=_tmp_hist_r[i].begin();
			iter != _tmp_hist_r[i].end(); ++iter ) {
				if ( *iter != 0.0 )
					cerr << k << "\t" << *iter << endl;
				++k;
		}
	}
#endif
	//------------------

	_shadow_first_column_results();
}

template <typename HIST_FUN,  typename BIN_MAP>
inline void piHistSlider<HIST_FUN,BIN_MAP>::_inc_fun_refresh( result_t& r, result_t* dim_r, InternalHist& hist, const ColumnHist& inc  ) {
	for ( typename ColumnHist::const_iterator iter = inc.begin();
			iter != inc.end(); ++iter ) {
		size_t k = iter->idx;
		r -= dim_r[k];
		r += ( dim_r[k] = _hist_fun->operator()(k, (hist[k]+=iter->val) ) );
	}
}


template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_add_bodyInc2hist( const PArrPPair& bodyInc, iPoint p, ColumnHist& h ) {
	const PixelArr& r1 = *bodyInc.e1,
		r2 = *bodyInc.e2;
	for ( PixelArr::const_iterator iter = r1.begin(); 
		iter != r1.end(); ++iter ) {
			iPoint pp(*iter + p);
			for (size_t v = 0; v < _binmap_num; ++v) {
				bin_t b = _binmaps[v]( pp );
				weight_t	wi = _binmaps[v][pp];
				h.inc_s( b, -wi );
			}
	}
	for ( PixelArr::const_iterator iter = r2.begin(); 
		iter != r2.end(); ++iter ) {
			iPoint pp(*iter + p);
			for (size_t v = 0; v < _binmap_num; ++v) {
				bin_t b = _binmaps[v]( pp );
				weight_t	wi = _binmaps[v][pp];
				h.inc_s( b, wi );
			}
	}
}

template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_row_slide() {

	_current.x = 0;
	_currentRP.x = 0;
	++_currentRP.y;
	_current.y += _windows.GetYStride();

	ColumnHist* oldBoundary = _v_boundary[_rbIdx],
			  * newBoundary = _v_boundary[1-_rbIdx];
	_rbIdx = 1-_rbIdx;

	_avail_result_num = 0;
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		TEST_AVAILABILITY(_current,i,_avail_result_num);
		
		const swPolygonGroup& w = _windows.win[i];
	
		_buffer_hist.clear();
		//Interior Intersect
		_add_bodyInc2hist( w.bodyResidual_V, _current, _buffer_hist );

		//Boundary residual
		_buffer_hist -= oldBoundary[i];
		ColumnHist& nBoundary = newBoundary[i];
		nBoundary.clear();
		_add_boundary2hist( w.boundary , _current, nBoundary );
		_buffer_hist += nBoundary;

		//Incremental estimate
		_inc_fun_refresh( _fun_results_r[i], _dim_results_r[i], _tmp_hist_r[i], _buffer_hist );
	}
	_shadow_first_column_results();
}

template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_column_slide_edgeInc( iPoint ps, const swPolygonGroup& w ) {
	_buffer_hist.clear();
	for ( swPolygonGroup::EdgeArr::const_iterator iter = w.edges.begin();
		iter !=w.edges.end(); ++iter ) {
			if (!iter->line.vec.y)
				continue;
			iPoint p(ps + iter->line.sp);
			const ColumnHist& e = _get_seg( iter->seg, p );
			if (iter->is_right) {
				_buffer_hist += e;
			} else {
				_buffer_hist -= e;
			}
	}
}

template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_column_slide_naiveInc( iPoint ps, const swPolygonGroup& w, ColumnHist& oBoundary, ColumnHist& nBoundary ) {
	_buffer_hist.clear();

	//Interior intersect
	_add_bodyInc2hist( w.bodyResidual_H, ps, _buffer_hist );
	//Boundary residual
	_buffer_hist -= oBoundary;
	nBoundary.clear();
	_add_boundary2hist( w.boundary , ps, nBoundary );
	_buffer_hist += nBoundary;

}

template <typename HIST_FUN,  typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_column_slide() {

	_current.x += _windows.GetXStride();
	++_currentRP.x;

	ColumnHist* oldBoundary = _h_boundary[_rb_h_Idx],
		* newBoundary = _h_boundary[1-_rb_h_Idx];
	_rb_h_Idx = 1-_rb_h_Idx;

	_avail_result_num = 0;
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		TEST_AVAILABILITY(_current,i,_avail_result_num);


		const swPolygonGroup& w = _windows.win[i];


		if (w.is_naiveInc) {
			_column_slide_naiveInc( _current, w, oldBoundary[i], newBoundary[i]);
		} else {
			_column_slide_edgeInc( _current, w );
		}

		//Incremental estimate
		_inc_fun_refresh( _fun_results[i], _dim_results[i], _tmp_hist[i], _buffer_hist );
	}
}




template <typename HIST_FUN,  typename BIN_MAP>
int	piHistSlider<HIST_FUN,BIN_MAP>::stepSlide() {
	(this->*_column_slide_c)();
	if (!_avail_result_num)
		_row_slide();
	return _avail_result_num;
}

//----------------------------------------------------------------------------------------------
template <typename HIST_FUN, typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_generateResultAllocation ( ResultOutputPVec& rp ) {
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		ResultOutput* ro = new ResultOutput( uISize(_result_size[i]), _hist_fun->zero(), typename ResultOutput::WithFixedConstructor() );
		rp.push_back( ro );
	}
}

template <typename HIST_FUN, typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_generateResultAllocation ( ResultOutputPVec& rp, int ) {
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		ResultOutput* ro = new ResultOutput( uISize(_result_size[i]) );
		rp.push_back( ro );
	}
}


template <typename HIST_FUN, typename BIN_MAP>
inline void piHistSlider<HIST_FUN,BIN_MAP>::generateFastResults( ResultOutputPVec& rp ) {
	_generateResultAllocation( rp );
	_generateFastResults(rp);
}

template <typename HIST_FUN, typename BIN_MAP>
inline void  piHistSlider<HIST_FUN,BIN_MAP>::generateFastResults( ResultOutputPVec& rp, int ) {
	_generateResultAllocation( rp, 0 );
	_generateFastResults(rp);
}

template <typename HIST_FUN, typename BIN_MAP>
void piHistSlider<HIST_FUN,BIN_MAP>::_generateFastResults(ResultOutputPVec& r) {
	iPoint p(0,0);
	_initial_slide();
	while ( _avail_result_num ) {
		p.y = _current.y;
		for ( size_t i=0; i<_windows.win.size(); ++i ) {
			if ( p.y> _slide_range[i].height )
				continue;
			//Copy First Column Result
			result_t* cResults = r[i][_currentRP.y];
			*cResults = _fun_results_r[i];

			//Generate Row Results
			const swPolygonGroup& w = _windows.win[i];
			result_t&		fun_r = _fun_results[i];
			result_t*		dim_r = _dim_results[i];
			InternalHist&	th = _tmp_hist[i];
			int X = _slide_range[i].width+1;
			if (w.is_naiveInc) {
				ColumnHist* oBoundary = _h_boundary[_rb_h_Idx]+i,
					* nBoundary = _h_boundary[1-_rb_h_Idx]+i;
				for ( p.x = _windows.GetXStride(); p.x<X; p.x+=_windows.GetXStride() ) {
					//_rb_h_Idx = 1-_rb_h_Idx;
					_column_slide_naiveInc( p, w, *oBoundary, *nBoundary);
					_inc_fun_refresh( fun_r, dim_r, th, _buffer_hist );
					*(++cResults) = fun_r;
					std::swap(oBoundary,nBoundary);
				}
			} else {
				for ( p.x = _windows.GetXStride(); p.x<X; p.x+=_windows.GetXStride() ) {
					_column_slide_edgeInc( p, w );	//Buffer hist is wrong for multi-window*****************
					//----DEBUG
#if 0
					if (p.x == _windows.GetXStride()) {
						cerr << "No.\t" << i << endl;
						cerr << "fun_r:\t" << fun_r << endl;
						cerr << "dim_r:" << endl;
						for ( bin_t b=0; b<_binmap.binnum(); ++b ) {
								if ( dim_r[b] != 0.0 )
									cerr << b << "\t" << dim_r[b] << endl;
						}
						cerr << "th:" << endl;
						int k=0;
						for ( InternalHist::iterator iter=th.begin();
							iter != th.end(); ++iter ) {
								if ( *iter != 0.0 )
									cerr << k << "\t" << *iter << endl;
								++k;
						}
						cerr << "buffer:" << endl;
						for ( ColumnHist::iterator iter=_buffer_hist.begin();
							iter != _buffer_hist.end(); ++iter ) {
								cerr << iter->idx << "\t" << iter->val << endl;
						}
					}
#endif
					//---------
					_inc_fun_refresh( fun_r, dim_r, th, _buffer_hist );
					*(++cResults) = fun_r;
				}
			}
		}
		_row_slide();
	}
}

template <typename HIST_FUN, typename BIN_MAP>
inline void  piHistSlider<HIST_FUN,BIN_MAP>::generateNaiveResults( ResultOutputPVec& rp ) {
	_generateResultAllocation( rp );
	_generateNaiveResults(rp);
}

template <typename HIST_FUN, typename BIN_MAP>
inline void  piHistSlider<HIST_FUN,BIN_MAP>::generateNaiveResults( ResultOutputPVec& rp, int ) {
	_generateResultAllocation( rp, 0 );
	_generateNaiveResults(rp);
}

template <typename HIST_FUN, typename BIN_MAP>
void  piHistSlider<HIST_FUN,BIN_MAP>::_generateNaiveResults( ResultOutputPVec& r ) {
	int ry = 0, rx;
	for ( int y=0; y<=_max_slide_range.height; y+=_windows.GetYStride() ) {
		rx = 0;
		for ( int x=0; x<=_max_slide_range.width; x+=_windows.GetXStride() ) {
			this->funAt( x, y );
			for ( size_t i=0; i<_windows.win.size(); ++i )
				if ( this->ResultMask[i] )
					r[i][ry][rx] = this->ResultsFixed[i];
			++rx;
		}
		++ry;
	}
}

/* -------------------------------------------------------------------------------
FOR Integral Image
----------------------------------------------------------------------------------*/

template <typename HIST_FUN, typename BIN_MAP>
void  piHistSlider<HIST_FUN,BIN_MAP>::InitializeIntegralImage() {
	//Release Other Buffers
	_arr_buffer_holder.reset();
	//Find maximum uvec Y
	int maxUnitY = 0;
	_intZeroSlope	= -1;
	_intInfSlope	= -1;
	for ( size_t i=0; i<_windows.unit.size(); ++i ) {
		iPoint uvec = _windows.unit[i].vec;
		PI_CONDITIONAL_REFRESH(maxUnitY, uvec.y, <);
		if ( uvec.y == 0 )
			_intZeroSlope = i;
		if ( uvec.x == 0 )
			_intInfSlope  = i;
	}

	const piSlidingWindows& sw = _windows.OldSlidingWins();
	int X = _binmap.size().width+1, 
		Y = (sw.Shapes.RightBottom-sw.Shapes.LeftTop).y+maxUnitY+1+1;	//I don't know whether the second +1 is necessary or not
	int D = _windows.unit.size();
	int B = _binmap.binnum();
	if ( _intZeroSlope<0 ) {
		_intZeroSlope = D;
		++D;
	}
	if ( _intInfSlope<0 ) {
		_intInfSlope = D;
		++D;
	}

	IntegralType2* integralimg = 
		new IntegralType2( uISize(X,Y), uISize(B,D), typename IntegralType2::WithFixedConstructor() );
	_imap = IMapPtr(integralimg);
	//initialize
	for ( int x=0; x<X; ++x ) {
		for ( int d=0; d<D; ++d ) {
			for ( int b=0; b<B; ++b) {
				(*_imap)[0][x][d][b] = 0;
			}
		}
	}
	_intAtY = 0;
	_intSize.width  = X;
	_intSize.height = Y;
	_intD	= D;

	//initial integral image
	_current = iPoint(-_windows.GetXStride(),0);
	_currentRP = iPoint(-1,0);
	const piPolyGroupArr& shapes = _windows.OldSlidingWins().Shapes;
	_makeIntegral( shapes.RightBottom.y );
}



template <typename HIST_FUN, typename BIN_MAP>
void  piHistSlider<HIST_FUN,BIN_MAP>::_makeIntegral( int _Y ) {
	using namespace std;
	int Y = min(_Y, (int)_binmap.size().height);
	if ( _intAtY>=Y )
		return;
	IntegralType2& imap = *_imap;
	for ( int y=_intAtY+1; y<=Y; ++y ) {
		int cy = y % _intSize.height;
		IntegralType1* mrow = imap[cy];
		for ( int x=0; x<_intSize.width; ++x ) {
			IntegralType1& mdim = mrow[x];
			iPoint pB = iPoint(x-1,y-1); 
			for (size_t v = 0; v < _binmap_num; ++v) {
				int nowB  = (x)?((int)_binmaps[v](pB)):(-1);
				//Special integral image
				{
					IntegralType1& odim = imap[(y-1)% _intSize.height][x];
					for ( size_t b=0; b<_binmap.binnum(); ++b )
						mdim[_intInfSlope][b] = odim[_intInfSlope][b];
					if (nowB>=0) {
						weight_t wi = _binmaps[v][pB];
						mdim[_intInfSlope][nowB] += wi;
						IntegralType1& odim2 = imap[cy][x-1];
						for ( size_t b=0; b<_binmap.binnum(); ++b )
							mdim[_intZeroSlope][b] = odim2[_intZeroSlope][b]+mdim[_intInfSlope][b];
					} else {
						for ( size_t b=0; b<_binmap.binnum(); ++b )
							mdim[_intZeroSlope][b] = 0;
					}
				}
			}
			//Ordinary integral image
			for ( int d=0; d<_intD; ++d ) {
				if ( d==_intInfSlope || d==_intZeroSlope )
					continue;
				const swUnitStruct& unit=_windows.unit[d];
				const piUnitRasterizedLineSeg& urls = *(unit.rdata);
				iPoint p1(x,y);
				iPoint p0( p1-unit.vec );
				count_t* SI1 = mdim[d];
				if ( (unit.vec.x>0 && x==0) || (unit.vec.x<0 && x==(int)_binmap.size().width) ) {
					for ( size_t b=0; b<_binmap.binnum(); ++b)
						SI1[b] = 0;
				}
				else if ( p0.x>=0 && p0.y>=0 && IsPointLessThanRange(p0, iISize(_binmap.size())) ) {
					//p0 is inside the image
					//Pre
					int cy0 = p0.y % _intSize.height;
					count_t* NI1 = imap[cy0][x][_intZeroSlope];
					count_t* NI0 = imap[cy0][p0.x][_intZeroSlope];
					count_t* SI0 = imap[cy0][p0.x][d];
					for ( size_t b=0; b<_binmap.binnum(); ++b)
						SI1[b] = abs(NI1[b]-NI0[b])+SI0[b];
					//New
					if ( unit.vec.x>0 ) {
						for ( size_t i = 0; i<urls.PixelCount; ++i ) {
							iPoint p = p0+urls.Pixels[i];	// no need to -iPoint(1,1);
							for (size_t v = 0; v < _binmap_num; ++v) {
								bin_t nowB = _binmaps[v](p);
								weight_t wi= _binmaps[v][p];
								SI1[nowB] += urls.RightPortions[i]*wi;
							}
						}
					} else {
						for ( size_t i = 0; i<urls.PixelCount; ++i ) {
							iPoint p = p0+urls.Pixels[i];	//-iPoint(1,1);
							for (size_t v = 0; v < _binmap_num; ++v) {
								bin_t nowB = _binmaps[v](p);
								weight_t wi= _binmaps[v][p];
								SI1[nowB] += urls.LeftPortions[i]*wi;
							}
						}
					}
					for ( size_t i = 0; i<urls.InnerYs.size(); ++i ) {
						int xp = (unit.vec.x>0)?(p0.x+(int)i):(p0.x-(int)i-1);	//Please Recheck
						int YP = p0.y+urls.InnerYs[i];
						for ( int yp=p0.y; yp<=YP; ++yp ) {
							for (size_t v = 0; v < _binmap_num; ++v) {
								bin_t nowB = _binmaps[v](xp,yp);
								weight_t wi= _binmaps[v][yp][xp];
								SI1[nowB] += wi;
							}
						}
					}
				} else {

					//p0 is outside the image
					int startIdx;
					iPoint pp;
					if ( unit.vec.x>0 ) {
						for ( startIdx=0; startIdx<(int)urls.PixelCount; ++startIdx ) {
							pp = urls.Pixels[startIdx]+p0;
							if ( pp.x>=0 && pp.y>=0 )
								break;
						}
					} else {
						for ( startIdx=0; startIdx<(int)urls.PixelCount; ++startIdx ) {
							pp = urls.Pixels[startIdx]+p0;
							++pp.x;
							if ( pp.x<=(int)_binmap.size().width && pp.y>=0 )
								break;
						}
					}
					
					//Pre
					if ( pp.y>0 ) {
						int cy0 = pp.y % _intSize.height;;
						count_t* NI1 = imap[cy0][x][_intZeroSlope];
						count_t* NI0 = imap[cy0][pp.x][_intZeroSlope];
						for ( size_t b=0; b<_binmap.binnum(); ++b)
							SI1[b] = abs(NI1[b]-NI0[b]);
					} else {
						for ( size_t b=0; b<_binmap.binnum(); ++b)
							SI1[b] = 0;
					}
					//New
					if ( unit.vec.x>0 ) {
						for ( size_t i = startIdx; i<urls.PixelCount; ++i ) {
							iPoint p = p0+urls.Pixels[i];	// no need to -iPoint(1,1);
							for (size_t v = 0; v < _binmap_num; ++v) {
								bin_t nowB = _binmaps[v](p);
								weight_t wi= _binmaps[v][p];
								SI1[nowB] += urls.RightPortions[i]*wi;
							}
						}
					} else {
						for ( size_t i = startIdx; i<urls.PixelCount; ++i ) {
							iPoint p = p0+urls.Pixels[i];	//-iPoint(1,1);
							for (size_t v = 0; v < _binmap_num; ++v) {
								bin_t nowB = _binmaps[v](p);
								weight_t wi= _binmaps[v][p];
								SI1[nowB] += urls.LeftPortions[i]*wi;
							}
						}
					}
					for ( size_t i = abs(pp.x-p0.x); i<urls.InnerYs.size(); ++i ) {
						int xp = (unit.vec.x>0)?(p0.x+(int)i):(p0.x-(int)i-1);	//Please Recheck
						int YP = p0.y+urls.InnerYs[i];
						for ( int yp=pp.y; yp<=YP; ++yp ) {
							for (size_t v = 0; v < _binmap_num; ++v) {
								bin_t nowB = _binmaps[v](xp,yp);
								weight_t wi= _binmaps[v][yp][xp];
								SI1[nowB] += wi;
							}
						}
					}
					
				}
			}
		}
	}
	_intAtY = Y;
}

template <typename HIST_FUN, typename BIN_MAP>
void  piHistSlider<HIST_FUN,BIN_MAP>::_intHistInternal( iPoint p ) {
	using namespace std;
	IntegralType2& imap = *_imap;
	_avail_result_num = 0;
	for ( size_t i=0; i<_windows.win.size(); ++i ) {
		TEST_AVAILABILITY(_current,i,_avail_result_num);
		//test available
		const piPolygonGroup& shape = _windows.OldSlidingWins().Shapes[i];
		const piSlidingWindows::SegIdxVecPVec& sivp =  _windows.OldSlidingWins().SegIdx[i];
		InternalHist& h = _fixed_hist[i];
		h.set( count_t(0) );
		for ( size_t j=0; j<shape.size(); ++j ) {
			const piPolygon& poly = shape[j];
			const piSlidingWindows::SegIdxVec& siv = sivp[j];
			for (size_t k=0; k<poly.Vertice.size(); ++k ) {
				int ugid = _windows.OldSlidingWins().Seg2Unit[siv[k]];
				int k2 = (k+1) % poly.Vertice.size();
				iPoint p1 = poly.Vertice[k]+p, p2 = poly.Vertice[k2]+p;
				count_t* m1 = imap[p1.y%_intSize.height][p1.x][ugid];
				count_t* m2 = imap[p2.y%_intSize.height][p2.x][ugid];
				if (p1.x>p2.x) {
					for ( size_t b = 0; b<_binmap.binnum(); ++b )
						h[b] += abs(m1[b]-m2[b]);
				} 
				else if (p1.x<p2.x) {
					for ( size_t b = 0; b<_binmap.binnum(); ++b )
						h[b] -= abs(m1[b]-m2[b]);
				}
			}
		}
	}
}

template <typename HIST_FUN, typename BIN_MAP>
int piHistSlider<HIST_FUN,BIN_MAP>::IntegralSlide() {
	_current.x += _windows.GetXStride();
	++_currentRP.x;
	_intHistInternal(_current);
	if (!_avail_result_num) {
		_current.x = 0;
		_currentRP.x = 0;
		_current.y += _windows.GetYStride();
		++_currentRP.y;
		_makeIntegral( _intAtY+_windows.GetYStride() );
		_intHistInternal(_current);
	}
	_full_cal_fun( _fixed_hist, _fixed_results );
	return _avail_result_num;
}


template <typename HIST_FUN, typename BIN_MAP>
inline void  piHistSlider<HIST_FUN,BIN_MAP>::generateIntegralResults( ResultOutputPVec& rp ) {
	_generateResultAllocation( rp );
	_generateIntegralResults(rp);
}

template <typename HIST_FUN, typename BIN_MAP>
inline void  piHistSlider<HIST_FUN,BIN_MAP>::generateIntegralResults( ResultOutputPVec& rp, int ) {
	_generateResultAllocation( rp, 0 );
	_generateIntegralResults(rp);
}

template <typename HIST_FUN, typename BIN_MAP>
void  piHistSlider<HIST_FUN,BIN_MAP>::_generateIntegralResults( ResultOutputPVec& r ) {
	InitializeIntegralImage();
	while ( IntegralSlide() ) {
		for ( size_t i=0; i<_windows.win.size(); ++i )
			if ( this->ResultMask[i] )
				r[i][_currentRP.y][_currentRP.x] = this->ResultsFixed[i];
	}
}

#undef TEST_AVAILABILITY

