/*
 * sliding_hist.cpp
 *
 *  Created on: 2010-10-26
 *      Author: zhangyuting
 */

#include "zutils.h"
#include "sliding_hist.h"
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <memory>
#include <cmath>
#include <set>
#include <boost/math/common_factor_rt.hpp>
using namespace std;

/* ----------- Polygon group ------------- */

void piPolygonGroup::transit(iPoint tvec) {
	_lefttop = _lefttop + tvec;
	_rightbottom = _rightbottom + tvec;
	for ( size_t i=0; i<size() ; ++i )
		(*this)[i].transit(tvec);
}

void piPolygonGroup::locate(iPoint loc) {
	transit(loc-_lefttop);
}

void piPolygonGroup::locateOrig() {
	transit(iPoint(0,0)-_lefttop);
}

/* ----------- Sliding window ------------ */

#define PI_SLIDINGWINDOWS_INI_POSTFIX \
		Shapes(_shapes), WinNum(_winNum), WinSize(_winSize), PolyNum(_polyNum), \
		InteriorPixels(_interiorPixels), BodyIntersect(_bodyIntersect), BodyIntersectC(_bodyIntersectC), \
		BoundaryPixels (_boundaryPixels), BoundaryPortions(_boundaryPortions),	\
		SegIdx(_segIdx), Seg2Unit(_seg2unit), SegFrac(_segFrac), Unit2Seg(_unit2seg), \
		UnitSegPixels(_unitSegPixels), UnitSegPixelPortions(_unitSegPixelPortions), \
		SegYLowest(_segYlowest), SegYSpan(_segYspan), UnitSegYSpan(_unitYspan), \
		SegXSpan(_segXspan), UnitSegXSpan(_unitXspan), Idx2Seg(_idx2seg), Stride(_stride)

piSlidingWindows::piSlidingWindows( const piPolyGroupArr& shapes, const iPoint stride ) :
	_shapes(shapes), _stride(stride), _winNum(shapes.size()), _winSize(_winNum), _polyNum(_winNum),
	_segIdx(_winNum), 
	PI_SLIDINGWINDOWS_INI_POSTFIX {
	_initialize( shapes, stride );
}

#undef PI_SLIDINGWINDOWS_INI_POSTFIX


struct less4iPoint : binary_function <iPoint,iPoint,bool> {
	bool operator() ( const iPoint& a, const iPoint& b ) const {
		return (a.y<b.y)||(a.y==b.y&&a.x<b.x);
	}
};

struct iPointWithInt {
	iPoint	p;
	int		n;
public:
	iPointWithInt() {}
	iPointWithInt( iPoint _p, int _n ) : p(_p), n(_n) {}
};

struct less4iPointWithInt : binary_function <iPointWithInt,iPointWithInt,bool> {
	bool operator() ( const iPointWithInt& a, const iPointWithInt& b ) const {
		less4iPoint lp;
		return lp(a.p,b.p) || (a.p==b.p&&a.n<b.n);
	}
};


inline lineseg_t Find_Seg_Rectangle( const vector<iPoint>& lsp ) {
	typedef vector<iPoint> PosVec;
	iPoint lt( INT_MAX, INT_MAX ), rb( INT_MIN, INT_MIN );
	for ( PosVec::const_iterator iter = lsp.begin();
		iter != lsp.end(); ++iter ) {
			PI_CONDITIONAL_REFRESH( lt.x, iter->x, > );
			PI_CONDITIONAL_REFRESH( lt.y, iter->y, > );
			PI_CONDITIONAL_REFRESH( rb.x, iter->x, < );
			PI_CONDITIONAL_REFRESH( rb.y, iter->y, < );
	}
	return lineseg_t(lt, rb-lt);
}

void piSlidingWindows::_initialize( const piPolyGroupArr& shapes, const iPoint stride  ) {

	int strideGCD = boost::math::gcd(_stride.x,_stride.y);

	//Interior & Interior intersection
	//Boundary portions
	{
		for ( size_t i=0; i<_winNum; ++i ) {
			const piPolygonGroup& pg = shapes[i];
			iISize shapeSize(iISize(0,0)+pg.RightBottom);//-pg.LeftTop);
			_winSize[i] = shapeSize;
			_polyNum[i] = pg.size();


			//Definition for interior
			std::list<iPoint> interiorList;
			ptrvec<piPixelIntersection> bodyIntersectArr, bodyIntersectArrC;
			//Definitions for Boundary
			PixelVec&	bPix = *( new PixelVec );
			PortionVec&	bPor = *( new PortionVec );
			_boundaryPixels.push_back( &bPix );
			_boundaryPortions.push_back( &bPor );

			for ( size_t j=0; j<pg.size(); ++j ) {
				const piPolygon& poly = pg[j];

				//Mask
				std::list<iPoint> interiorListSingle;
				const piPolygon::MaskType& mask = poly.GetMask();
				for ( size_t y=0; y<mask.size().height; ++y ) {
					for ( size_t x=0; x<mask.size().width; ++x ) {
						if (mask[y][x]==piPolygon::InsidePolygon)
							interiorListSingle.push_back(iPoint(x,y)+poly.LeftTop);
					}
				}
				interiorList.insert( interiorList.end(), interiorListSingle.begin(), interiorListSingle.end() );

				//Interior Intersection
				{
					PixelArr* interiorArrSinglePtr = new PixelArr(interiorListSingle.begin(),interiorListSingle.end(),PixelArr::FromIterator());
					{
						PixelArr interiorArr1(*interiorArrSinglePtr);
						foreach( iPoint& p, interiorArr1 )
							p.y -= _stride.y;
						bodyIntersectArr.push_back(new piPixelIntersection(interiorArr1,*interiorArrSinglePtr));
					}
					{
						PixelArr interiorArr1(*interiorArrSinglePtr);
						foreach( iPoint& p, interiorArr1 )
							p.x -= _stride.x;
						bodyIntersectArrC.push_back(new piPixelIntersection(interiorArr1,*interiorArrSinglePtr));
					}
				}

				//Boundary ---------------------------------------------------------------
				//Definitions for Boundary-------------
				typedef std::map<iPoint, PortionType, less4iPoint>	P2PorType;
				typedef std::map<iPoint, IdxType, less4iPoint>		P2IdxType;
				typedef std::vector<PortionType>					PortionVec;
				typedef std::vector<bool>							BoolVec;
				P2PorType	p2por, p2topPor;
				P2IdxType	idx4conflicted;
				ptrvec< PortionVec >	conflicted_portions;
				PortionVec	min_portions;
				BoolVec		min_is_top;
				//Boundary portions --------------------
				for ( size_t k=0; k<poly.Boundary.size(); ++k ) {
					const piRasterizedLineSeg& rls = poly.Boundary[k];
					if (!rls.SegUnit.Vec.x)
						continue;
					bool	is_topEdge = (poly.TopPPR[k] == piPolygon::OutsidePolygon);
					for ( size_t m=0; m<rls.PixelCount; ++m ) {
						iPoint p = rls.Pixels[m];
						P2PorType::iterator	iter = p2topPor.find(p);
						if ( iter == p2topPor.end() ) {
							if (is_topEdge)
								p2por[p] = rls.BottomPortions[m];
							else
								p2por[p] = rls.TopPortions[m];
							p2topPor[p] = rls.TopPortions[m];
						} else {
							//Handle the case that multiple lines pass one pixel
							P2IdxType::iterator	ci = idx4conflicted.find(p);
							if ( ci == idx4conflicted.end() ) {
								PortionVec& pv = *(new PortionVec);
								conflicted_portions.push_back(&pv);

								idx4conflicted[p] = conflicted_portions.size()-1;
								ci = idx4conflicted.find(p);
								pv.push_back( iter->second );
								min_portions.push_back( iter->second );
								min_is_top.push_back( is_topEdge );
							}
							IdxType n = ci->second;
							PortionType tp = rls.TopPortions[m];
							conflicted_portions[n].push_back(tp);
							if ( tp < min_portions[n] ) {
								min_portions[n] = tp;
								min_is_top[n] = ( is_topEdge );
							}
						}
					}
				}
				for ( P2IdxType::iterator ci = idx4conflicted.begin();
						ci!=idx4conflicted.end(); ++ci ) {
					iPoint  p = ci->first;
					IdxType n = ci->second;
					PortionVec& pv = conflicted_portions[n];
					PortionType por = PortionType(0);
					pv.push_back(PortionType(0));
					pv.push_back(PortionType(1));
					sort(pv.begin(),pv.end());
					for ( IdxType m = 1; m<pv.size(); m+=2 ) {
						por += pv[m]-pv[m-1];
					}
					if ( !min_is_top[n] )
						por = 1-por;
					p2por[p] = por;
				}
				//Write into Boundary vector --------------------
				for ( P2PorType::iterator iter = p2por.begin();
						iter!=p2por.end(); ++iter ) {
					bPix.push_back(iter->first);
					bPor.push_back(iter->second);
				}

			}
			//Body intersect
			{
				PixelArr* interiorArrPtr = new PixelArr(interiorList.begin(),interiorList.end(),PixelArr::FromIterator());
				_interiorPixels.push_back( interiorArrPtr );
				{
					array1d<piPixelIntersection> ia( bodyIntersectArr.size(), 0, array1d<piPixelIntersection>::WithoutConstruction() );
					for ( size_t k=0; k<ia.size(); ++k )
						new (&(ia[k])) piPixelIntersection(bodyIntersectArr[k]);
					_bodyIntersect.push_back(new piPixelIntersection(ia));
				}
				{
					array1d<piPixelIntersection> ia( bodyIntersectArrC.size(), 0, array1d<piPixelIntersection>::WithoutConstruction() );
					for ( size_t k=0; k<ia.size(); ++k )
						new (&(ia[k])) piPixelIntersection(bodyIntersectArrC[k]);
					_bodyIntersectC.push_back(new piPixelIntersection(ia));
				}
			}
		}
	}

	//LineSeg
	{
		//Grouping LineSeg & UnitLineSeg
		map<iPointWithInt, IdxType, less4iPointWithInt> UnitSeg2Idx;
		map<iPointWithInt, IdxType, less4iPointWithInt> Seg2Idx;
		typedef vector<iPoint> PosVec;
		ptrvec<PosVec> LineSegPos;
		ptrvec<PosVec> LineSegStartPoint;
		ptrvec<PosVec> LineSegEndPoint;
		ptrvec<PosVec> UnitLineSegPos;
		ptrvec<SegIdxVec>	UnitIdx2SegIdx;
		IdxType latest = 0, ulatest = 0;
		for ( size_t i=0; i<_winNum; ++i ) {
			const piPolygonGroup& pg = shapes[i];
			SegIdxVecPVec& segi_vpv= _segIdx[i];
			for ( size_t j=0; j<pg.size(); ++j ) {
				const piPolygon& poly = pg[j];
				SegIdxVec& segi_v = *(new SegIdxVec);
				segi_vpv.push_back(&segi_v);
				for ( size_t k=0; k<poly.NormEdges.size(); ++k ) {
					lineseg_t ne = poly.NormEdges[k];
					ne.sp = ne.sp + poly.LeftTop;
					IdxType gid, ugid;
					bool is_new_unit = false;
					const piUnitRasterizedLineSeg& urls = poly.Boundary[k].SegUnit;
					const piRasterizedLineSeg& rls = poly.Boundary[k];
					iPointWithInt upwi( urls.Vec, (ne.sp.y % _stride.y)*_stride.x + (ne.sp.x % _stride.x) );
					if ( UnitSeg2Idx.find(upwi)==UnitSeg2Idx.end() ) {
						UnitSeg2Idx[upwi] = ulatest;
						ugid = ulatest;
						UnitLineSegPos.push_back(new PosVec);
						UnitIdx2SegIdx.push_back(new SegIdxVec);
						++ulatest;
						is_new_unit = true;
					} else {
						ugid = UnitSeg2Idx[upwi];
					}
					UnitLineSegPos[ugid].push_back(ne.sp+ne.vec/_stride);	//Problematic*****&&UnitSeg Sharing has problem, too
					UnitLineSegPos[ugid].push_back(ne.sp);

					iPointWithInt pwi( ne.vec, ne.sp.y * _stride.x + (ne.sp.x % _stride.x) );
					if ( is_new_unit || Seg2Idx.find(pwi)==Seg2Idx.end() ) {
						Seg2Idx[pwi]=latest;
						gid = latest;
						LineSegPos.push_back(new PosVec);
						LineSegStartPoint.push_back(new PosVec);
						LineSegEndPoint.push_back(new PosVec);
						++latest;
						_idx2seg.push_back(&(poly.Boundary[k]));
						_seg2unit.push_back(ugid);
						if (is_new_unit) 
							_unit2seg.push_back(gid);
						//Frac
						_segFrac.push_back( urls.Frac );
						_segYspan.push_back( (urls.Vec.x)?(urls.Vec.y*strideGCD/_stride.y):(1) );
					} else {
						gid = Seg2Idx[pwi];
					}
					UnitIdx2SegIdx[ugid].push_back(gid);
					LineSegPos[gid].push_back( ne.sp+rls.Line.vec );
					LineSegPos[gid].push_back( ne.sp );
					LineSegStartPoint[gid].push_back( ne.sp );
					LineSegEndPoint[gid].push_back( poly.RightBottom - ne.sp );
					segi_v.push_back(gid);
				}
			}

		}

		// y-span is grid-based
		// x-span is pixel-based

		//Seg Span
		for ( size_t i=0; i<LineSegPos.size(); ++i ) {
			lineseg_t rec;
			rec = Find_Seg_Rectangle(LineSegPos[i]);
			_segYlowest.push_back( rec.sp.y + 
				((_idx2seg[i]->Line.vec.x)?(_idx2seg[i]->SegUnit.Vec.y*strideGCD):(_stride.y)) );

			iSpan xspan;
			rec = Find_Seg_Rectangle(LineSegStartPoint[i]);
			xspan.e1 = rec.sp.x;
			rec = Find_Seg_Rectangle(LineSegEndPoint[i]);
			xspan.e2 = rec.sp.x;
			_segXspan.push_back( xspan );
		}

		//UnitSeg Span
		for ( size_t i=0; i<UnitLineSegPos.size(); ++i ) {
			lineseg_t rec = Find_Seg_Rectangle(UnitLineSegPos[i]);
			_unitYspan.push_back( rec.vec.y+1 );
			iSpan xspan(INT_MAX,INT_MAX);
			int vecX = 0, uvecX = 0;
			for (size_t j=0; j<UnitIdx2SegIdx[i].size(); ++j) {
				IdxType sidx= UnitIdx2SegIdx[i][j];
				iSpan sx = _segXspan[sidx];
				PI_CONDITIONAL_REFRESH( xspan.e1, sx.e1, > );
				PI_CONDITIONAL_REFRESH( xspan.e2, sx.e2, > );
				const piRasterizedLineSeg& rls = *_idx2seg[sidx];
				if ( abs(vecX)<abs(rls.Line.vec.x) ) {
					vecX  = rls.Line.vec.x;
					uvecX = rls.SegUnit.Vec.x;	//redundant
				}
			}
			if (vecX>=0)
				xspan.e2 -= vecX - uvecX * _stride.x;
			else
				xspan.e1 += vecX - uvecX * _stride.x;

			_unitXspan.push_back( xspan );
		}

		//Colunm-wise increase
		for ( size_t i=0; i<_unit2seg.size(); ++i ) {
			IdxType seg_idx=_unit2seg[i];
			const piUnitRasterizedLineSeg& urls = _idx2seg[seg_idx]->SegUnit;

			//find incremental region
			//for sparse histogram, we need to add the interior pixels, and use incremental refreshing
			const piRasterizedLineSeg::PixelArr& r = urls.Pixels;
			PixelArr l(r);
			for ( size_t j=0; j<r.size(); ++j )
				l[j].x-=stride.x;

			set<iPoint,less4iPoint> mid_region;
			{
				for ( int k=1; k<stride.x; ++k ) {
					for ( size_t j=0; j<r.size(); ++j )
						mid_region.insert(r[j]-iPoint(k,0));
				}
			}
			PixelArr mid(mid_region.size());
			copy(mid_region.begin(),mid_region.end(),mid.begin());


			if (!urls.Vec.x) {
				if ( stride.x==1 ) {
					_unitSegPixels.push_back(new PixelArr(l));
					_unitSegPixelPortions.push_back(new PortionArr(urls.RightPortions));
				} else {
					PixelArr* pxArr= new PixelArr(l.size()+mid.size());
					copy(l.begin(),l.end(),pxArr->begin());
					copy(mid.begin(),mid.end(),pxArr->begin()+l.size());
					PortionArr* porArr = new PortionArr(l.size()+mid.size());
					copy(urls.RightPortions.begin(),urls.RightPortions.end(),porArr->begin());
					fill(porArr->begin()+l.size(),porArr->end(),piFloat(1.0));
					_unitSegPixels.push_back(pxArr);
					_unitSegPixelPortions.push_back(porArr);
				}
				continue;
			}

			map<iPoint, PortionType, less4iPoint> l_rp, r_lp;

			for ( size_t j=0; j<r.size(); ++j ) {
				l_rp[l[j]] = urls.RightPortions[j];
				r_lp[r[j]] = urls.LeftPortions[j];
			}

			piPixelIntersection it( l, r );

			piPixelIntersection mid_l_it( mid, l );
			piPixelIntersection mid_lr_it( mid_l_it.Residual1, r );

			size_t tn = it.Residual1.size() + it.Intersection.size() + it.Residual2.size() +
					mid_lr_it.Residual1.size();
			PixelArr& incarr = *( new PixelArr( tn ) );
			PortionArr& incpor = *( new PortionArr( tn ) );
			_unitSegPixels.push_back(&incarr);
			_unitSegPixelPortions.push_back(&incpor);

			size_t k=0;
			for ( size_t j=0; j<it.Residual1.size(); ++j ) {
				iPoint p = it.Residual1[j];
				incarr[k]=p;
				incpor[k]=l_rp[p];
				++k;
			}
			for ( size_t j=0; j<it.Intersection.size(); ++j ) {
				iPoint p = it.Intersection[j];
				incarr[k]=p;
				incpor[k]=(r_lp[p]+l_rp[p])-piFloat(1.0);
				++k;
			}
			for ( size_t j=0; j<it.Residual2.size(); ++j ) {
				iPoint p = it.Residual2[j];
				incarr[k]=p;
				incpor[k]=r_lp[p];
				++k;
			}
			for ( size_t j=0; j<mid_lr_it.Residual1.size(); ++j ) {
				iPoint p = mid_lr_it.Residual1[j];
				incarr[k]=p;
				incpor[k]=piFloat(1.0);
				++k;
			}
		}
	}
}

//------------------------------------------------------------------------------------------------------------

void swBetterSlidingWindows::_initialize() {
	piSlidingWindows& psw =*(_psw);

	typedef piSlidingWindows::PixelArr		PixelArr;
	typedef piSlidingWindows::PortionArr	PortionArr;
	typedef piSlidingWindows::PixelVec		PixelVec;
	typedef piSlidingWindows::PortionVec	PortionVec;

	this->unitDupFrac4sparse_buffer = boost::math::gcd( GetXStride(), GetYStride() );

	//Unit
	for (size_t i=0; i<unit.size(); ++i) {
		swUnitStruct&	u = unit[i];
		u.Yspan	= psw.UnitSegYSpan[i];
		u.Xispan = psw.UnitSegXSpan[i];

		const PixelArr& posA = psw.UnitSegPixels[i];
		const PortionArr& porA = psw.UnitSegPixelPortions[i];
		array1d<swUnitPixel> sup(posA.size());
		for ( size_t j=0; j<sup.size(); ++j ) {
			sup[j].pos = posA[j];
			sup[j].por = porA[j];
		}
		u.pixels.Set2(sup);
		const piRasterizedLineSeg* rls = psw.Idx2Seg[ psw.Unit2Seg[i] ];
		u.vec = rls->SegUnit.Vec;
		u.rdata = &(rls->SegUnit);
	}

	//Seg
	for (size_t i=0; i<seg.size(); ++i) {
		swSegStruct&	s = seg[i];
		s.frac   = psw.SegFrac[i];
		s.fracY4sparse_buffer = s.frac/this->unitDupFrac4sparse_buffer;
		s.lowY   = psw.SegYLowest[i];
		s.Yspan  = psw.SegYSpan[i];
		s.Xispan = psw.SegXSpan[i];
		s.unit   = &(unit[psw.Seg2Unit[i]]);
		s.vec	 = psw.Idx2Seg[i]->Line.vec;
		s.inc_vec = s.unit->vec * unitDupFrac4sparse_buffer;

		int cf = boost::math::gcd( GetXStride(), GetYStride() );
		if ( _conf.stride.x>1 || _conf.stride.y>1 ) {	//Sparse grid
			if ( !s.vec.x ) {
				s.refreshType = swSegStruct::SS_Vertical;
			} else if ( (s.vec.y % GetYStride()) || (s.vec.x % GetXStride())	//Can't be accelerated
				|| s.fracY4sparse_buffer <= 1 ) {
				s.refreshType = swSegStruct::SS_UnitSeg;
			} else {
				if ( s.unit->pixels.size() * GetYStride() < _conf.buffer_pixel_threshold ) {
					s.refreshType = swSegStruct::SS_NoBuffer;
				} else {
					s.refreshType = swSegStruct::SS_Buffer;
					//s.refreshType = swSegStruct::SS_NoBuffer;
				}
				//s.refreshType = swSegStruct::SS_UnitSeg;
			}
			//s.refreshType = swSegStruct::SS_UnitSeg;
			
		} else {	//Dense grid
			if ( !s.vec.x ) {
				s.refreshType = swSegStruct::SS_Vertical;
			} else if ( s.frac == 1 ) {
				s.refreshType = swSegStruct::SS_UnitSeg;
			} else if ( s.unit->pixels.size() < _conf.buffer_pixel_threshold ) {
				s.refreshType = swSegStruct::SS_NoBuffer;
			} else {
				s.refreshType = swSegStruct::SS_Buffer;
			}
		}
		s.unit->seg.push_back(&s);
	}

	//Win
	for (size_t i=0; i<win.size(); ++i) {
		swPolygonGroup& w = win[i];

		//interior
		w.size = psw.WinSize[i];
		w.interior.Set2( psw.InteriorPixels[i] );
		w.interior_v_r1.Set2( psw.BodyIntersect[i].Residual1 );
		w.interior_v_r2.Set2( psw.BodyIntersect[i].Residual2 );
		w.interior_h_r1.Set2( psw.BodyIntersectC[i].Residual1 );
		w.interior_h_r2.Set2( psw.BodyIntersectC[i].Residual2 );

		{
			//Run length encode the interior
			if ( w.interior.size() ) {
				swPolygonGroup::ContinuousLine cl;
				PixelArr::const_iterator iter = w.interior.begin();
				cl.num = iter->y;
				cl.span.e1 = iter->x;
				while ((++iter)!=w.interior.end()) {
					if (iter->x-(iter-1)->x!=1 || iter->y!=cl.num) {
						cl.span.e2 = (iter-1)->x;
						w.interiorRL.push_back(cl);
						cl.num = iter->y;
						cl.span.e1 = iter->x;
					}
				}
				cl.span.e2 = (iter-1)->x;
				w.interiorRL.push_back(cl);
			}
		}

		//boundary
		const PixelVec&		bdp = psw.BoundaryPixels[i];
		const PortionVec&	ptp = psw.BoundaryPortions[i];
		swPolygonGroup::UPixelArr upa( bdp.size() );
		for ( size_t j=0; j<upa.size(); ++j ) {
			upa[j].pos = bdp[j];
			upa[j].por = ptp[j];
		}
		w.boundary.Set2(upa);

		//edges
		size_t nEdgeIncPixel = 0;
		const piPolygonGroup& pg = psw.Shapes[i];
		size_t edgeNum = 0;
		for (size_t j=0; j<pg.size(); ++j)
			edgeNum += pg[j].NormEdges.size();

		swPolygonGroup::EdgeArr ea(edgeNum);
		size_t m = 0;
		for (size_t j=0; j<pg.size(); ++j) {
			const piPolygon& poly = pg[j];
			for ( size_t k=0; k<poly.NormEdges.size(); ++k ) {
				ea[m].line = poly.NormEdges[k];
				ea[m].line.sp = ea[m].line.sp + poly.LeftTop;
				ea[m].is_right = (poly.RightPPR[k]==piPolygon::OutsidePolygon);
				size_t sidx = psw.SegIdx[i][j][k];
				swSegStruct* s = &(seg[ sidx ]);
				ea[m].seg  = s;
				++m;
				nEdgeIncPixel += s->unit->pixels.size()*s->frac;
			}
		}
		w.edges.Set2(ea);

		//inc method
		size_t nNaiveIncPixels = w.interior_v_r1.size() + w.interior_v_r2.size() + w.boundary.size();
		w._is_naiveInc = ( nEdgeIncPixel >= nNaiveIncPixels );
	}
	ChangeIncMethod(_conf.inc_method);

}

#define SW_BETTER_SLIDING_WINDOWS_CONSTRUCT_INI \
	_psw(new piSlidingWindows(_shapes, _conf.stride)),		\
	unit( _psw->Unit2Seg.size() ), seg( _psw->Seg2Unit.size() ),	\
	win( _psw->WinNum )

swBetterSlidingWindows::swBetterSlidingWindows( const piPolyGroupArr& shapes, const Conf& conf ) : 
	_shapes(shapes), _conf(conf), SW_BETTER_SLIDING_WINDOWS_CONSTRUCT_INI {
	_initialize();
}

swBetterSlidingWindows::swBetterSlidingWindows( const swBetterSlidingWindows& _orig ) :
	_shapes(_orig._shapes), _conf(_orig._conf), SW_BETTER_SLIDING_WINDOWS_CONSTRUCT_INI {
	_initialize();
}

swBetterSlidingWindows::~swBetterSlidingWindows() {
	//Release original
	delete _psw;
}

void swBetterSlidingWindows::ChangeIncMethod( int incMethod ) {
	_conf.inc_method = incMethod;
	for (size_t i=0; i<win.size(); ++i) {
		swPolygonGroup& w = win[i];
		if (incMethod == Conf::Inc_Edge ) {
			w.is_naiveInc = false;
		} else if (incMethod == Conf::Inc_Naive) {
			w.is_naiveInc = true;
		} else {	//Auto
			w.is_naiveInc = w._is_naiveInc;
		}
	}
}

