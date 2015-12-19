/*
 * polyint.cpp
 *
 *  Created on: 2010-10-21
 *      Author: zhangyuting
 */

#include "polyint.h"
#include <limits>
#include <climits>
#include <cmath>
#include <boost/rational.hpp>
#include <vector>
using namespace std;

#define Z_CONSTRUCTOR_ORIG	_orig

/*------------------ piPolygon ------------------*/

#define PI_POLYGON_INI_POSTFIX	 \
		Vertice(_vertice), NormVertice(_norm_vertice), \
		NormEdges(_norm_edges), RightPPR(_rightPPR), TopPPR(_topPPR), \
		LeftTop(_lefttop), RightBottom(_rightbottom), \
		Boundary(_boundary), Size(_size)

piPolygon::piPolygon(const VertexArr& vertice) :
	_vertice(vertice), _norm_vertice(vertice.size()), _norm_edges(vertice.size()),
	_rightPPR(vertice.size()), _topPPR(vertice.size()),
	_lefttop(INT_MAX, INT_MAX),
	_rightbottom(INT_MIN, INT_MIN),
	_boundary(), PI_POLYGON_INI_POSTFIX
{
	_initialize();
}

piPolygon::piPolygon(const MyType& _orig) :
	Z_FCC(_vertice), Z_FCC(_norm_vertice), Z_FCC(_norm_edges),
	Z_FCC(_rightPPR), Z_FCC(_topPPR), Z_FCC(_lefttop), Z_FCC(_rightbottom),
	Z_FCC(_boundary), Z_FCC(_size), PI_POLYGON_INI_POSTFIX
{
	_mask = MaskPtr(new MaskType(_orig.GetMask()) );
}

#undef PI_POLYGON_INI_POSTFIX

void piPolygon::_initialize() {
	//NormVertice
	for ( size_t i=0; i<_vertice.size(); ++i ) {
		iPoint  p= _vertice[i];
		if ( p.x< _lefttop.x )
			_lefttop.x = p.x;
		if ( p.x> _rightbottom.x )
			_rightbottom.x = p.x;
		if ( p.y< _lefttop.y )
			_lefttop.y = p.y;
		if ( p.y> _rightbottom.y )
			_rightbottom.y = p.y;
	}
	_size = iISize(0,0) + (_rightbottom - _lefttop);
	
	for ( size_t i=0; i<_vertice.size(); ++i )
		_norm_vertice[i] = _vertice[i] - _lefttop;
	VertexArr nv(_norm_vertice.size()+1,_norm_vertice);
	nv[_norm_vertice.size()] = _norm_vertice[0];

	//NormEdges & PPR
	for ( size_t i=0; i<_norm_edges.size(); ++i ) {
		lineseg_t& ne =  _norm_edges[i];
		iPoint vec = nv[i+1] - nv[i];
		bool RightOutside = (vec.y>0) || (vec.y==0 && vec.x > 0);
		bool TopOutside   = (vec.x>0) || (vec.x==0 && vec.y > 0);
		_topPPR[i] = (TopOutside)?(OutsidePolygon):(InsidePolygon);
		if ( RightOutside ) {
			_rightPPR[i] = OutsidePolygon;
			ne.sp  = nv[i];
			ne.vec = vec;
		} else {
			_rightPPR[i] = InsidePolygon;
			ne.sp  = nv[i+1];
			ne.vec = iPoint(0,0)-vec;
		}
	}

	//Boundary
	for ( size_t i=0; i<_vertice.size() ; ++i ) {
		lineseg_t lns;
		lns = _norm_edges[i];
		lns.sp = lns.sp + _lefttop;
		piRasterizedLineSeg* rp = new piRasterizedLineSeg(lns);
		_boundary.push_back( rp );
	}
	
	//Mask
	MaskType* mptr = new MaskType(_size);
	MaskType& mask = *mptr;
	_mask = MaskPtr(mptr);
	for ( size_t y=0; y<_size.height; ++y ) {
		for  ( size_t x=0; x<_size.width; ++x ) {
			mask[y][x] = (IsPointInPolygon( iPoint(x,y) ))?(InsidePolygon):(OutsidePolygon);
		}
	}

	for ( size_t i=0; i<_vertice.size() ; ++i ) {
		const piRasterizedLineSeg& rls = _boundary[i];
		if (!rls.Line.vec.x) continue;
		const piRasterizedLineSeg::PixelArr& p = rls.Pixels;
		for ( piRasterizedLineSeg::PixelArr::const_iterator iter = p.begin();
				iter != p.end(); ++iter ) {
					iPoint p = *iter - _lefttop;
					mask[p.y][p.x] = OnPolygon;
		}
	}
	
}

bool piPolygon::IsPointInPolygon( iPoint p ) const {	//Normlized point according to lefttop
	piFloat x = p.x+0.5f, y = p.y+0.5f;
	size_t n = 0;
	for ( size_t i=0; i < _norm_edges.size(); ++i ) {
		const lineseg_t& l = _norm_edges[i];
		iPoint p1 = l.sp, p2 = l.sp + l.vec;
		if ( (x>p1.x&&x<p2.x) || (x>p2.x&&x<p1.x) ) {
			piFloat u = _boundary[i].Slope * piFloat(x - p1.x) + (piFloat)p1.y;
			if (u<y)
				++n;
		}
	}
	return bool(n&1u);
}

void piPolygon::transit ( iPoint tvec ) {
	_lefttop = _lefttop + tvec;
	_rightbottom = _rightbottom + tvec;
	_vertice += tvec;
	BATCH_POST_OP4VECTOR(_boundary, .transit(tvec) );
}

void piPolygon::locate ( iPoint loc ) {
	transit( loc - _lefttop );
}

/*------------------ piUnitRasterizedLineSeg ------------------*/

#define PI_UNIT_LINESEG_INI_POSTFIX	 \
		Vec(_vec), Pixels(_pixels), LeftPortions(_left_portions), RightPortions(_right_portions), \
		TopPortions(_top_portions), BottomPortions(_bottom_portions), InnerYs(_inner_ys), TopYs(_top_ys)

piUnitRasterizedLineSeg::piUnitRasterizedLineSeg(iPoint vec) :
	Frac(boost::gcd(abs(vec.x), abs(vec.y))), _vec( vec / Frac ),
	Slope((_vec.x)?(piFloat(_vec.y)/piFloat(_vec.x)):(numeric_limits<piFloat>::infinity())),
	PixelCount((_vec.x)?(abs(_vec.x)+abs(_vec.y)-1):(1)), _pixels(PixelCount),
	_left_portions(PixelCount), _right_portions(PixelCount), _top_portions(PixelCount), _bottom_portions(PixelCount),
	_inner_ys(abs(_vec.x)), _top_ys(abs(_vec.x)),
	PI_UNIT_LINESEG_INI_POSTFIX {
	_initialize();
}

piUnitRasterizedLineSeg::piUnitRasterizedLineSeg(const MyType& _orig) :
	Z_FCC(Frac), Z_FCC( _vec ), Z_FCC(Slope), Z_FCC(PixelCount),
	Z_FCC(_pixels), Z_FCC(_left_portions), Z_FCC(_right_portions), Z_FCC(_top_portions), Z_FCC(_bottom_portions),
	Z_FCC(_inner_ys), Z_FCC(_top_ys),
	PI_UNIT_LINESEG_INI_POSTFIX {
}

void piUnitRasterizedLineSeg::_initialize() {
	if (!PixelCount) return;

	if (!_vec.x) {
		_pixels[0] = iPoint(0,0);
		_right_portions[0]=1;
		_left_portions[0] =0;
		_top_portions[0]  =0;
		_bottom_portions[0]=0;
		return;
	}

	enum {
		XChange = 1,
		YChange = 0
	};

	piFloat PSlope = abs(Slope);


	int C0, C1;
	piFloat Q0, Q1;

	iPoint p(0,0);
	C1 = XChange;
	Q1 = 0;
	_inner_ys[0] = -1;
	_top_ys[0]   = 0;


	for (size_t i=0; i<PixelCount; ++i) {
		_pixels[i] = p;
		C0 = C1; Q0 = Q1;
		if ((p.x+1) * PSlope>p.y+1+numeric_limits<piFloat>::epsilon()) {
			++p.y;
			Q1 = p.y / PSlope;
			Q1 = Q1 - ceil(Q1) + 1;
			C1 = YChange;
			if ( C0 == XChange )
				_right_portions[i] = 1.0f-((1-Q0)*Q1)*0.5f; 	// XY: 1 - Upper triangle
			else
				_right_portions[i] = 1.0f-(Q0+Q1)*0.5f;	// YY: right trapezoidal
			++_top_ys[p.x];
		} else {
			++p.x;
			Q1 = p.x * PSlope;
			Q1 = Q1 - ceil(Q1) + 1;
			C1 = XChange;
			if ( C0 == XChange )
				_right_portions[i] = (Q0+Q1)*0.5f;	// XX: lower trapezoidal
			else
				_right_portions[i] = ((1-Q0)*Q1)*0.5f;	// YX: lower triangle
			if (p.x<(int)_top_ys.size()) {
				_top_ys[p.x] = p.y;
				_inner_ys[p.x] = p.y-1;
			}
		}
		_top_portions[i] = _right_portions[i];
	}


	if (_vec.x<0) {
		for (size_t j=0; j<PixelCount ; ++j) {
			_pixels[j].x = -(_pixels[j].x+1);
			_right_portions[j] = 1-_right_portions[j];
		}
	}
	if (_vec.y<0) {
		for (size_t j=0; j<PixelCount ; ++j) {
			_pixels[j].y = -(_pixels[j].y+1);
			_top_portions[j] = 1-_top_portions[j];
		}
	}
	for (size_t j=0; j<PixelCount ; ++j) {
		_left_portions[j]   = 1-_right_portions[j];
		_bottom_portions[j] = 1-_top_portions[j];
	}
}

#undef PI_UNIT_LINESEG_INI_POSTFIX

#define PI_LINESEG_INI_POSTFIX \
	Pixels(_pixels), LeftPortions(_left_portions), RightPortions(_right_portions), \
	TopPortions(_top_portions), BottomPortions(_bottom_portions), InnerYs(_inner_ys), TopYs(_top_ys), SegUnit(_unitseg)

piRasterizedLineSeg::piRasterizedLineSeg(LineSeg lnseg) :
	Line(lnseg), _unitseg(lnseg.vec), Slope(_unitseg.Slope), PixelCount(_unitseg.PixelCount*_unitseg.Frac),
	_pixels(PixelCount), _left_portions(PixelCount), _right_portions(PixelCount), _top_portions(PixelCount), _bottom_portions(PixelCount),
	_inner_ys(abs(lnseg.vec.x)), _top_ys(abs(lnseg.vec.x)),
	PI_LINESEG_INI_POSTFIX {

	arrcpy( _pixels, _unitseg.Pixels, _unitseg.PixelCount, _unitseg.Frac, _unitseg.Vec, Line.sp );
	arrcpy( _left_portions,		_unitseg.LeftPortions, 		_unitseg.PixelCount, _unitseg.Frac, 0.0f, 0.0f );
	arrcpy( _right_portions,	_unitseg.RightPortions,		_unitseg.PixelCount, _unitseg.Frac, 0.0f, 0.0f );
	arrcpy( _top_portions,		_unitseg.TopPortions, 		_unitseg.PixelCount, _unitseg.Frac, 0.0f, 0.0f );
	arrcpy( _bottom_portions,	_unitseg.BottomPortions,	_unitseg.PixelCount, _unitseg.Frac, 0.0f, 0.0f );
	//arrcpy( _inner_ys, _unitseg.InnerYs, _unitseg.InnerYs.size(), _unitseg.Frac, _unitseg.Vec.y, Line.sp.y );
	//arrcpy( _top_ys, _unitseg.TopYs, _unitseg.TopYs.size(), _unitseg.Frac, _unitseg.Vec.y, Line.sp.y );
	arrcpy( _inner_ys, _unitseg.InnerYs, _unitseg.InnerYs.size(), _unitseg.Frac, 0, Line.sp.y );
	arrcpy( _top_ys, _unitseg.TopYs, _unitseg.TopYs.size(), _unitseg.Frac, 0, Line.sp.y );

}

piRasterizedLineSeg::piRasterizedLineSeg(const MyType& _orig) :
		Z_FCC(Line), Z_FCC(_unitseg), Z_FCC(Slope), Z_FCC(PixelCount),
		Z_FCC(_pixels), Z_FCC(_left_portions), Z_FCC(_right_portions), Z_FCC(_top_portions), Z_FCC(_bottom_portions),
		Z_FCC(_inner_ys), Z_FCC(_top_ys), PI_LINESEG_INI_POSTFIX {
}


#undef PI_LINESEG_INI_POSTFIX

#define PI_PIXEL_INTERSECTION_INI_POSTFIX \
	Intersection(_inter->_in), Residual1(_inter->_r1), Residual2(_inter->_r2)


piPixelIntersection::piPixelIntersection( const PixelArr& p1, const PixelArr& p2 ) :
	_inter(_internal_t_factory(p1,p2)), PI_PIXEL_INTERSECTION_INI_POSTFIX {}

piPixelIntersection::piPixelIntersection( const MyType& _orig ) :
	_inter(new _internal_t(*(_orig._inter))),  PI_PIXEL_INTERSECTION_INI_POSTFIX {}

piPixelIntersection::piPixelIntersection( const array1d< MyType >& _origArr )  :
	_inter(_internal_t_factory(_origArr)),  PI_PIXEL_INTERSECTION_INI_POSTFIX {}

piPixelIntersection::~piPixelIntersection() {
	delete _inter;
}

void piRasterizedLineSeg::transit( iPoint tvec ) {
	Line.sp = Line.sp + tvec;
	_pixels += tvec;
	//_inner_ys += tvec.y;
	//_top_ys +=tvec.y;
}

void piRasterizedLineSeg::locate( iPoint loc ) {
	transit( loc - Line.sp );
}





bool iPoint_cmp( iPoint a, iPoint b ) {
	return (a.y<b.y)||(a.y==b.y&&a.x<b.x);
}

piPixelIntersection::_internal_t* piPixelIntersection::_internal_t_factory( const PixelArr& p1, const PixelArr& p2 ) {
	std::vector<iPoint> vec1(p1.begin(), p1.end()), vec2(p2.begin(), p2.end());
	sort(vec1.begin(),vec1.end(),iPoint_cmp);
	sort(vec2.begin(),vec2.end(),iPoint_cmp);
	PixelArr arr1(vec1.begin(),vec1.end(),PixelArr::FromIterator()), 
		arr2(vec2.begin(),vec2.end(),PixelArr::FromIterator());

	typedef vector<iPoint>	PixelVec;
	PixelVec in, r1, r2;
	PixelArr::iterator iter1 = arr1.begin(), iter2 = arr2.begin();
	for (;;) {
		if (iter1 == arr1.end()) {
			if (iter2 == arr2.end())
				break;
			r2.push_back(*iter2);
			++iter2;
		}
		else if ( iter2 == arr2.end() ) {
			r1.push_back(*iter1);
			++iter1;
		}
		else if (*iter1==*iter2) {
			in.push_back(*iter1);
			++iter1;++iter2;
		}
		else if (iPoint_cmp(*iter1,*iter2)) {
			r1.push_back(*iter1);
			++iter1;
		}
		else {
			r2.push_back(*iter2);
			++iter2;
		}
	}

	_internal_t* a = new _internal_t(in.size(),r1.size(),r2.size());
	arrcpy( a->_in, in, in.size() );
	arrcpy( a->_r1, r1, r1.size() );
	arrcpy( a->_r2, r2, r2.size() );

	return a;
}

piPixelIntersection::_internal_t* piPixelIntersection::_internal_t_factory( const array1d< MyType >& _origArr ) {
	size_t in_num = 0, r1_num = 0, r2_num = 0;
	for (array1d< MyType >::const_iterator iter = _origArr.begin();
		iter != _origArr.end(); ++iter ) {
			in_num += iter->Intersection.size();
			r1_num += iter->Residual1.size();
			r2_num += iter->Residual1.size();
	}
	_internal_t* a = new _internal_t(in_num,r1_num,r2_num);
	PixelArr::iterator in_iter = a->_in.begin(),
		r1_iter = a->_r1.begin(),
		r2_iter = a->_r2.begin();
	for (array1d< MyType >::const_iterator iter = _origArr.begin();
		iter != _origArr.end(); ++iter ) {
			copy( iter->Intersection.begin(), iter->Intersection.end(), in_iter);
			copy( iter->Residual1.begin(), iter->Residual1.end(), r1_iter);
			copy( iter->Residual2.begin(), iter->Residual2.end(), r2_iter);
			in_iter += iter->Intersection.size();
			r1_iter += iter->Residual1.size();
			r2_iter += iter->Residual2.size();
	}
	return a;
}

#undef PI_PIXEL_INTERSECTION_INI_POSTFIX

#undef Z_CONSTRUCTOR_ORIG
