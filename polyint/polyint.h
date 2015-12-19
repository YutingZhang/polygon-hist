/*
 * polyint.h
 *
 *  Created on: 2010-10-21
 *      Author: zhangyuting
 */

#ifndef POLYINT_H_
#define POLYINT_H_

#include "config.h"

#include "zutils.h"
#include <memory>

NEW_PAIR_TEMPLATE(piLineSegT, sp, vec);	//sp - start point, vec - vector

typedef piLineSegT<iPoint> lineseg_t;

class piRasterizedLineSeg;

//Only for simple polygons with integer coordinates,
//and CLOCKWISE vertice
class piPolygon {
public:
	typedef piPolygon		MyType;

	enum PPRelation {	//Point-Polygon relation
		OutsidePolygon	=  -1,
		InsidePolygon	=  1,
		OnPolygon		=  0,
	};

	typedef iPoint	vertex_t;

	typedef array1d<vertex_t>	VertexArr;
	typedef array1d<lineseg_t>	EdgeArr;
	typedef array1d<int> 		PPRArr;
	typedef ptrvec<piRasterizedLineSeg>	BoundaryArr;
	typedef array2d<signed char>	MaskType;
	typedef std::auto_ptr<MaskType>	MaskPtr;
private:
	VertexArr	_vertice;
	VertexArr	_norm_vertice;
	EdgeArr		_norm_edges;
	PPRArr		_rightPPR;
	PPRArr		_topPPR;
	iPoint		_lefttop;
	iPoint		_rightbottom;
	BoundaryArr	_boundary;
	uISize		_size;
	MaskPtr		_mask;
public:
	const VertexArr&	Vertice;		//original vertice
	const VertexArr&	NormVertice;	//applied transition (-Left,-Top)
	const EdgeArr&		NormEdges;		//Edges on NormVertice, and from top to bottom
	const PPRArr&		RightPPR;		//The right-side region properties: in, out, on; if slope==0, then top region
	const PPRArr&		TopPPR;		//The right-side region properties: in, out, on; if slope==0, then top region
	const iPoint&		LeftTop;		//The LeftTop coordinate
	const iPoint&		RightBottom;	//The LeftTop coordinate
	const BoundaryArr&	Boundary;		//The rasterized edge
	const uISize&		Size;			//The ploygon size
	piPolygon( const VertexArr& vertice );
	piPolygon( const MyType& _orig );
	~piPolygon() {}
	void transit( iPoint tvec );
	void locate( iPoint loc );
	const MaskType& GetMask() const { return *_mask; }
	bool   IsPointInPolygon( iPoint p ) const;	//Does NOT work when p is in the boundary region
private:
	const MyType& operator = ( const MyType& ) { return (*this); }
	void _initialize();
};

class piUnitRasterizedLineSeg {
public:
	typedef piUnitRasterizedLineSeg MyType;
	typedef array1d<iPoint> PixelArr;
	typedef piFloat					PortionType;
	typedef array1d<PortionType> 	PortionArr;
	typedef array1d<int>	HeightArr;
public:
	const int	Frac;
private:
	iPoint		_vec;
public:
	const piFloat 	Slope;
	const size_t	PixelCount;
private:
	PixelArr	_pixels;
	PortionArr	_left_portions;
	PortionArr	_right_portions;
	PortionArr	_top_portions;
	PortionArr	_bottom_portions;
	HeightArr	_inner_ys;
	HeightArr	_top_ys;
public:
	const iPoint&		Vec;
	const PixelArr&		Pixels;
	const PortionArr&	LeftPortions;	//Pixel Portion of the left-side
	const PortionArr&	RightPortions;	//Pixel Portion of the right-side
	const PortionArr&	TopPortions;	//Pixel Portion of the top-side
	const PortionArr&	BottomPortions;	//Pixel Portion of the bottom-side
	const HeightArr&	InnerYs;	//Interior Y coordinates, maybe -1, as to x axis
	const HeightArr&	TopYs;		//Boundary Y coordinates

	piUnitRasterizedLineSeg(iPoint vec);
	piUnitRasterizedLineSeg(const MyType& _orig);
	~piUnitRasterizedLineSeg(){}
private:
	const MyType& operator = ( const MyType& ) { return (*this); }
	void _initialize();
};

class piRasterizedLineSeg {
public:
	typedef piRasterizedLineSeg				MyType;
	typedef piUnitRasterizedLineSeg::PixelArr	PixelArr;
	typedef piUnitRasterizedLineSeg::PortionType	PortionType;
	typedef piUnitRasterizedLineSeg::PortionArr 	PortionArr;
	typedef piUnitRasterizedLineSeg::HeightArr 	HeightArr;
	typedef lineseg_t		LineSeg;
public:
	LineSeg		Line;
private:
	piUnitRasterizedLineSeg _unitseg;
public:
	const piFloat 	Slope;
	const size_t	PixelCount;
private:
	PixelArr	_pixels;
	PortionArr	_left_portions;
	PortionArr	_right_portions;
	PortionArr	_top_portions;
	PortionArr	_bottom_portions;
	HeightArr	_inner_ys;
	HeightArr	_top_ys;
public:
	const PixelArr&		Pixels;
	const PortionArr&	LeftPortions;	//Pixel Portion of the left-side
	const PortionArr&	RightPortions;	//Pixel Portion of the right-side
	const PortionArr&	TopPortions;	//Pixel Portion of the top-side
	const PortionArr&	BottomPortions;	//Pixel Portion of the bottom-side
	const HeightArr&	InnerYs;
	const HeightArr&	TopYs;
	const piUnitRasterizedLineSeg& SegUnit;

	piRasterizedLineSeg(LineSeg lnseg);
	piRasterizedLineSeg(const MyType& _orig);
	~piRasterizedLineSeg(){}
	void transit( iPoint tvec );
	void locate( iPoint loc );
private:
	const MyType& operator = ( const MyType& ) { return (*this); }
};

class piPixelIntersection {
public:
	typedef piPixelIntersection MyType;
	typedef array1d<iPoint> PixelArr;
private:
	typedef struct piPixelIntersection_internal {
		PixelArr _in;	//intersection
		PixelArr _r1;	//residual 1
		PixelArr _r2;	//residual 2
		piPixelIntersection_internal( size_t in_num, size_t r1_num, size_t r2_num ) :
			_in(in_num), _r1(r1_num), _r2(r2_num) { }
		piPixelIntersection_internal( const piPixelIntersection_internal& _orig ) :
			_in(_orig._in), _r1(_orig._r1), _r2(_orig._r2) { }
	} _internal_t;
	_internal_t* _inter;
public:
	const PixelArr& Intersection;
	const PixelArr& Residual1;
	const PixelArr& Residual2;
	piPixelIntersection( const PixelArr& p1, const PixelArr& p2 );
	piPixelIntersection( const MyType& _orig );
	piPixelIntersection( const array1d< MyType >& _origArr );	//Might cause consistency issue
	~piPixelIntersection();
private:
	_internal_t* _internal_t_factory( const PixelArr& p1, const PixelArr& p2 );
	_internal_t* _internal_t_factory( const array1d< MyType >& _origArr );
};

#endif /* POLYINT_H_ */
