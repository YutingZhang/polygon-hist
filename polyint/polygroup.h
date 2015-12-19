/*
 * polygroup.h
 *
 *  Created on: 2010-11-8
 *      Author: zhangyuting
 */

#ifndef POLYGROUP_H_
#define	POLYGROUP_H_

#include "polyint.h"
#include <climits>


/* ------------------ Polygon structure -------------- */

//class piPolygon;

template<typename T>
class piGeometryGroup : public ptrvec<T> {
public:
	typedef piGeometryGroup MyType;
	typedef ptrvec<T> PolyVec;
protected:
	iPoint _lefttop;
	iPoint _rightbottom;
public:
	const iPoint& LeftTop;
	const iPoint& RightBottom;
	piGeometryGroup() : _lefttop( INT_MAX, INT_MAX ),
		_rightbottom( INT_MIN , INT_MIN ),
		LeftTop(_lefttop), RightBottom(_rightbottom) {}
	piGeometryGroup( const MyType& r ) : PolyVec(r), _lefttop(r._lefttop), _rightbottom(r._rightbottom),
		LeftTop(_lefttop), RightBottom(_rightbottom) {}
	virtual ~piGeometryGroup() {}
	virtual void push_back(T* ptr);
};

template<typename T>
void piGeometryGroup<T>::push_back(T* ptr) {

	PolyVec::push_back( ptr );

	PI_CONDITIONAL_REFRESH( _lefttop.x, ptr->LeftTop.x, > );
	PI_CONDITIONAL_REFRESH( _lefttop.y, ptr->LeftTop.y, > );
	PI_CONDITIONAL_REFRESH( _rightbottom.x, ptr->RightBottom.x, < );
	PI_CONDITIONAL_REFRESH( _rightbottom.y, ptr->RightBottom.y, < );

}

//A shape combined by multiple polygons
class piPolygonGroup : public piGeometryGroup<piPolygon> {
private:
	typedef piGeometryGroup<piPolygon> MyBase;
public:
	piPolygonGroup() {}
	virtual ~piPolygonGroup() {}
	//do NOT call transit or locate before you push back all the polygons
	//  unless you really know will happen
	void transit(iPoint tvec);
	void locate(iPoint loc);
	void locateOrig();
};
typedef piGeometryGroup<piPolygonGroup>		piPolyGroupArr;	//Separated shapes


#endif
