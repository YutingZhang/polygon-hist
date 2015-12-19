/*
 * zutil.h
 *
 *  Created on: 2010-10-21
 *      Author: zhangyuting
 */

#ifndef ZUTIL_H_
#define ZUTIL_H_

//--------------_DEBUG
#include <typeinfo>
#include <iostream>
//--------------

#include "config.h"

#include <cctype>
#include <algorithm>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/noncopyable.hpp>

namespace zUtilities_PolyHistInternal {

/* ------- Basic Operations ---------*/
#define foreach BOOST_FOREACH

#define BATCH_POST_OP4VECTOR( VECTOR, OPERATION ) \
	for ( size_t _i = 0; _i<VECTOR.size(); ++_i )	\
		(VECTOR[_i]) OPERATION;	\

#define Z_CC( ORIG, PARA )	PARA(ORIG.PARA)
#define Z_FCC( PARA )	Z_CC(Z_CONSTRUCTOR_ORIG,PARA)

#define PI_CONDITIONAL_REFRESH( DST, SRC, RELATION ) \
	if ( DST RELATION SRC ) \
		DST = SRC;

/* ----------------- Pair ----------------- */

//NEW_PAIR_TEMPLATE(Class_name, Element_name1, Element_name2)

#define NEW_PAIR_TEMPLATE(CNAME, ELE1, ELE2) \
	template <typename T>	\
	struct CNAME {	\
		T ELE1;	\
		T ELE2;	\
	public:	\
		CNAME() {}	\
		explicit CNAME(T _ ## ELE) : \
			ELE1(_ ## ELE), ELE2(_ ## ELE) {}	\
		explicit CNAME(T _ ## ELE1, T _ ## ELE2) : \
			ELE1(_ ## ELE1), ELE2(_ ## ELE2) {}	\
		template<typename TR>	\
		explicit CNAME(const CNAME<TR>& _ ## CNAME  ## _i) : \
			ELE1(_ ## CNAME  ## _i.ELE1), ELE2(_ ## CNAME  ## _i.ELE2) {} \
		template<typename TR>\
		const CNAME<T>& operator = (const CNAME<TR>& _ ## CNAME  ## _i) { \
			ELE1 = _ ## CNAME  ## _i.ELE1; ELE2 = _ ## CNAME  ## _i.ELE2; \
			return *this; } \
		bool operator == (const CNAME<T>& _ ## CNAME  ## _i) const { \
			return (ELE1 == _ ## CNAME ## _i.ELE1) && (ELE2 == _ ## CNAME ## _i.ELE2); } \
		~CNAME() {}	\
	};


#define NEW_PAIR_BI_OP2(CNAME0, CNAME1, ELE1_1, ELE1_2, CNAME2, ELE2_1, ELE2_2, OP) \
	template<typename T, typename T2> \
	const CNAME0<T> operator OP (const CNAME1<T>& _left, const CNAME2<T2>& _right) { \
		return CNAME0<T>(_left.ELE1_1 OP _right.ELE2_1, _left.ELE1_2 OP _right.ELE2_2);	\
	}

#define NEW_PAIR_BI_OP(CNAME, ELE1, ELE2, OP) \
	NEW_PAIR_BI_OP2(CNAME, CNAME, ELE1, ELE2, CNAME, ELE1, ELE2, OP)


#define NEW_PAIR_PRE_UNI_BI_OP2(CNAME0, CNAME1, ELE1_1, ELE1_2, OP) \
	template<typename T> \
	const CNAME0<T> operator OP (const T& _left, const CNAME1<T>& _right) { \
		return CNAME0<T>(_left OP _right.ELE1_1, _left OP _right.ELE1_2);	\
	}

#define NEW_PAIR_PRE_UNI_BI_OP(CNAME, ELE1, ELE2, OP) \
	NEW_PAIR_PRE_UNI_BI_OP2(CNAME, CNAME, ELE1, ELE2, OP)

#define NEW_PAIR_POST_UNI_BI_OP2(CNAME0, CNAME1, ELE1_1, ELE1_2, OP) \
	template<typename T> \
	const CNAME0<T> operator OP (const CNAME1<T>& _left, const T& _right) { \
		return CNAME0<T>(_left.ELE1_1 OP _right, _left.ELE1_2 OP _right);	\
	}

#define NEW_PAIR_POST_UNI_BI_OP(CNAME, ELE1, ELE2, OP) \
	NEW_PAIR_POST_UNI_BI_OP2(CNAME, CNAME, ELE1, ELE2, OP)

#define NEW_PAIR_UNI_BI_OP2(CNAME0, CNAME1, ELE1_1, ELE1_2, OP) \
	NEW_PAIR_PRE_UNI_BI_OP2(CNAME0, CNAME1, ELE1_1, ELE1_2, OP) \
	NEW_PAIR_POST_UNI_BI_OP2(CNAME0, CNAME1, ELE1_1, ELE1_2, OP)

#define NEW_PAIR_UNI_BI_OP(CNAME, ELE1, ELE2, OP) \
	NEW_PAIR_UNI_BI_OP2(CNAME, CNAME, ELE1, ELE2, OP)


NEW_PAIR_TEMPLATE(zPointT, x, y);			//Point
NEW_PAIR_BI_OP(zPointT, x, y, -);
NEW_PAIR_BI_OP(zPointT, x, y, +);
NEW_PAIR_BI_OP(zPointT, x, y, *);
NEW_PAIR_BI_OP(zPointT, x, y, /);
NEW_PAIR_UNI_BI_OP(zPointT, x, y, +);
NEW_PAIR_UNI_BI_OP(zPointT, x, y, -);
NEW_PAIR_UNI_BI_OP(zPointT, x, y, *);
NEW_PAIR_UNI_BI_OP(zPointT, x, y, /);

NEW_PAIR_TEMPLATE(zISizeT, width, height);	//Image Size
NEW_PAIR_BI_OP2(zISizeT, zISizeT, width, height, zPointT, x, y, -);
NEW_PAIR_BI_OP2(zISizeT, zISizeT, width, height, zPointT, x, y, +);
NEW_PAIR_BI_OP(zISizeT, width, height, -);
NEW_PAIR_BI_OP(zISizeT, width, height, +);
NEW_PAIR_UNI_BI_OP(zISizeT, width, height, -);
NEW_PAIR_UNI_BI_OP(zISizeT, width, height, +);

NEW_PAIR_TEMPLATE(zMSizeT, rows, cols);		//Matrix Size
NEW_PAIR_TEMPLATE(zPairT, e1, e2);			//General pair

/* ----------------- Array ----------------- */

template<typename T1, typename T2>
void zArrCopy(T1& dst, const T2& src, size_t count) {
	for (size_t i = 0; i<count; i++)
		dst[i] = src[i];
}

template<typename T1, typename T2, typename T3, typename T4>
void zArrCopy(T1& dst, const T2& src, size_t count, size_t repeat, const T3 incremental, const T4 delta) {
	size_t n = 0;
	for (size_t k=0; k<repeat; ++k) {
		T3 inc = delta + incremental * (int)k;
		for (size_t i = 0; i<count; i++) {
			dst[n++] = src[i]+inc;
		}
	}
}

template <typename T>
class zArray1D {
	size_t	_size;
	T*		_arr;
	bool	_is_replacement;
public:
	class FromIterator{};
	class WithFixedConstructor{};
	class WithoutConstruction{};	//Dangrous!!!!!!
public:
	typedef T	value_type;
	typedef zArray1D<T> MyType;
	//explicit zArray1D() {}
	explicit zArray1D(size_t size) throw() : _size(size), _arr(new T[size]), _is_replacement(false) {}
	zArray1D(const MyType&  orig) : Z_CC(orig,_size), _arr(new T[orig._size]), _is_replacement(false) {
		zArrCopy(_arr, orig._arr, orig._size);
	}
	template<typename T2>	//Should be extended for element without default constructor
	explicit zArray1D(size_t size, const T2&  orig, size_t _len = 0) throw() :
		_size(size), _arr(new T[size]), _is_replacement(false) {
		zArrCopy(_arr, orig, ( (_len) ? _len : std::min(size,orig._size)));
	}

	template<typename T2>	//Should be extended for element without default constructor
	explicit zArray1D( T2 _begin, T2 _end, const FromIterator ) throw() : _is_replacement(false) {
		_size = std::distance( _begin, _end );
		_arr = new T[_size];
		std::copy(_begin,_end,_arr);
	}

	template<typename T2>
	explicit zArray1D( size_t size, const T2& _ini_para, const WithFixedConstructor ) throw() :
		_size(size), _arr((T*)(void*)new char[sizeof(T)*size]),  _is_replacement(true) {
		for ( iterator iter=begin(); iter!=end() ; ++iter ) {
			new (iter) T(_ini_para);
		}
	}

	explicit zArray1D( size_t size, const int, const WithoutConstruction ) throw() :
		_size(size), _arr((T*)(void*)new char[sizeof(T)*size]),  _is_replacement(true) {
	}

	operator T*() {return _arr;}
	operator const T*() const {return _arr;}
	size_t size() const {return _size;}
private:
	void _release() {
		if (_is_replacement) {
			for ( iterator iter=begin(); iter!=end() ; ++iter ) {
				iter->~T();
			}
			delete[] (char*)(void*)_arr;
		} else {
			delete[] _arr;
		}
	}
public:
	~zArray1D() {
		_release();
	}

	//--------------------------
	template<typename T2>
	const MyType& operator += ( const T2& r ) {
		for ( iterator iter=begin(); iter!=end(); ++iter )
			*iter = *iter + r;
		return *this;
	}
public:
	typedef T* iterator;
	typedef const T* const_iterator;
	iterator begin() {return _arr;}
	const_iterator begin() const {return _arr;}
	iterator end() {return _arr+_size;}
	const_iterator end() const {return _arr+_size;}
	template<typename VT>
	void set( const VT& val ) {
		for(iterator iter=begin(); iter!=end(); ++iter )
			*iter = val;
	}
	const MyType& operator = ( const MyType& r ) { 
		using namespace std;
		zArrCopy(*this, r, min(_size,r._size));
		return (*this);
	}
	void Set2( const MyType& r ) {	//Similar to copy construction
		_release();
		_is_replacement = false;
		_size = r._size;
		_arr  = new T[_size];
		zArrCopy(_arr, r._arr, _size);
	}
};


typedef zISizeT<size_t>	uISize;
typedef zPointT<size_t>	uPoint;
typedef zPointT<int>	iPoint;

template <typename T>
class zArray2D {
public:
	typedef T	value_type;
	typedef zArray2D<T>		MyType;
	typedef uISize			ArraySize;
	typedef uPoint			ArrayIDX;
	typedef iPoint			Point;
private:
	ArraySize _size;
	T* _arr;
	void* _raw;
	size_t _widthStep;
	bool _is_replacement;
public:
	class WithFixedConstructor{};
public:
	//explicit zArray2D() {}
	explicit zArray2D( ArraySize size )  throw()
		: _size(size), _is_replacement(false)
	{
			_arr = new T[_size.height*_size.width];
			_raw = _arr;
			_widthStep = sizeof(T)*_size.width;
	}
	explicit zArray2D( void* data, ArraySize size, size_t widthStep = 0 )  throw()
		: _size(size), _arr(NULL), _raw(data), _is_replacement(false) {
			_widthStep = (widthStep)?(widthStep):(sizeof(T)*_size.width);
	}

	template<typename T2>
	explicit zArray2D( ArraySize size, const T2& _ini_para, const WithFixedConstructor, ArrayIDX origpoint = ArrayIDX(0,0) )  throw()
		:_size(size), _is_replacement(true) {
		size_t n=_size.height*_size.width;
		_arr = (T*)(void*)new char[sizeof(T)*n];
		_raw = _arr;
		_widthStep = sizeof(T)*_size.width;
		for ( size_t i=0; i<n ; ++i )
			new (_arr+i) T(_ini_para);
	}

	zArray2D( const MyType& orig) throw() : _size(orig._size), _is_replacement(false)
	{
			_arr = new T[_size.height*_size.width];
			_raw = _arr;
			_widthStep = sizeof(T)*_size.width;
			T* ptr = _arr;
			for ( size_t y=0; y<_size.height; ++y ) {
				zArrCopy( ptr, orig[y], _size.width);
				ptr+=_size.width;
			}
	}

	~zArray2D() {
		if (_arr) {
			if (_is_replacement) {
				size_t n=_size.height*_size.width;
				for ( size_t i=0; i<n ; ++i )
					_arr[i].~T();
				delete[] (char*)(void*)_arr;
			} else {
				delete[] _arr;
			}
		}
	}

	const T* operator [] (int row) const {
		return ((const T*)((const char*)_raw+_widthStep*(unsigned)row));
	}
	T* operator [] (int row) {
		return const_cast<T*>(((const MyType*)(this))->operator [] (row));
	}
	const T& operator () (int x, int y) const { return (*this)[y][x]; }
	T& operator () (int x, int y) { return (*this)[y][x]; }

	ArraySize size() const { return _size; }
private:
	const MyType& operator = (const MyType&) {}
};

template<typename T>
class zPtrVector {
public:
	typedef zPtrVector<T> MyType;
private:
	typedef std::vector<T*> PtrVec;
	PtrVec _vec;
public:
	zPtrVector() {}
	zPtrVector( const MyType& orig ) {
		for ( size_t i=0; i<orig._vec.size(); ++i ) {
			_vec.push_back(new T(*(orig._vec[i])));
		}
	}
	virtual void push_back(T* ptr) {
		_vec.push_back(ptr);
	}
	size_t size() const { return _vec.size(); }
	T& operator [] ( size_t idx ) { return *(_vec[idx]); }
	const T& operator [] ( size_t idx ) const { return *(_vec[idx]); }
	void clear() {
		typename PtrVec::iterator iter;
		for ( iter = _vec.begin(); iter != _vec.end(); ++iter )
			delete *iter;
	}
	virtual ~zPtrVector() {
		clear();
	}
};

/* -------------- Namespace END -------------- */

}


#endif /* ZARRAY_H_ */
