/*
 * hashspecial.h
 *
 *  Created on: 2010-10-25
 *      Author: zhangyuting
 */

#ifndef HASHSPECIAL_H_
#define HASHSPECIAL_H_

#include "config.h"

#include "zutils.h"
#include "zero_scheme.h"
#include <memory>

template <typename T>
class piPreallocatedList_Iter;
template <typename T>
class piPreallocatedList_Const_Iter;

template<typename T>
struct piPreallocatedList_ElementHull {
	typedef piPreallocatedList_ElementHull<T> MyType;
	T  val;
	MyType* next;
};

template <typename T>
class piPreallocatedList {
public:
	typedef piPreallocatedList<T>		MyType;
	typedef T							value_type;
	typedef piPreallocatedList_Iter<T>			iterator;
	typedef piPreallocatedList_Const_Iter<T>	const_iterator;
private:
	typedef piPreallocatedList_ElementHull<T>	ElementHull;
	typedef array1d<ElementHull>				PoolType;
private:
	PoolType _pool;
	ElementHull* _front;
	ElementHull* _end;
	size_t _size;
	size_t _max_size;
public:
	explicit piPreallocatedList(size_t pool_size);
	explicit piPreallocatedList(const MyType& _orig );
	~piPreallocatedList() {}
	void clear();
	void erase( iterator iter );
	iterator push_back( const T& val );
	const T& operator [] ( size_t idx) const { return _pool[idx].val; } //Not safe
	size_t size() const { return _size; }
	size_t max_size() const { return _max_size; };
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
private:
	void _initialize();
};

template <typename T>
class piPreallocatedList_Iter {
public:
	typedef piPreallocatedList_Iter<T>		MyType;
	typedef piPreallocatedList_Const_Iter<T>		ConstMyType;
	typedef piPreallocatedList<T>			ListType;
	typedef piPreallocatedList_ElementHull<T> ElementHull;
private:
	ElementHull* _eh;
public:
	friend class piPreallocatedList<T>;
	piPreallocatedList_Iter() {}
	piPreallocatedList_Iter(ElementHull* eh) : _eh(eh) {}
	piPreallocatedList_Iter(ElementHull& eh) : _eh(&eh) {}
	piPreallocatedList_Iter(const MyType& _orig) : _eh(_orig._eh) {}
	const MyType& operator++() { _eh = _eh->next ; return *this; }
	const MyType& operator++( int ) { MyType before(*this); _eh = _eh->next ; return before; }
	T& operator*() const { return (_eh->next)->val; }
	T* operator->() const {
		return &((_eh->next)->val);
	}
	bool is_null() const { return !_eh; }
	const MyType& operator = ( const MyType& r ) {_eh=r._eh; return *this;}
	const MyType& operator = ( ElementHull* eh ) {_eh=eh; return *this;}
	bool operator == ( const MyType& r ) {return (_eh==r._eh);}
	bool operator == ( const ConstMyType& r ) {return (_eh==r._eh);}
	bool operator != ( const MyType& r ) {return (_eh!=r._eh);}
	bool operator != ( const ConstMyType& r ) {return (_eh!=r._eh);}
	void vecAdvance() { ++_eh; }	//Dangerous
	ElementHull* prior() const { return _eh;}
	friend class piPreallocatedList_Const_Iter<T>;
};

template <typename T>
class piPreallocatedList_Const_Iter {
public:
	typedef piPreallocatedList_Const_Iter<T>		MyType;
	typedef piPreallocatedList_Iter<T>		NonconstMyType;
	typedef piPreallocatedList<T>			ListType;
	typedef piPreallocatedList_ElementHull<T> ElementHull;
private:
	ElementHull* _eh;
public:
	explicit piPreallocatedList_Const_Iter() {}
	explicit piPreallocatedList_Const_Iter(ElementHull* eh) : _eh(eh) {}
	explicit piPreallocatedList_Const_Iter(ElementHull& eh) : _eh(&eh) {}
	piPreallocatedList_Const_Iter(const MyType& _orig) : _eh(_orig._eh) {}
	piPreallocatedList_Const_Iter(const NonconstMyType& _orig) : _eh(_orig._eh) {}
	const MyType& operator++() { _eh = _eh->next ; return (*this); }
	const MyType& operator++( int ) { MyType before(*this); _eh = _eh->next ; return before; }
	const T& operator*() const { return (_eh->next)->val; }
	const T* operator->() const {
		return &((_eh->next)->val);
	}
	bool is_null() const { return !_eh; }
	const MyType& operator = ( const MyType& r ) {_eh=r._eh; return *this;}
	const MyType& operator = ( const NonconstMyType& r ) {_eh=r._eh; return *this;}
	const MyType& operator = ( ElementHull* eh ) {_eh=eh; return *this;}
	bool operator == ( const MyType& r ) {return (_eh==r._eh);}
	bool operator == ( const NonconstMyType& r ) {return (_eh==r._eh);}
	bool operator != ( const MyType& r ) {return (_eh!=r._eh);}
	bool operator != ( const NonconstMyType& r ) {return (_eh!=r._eh);}
	void vecAdvance() { ++_eh; }	//Dangerous
	ElementHull* prior() const { return _eh;}
};


/* ---------------------------------------------------------------------------- */

template <typename T, class ZERO_SCHEME = piZeroScheme::Basic<T> >
class piHashSpecialBase {
public:
	typedef piHashSpecialBase<T,ZERO_SCHEME> MyType;
	struct ListElement {
		T val;
		size_t idx;
	};
public:
	typedef piPreallocatedList<ListElement>		ListType;
	typedef typename ListType::iterator			iterator;
	typedef typename ListType::const_iterator	const_iterator;
	typedef array1d<iterator>			TableType;
protected:
	ZERO_SCHEME _zero_scheme;
	TableType& _tab;
	ListType& _l;
public:
	void _erase( iterator& iter );
public:
	//This constructor is not safe, in that tab and l may be inconsisit with each other.
	explicit piHashSpecialBase( TableType& tab, ListType& l, const ZERO_SCHEME zero_scheme = ZERO_SCHEME() ) 
		: _zero_scheme(zero_scheme), _tab(tab), _l(l) {}
	piHashSpecialBase( const MyType& _orig );
	void clear ();
	void clearTab ();	//Dangrous!!!!
	const T& ReadTable(size_t idx) const;
	size_t size() const {return _tab.size();}	//table size
	size_t nonzero_size() const { return _l.size(); }	//list size
	size_t max_size() const { return _l.max_size(); }	//maximum list size
	void WriteTable(size_t idx, const T& value);
	void inc_s( size_t idx, const T& val );	//safe, val can be negative
	void inc( size_t idx, const T& val );	//unsafe, val must be postive
	void dec( size_t idx, const T& val );	//unsafe, val must be postive
	const T& operator [] ( size_t idx ) const { return ReadTable(idx); }
	iterator begin() {return _l.begin();}
	const_iterator begin() const {return _l.begin();}
	iterator end() {return _l.end();}
	const_iterator end() const {return _l.end();}
	void insert(const_iterator _start, const_iterator _end);
	const MyType& operator = ( const MyType& r );
	bool operator == ( const MyType& r ) const;
	template< class TA >
	bool is_equal2arr( const TA& arr ) const;
};

template <typename T, class ZERO_SCHEME = piZeroScheme::Basic<T> >
class piHashSpecial : public piHashSpecialBase<T,ZERO_SCHEME> {
public:
	typedef piHashSpecial<T,ZERO_SCHEME>	MyType;
	typedef piHashSpecialBase<T,ZERO_SCHEME> BaseType;
	typedef typename BaseType::TableType	TableType;
	typedef typename BaseType::ListType		ListType;
	typedef typename BaseType::iterator    		iterator;
	typedef typename BaseType::const_iterator	const_iterator;
private:
	using BaseType::_tab;
	using BaseType::_l;
	typedef std::auto_ptr<TableType>	TablePtr;
	typedef std::auto_ptr<ListType>		ListPtr;
	TablePtr p_tab;
	ListPtr  p_l;
private:
	TableType*	_create_table( size_t size ) { return (new TableType(size));}
	ListType*	_create_list( size_t size ) { return (new ListType(size));}
public:
	explicit piHashSpecial( size_t table_size, size_t list_size = 0, const ZERO_SCHEME zero_scheme = ZERO_SCHEME() );
	piHashSpecial( const MyType& _orig );
	const MyType& operator = ( const MyType& r );
	bool operator == ( const MyType& r ) const;
private:
	void clearTab() {}
};

#include "hashspecial_t.hpp"

#endif
