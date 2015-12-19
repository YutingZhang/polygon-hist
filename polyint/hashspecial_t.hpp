/*
 * hashspecial_t.hpp
 *
 *  Created on: 2010-10-25
 *      Author: zhangyuting
 */

#include <cstring>

#define Z_CONSTRUCTOR_ORIG	_orig

template <typename T>
inline piPreallocatedList<T>::piPreallocatedList(size_t pool_size) :
	_pool(pool_size+1), _size(0), _max_size(pool_size) {
	_initialize();
}



template <typename T>
inline piPreallocatedList<T>::piPreallocatedList(const MyType& _orig ) :
	_pool(_pool.size()), Z_FCC(_size), Z_FCC(_max_size) {
	_initialize();
	for ( const_iterator iter=_orig.begin();
			iter != _orig.end(); ++iter ) {
		push_back(*iter);
	}
}

template <typename T>
inline void piPreallocatedList<T>::_initialize() {
	_front = _pool.begin();
	_end   = _pool.begin();
	ElementHull* ptr = _front;
	for ( size_t i=0; i<_max_size; ++i ) {
		ElementHull* next = ptr+1;
		ptr->next = next;
		ptr = next;
	}
	ptr->next = _front;
}

template <typename T>
inline void piPreallocatedList<T>::clear() {
	_size = 0;
	_end = _front;
}

template <typename T>
inline piPreallocatedList_Iter<T> piPreallocatedList<T>::push_back( const T& val ) {
	iterator back(_end);
	_end = _end->next;
	_end->val = val;
	++_size;
	return back;
}

template <typename T>
inline void	piPreallocatedList<T>::erase( iterator iter ) {
	ElementHull* t = iter._eh->next;
	if ( t == _end ) {
		_end = iter._eh;
	} else {
		iter._eh->next = t->next;
		t->next  = _end->next;
		_end->next = t;
	}

	--_size;
}

template <typename T>
inline piPreallocatedList_Iter<T>	 piPreallocatedList<T>::begin() {
	return iterator(_front);
}

template <typename T>
inline piPreallocatedList_Const_Iter<T>	 piPreallocatedList<T>::begin() const {
	return const_iterator(_front);
}

template <typename T>
inline piPreallocatedList_Iter<T>	 piPreallocatedList<T>::end() {
	return iterator(_end);
}

template <typename T>
inline piPreallocatedList_Const_Iter<T> piPreallocatedList<T>::end() const {
	return const_iterator(_end);
}


//--------------------------------------------------------------------------------


template <typename T, class ZERO_SCHEME >
piHashSpecialBase<T,ZERO_SCHEME>::piHashSpecialBase( const MyType& _orig ) :
	Z_FCC(_zero_scheme), Z_FCC(_tab), Z_FCC(_l) {}

template <typename T, class ZERO_SCHEME >
void piHashSpecialBase<T,ZERO_SCHEME>::clear () {
	clearTab();
	_l.clear();
}

template <typename T, class ZERO_SCHEME >
void piHashSpecialBase<T,ZERO_SCHEME>::clearTab () {
	iterator nilIter(NULL);
	if ( _l.size() < (_tab.size()>>1) ) {	//divide by 2
		for ( iterator iter=begin(); iter!=end(); ++iter )
			_tab[iter->idx] = nilIter;
	} else {
		_tab.set(nilIter);
	}
}

template <typename T, class ZERO_SCHEME>
inline void piHashSpecialBase<T,ZERO_SCHEME>::_erase( iterator& iter ) {
	iterator iter2 = iter;
	iter = NULL;
	_l.erase(iter2);
	iterator iterp(iter2.prior());
	if (iterp!=end())
		_tab[iterp->idx] = iterp;
}


template <typename T, class ZERO_SCHEME>
inline void piHashSpecialBase<T,ZERO_SCHEME>::inc( size_t idx, const T& val )  {
#if 0
	iterator& iter = _tab[idx];
	if ( iter.is_null() ) {
		ListElement e;
		e.idx = idx;
		e.val = val;
		iter = _l.push_back(e);
	}
	else {
		iter->val += val;
	}
#else
	inc_s(idx,val);
#endif
}

template <typename T, class ZERO_SCHEME>
inline void piHashSpecialBase<T,ZERO_SCHEME>::inc_s( size_t idx, const T& val )  {
	iterator& iter = _tab[idx];
	if ( iter.is_null() ) {
		ListElement e;
		e.idx = idx;
		e.val = val;
		iter = _l.push_back(e);
	}
	else {
		iter->val += val;
		if ( _zero_scheme.is_zero(iter->val) )
			_erase(iter);
	}
}

template <typename T, class ZERO_SCHEME>
inline void piHashSpecialBase<T,ZERO_SCHEME>::dec( size_t idx, const T& val )  {
#if 0
	iterator& iter = _tab[idx];
	iter->val -= val;
	if ( _zero_scheme.is_zero_p(iter->val) )	//Only concerns about positive
		_erase(iter);
#else
	inc_s(idx,-val);
#endif
}

template <typename T, class ZERO_SCHEME>
void piHashSpecialBase<T,ZERO_SCHEME>::WriteTable( size_t idx, const T& val )  {
	iterator& iter = _tab[idx];
	if ( iter.is_null() ) {
		if ( !_zero_scheme.is_zero(val) ) {
			ListElement e;
			e.idx = idx;
			e.val = val;
			iter = _l.push_back(e);
		}
	}
	else {
		if ( _zero_scheme.is_zero(val) ) {
			_erase(iter);
		} else {
			iter->val = val;
		}
	}
}

template <typename T, class ZERO_SCHEME>
const T& piHashSpecialBase<T,ZERO_SCHEME>::ReadTable( size_t idx ) const {
	const iterator& iter = _tab[idx];
	if ( iter.is_null() )
		return _zero_scheme.zero();
	return iter->val;
}

template <typename T, class ZERO_SCHEME>
void piHashSpecialBase<T,ZERO_SCHEME>::insert(const_iterator _start, const_iterator _end) {
	for ( const_iterator iter=_start;
			iter != _end; ++iter ) {
		WriteTable(iter->idx, iter->val);
	}
}

template <typename T, class ZERO_SCHEME>
const piHashSpecialBase<T,ZERO_SCHEME>& piHashSpecialBase<T,ZERO_SCHEME>::operator = ( const MyType& r ) {
	clear();
	insert( r.begin(), r.end() );
	return *this;
}

template <typename T, class ZERO_SCHEME>
bool piHashSpecialBase<T,ZERO_SCHEME>::operator == ( const MyType& r ) const {
	if (r.nonzero_size()!=nonzero_size())
		return false;
	for ( typename MyType::const_iterator iter = r.begin();
			iter!=r.end(); ++iter ) {
		if ( !_zero_scheme.is_equal(iter->val, ReadTable(iter->idx)) )
			return false;
	}
	return true;
}


template <typename T, class ZERO_SCHEME>
template< class TA >
bool piHashSpecialBase<T,ZERO_SCHEME>::is_equal2arr( const TA& arr ) const {
	for ( size_t i=0; i<size(); ++i ) {
		if ( !_zero_scheme.is_equal(arr[i], ReadTable(i)) )
			return false;
	}
	return true;
}


//-----------------------------------
template <typename T, class ZERO_SCHEME >
piHashSpecial<T,ZERO_SCHEME>::piHashSpecial( size_t table_size, size_t list_size, const ZERO_SCHEME zero_scheme) :
	BaseType(*_create_table(table_size), *_create_list((list_size)?(list_size):(table_size)), zero_scheme),
		p_tab(&_tab), p_l(&_l) {
	_tab.set(iterator(NULL));
}

template <typename T, class ZERO_SCHEME >
piHashSpecial<T,ZERO_SCHEME>::piHashSpecial( const MyType& _orig ) 
	: BaseType(*_create_table(_orig.size()), *_create_list(_orig.max_size()), _orig._zero_scheme), 
		p_tab(&_tab), p_l(&_l) {
	_tab.set(iterator(NULL));
	this->insert( _orig.begin(), _orig.end() );
}

template <typename T, class ZERO_SCHEME>
const piHashSpecial<T,ZERO_SCHEME>& piHashSpecial<T,ZERO_SCHEME>::operator = ( const MyType& r ) {
	BaseType::operator = ( r );
	return *this;
}

template <typename T, class ZERO_SCHEME>
bool piHashSpecial<T,ZERO_SCHEME>::operator == ( const MyType& r ) const {
	return BaseType::operator ==(r);
}


#undef Z_CONSTRUCTOR_ORIG
