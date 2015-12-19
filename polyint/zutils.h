/*
 * zutils.h
 *
 *  Created on: 2010-10-21
 *      Author: zhangyuting
 */

#ifndef ZUTILS_H_
#define ZUTILS_H_

#include "zutil.h"


typedef zUtilities_PolyHistInternal::zPointT<int> iPoint;
typedef zUtilities_PolyHistInternal::zISizeT<int> iISize;
typedef zUtilities_PolyHistInternal::zMSizeT<int> iMSize;

#define pair_t zUtilities_PolyHistInternal::zPairT
typedef pair_t<int>	iSpan;

#define	array1d	zUtilities_PolyHistInternal::zArray1D
#define	array2d	zUtilities_PolyHistInternal::zArray2D
using zUtilities_PolyHistInternal::uISize;
using zUtilities_PolyHistInternal::uPoint;

#define arrcpy	zUtilities_PolyHistInternal::zArrCopy

//#define array1d_rotref_const zUtilities_PolyHistInternal::zArray1D_RotRef_const
//#define array1d_rotref_delta zUtilities_PolyHistInternal::zArray1D_RotRef_Delta

#define ptrvec zUtilities_PolyHistInternal::zPtrVector

#endif /* ZUTILS_H_ */
