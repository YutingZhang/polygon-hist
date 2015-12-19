/*
 * zero_scheme.h
 *
 *  Created on: 2010-10-26
 *      Author: zhangyuting
 */

#ifndef ZERO_SCHEME_H_
#define ZERO_SCHEME_H_

#include "config.h"

#include<limits>

class piZeroScheme {
public:
	template <typename T>
	class Basic {
		const T _zero;
	public:
		Basic() : _zero(T(0)) {}
		Basic(const Basic& _orig) : _zero(T(0)) {}
		bool is_zero(const T& a) const {
			return (a==T(0));
		}
		bool is_zero_p(const T& a) const {	//only concerns positive
			return (a==T(0));
		};
		bool is_zero_n(const T& a) const {	//only concerns negative
			return (a==T(0));
		};
		bool is_equal(const T& a, const T& b) const {
			return (a==b);
		}
		const T& zero() const {
			return _zero;
		}
	};


	template <typename T>
	class Epsilon {
		const T _zero;
		const T _epsilon;
	public:
		Epsilon() : _zero(T(0)), _epsilon(piFloatEpsilon) {}
		Epsilon(T epsilon) : _zero(T(0)), _epsilon(epsilon) {}
		Epsilon(const Epsilon& _orig) : _zero(T(0)), _epsilon(_orig._epsilon) {}
		bool is_zero(const T& a) const {
			return (a<T(0)+_epsilon && a>T(0)-_epsilon);
		};
		bool is_zero_p(const T& a) const {	//only concerns positive
			return (a<T(0)+_epsilon);
		};
		bool is_zero_n(const T& a) const {	//only concerns negative
			return (a>T(0)-_epsilon);
		};
		bool is_equal(const T& a, const T& b) const {
			return is_zero(a-b);
		}
		const T& zero() const {
			return _zero;
		}
	};
};

#endif

