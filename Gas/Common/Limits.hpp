/*
 * Copyright (c) 2008, Davide Ferrarese & Mattia Penati
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the <ORGANIZATION> nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _GAS_LIMITS_H_
#define _GAS_LIMITS_H_

#include "Math.hpp"

namespace Common {
	/** Limits class **/
	template<typename T>
	struct Limits {
		static T const Epsilon = 0;
		static T const Max = 0;
		static T const Min = 0;
		static inline bool Equal(T const &a, T const &b) { return Common::Math::Abs(a - b) == 0; }
		static inline bool IsZero(T const &a) { return Common::Limits<T>::Equal(a, 0); }
	};
	/** Limits class specialization: float **/
	template<>
	struct Limits<float> {
		static float const Epsilon = __FLT_EPSILON__;
		static float const Max = __FLT_MAX__;
		static float const Min = __FLT_MIN__;
		static inline bool Equal(float const &a, float const &b) { return Common::Math::Abs(a - b) <= 8 * Limits<float>::Epsilon; }
		static inline bool IsZero(float const &a) { return Common::Limits<float>::Equal(a, 0.); }
	};
	/** Limits class specialization: double **/
	template<>
	struct Limits<double> {
		static double const Epsilon = __DBL_EPSILON__;
		static double const Max = __DBL_MAX__;
		static double const Min = __DBL_MIN__;
		static inline bool Equal(double const &a, double const &b) { return Common::Math::Abs(a - b) <= 8 * Limits<double>::Epsilon; }
		static inline bool IsZero(double const &a) { return Common::Limits<double>::Equal(a, 0); }
	};
	/** Limits class specialization: long double **/
	template<>
	struct Limits<long double> {
		static long double const Epsilon = __LDBL_EPSILON__;
		static long double const Max = __LDBL_MAX__;
		static long double const Min = __LDBL_MIN__;
		static inline bool Equal(long double const &a, long double const &b) { return Common::Math::Abs(a - b) <= 8 * Limits<long double>::Epsilon; }
		static inline bool IsZero(long double const &a) { return Common::Limits<long double>::Equal(a, 0); }
	};
}

#endif
