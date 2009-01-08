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

#ifndef _GAS_MATH_H_ /* BEGIN _GAS_MATH_H_ */
#define _GAS_MATH_H_

#include <cmath>
#include <complex>

/** @namespace Math
 *  @brief Contains some useful mathematical functions to use in your program **/
namespace Common { namespace Math {
	/** Identity **/
	template<typename T> inline T Id(T const &x) { return x; }

	/** Four basic operations **/
	template<typename T> inline T Sum(T const &a, T const &b) { return a + b; }
	template<typename T> inline T Sub(T const &a, T const &b) { return a - b; }
	template<typename T> inline T Mul(T const &a, T const &b) { return a * b; }
	template<typename T> inline T Div(T const &a, T const &b) { return a / b; }

	/** Division mod **/
	template<typename T> inline T Mod(T const &a, T const &b) { return a % b; }

	/** Max **/
	template<typename T> inline T Max(T const &a, T const &b) { return ((a < b) ? b : a); }
	/** Min **/
	template<typename T> inline T Min(T const &a, T const &b) { return ((a < b) ? a : b); }

	/** Absolute value **/
	template<typename T> inline T Abs(T const &x) { return std::abs(x); };

	/** Complex Conjugate **/
	template<typename T> inline T Conj(T const &x) { return x; };
	template<typename T> inline std::complex<T> Conj(std::complex<T> const &x) { return std::conj(x); };

	/** Square Root **/
	template<typename T> inline T Sqrt(T const &x) { return std::sqrt(x); }
}}

#endif /* END _GAS_MATH_H_ */
