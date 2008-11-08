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

#ifndef _GAS_INTEGRATION1D_H_
#define _GAS_INTEGRATION1D_H_

#include <cstdlib>
#include <iostream>
#include "Gas/Common/Common.h"

namespace Integration {
	namespace Integrator {
		/** Virtual class for 1-dimensional numerical integration **/
		template<typename T, T F(T)>
		class Integrator1D {
			public:
				virtual T Integrate(T const, T const) = 0;
		};

		/** The method of rectangles or middle points for numerical integration **/
		template<typename T, T F(T)>
		class Rectangle: virtual public Integrator1D<T, F> {
			private:
				size_t N;
			public:
				Rectangle(): N(30) {};
				Rectangle(size_t n): N(n) {};
				T Integrate(T const a, T const b) {
					if (a == b) return 0;
					T h = (b - a) / N;
					T I = 0;
					range(i, 0, N) I += F(a + (i*h) + (h/2));
					return h * I;
				}
		};
		
		/** The method of trapeziums for numerical integration **/
		template<typename T, T F(T)>
		class Trapezium: virtual public Integrator1D<T, F> {
			private:
				size_t N;
			public:
				Trapezium(): N(30) {};
				Trapezium(size_t n): N(n) {};
				T Integrate(T const a, T const b) {
					if (a == b) return 0;
					T h = (b - a) / N;
					T I = 0;
					range(i, 0, N) I += F(a + (i * h)) + F(a + ((i + 1) * h));
					return h * I / 2;
				}
		};

		/** The method of Cavalieri-Simpson for numerical integration **/
		template<typename T, T F(T)>
		class Simpson: virtual public Integrator1D<T, F> {
			private:
				size_t N;
			public:
				Simpson(): N(30) {};
				Simpson(size_t n): N(n) {};
				T Integrate(T const a, T const b) {
					if (a == b) return 0;
					T h = (b - a) / N;
					T I = 0;
					range(i, 0, N) I += F(a + (i * h)) + 4 * F(a + (h / 2)+(i * h))+F(a + ((i + 1) * h));
					return h * I / 6;
				}
		};
	}

	template<typename T, T F(T)>
	T Integrate(T, T, Integrator::Integrator1D<T, F> &);

	template<typename T, T F(T)>
	T Integrate(T, T);
}

template<typename T, T F(T)>
T Integration::Integrate(T a, T b, Integrator::Integrator1D<T, F> &i) {
	if (a > b) std::swap(a, b);
	return i.Integrate(a, b);
}

template<typename T, T F(T)>
T Integration::Integrate(T a, T b) {
	return Integration::Integrate(a, b, Integration::Integrator::Rectangle<T, F>());
}

#endif
