/*
 * Copyright (c) 2009, Politecnico di Milano
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
 * 3. Neither the name of the Politecnico di Milano nor the names of its
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

/*!
 * @file P1_triangle.h
 * @brief The Lagrange base of order 1 on the triangle
 */

#ifndef GAS_FUNCTIONAL_BASE_P1_TRIANGLE_H
#define GAS_FUNCTIONAL_BASE_P1_TRIANGLE_H

#include "base"
#include "../../geometry/unit/unit"
#include "../../gas"

namespace gas { namespace functional { namespace base {

/*!
 * @brief The Lagrange base of order 1 on the triangle
 * @see gas::geometry::unit::triangle
 *
 * The base is defined by \f$\varphi_i:\hat{\mathcal{T}}\to[0,1]\f$, where
 * \f[
 * \varphi_0(X,Y) = 1-X-Y \qquad
 * \varphi_1(X,Y) = X \qquad
 * \varphi_2(X,Y) = Y
 * \f]
 */
template <>
class P1<gas::geometry::unit::triangle> {

public:
	/*! @brief The basic shape on which is defined */
	typedef gas::geometry::unit::triangle unit_t;

	/*! @brief The number of function */
	static unsigned int const n = 3u;

public:
	/*!
	 * @brief The base function
	 * @param i The index of basis function
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The evaluation of i-th base function in (X,Y)
	 */
	static inline double b (int const i, double const X, double const Y) {
		GAS_ASSERT(unit_t::in(X, Y)); // The point must be in the triangle
		GAS_ASSERT(i < n);            // A valid index
		switch (i) {
		case 0: return (1. - X - Y);
		case 1: return X;
		case 2: return Y;
		}
		return 0.;
	}

	/*!
	 * @brief The derivative along the first coordinate of base function
	 * @param i The index of basis function
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The evaluation of i-th base function in (X,Y)
	 */
	static inline double dbdX (int const i, double const X, double const Y) {
		GAS_ASSERT(unit_t::in(X, Y)); // The point must be in the triangle
		GAS_ASSERT(i < n);            // A valid index
		switch (i) {
		case 0: return -1.;
		case 1: return 1.;
		case 2: return 0.;
		}
		return 0.;
	}

	/*!
	 * @brief The derivative along the second coordinate of base function
	 * @param i The index of basis function
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The evaluation of i-th base function in (X,Y)
	 */
	static inline double dbdY (int const i, double const X, double const Y) {
		GAS_ASSERT(unit_t::in(X, Y)); // The point must be in the triangle
		GAS_ASSERT(i < n);            // A valid index
		switch (i) {
		case 0: return -1.;
		case 1: return 0.;
		case 2: return 1.;
		}
		return 0.;
	}

};

} } }

#endif // GAS_FUNCTIONAL_BASE_P1_TRIANGLE_H
