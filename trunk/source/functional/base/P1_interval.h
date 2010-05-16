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
 * @file P1_interval.h
 * @brief The Lagrange base of order 1 on the interval
 */

#ifndef GAS_FUNCTIONAL_BASE_P1_INTERVAL_H
#define GAS_FUNCTIONAL_BASE_P1_INTERVAL_H

#include "../../geometry/unit/unit"
#include "../../gas"

namespace gas { namespace functional { namespace base {

/*!
 * @brief The Lagrange base of order 1 on the interval
 * @see gas::geometry::unit::interval
 *
 * The base is defined by \f$\varphi_i:[-1,1]\to[0,1]\f$, where
 * \f[
 * \varphi_0(X) = \frac{1-X}{2} \qquad
 * \varphi_1(X) = \frac{1+X}{2}
 * \f]
 */
template <>
class P1<gas::geometry::unit::interval> {

public:
	/*! @brief The basic shape on which is defined */
	typedef gas::geometry::unit::interval unit_t;

	/*! @brief The number of function */
	static unsigned int const n = 2u;

public:
	/*!
	 * @brief The base function
	 * @param i The index of basis function
	 * @param X The coordinate
	 * @return The evaluation of i-th base function in X
	 */
	static inline double b (int const i, double const X) {
		GAS_ASSERT(unit_t::in(X)); // The point must be in the interval
		GAS_ASSERT(i < n);         // A valid index
		switch (i) {
		case 0: return (1. - X) * 0.5;
		case 1: return (1. + X) * 0.5;
		}
		return 0.;
	}

	/*!
	 * @brief The derivative of base function
	 * @param i The index of basis function
	 * @param X The coordinate
	 * @return The evaluation of i-th base function in X
	 */
	static inline double dbdX (int const i, double const X) {
		GAS_ASSERT(unit_t::in(X)); // The point must be in the interval
		GAS_ASSERT(i < n);         // A valid index
		switch (i) {
		case 0: return -0.5;
		case 1: return +0.5;
		}
		return 0.;
	}

};

} } }

#endif // GAS_FUNCTIONAL_BASE_P1_INTERVAL_H
