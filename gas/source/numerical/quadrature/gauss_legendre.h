/*
 * Copyright (c) 2008, Politecnico di Milano
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
 * @file gauss_legendre.h
 * @brief The Guass-Legendre quadrature formulae
 */

#ifndef _gas_numerical_quadrature_gauss_legendre_
#define _gas_numerical_quadrature_gauss_legendre_

#include "../../geometry/unit/interval.h"
#include "formula.h"
#include "method.h"

namespace gas { namespace numerical { namespace quadrature {

/*!
 * @brief The Gauss-Legendre quadrature formula
 * @param unit_ The unit geometry on which is defined
 * @param nodes_ The number of nodes
 */
template <typename unit_, unsigned int nodes_>
class gauss_legendre;

/*!
 * @brief The Gauss-Legendre formula for interval
 * @param nodes_ The number of nodes
 */
template <unsigned int nodes_>
struct gauss_legendre<gas::geometry::unit::interval, nodes_>:
	public method_1<gas::geometry::unit::interval, nodes_, gauss_legendre<gas::geometry::unit::interval, nodes_> > {

	/*! @brief The coordinates of nodes */
	static double const x_[nodes_];

	/*! @brief The weights of nodes */
	static double const w_[nodes_];

	template <typename method_, typename map_>
	friend class formula;

};

/*! @brief Coordinates for 2 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 2u>::x_[2] = {
	-0.5773502691896257645091488, +0.5773502691896257645091488
};

/*! @brief Weights for 2 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 2u>::w_[2] = {
	1.0000000000000000000000000, 1.0000000000000000000000000
};


/*! @brief Coordinates for 3 nodes formula */
template<>
double const gauss_legendre<gas::geometry::unit::interval, 3u>::x_[3] = {
	-0.7745966692414833770358531, 0.0000000000000000000000000, +0.7745966692414833770358531
};

/*! @brief Weights for 3 nodes formula */
template<>
double const gauss_legendre<gas::geometry::unit::interval, 3u>::w_[3] = {
	0.5555555555555555555555556, 0.8888888888888888888888889, 0.5555555555555555555555556
};

} } }

#endif // _gas_numerical_quadrature_gauss_legendre_
