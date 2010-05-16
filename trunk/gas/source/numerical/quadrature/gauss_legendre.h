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
 * @file gauss_legendre.h
 * @brief The Guass-Legendre quadrature formulae
 */

#ifndef GAS_NUMERICAL_QUADRATURE_GAUSS_LEGENDRE_H
#define GAS_NUMERICAL_QUADRATURE_GAUSS_LEGENDRE_H

#include "../../geometry/unit/unit"
#include "quadrature"

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
class gauss_legendre<gas::geometry::unit::interval, nodes_>:
	public method_1<gas::geometry::unit::interval, nodes_, gauss_legendre<gas::geometry::unit::interval, nodes_> > {

private:
	/*! @brief The coordinates of nodes */
	static double const x_[nodes_];

	/*! @brief The weights of nodes */
	static double const w_[nodes_];

	/*! @brief The basic shape on which is defined */
	typedef gas::geometry::unit::interval unit_t;

	/*! @brief The method */
	typedef gauss_legendre<gas::geometry::unit::interval, nodes_> method_t;

	/*! @brief The number of nodes */
	static int const n_ = nodes_;

	/*! @brief Degree of exactness */
	static int const degree_ = 2u * nodes_ - 1u;

	template <typename type__> friend class info;

	friend class method_1<gas::geometry::unit::interval, nodes_,
		gauss_legendre<gas::geometry::unit::interval, nodes_> >;

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


/*! @brief Coordinates for 4 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 4u>::x_[4] = {
	-0.8611363115940525752239465, -0.3399810435848562648026658, +0.3399810435848562648026658,
	+0.8611363115940525752239465
};

/*! @brief Weights for 4 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 4u>::w_[4] = {
	0.3478548451374538573730639, 0.6521451548625461426269361, 0.6521451548625461426269361,
	0.3478548451374538573730639
};


/*! @brief Coordinates for 5 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 5u>::x_[5] = {
	-0.9061798459386639927976269, -0.5384693101056830910363144, 0.0000000000000000000000000,
	+0.5384693101056830910363144, +0.9061798459386639927976269
};

/*! @brief Weights for 5 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 5u>::w_[5] = {
	0.2369268850561890875142640, 0.4786286704993664680412915, 0.5688888888888888888888889,
	0.4786286704993664680412915, 0.2369268850561890875142640
};


/*! @brief Coordinates for 6 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 6u>::x_[6] = {
	-0.9324695142031520278123016, -0.6612093864662645136613996, -0.2386191860831969086305017,
	+0.2386191860831969086305017, +0.6612093864662645136613996, +0.9324695142031520278123016
};

/*! @brief Weights for 6 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 6u>::w_[6] = {
	0.1713244923791703450402961, 0.3607615730481386075698335, 0.4679139345726910473898703,
	0.4679139345726910473898703, 0.3607615730481386075698335, 0.1713244923791703450402961
};


/*! @brief Coordinates for 7 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 7u>::x_[7] = {
	-0.9491079123427585245261897, -0.7415311855993944398638648, -0.4058451513773971669066064,
	 0.0000000000000000000000000, +0.4058451513773971669066064, +0.7415311855993944398638648,
	+0.9491079123427585245261897
};

/*! @brief Weights for 7 nodes formula */
template <>
double const gauss_legendre<gas::geometry::unit::interval, 7u>::w_[7] = {
	0.1294849661688696932706114, 0.2797053914892766679014678, 0.3818300505051189449503698,
	0.4179591836734693877551020, 0.3818300505051189449503698, 0.2797053914892766679014678,
	0.1294849661688696932706114
};

} } }

#endif // GAS_NUMERICAL_QUADRATURE_GAUSS_LEGENDRE_H
