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
 * @file newton_cotes.h
 * @brief The Newton-Cotes quadrature formulae
 */

#ifndef _gas_numerical_quadrature_newton_cotes_
#define _gas_numerical_quadrature_newton_cotes_

#include "../../geometry/unit/unit"
#include "quadrature"

namespace gas { namespace numerical { namespace quadrature {

/*!
 * @brief The Newton-Cotes quadrature formula
 * @param unit_ The unit geometry on which is defined
 * @param nodes_ The number of nodes
 */
template <typename unit_, unsigned int nodes_>
class newton_cotes;

/*!
 * @brief The Gauss-Legendre formula for triangle
 * @param nodes_ The number of nodes
 */
template <unsigned int nodes_>
class newton_cotes<gas::geometry::unit::triangle, nodes_>:
	public method_2<gas::geometry::unit::triangle, nodes_, newton_cotes<gas::geometry::unit::triangle, nodes_> > {

private:
	/*! @brief The coordinates of nodes */
	static double const x_[nodes_];

	/*! @brief The coordinates of nodes */
	static double const y_[nodes_];

	/*! @brief The weights of nodes */
	static double const w_[nodes_];

	/*! @brief The basic shape on which is defined */
	typedef gas::geometry::unit::triangle unit_t;

	/*! @brief The method */
	typedef newton_cotes<gas::geometry::unit::triangle, nodes_> method_t;

	/*! @brief The number of nodes */
	static unsigned int const n_ = nodes_;

	/*! @brief Degree of exactness */
	static unsigned int const degree_;

	template <typename type__> friend class info;

	template<typename method_, typename map_>
	friend class formula;

	friend class method_2<gas::geometry::unit::triangle, nodes_,
		newton_cotes<gas::geometry::unit::triangle, nodes_> > ;

};


/*! @brief Coordinates for 3 nodes formula */
template <>
double const newton_cotes<gas::geometry::unit::triangle, 3u>::x_[3] = {
	0.1666666666667, 0.6666666666667, 0.1666666666667
};

/*! @brief Coordinates for 3 nodes formula */
template <>
double const newton_cotes<gas::geometry::unit::triangle, 3u>::y_[3] = {
	0.6666666666667, 0.1666666666667, 0.1666666666667
};

/*! @brief Weights for 3 nodes formula */
template <>
double const newton_cotes<gas::geometry::unit::triangle, 3u>::w_[3] = {
	0.1666666666667, 0.1666666666667, 0.1666666666667
};

/*! @brief Degree of exactness */
template <>
unsigned int const newton_cotes<gas::geometry::unit::triangle, 3u>::degree_ = 2u;


/*! @brief Coordinates for 6 nodes formula */
template <>
double const newton_cotes<gas::geometry::unit::triangle, 6u>::x_[6] = {
	0.0915762135098, 0.8168475729805, 0.0915762135098,
	0.1081030181681, 0.4459484909160, 0.4459484909160
};

/*! @brief Coordinates for 6 nodes formula */
template <>
double const newton_cotes<gas::geometry::unit::triangle, 6u>::y_[6] = {
	0.0915762135098, 0.0915762135098, 0.8168475729805,
	0.4459484909160, 0.1081030181681, 0.4459484909160
};

/*! @brief Coordinates for 6 nodes formula */
template <>
double const newton_cotes<gas::geometry::unit::triangle, 6u>::w_[6] = {
	0.0549758718277, 0.0549758718277, 0.0549758718277,
	0.1116907948390, 0.1116907948390, 0.1116907948390
};

/*! @brief Degree of exactness */
template <>
unsigned int const newton_cotes<gas::geometry::unit::triangle, 6u>::degree_ = 4u;

} } }

#endif // _gas_numerical_quadrature_newton_cotes_
