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
 * @file method.h
 * @brief Some generical implementation of quadrature formulae
 */

#ifndef _gas_numerical_quadrature_method_method_
#define _gas_numerical_quadrature_method_method_

#include <cmath>
#include "../../gas/macro.h"

namespace gas { namespace numerical { namespace quadrature {

/*!
 * @brief A general implementation, with a list of all nodes and weights
 * @param data_ The struct that contains the nodes and weights
 * @param nodes_ The number of nodes
 */
template <typename unit_, unsigned int nodes_, typename data_>
class method_1 {

public:
	/*! @brief The self type */
	typedef method_1<unit_, nodes_, data_> self_t;

	/*! @brief The unit geometry */
	typedef unit_ unit_t;

private:
	/*! @brief The dimension of unit_t */
	static unsigned int const d = unit_t::d;

	/*! @brief The number of nodes */
	static unsigned int const n = nodes_;

	/*! @brief The class containing all data */
	typedef data_ data_t;

public:
	/*!
	 * @brief Initialize the method, copying the value of nodes and weights
	 */
	method_1 () {
		gas_static_assert((d == 1), This_is_a_1d_quadrature_formula);
		rangeu(i, n) x_[i] = data_t::x_[i];
		rangeu(i, n) w_[i] = data_t::w_[i];
	}

private:
	/*!
	 * @brief Map the nodes and weight
	 * @param m The map used to recalculate nodes and weights
	 */
	template <typename map_>
	void map (map_ const & m) {
		rangeu(i, n)
			self_t::x_[i] = m.x(data_t::x[i]);
		rangeu(i, n)
			self_t::w_[i] = std::abs(m.DET(self_t::x_[i])) * data_t::w_[i];
	}

	/*!
	 * @brief Apply the quadrature formula with current nodes and weights
	 * @param f Any kind of object with the possibility to call f(...)
	 * @return The computed numerical integral
	 */
	template <typename function_>
	double apply (function_ const & f) {
		double r(0.);
		rangeu(i, n)
			r += f(self_t::x_[i]) * self_t::w_[i];
		return r;
	}

private:
	/*! @brief The list of nodes of quadrature formula */
	double x_[n];

	/*! @brief The list of weights of quadrature formula */
	double w_[n];

	template <typename method_, typename map_>
	friend class formula;

};

} } }

#endif // _gas_numerical_quadrature_method_method_
