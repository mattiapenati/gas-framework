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
 * @file formula.h
 * @brief The container for quadrature formula
 */

#ifndef _gas_numerical_quadrature_formula_
#define _gas_numerical_quadrature_formula_

#include "../../gas/assertion.h"
#include "../../gas/type.h"
#include "../../geometry/map/affine.h"

namespace gas { namespace numerical { namespace quadrature {

/*!
 * @brief A quadrature formulae
 * @param method_ The quadrature method to use
 * @param map_ The map to transform the coordinates of points
 */
template <typename method_, typename map_ = gas::geometry::map::affine<typename method_::unit_t> >
class formula {

public:
	/*! @brief The self type */
	typedef formula<method_, map_> self_t;

private:
	/*! @brief The quadrature method to use */
	typedef method_ method_t;

	/*! @brief The map to transform the point */
	typedef map_ map_t;

public:
	/*!
	 * @brief The constructor
	 */
	inline formula(): m_() {
	}

	/*!
	 * @brief Set the domain of integration
	 * @param domain The domain of integration (a suitable object to create the
	 *        map)
	 * @return A reference to the current object
	 */
	template <typename domain_>
	inline self_t & operator() (domain_ const & domain) {
		typedef typename method_t::unit_t method_unit_t; // geometry of method
		typedef typename map_t::unit_t map_unit_t;       // geometry of map
		gas_static_assert(
			(gas::same_type<method_unit_t, map_unit_t>::value),
			Method_and_map_must_act_on_the_same_unit_geometry
		); // This static assertion to check that the method and the map act on
		   // the same geometry
		m_.map(map_(domain));
		return *this;
	}

	/*!
	 * @brief Apply the method to a function
	 * @param f Any kind of object with the possibility to call f(...)
	 */
	template <typename function_>
	inline double integrate (function_ const & f) {
		return m_.apply(f);
	}

private:
	method_ m_;

};

} } }

#endif // _gas_numerical_quadrature_formula_
