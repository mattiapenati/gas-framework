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

#ifndef _gas_functional_element_
#define _gas_functional_element_

#include "../gas/assertion.h"
#include "../gas/type.h"

namespace gas { namespace functional {

/*!
 * @brief A base function
 * @param d_ The dimension of geometry
 * @param element_ The element to use
 */
template <unsigned int d_, typename element_>
class base_function;

/*!
 * @brief A local element, link a base with a map
 * @param base_ The base to use
 * @param map_ The map to transform the coordinates
 */
template <typename base_, typename map_>
class element {

public:
	/*! @brief The self type */
	typedef element<base_, map_> self_t;

	/*! @brief The base */
	typedef base_ base_t;

	/*! @brief The map */
	typedef map_ map_t;

	/*! @brief The dimensione of domain */
	static unsigned int const d = map_t::unit_t::d;

public:
	/*!
	 * @brief The constructor
	 * @param domain A suitable object from which construct a new map
	 */
	template <typename domain_>
	element (domain_ const & domain): m_(domain) {
		gas_static_assert(
			(gas::same_type<typename base_t::unit_t, typename map_t::unit_t>::value),
			Base_and_map_must_be_defined_on_the_same_geometry
		);
	}

private:
	/*! @brief The map */
	map_t const m_;

	template <unsigned int d__, typename element__>
	friend class base_function;

};

/*!
 * @brief A base function (1d domain)
 */
template <typename element_>
class base_function<1u, element_> {

public:
	/*! @brief The self type */
	typedef base_function<1u, element_> self_t;

private:
	/*! @brief The base */
	typedef element_ element_t;

	/*! @brief The base */
	typedef typename element_t::base_t base_t;

private:
	/*!
	 * @brief The constructor
	 * @param i The index of base
	 */
	inline base_function (unsigned int const & i, element_ const & element): i_(i), el_(element) {
	}

public:
	inline double operator() (double const & x) {
		double const _X(el_.m_.X(x));
		return base_t::b(i_, _X);
	}

private:
	unsigned int const i_;
	element_ const & el_;

	template <typename base__, typename map__>
	friend class element;

};

} }

#endif // _gas_functional_element_
