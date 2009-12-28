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
 * @file element.h
 * @brief The element definition
 */

#ifndef _gas_functional_element_
#define _gas_functional_element_

#include "../gas"
#include "../geometry/map/map"
#include "derivative.h"
#include "function.h"
#include "base/base"

namespace gas { namespace functional {

template <unsigned int d_, typename element_> class base_function;
template <unsigned int d_, typename element_> class d_base_function;
template <unsigned int d_, typename element_> class dx_base_function;
template <unsigned int d_, typename element_> class dy_base_function;
template <unsigned int d_, typename element_> class dz_base_function;

/*!
 * @brief A local element, link a base with a map
 * @param base_ The base to use
 * @param map_ The map to transform the coordinates
 */
template <typename base_, typename map_ = gas::geometry::map::affine<typename gas::functional::base::info<base_>::unit_t> >
class element {

public:
	/*! @brief The self type */
	typedef element<base_, map_> self_t;

	/*! @brief The base */
	typedef base_ base_t;

	/*! @brief The map */
	typedef map_ map_t;

	/*! @brief The dimensione of domain */
	static unsigned int const d = gas::geometry::map::info<map_t>::d;

public:
	inline element (): m_() {}
	/*!
	 * @brief The constructor
	 * @param domain A suitable object from which construct a new map
	 */
	template <typename domain_>
	inline element (domain_ const & domain): m_(domain) {
		typedef typename gas::functional::base::info<base_t>::unit_t base_unit_t;
		typedef typename gas::geometry::map::info<map_t>::unit_t map_unit_t;
		gas_static_assert(
			(gas::same_type<base_unit_t, map_unit_t>::value),
			Base_and_map_must_be_defined_on_the_same_geometry
		);
	}

	/*!
	 * @brief The base function associated to the index i
	 * @param i The index of base function
	 * @return The base function
	 */
	inline base_function<d, self_t> b (unsigned int const & i) const {
		return base_function<d, self_t>(i, *this);
	}

private:
	/*! @brief The map */
	map_t const m_;

	template <unsigned int d__, typename element__> friend class base_function;
	template <unsigned int d__, typename element__> friend class d_base_function;
	template <unsigned int d__, typename element__> friend class dx_base_function;
	template <unsigned int d__, typename element__> friend class dy_base_function;
	template <unsigned int d__, typename element__> friend class dz_base_function;

};

/*!
 * @brief A base function (1d domain)
 */
template <typename element_>
class base_function<1u, element_>: public function<1u, base_function<1u, element_> > {

/* TODO derivate */

private:
	/*! @brief The base */
	typedef typename element_::base_t base_t;

private:
	/*!
	 * @brief The constructor
	 * @param i The index of base
	 */
	inline base_function (unsigned int const & i, element_ const & element): i_(i), el_(element) {
	}

public:
	inline double operator() (double const & x) const {
		double const _X(el_.m_.X(x));
		return base_t::b(i_, _X);
	}

private:
	unsigned int const i_;
	element_ const & el_;

	template <typename base__, typename map__> friend class element;
	friend class d_base_function<2u, element_>;

};

/*!
 * @brief A base function (2d domain)
 */
template <typename element_>
class base_function<2u, element_>: public function<2u, base_function<2u, element_> > {

public:
	typedef dx_base_function<2u, element_> dx_t;
	typedef dy_base_function<2u, element_> dy_t;

private:
	/*! @brief The base */
	typedef typename element_::base_t base_t;

private:
	/*!
	 * @brief The constructor
	 * @param i The index of base
	 */
	inline base_function (unsigned int const & i, element_ const & element): i_(i), el_(element) {
	}

public:
	inline double operator() (double const & x, double const & y) const {
		double const _X(el_.m_.X(x, y));
		double const _Y(el_.m_.Y(x, y));
		return base_t::b(i_, _X, _Y);
	}

private:
	unsigned int const i_;
	element_ const & el_;

	template <typename base__, typename map__> friend class element;
	friend class d_base_function<2u, element_>;
	friend class dx_base_function<2u, element_>;
	friend class dy_base_function<2u, element_>;

};

/*!
 * @brief A base function (3d domain)
 */
template <typename element_>
class base_function<3u, element_>: public function<3u, base_function<3u, element_> > {

/* TODO derivate */

private:
	/*! @brief The base */
	typedef typename element_::base_t base_t;

private:
	/*!
	 * @brief The constructor
	 * @param i The index of base
	 */
	inline base_function (unsigned int const & i, element_ const & element): i_(i), el_(element) {
	}

public:
	inline double operator() (double const & x, double const & y, double const & z) const {
		double const _X(el_.m_.X(x, y, z));
		double const _Y(el_.m_.Y(x, y, z));
		double const _Z(el_.m_.Z(x, y, z));
		return base_t::b(i_, _X, _Y, _Z);
	}

private:
	unsigned int const i_;
	element_ const & el_;

	template <typename base__, typename map__>
	friend class element;

};

/*!
 * @brief The derivative of base function (2d domain)
 */
template <typename element_>
class dx_base_function<2u, element_>: public function<2u, dx_base_function<2u, element_> > {

private:
	/*! @brief The base */
	typedef typename element_::base_t base_t;

private:
	dx_base_function (base_function<2u, element_> const & base): i_(base.i_), el_(base.el_) {
	}

public:
	inline double operator() (double const & x, double const & y) const {
		double const _X(el_.m_.X(x, y));
		double const _Y(el_.m_.Y(x, y));
		return base_t::dbdX(i_, _X, _Y) * el_.m_.dXdx(x, y) + base_t::dbdY(i_, _X ,_Y) * el_.m_.dYdx(x, y);
	}

	unsigned int const & i_;
	element_ const & el_;

	template <typename function__>
	friend typename function__::dx_t dx (function__ const & f);

};

/*!
 * @brief The derivative of base function (2d domain)
 */
template <typename element_>
class dy_base_function<2u, element_>: public function<2u, dy_base_function<2u, element_> > {

private:
	/*! @brief The base */
	typedef typename element_::base_t base_t;

private:
	dy_base_function (base_function<2u, element_> const & base): i_(base.i_), el_(base.el_) {
	}

public:
	inline double operator() (double const & x, double const & y) const {
		double const _X(el_.m_.X(x, y));
		double const _Y(el_.m_.Y(x, y));
		return base_t::dbdX(i_, _X, _Y) * el_.m_.dXdy(x, y) + base_t::dbdY(i_, _X ,_Y) * el_.m_.dYdy(x, y);
	}

	unsigned int const & i_;
	element_ const & el_;

	template <typename function__>
	friend typename function__::dy_t dy (function__ const & f);

};

} }

#endif // _gas_functional_element_
