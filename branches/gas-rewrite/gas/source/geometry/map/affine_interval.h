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
 * @file affine_interval.h
 * @brief The specialization of affine map for interval
 */

#ifndef _gas_geometry_map_affine_interval_
#define _gas_geometry_map_affine_interval_

#include "../unit/interval.h"
#include "../../gas/assertion.h"

namespace gas { namespace geometry { namespace map {

/*!
 * @brief The specialization of affine map for interval
 *
 * The transformation is \f$\varphi:X\in(-1,1)\to x\in(a,b)\f$, where
 * \f[
 *  x=\varphi(X)=\alpha X + \beta \qquad\mbox{and}\qquad
 *  X=\varphi^{-1}(x)=\frac{x-\beta}{\alpha}.
 * \f]
 */
template <>
class affine<gas::geometry::unit::interval> {

public:
	/*! @brief The self type */
	typedef affine<gas::geometry::unit::interval> self_t;

	/*! @brief The basic shape on which is defined */
	typedef gas::geometry::unit::interval unit_t;

public:
	/*!
	 * @brief Default constructor, creates the identity map
	 */
	inline affine (): alfa_(1.), beta_(0.) {
	}

	/*!
	 * @brief Copy constructor
	 * @param map An other affine map
	 */
	inline affine (self_t const & map): alfa_(map.alfa_), beta_(map.beta_) {
	}

	/*!
	 * @brief The constructor from a general interval
	 * @param interval An object that represent an interval
	 * @pre interval.a() < interval.b()
	 *
	 * This method could be specialized for a specific type, otherwise the
	 * <tt>interval_</tt> type have to implements the methods <tt>a()</tt> and
	 * <tt>b()</tt>.
	 */
	template <typename interval_>
	inline affine (interval_ const & interval) {
		double const a(interval.a());
		double const b(interval.b());
		gas_pre(a < b);
		alfa_ = (b - a) / 2.;
		beta_ = (b + a) / 2.;
	}

	/*!
	 * @brief Copy operator
	 * @param map An other affine map
	 * @return A reference to the current object
	 */
	inline self_t & operator= (self_t const & map) {
		alfa_ = map.alfa_;
		beta_ = map.beta_;
		return *this;
	}

	/*!
	 * @brief The transformation \f$\varphi\f$
	 * @param X The coordinate
	 * @return The new coordinate \f$x=\varphi(X)\f$
	 */
	inline double x (double const & X) const {
		gas_assert(unit_t::in(X)); // The point must be in the interval
		return alfa_ * X + beta_;
	}

	/*!
	 * @brief The transformation \f$\varphi^{-1}\f$
	 * @param x The coordinate
	 * @return The new coordinate \f$X=\varphi^{-1}(x)\f$
	 */
	inline double X (double const & x) const {
		double const _X((x-beta_)/alfa_);
		gas_assert(unit_t::in(_X)); // The point must be in the interval
		return _X;
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The coordinate
	 * @return The value of derivative \f$\frac{dx}{dX}=\varphi'(X)\f$
	 */
	inline double dxdX (double const & X) const {
		gas_assert(unit_t::in(X)); // The point must be in the interval
		return alfa_;
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The coordinate
	 * @return The value of derivative \f$\frac{dX}{dx}=(\varphi^{-1})'(x)\f$
	 */
	inline double dXdx (double const & x) const {
		double const _X((x-beta_)/alfa_);
		gas_assert(unit_t::in(_X)); // The point must be in the interval
		return 1./alfa_;
	}

	/*!
	 * @brief Jacobian of transformation
	 * @param X The coordinate
	 * @return The value of Jacobian \f$\det\left[\frac{dx_i}{dX_j}\right](X)\f$
	 */
	inline double det (double const & X) const {
		gas_assert(unit_t::in(X)); // The point must be in the interval
		return alfa_;
	}

	/*!
	 * @brief Jacobian of transformation
	 * @param x The coordinate
	 * @return The value of Jacobian \f$\det\left[\frac{dX_i}{dx_j}\right](x)\f$
	 */
	inline double DET (double const & x) const {
		double const _X((x-beta_)/alfa_);
		gas_assert(unit_t::in(_X)); // The point must be in the interval
		return 1./alfa_;
	}

private:
	/*! @brief The first parameter of map */
	double alfa_;

	/*! @brief The second parameter of map */
	double beta_;

};

} } }

#endif // _gas_geometry_map_affine_interval_
