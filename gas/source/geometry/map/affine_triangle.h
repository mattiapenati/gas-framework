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
 * @file affine_triangle.h
 * @brief The specialization of affine map for triangle
 */

#ifndef _gas_geometry_map_affine_triangle_
#define _gas_geometry_map_affine_trinagle_

#include "../unit/triangle.h"
#include "../../gas/assertion.h"
#include "../../numerical/tiny/det.h"
#include "../../numerical/tiny/matrix.h"
#include "../../numerical/tiny/mul.h"
#include "../../numerical/tiny/utility.h"
#include "../../numerical/tiny/vector.h"

namespace gas { namespace geometry { namespace map {

/*!
 * @brief The specialization of affine map for triangle
 *
 * The transformation is \f$\varphi:\mathbf{X}\in\hat{\mathcal{T}}\to\mathbf{x}\in\mathcal{T}\f$, where
 * \f[
 *  \mathbf{x}=\varphi(\mathbf{X})=A \mathbf{X} + \mathbf{b} \qquad\mbox{and}\qquad
 *  \mathbf{X}=\varphi^{-1}(\mathbf{x})=\mathbf{A}^{-1}(\mathbf{x}-\mathbf{b}).
 * \f]
 */
template <>
class affine<gas::geometry::unit::triangle> {

public:
	/*! @brief The self type */
	typedef affine<gas::geometry::unit::triangle> self_t;

	/*! @brief The basic shape on which is defined */
	typedef gas::geometry::unit::triangle unit_t;

public:
public:
	/*!
	 * @brief Default constructor, creates the identity map
	 */
	inline affine (): A_(0.), b_(0.), detA_(1.) {
		A_(0,0) = A_(1,1) = 1.;
	}

	/*!
	 * @brief Copy constructor
	 * @param map An other affine map
	 */
	inline affine (self_t const & map): A_(map.A_), b_(map.b_), detA_(map.detA_) {
	}

	/*!
	 * @brief The constructor from a general triangle
	 * @param triangle An object that represent an interval
	 *
	 * This method could be specialized for a specific type, otherwise the
	 * <tt>interval_</tt> type have to implements the methods <tt>x(i)</tt> and
	 * <tt>y(i)</tt>.
	 */
	template <typename triangle_>
	inline affine (triangle_ const & triangle) {
		/* matrice A */
		A_(0,0) = triangle.x(1) - triangle.x(0);
		A_(0,1) = triangle.x(2) - triangle.x(0);
		A_(1,0) = triangle.y(1) - triangle.y(0);
		A_(1,1) = triangle.y(2) - triangle.y(0);

		/* vettore b */
		b_(0) = triangle.x(0);
		b_(1) = triangle.y(0);

		/* determinante di A */
		detA_ = gas::numerical::tiny::det(A_);
	}

	/*!
	 * @brief Copy operator
	 * @param map An other affine map
	 * @return A reference to the current object
	 */
	inline self_t & operator= (self_t const & map) {
		A_ = map.A_;
		b_ = map.b_;
		detA_ = map.detA_;
		return *this;
	}

	/*!
	 * @brief The transformation \f$\varphi\f$
	 * @param X First coordinate
	 * @param Y Second coordinate
	 * @return The new coordinate \f$x=\varphi(\mathbf{X})\f$
	 */
	inline double x (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return A_(0,0) * X + A_(0,1) * Y + b_(0);
	}

	/*!
	 * @brief The transformation \f$\varphi\f$
	 * @param X First coordinate
	 * @param Y Second coordinate
	 * @return The new coordinate \f$y=\varphi(\mathbf{X})\f$
	 */
	inline double y (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return A_(1,0) * X + A_(1,1) * Y + b_(1);
	}

	/*!
	 * @brief The transformation \f$\varphi^{-1}\f$
	 * @param x First coordinate
	 * @param y Second coordinate
	 * @return The new coordinate \f$X=\varphi^{-1}(\mathbf{x})\f$
	 */
	inline double X (double const & x, double const & y) const {
		/* TODO Sostituire con operazioni con matrici */
		double const X((A_(1,1) * (x - b_(0)) - A_(0,1) * (y - b_(1))) / detA_);
		double const Y((A_(0,0) * (y - b_(1)) - A_(1,0) * (x - b_(0))) / detA_);
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return X;
	}

	/*!
	 * @brief The transformation \f$\varphi^{-1}\f$
	 * @param x First coordinate
	 * @param y Second coordinate
	 * @return The new coordinate \f$Y=\varphi^{-1}(\mathbf{x})\f$
	 */
	inline double Y (double const & x, double const & y) const {
		/* TODO Sostituire con operazioni con matrici */
		double const X((A_(1,1) * (x - b_(0)) - A_(0,1) * (y - b_(1))) / detA_);
		double const Y((A_(0,0) * (y - b_(1)) - A_(1,0) * (x - b_(0))) / detA_);
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return Y;
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X First coordinate
	 * @param Y Second coordinate
	 * @return The value of derivative \f$\frac{dx}{dX}=\varphi'(\mathbf{X})\f$
	 */
	inline double dxdX (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return A_(0,0);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X First coordinate
	 * @param Y Second coordinate
	 * @return The value of derivative \f$\frac{dx}{dY}=\varphi'(\mathbf{X})\f$
	 */
	inline double dxdY (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return A_(0,1);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X First coordinate
	 * @param Y Second coordinate
	 * @return The value of derivative \f$\frac{dy}{dX}=\varphi'(\mathbf{X})\f$
	 */
	inline double dydX (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return A_(1,0);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X First coordinate
	 * @param Y Second coordinate
	 * @return The value of derivative \f$\frac{dy}{dY}=\varphi'(\mathbf{X})\f$
	 */
	inline double dydY (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return A_(1,1);
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x First coordinate
	 * @param y Second coordinate
	 * @return The value of derivative \f$\frac{dX}{dx}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dXdx (double const & x, double const & y) const {
		double const X((A_(1,1) * (x - b_(0)) - A_(0,1) * (y - b_(1))) / detA_);
		double const Y((A_(0,0) * (y - b_(1)) - A_(1,0) * (x - b_(0))) / detA_);
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return  A_(1,1)/detA_;
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x First coordinate
	 * @param y Second coordinate
	 * @return The value of derivative \f$\frac{dX}{dy}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dXdy (double const & x, double const & y) const {
		double const X((A_(1,1) * (x - b_(0)) - A_(0,1) * (y - b_(1))) / detA_);
		double const Y((A_(0,0) * (y - b_(1)) - A_(1,0) * (x - b_(0))) / detA_);
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return -A_(0,1)/detA_;
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x First coordinate
	 * @param y Second coordinate
	 * @return The value of derivative \f$\frac{dY}{dx}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dYdx (double const & x, double const & y) const {
		double const X((A_(1,1) * (x - b_(0)) - A_(0,1) * (y - b_(1))) / detA_);
		double const Y((A_(0,0) * (y - b_(1)) - A_(1,0) * (x - b_(0))) / detA_);
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return -A_(1,0)/detA_;
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x First coordinate
	 * @param y Second coordinate
	 * @return The value of derivative \f$\frac{dY}{dy}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dYdy (double const & x, double const & y) const {
		double const X((A_(1,1) * (x - b_(0)) - A_(0,1) * (y - b_(1))) / detA_);
		double const Y((A_(0,0) * (y - b_(1)) - A_(1,0) * (x - b_(0))) / detA_);
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return  A_(0,0)/detA_;
	}

	/*!
	 * @brief Jacobian of transformation
	 * @param X First coordinate
	 * @param Y Second coordinate
	 * @return The value of Jacobian \f$\det\left[\frac{dx_i}{dX_j}\right](\mathbf{X})\f$
	 */
	inline double det (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return detA_;
	}

	/*!
	 * @brief Jacobian of transformation
	 * @param x First coordinate
	 * @param y Second coordinate
	 * @return The value of Jacobian \f$\det\left[\frac{dX_i}{dx_j}\right](\mathbf{x})\f$
	 */
	inline double DET (double const & x, double const & y) const {
		double const X((A_(1,1) * (x - b_(0)) - A_(0,1) * (y - b_(1))) / detA_);
		double const Y((A_(0,0) * (y - b_(1)) - A_(1,0) * (x - b_(0))) / detA_);
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return 1./detA_;
	}

private:
	/*! @brief The matrix of affine transformation */
	gas::numerical::tiny::matrix<2u,2u> A_;
	/*! @brief The vector of affine transformation */
	gas::numerical::tiny::vector<2u> b_;
	/*! @brief Determinant of transformation */
	double detA_;

};

} } }

#endif // _gas_geometry_map_affine_triangle_