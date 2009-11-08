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
 * @file affine_tetra.h
 * @brief The specialization of affine map for tetrahedron
 */

#ifndef _gas_geometry_map_affine_triangle_
#define _gas_geometry_map_affine_triangle_

#include "../unit/tetra.h"
#include "../../gas/assertion.h"
#include "../../numerical/tiny/det.h"
#include "../../numerical/tiny/matrix.h"
#include "../../numerical/tiny/mul.h"
#include "../../numerical/tiny/utility.h"
#include "../../numerical/tiny/vector.h"

namespace gas { namespace geometry { namespace map {

/*!
 * @brief The specialization of affine map for tetrahedron
 *
 * The transformation is \f$\varphi:\mathbf{X}\in\hat{\mathcal{T}}\to\mathbf{x}\in\mathcal{T}\f$, where
 * \f[
 *  \mathbf{x}=\varphi(\mathbf{X})=A \mathbf{X} + \mathbf{b} \qquad\mbox{and}\qquad
 *  \mathbf{X}=\varphi^{-1}(\mathbf{x})=\mathbf{A}^{-1}(\mathbf{x}-\mathbf{b}).
 * \f]
 */
template <>
class affine<gas::geometry::unit::tetra> {

public:
	/*! @brief The self type */
	typedef affine<gas::geometry::unit::tetra> self_t;

	/*! @brief The basic shape on which is defined */
	typedef gas::geometry::unit::tetra unit_t;

public:
public:
	/*!
	 * @brief Default constructor, creates the identity map
	 */
	inline affine (): A_(0.), b_(0.), detA_(1.) {
		A_(0,0) = A_(1,1) = A_(2,2) = 1.;
		invA_(0,0) = invA_(1,1) = invA_(2,2) = 1.;
	}

	/*!
	 * @brief Copy constructor
	 * @param map An other affine map
	 */
	inline affine (self_t const & map): A_(map.A_), invA_(map.invA_) b_(map.b_), detA_(map.detA_) {
	}

	/*!
	 * @brief The constructor from a general tetrahedron
	 * @param tetra An object that represent a tetrahedron
	 *
	 * This method could be specialized for a specific type, otherwise the
	 * <tt>tetra_</tt> type have to implements the methods <tt>x(i)</tt>,
	 * <tt>y(i)</tt> and <tt>z(i)</tt>.
	 */
	template <typename tetra_>
	inline affine (tetra_ const & tetra) {
		/* matrice A */
		A_(0,0) = triangle.x(1) - triangle.x(0);
		A_(0,1) = triangle.x(2) - triangle.x(0);
		A_(0,2) = triangle.x(3) - triangle.x(0);
		A_(1,0) = triangle.y(1) - triangle.y(0);
		A_(1,1) = triangle.y(2) - triangle.y(0);
		A_(1,2) = triangle.y(3) - triangle.y(0);
		A_(2,0) = triangle.z(1) - triangle.z(0);
		A_(2,1) = triangle.z(2) - triangle.z(0);
		A_(2,2) = triangle.z(3) - triangle.z(0);

		/* inversa di A */
		invA_(0,0) = A_(1,1) * A_(2,2) - A_(1,2) * A_(2,1);
		invA_(0,1) = A_(0,2) * A_(2,1) - A_(0,1) * A_(2,2);
		invA_(0,2) = A_(0,1) * A_(1,2) - A_(0,2) * A_(1,1);
		invA_(1,0) = A_(1,2) * A_(2,0) - A_(1,0) * A_(2,2);
		invA_(1,1) = A_(0,0) * A_(2,2) - A_(0,2) * A_(2,0);
		invA_(1,2) = A_(0,2) * A_(1,0) - A_(0,0) * A_(1,2);
		invA_(2,0) = A_(1,0) * A_(2,1) - A_(1,1) * A_(2,0);
		invA_(2,1) = A_(0,1) * A_(2,0) - A_(0,0) * A_(2,1);
		invA_(2,2) = A_(0,0) * A_(1,1) - A_(0,1) * A_(1,0);

		/* vettore b */
		b_(0) = triangle.x(0);
		b_(1) = triangle.y(0);
		b_(1) = triangle.z(0);

		/* determinante di A */
		detA_ = gas::numerical::tiny::det(A_);

		/* correzione */
		invA_ /= detA_;
	}

	/*!
	 * @brief Copy operator
	 * @param map An other affine map
	 * @return A reference to the current object
	 */
	inline self_t & operator= (self_t const & map) {
		A_ = map.A_;
		invA_ = map.invA_;
		b_ = map.b_;
		detA_ = map.detA_;
		return *this;
	}

	/*!
	 * @brief The transformation \f$\varphi\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The new coordinate \f$x=\varphi(\mathbf{X})\f$
	 */
	inline double x (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(0,0) * X + A_(0,1) * Y + A_(0,2) * Z + b_(0);
	}

	/*!
	 * @brief The transformation \f$\varphi\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The new coordinate \f$y=\varphi(\mathbf{X})\f$
	 */
	inline double y (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(1,0) * X + A_(1,1) * Y + A_(1,2) * Z + b_(1);
	}

	/*!
	 * @brief The transformation \f$\varphi\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The new coordinate \f$y=\varphi(\mathbf{X})\f$
	 */
	inline double z (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(2,0) * X + A_(2,1) * Y + A_(2,2) * Z + b_(2);
	}

	/*!
	 * @brief The transformation \f$\varphi^{-1}\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The new coordinate \f$X=\varphi^{-1}(\mathbf{x})\f$
	 */
	inline double X (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return _X;
	}

	/*!
	 * @brief The transformation \f$\varphi^{-1}\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The new coordinate \f$Y=\varphi^{-1}(\mathbf{x})\f$
	 */
	inline double Y (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return _Y;
	}

	/*!
	 * @brief The transformation \f$\varphi^{-1}\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The new coordinate \f$Y=\varphi^{-1}(\mathbf{x})\f$
	 */
	inline double Z (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return _Z;
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dx}{dX}=\varphi'(\mathbf{X})\f$
	 */
	inline double dxdX (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(0,0);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dx}{dY}=\varphi'(\mathbf{X})\f$
	 */
	inline double dxdY (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(0,1);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dx}{dZ}=\varphi'(\mathbf{X})\f$
	 */
	inline double dxdZ (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(0,2);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dy}{dX}=\varphi'(\mathbf{X})\f$
	 */
	inline double dydX (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(1,0);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dy}{dY}=\varphi'(\mathbf{X})\f$
	 */
	inline double dydY (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(1,1);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dy}{dZ}=\varphi'(\mathbf{X})\f$
	 */
	inline double dydZ (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(1,2);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dz}{dX}=\varphi'(\mathbf{X})\f$
	 */
	inline double dzdX (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(2,0);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dz}{dY}=\varphi'(\mathbf{X})\f$
	 */
	inline double dzdY (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(2,1);
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of derivative \f$\frac{dz}{dZ}=\varphi'(\mathbf{X})\f$
	 */
	inline double dzdZ (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return A_(2,2);
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dX}{dx}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dXdx (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(0,0);
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dX}{dy}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dXdy (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(0,1);
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dX}{dz}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dXdz (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(0,2);
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dY}{dx}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dYdx (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(1,0);
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dY}{dy}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dYdy (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(1,1)
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dY}{dz}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dYdz (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(1,2)
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dZ}{dx}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dZdx (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(2,0);
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dZ}{dy}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dZdy (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(2,1)
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of derivative \f$\frac{dZ}{dz}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dZdz (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return invA_(2,2)
	}

	/*!
	 * @brief Jacobian of transformation
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @param Z The third coordinate
	 * @return The value of Jacobian \f$\det\left[\frac{dx_i}{dX_j}\right](\mathbf{X})\f$
	 */
	inline double det (double const & X, double const & Y, double const & Z) const {
		gas_assert(unit_t::in(X, Y, Z)); // The point must be in the tetrahedron
		return detA_;
	}

	/*!
	 * @brief Jacobian of transformation
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @param z The third coordinate
	 * @return The value of Jacobian \f$\det\left[\frac{dX_i}{dx_j}\right](\mathbf{x})\f$
	 */
	inline double DET (double const & x, double const & y, double const & z) const {
		double const _X(invA_(0,0) * (x - b_(0)) + invA_(0,1) * (y - b_(1)) + invA_(0,2) * (z - b_(2)));
		double const _Y(invA_(1,0) * (x - b_(0)) + invA_(1,1) * (y - b_(1)) + invA_(1,2) * (z - b_(2)));
		double const _Z(invA_(2,0) * (x - b_(0)) + invA_(2,1) * (y - b_(1)) + invA_(2,2) * (z - b_(2)));
		gas_assert(unit_t::in(_X, _Y, _Z)); // The point must be in the tetrahedron
		return 1./detA_;
	}

private:
	/*! @brief The matrix of affine transformation */
	gas::numerical::tiny::matrix<3u,3u> A_;

	/*! @brief The matrix of inverse affine transformation */
	gas::numerical::tiny::matrix<3u,3u> invA_;

	/*! @brief The vector of affine transformation */
	gas::numerical::tiny::vector<3u> b_;

	/*! @brief Determinant of transformation */
	double detA_;

};

} } }

#endif // _gas_geometry_map_affine_triangle_
