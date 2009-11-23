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
 * @file affine_triangle.h
 * @brief The specialization of affine map for triangle
 */

#ifndef _gas_geometry_map_affine_triangle_
#define _gas_geometry_map_affine_triangle_

#include "../../../unit/unit"
#include "../../../../gas"

namespace gas { namespace geometry { namespace map {

template <>
class affine<gas::geometry::unit::triangle> {

public:
	/* The basic shape on which is defined */
	typedef gas::geometry::unit::triangle unit_t;

private:
	/* The self type */
	typedef affine<gas::geometry::unit::triangle> self_t;

public:
	/* Default constructor, creates the identity map */
	inline affine (): detA_(1.) {
		Ax_[0] = Ay_[1] = invAx_[0] = invAy_[1] = 1.;
		Ax_[1] = Ax_[0] = invAx_[1] = invAx_[0] = 0.;
		b_[0] = b_[1] = 0.;
	}

	/*!
	 * @brief Copy constructor
	 * @param map An other affine map
	 */
	inline affine (self_t const & map): detA_(map.detA_) {
		Ax_[0] = map.Ax_[0];
		Ax_[1] = map.Ax_[1];

		Ay_[0] = map.Ay_[0];
		Ay_[1] = map.Ay_[1];

		invAx_[0] = map.invAx_[0];
		invAx_[1] = map.invAx_[1];

		invAy_[0] = map.invAy_[0];
		invAy_[1] = map.invAy_[1];

		b_[0] = map.b_[0];
		b_[1] = map.b_[1];
	}

	/*!
	 * @brief The constructor from a general triangle
	 * @param triangle An object that represent a triangle
	 *
	 * This method could be specialized for a specific type, otherwise the
	 * <tt>triangle_</tt> type have to implements the methods <tt>x(i)</tt> and
	 * <tt>y(i)</tt>.
	 */
	template <typename triangle_>
	inline affine (triangle_ const & triangle) {
		/* trasformazione diretta */
		Ax_[0] = triangle.x(1) - triangle.x(0);
		Ax_[1] = triangle.x(2) - triangle.x(0);

		Ay_[0] = triangle.y(1) - triangle.y(0);
		Ay_[1] = triangle.y(2) - triangle.y(0);

		/* determinante di A */
		detA_ = Ax_[0] * Ay_[1] - Ax_[1] * Ay_[0];

		/* trasformazione inversa */
		invAx_[0] = Ay_[1] / detA_;
		invAx_[1] = - Ax_[1] / detA_;

		invAy_[0] = - Ay_[0] / detA_;
		invAy_[1] = Ax_[0] / detA_;

		/* vettore b */
		b_[0] = triangle.x(0);
		b_[1] = triangle.y(0);
	}

	/*!
	 * @brief Copy operator
	 * @param map An other affine map
	 * @return A reference to the current object
	 */
	inline self_t & operator= (self_t const & map) {
		Ax_[0] = map.Ax_[0];
		Ax_[1] = map.Ax_[1];

		Ay_[0] = map.Ay_[0];
		Ay_[1] = map.Ay_[1];

		invAx_[0] = map.invAx_[0];
		invAx_[1] = map.invAx_[1];

		invAy_[0] = map.invAy_[0];
		invAy_[1] = map.invAy_[1];

		b_[0] = map.b_[0];
		b_[1] = map.b_[1];

		detA_ = map.detA_;
		return *this;
	}

	/*!
	 * @brief The transformation \f$\varphi\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The new coordinate \f$x=\varphi(\mathbf{X})\f$
	 */
	inline double x (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		__m128d _X = _mm_set_pd(X, Y); // loading coordinate
		__m128d _a = _mm_load_pd(Ax_); // loading coefficient

		__m128d _x = _mm_mul_pd(_a, _X); // multiplication

		double _v[2] __attribute__((aligned(16)));
		_mm_store_pd(_v, _x);

		return _v[0] + _v[1] + b_[0];
	}

	/*!
	 * @brief The transformation \f$\varphi\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The new coordinate \f$y=\varphi(\mathbf{X})\f$
	 */
	inline double y (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		__m128d _X = _mm_set_pd(X, Y); // loading coordinate
		__m128d _a = _mm_load_pd(Ay_); // loading coefficient

		__m128d _x = _mm_mul_pd(_a, _X); // multiplication

		double _v[2] __attribute__((aligned(16)));
		_mm_store_pd(_v, _x);

		return _v[0] + _v[1] + b_[1];
	}

	/*!
	 * @brief The transformation \f$\varphi^{-1}\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @return The new coordinate \f$X=\varphi^{-1}(\mathbf{x})\f$
	 */
	inline double X (double const & x, double const & y) const {
		__m128d _x = _mm_set_pd(y, x); // loading coordinate
		__m128d _b = _mm_load_pd(b_);  // loading vector

		__m128d _ax = _mm_load_pd(invAx_); // loading coefficient
		__m128d _ay = _mm_load_pd(invAy_); // loading coefficient

		_x = _mm_sub_pd(_x, _b);

		__m128d _X = _mm_mul_pd(_ax, _x);
		__m128d _Y = _mm_mul_pd(_ay, _x);

		double _vx[2] __attribute__((aligned(16)));
		double _vy[2] __attribute__((aligned(16)));
		_mm_store_pd(_vx, _X);
		_mm_store_pd(_vy, _Y);

		double const X(_vx[0] + _vx[1]);
		double const Y(_vy[0] + _vy[1]);

		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		return X;
	}

	/*!
	 * @brief The transformation \f$\varphi^{-1}\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @return The new coordinate \f$Y=\varphi^{-1}(\mathbf{x})\f$
	 */
	inline double Y (double const & x, double const & y) const {
		__m128d _x = _mm_set_pd(y, x); // loading coordinate
		__m128d _b = _mm_load_pd(b_);  // loading vector

		__m128d _ax = _mm_load_pd(invAx_); // loading coefficient
		__m128d _ay = _mm_load_pd(invAy_); // loading coefficient

		_x = _mm_sub_pd(_x, _b);

		__m128d _X = _mm_mul_pd(_ax, _x);
		__m128d _Y = _mm_mul_pd(_ay, _x);

		double _vx[2] __attribute__((aligned(16)));
		double _vy[2] __attribute__((aligned(16)));
		_mm_store_pd(_vx, _X);
		_mm_store_pd(_vy, _Y);

		double const X(_vx[0] + _vx[1]);
		double const Y(_vy[0] + _vy[1]);

		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		return Y;
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The value of derivative \f$\frac{dx}{dX}=\varphi'(\mathbf{X})\f$
	 */
	inline double dxdX (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return Ax_[0];
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The value of derivative \f$\frac{dx}{dY}=\varphi'(\mathbf{X})\f$
	 */
	inline double dxdY (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return Ax_[1];
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The value of derivative \f$\frac{dy}{dX}=\varphi'(\mathbf{X})\f$
	 */
	inline double dydX (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return Ay_[0];
	}

	/*!
	 * @brief The derivative of transformation \f$\varphi'\f$
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The value of derivative \f$\frac{dy}{dY}=\varphi'(\mathbf{X})\f$
	 */
	inline double dydY (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return Ay_[1];
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @return The value of derivative \f$\frac{dX}{dx}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dXdx (double const & x, double const & y) const {
		__m128d _x = _mm_set_pd(y, x); // loading coordinate
		__m128d _b = _mm_load_pd(b_);  // loading vector

		__m128d _ax = _mm_load_pd(invAx_); // loading coefficient
		__m128d _ay = _mm_load_pd(invAy_); // loading coefficient

		_x = _mm_sub_pd(_x, _b);

		__m128d _X = _mm_mul_pd(_ax, _x);
		__m128d _Y = _mm_mul_pd(_ay, _x);

		double _vx[2] __attribute__((aligned(16)));
		double _vy[2] __attribute__((aligned(16)));
		_mm_store_pd(_vx, _X);
		_mm_store_pd(_vy, _Y);

		double const X(_vx[0] + _vx[1]);
		double const Y(_vy[0] + _vy[1]);

		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		return invAx_[0];
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @return The value of derivative \f$\frac{dX}{dy}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dXdy (double const & x, double const & y) const {
		__m128d _x = _mm_set_pd(y, x); // loading coordinate
		__m128d _b = _mm_load_pd(b_);  // loading vector

		__m128d _ax = _mm_load_pd(invAx_); // loading coefficient
		__m128d _ay = _mm_load_pd(invAy_); // loading coefficient

		_x = _mm_sub_pd(_x, _b);

		__m128d _X = _mm_mul_pd(_ax, _x);
		__m128d _Y = _mm_mul_pd(_ay, _x);

		double _vx[2] __attribute__((aligned(16)));
		double _vy[2] __attribute__((aligned(16)));
		_mm_store_pd(_vx, _X);
		_mm_store_pd(_vy, _Y);

		double const X(_vx[0] + _vx[1]);
		double const Y(_vy[0] + _vy[1]);

		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		return invAx_[1];
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @return The value of derivative \f$\frac{dY}{dx}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dYdx (double const & x, double const & y) const {
		__m128d _x = _mm_set_pd(y, x); // loading coordinate
		__m128d _b = _mm_load_pd(b_);  // loading vector

		__m128d _ax = _mm_load_pd(invAx_); // loading coefficient
		__m128d _ay = _mm_load_pd(invAy_); // loading coefficient

		_x = _mm_sub_pd(_x, _b);

		__m128d _X = _mm_mul_pd(_ax, _x);
		__m128d _Y = _mm_mul_pd(_ay, _x);

		double _vx[2] __attribute__((aligned(16)));
		double _vy[2] __attribute__((aligned(16)));
		_mm_store_pd(_vx, _X);
		_mm_store_pd(_vy, _Y);

		double const X(_vx[0] + _vx[1]);
		double const Y(_vy[0] + _vy[1]);

		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		return invAy_[0];
	}

	/*!
	 * @brief The derivative of transformation \f$(\varphi^{-1})'\f$
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @return The value of derivative \f$\frac{dY}{dy}=(\varphi^{-1})'(\mathbf{x})\f$
	 */
	inline double dYdy (double const & x, double const & y) const {
		__m128d _x = _mm_set_pd(y, x); // loading coordinate
		__m128d _b = _mm_load_pd(b_);  // loading vector

		__m128d _ax = _mm_load_pd(invAx_); // loading coefficient
		__m128d _ay = _mm_load_pd(invAy_); // loading coefficient

		_x = _mm_sub_pd(_x, _b);

		__m128d _X = _mm_mul_pd(_ax, _x);
		__m128d _Y = _mm_mul_pd(_ay, _x);

		double _vx[2] __attribute__((aligned(16)));
		double _vy[2] __attribute__((aligned(16)));
		_mm_store_pd(_vx, _X);
		_mm_store_pd(_vy, _Y);

		double const X(_vx[0] + _vx[1]);
		double const Y(_vy[0] + _vy[1]);

		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		return invAy_[1];
	}

	/*!
	 * @brief Jacobian of transformation
	 * @param X The first coordinate
	 * @param Y The second coordinate
	 * @return The value of Jacobian \f$\det\left[\frac{dx_i}{dX_j}\right](\mathbf{X})\f$
	 */
	inline double det (double const & X, double const & Y) const {
		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle
		return detA_;
	}

	/*!
	 * @brief Jacobian of transformation
	 * @param x The first coordinate
	 * @param y The second coordinate
	 * @return The value of Jacobian \f$\det\left[\frac{dX_i}{dx_j}\right](\mathbf{x})\f$
	 */
	inline double DET (double const & x, double const & y) const {
		__m128d _x = _mm_set_pd(y, x); // loading coordinate
		__m128d _b = _mm_load_pd(b_);  // loading vector

		__m128d _ax = _mm_load_pd(invAx_); // loading coefficient
		__m128d _ay = _mm_load_pd(invAy_); // loading coefficient

		_x = _mm_sub_pd(_x, _b);

		__m128d _X = _mm_mul_pd(_ax, _x);
		__m128d _Y = _mm_mul_pd(_ay, _x);

		double _vx[2] __attribute__((aligned(16)));
		double _vy[2] __attribute__((aligned(16)));
		_mm_store_pd(_vx, _X);
		_mm_store_pd(_vy, _Y);

		double const X(_vx[0] + _vx[1]);
		double const Y(_vy[0] + _vy[1]);

		gas_assert(unit_t::in(X, Y)); // The point must be in the triangle

		return 1./detA_;
	}

private:
	/* the matrix of transformation */
	double Ax_[2] __attribute__((aligned(16)));
	double Ay_[2] __attribute__((aligned(16)));

	/* the matrix of inverse transformation */
	double invAx_[2] __attribute__((aligned(16)));
	double invAy_[2] __attribute__((aligned(16)));

	/* the vector of affine transformation */
	double b_[2];

	/* the determinant of transformation */
	double detA_;

};

} } }

#endif // _gas_geometry_map_affine_triangle_
