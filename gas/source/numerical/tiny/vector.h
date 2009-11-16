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
 * @file vector.h
 * @brief Vector of fixed size at compile time
 */

#ifndef _gas_numerical_tiny_vector_
#define _gas_numerical_tiny_vector_

#include "../../gas"

namespace gas { namespace numerical { namespace tiny {

template <unsigned int size_>
class vector;

template <typename left_, typename right_, unsigned int size_, typename operator_>
class vector_binexp;

/*!
 * @brief Vector class with fixed size at compile time
 * @param size_ The dimension of vector
 */
template <unsigned int size_>
class vector {

public:
	/*! @brief The self type */
	typedef vector<size_> self_t;

	/*! @brief The size of vector */
	static const unsigned int size = size_;

public:
	/*!
	 *@brief The constructor by value, all components with the same value
	 *@param s The value for all components
	 */
	inline vector (double const & s = 0.) {
		gas_rangeu(i, size)
			v_[i] = s;
	}

	/*!
	 * @brief The copy constructor
	 * @param v An other vector
	 */
	inline vector (self_t const & v) {
		gas_rangeu(i, size)
			v_[i] = v(i);
	}

	/*!
	 * @brief The expression constructor
	 * @param exp An expression
	 */
	template <typename left_, typename right_, typename operator_>
	inline vector (vector_binexp<left_, right_, size, operator_> const & exp) {
		gas_rangeu(i, size)
			v_[i] = exp(i);
	}

	/*!
	 * @brief Access to a specif component
	 * @param i The index of the component
	 * @return The reference to the i-th component
	 */
	inline double & operator() (unsigned int const & i) {
		gas_pre(i < size);
		return v_[i];
	}

	/*!
	 * @brief Access to a specif component (const version)
	 * @param i The index of the component
	 * @return The const reference to the i-th component
	 */
	inline double const & operator() (unsigned int const & i) const {
		gas_pre(i < size);
		return v_[i];
	}

	/*!
	 * @brief The assign operator
	 * @param s The value for all components
	 * @return A reference to the current object
	 */
	inline self_t & operator= (double const & s) {
		gas_rangeu(i, size)
			v_[i] = s;
		return *this;
	}

	/*!
	 * @brief The copy operator
	 * @param v An other vector
	 * @return A reference to the current object
	 */
	inline self_t & operator= (self_t const & v) {
		gas_rangeu(i, size)
			v_[i] = v(i);
		return *this;
	}

	/*!
	 * @brief The expression copy operator
	 * @param exp An expression
	 * @return A reference to the current object
	 */
	template <typename left_, typename right_, typename operator_>
	inline self_t & operator= (vector_binexp<left_, right_, size, operator_> const & exp) {
		gas_rangeu(i, size)
			v_[i] = exp(i);
		return *this;
	}

	/*!
	 * @brief Self addition
	 * @param v A vector
	 * @return A reference to the current object
	 */
	inline self_t & operator+= (self_t const & v) {
		gas_rangeu(i, size)
			v_[i] += v(i);
		return *this;
	}

	/*!
	 * @brief Self addition (expression version)
	 * @param exp An expression
	 * @return A reference to the current object
	 */
	template <typename left_, typename right_, typename operator_>
	inline self_t & operator+= (vector_binexp<left_, right_, size, operator_> const & exp) {
		gas_rangeu(i, size)
			v_[i] += exp(i);
		return *this;
	}

	/*!
	 * @brief Self subtraction
	 * @param v A vector
	 * @return A reference to the current object
	 */
	inline self_t & operator-= (self_t const & v) {
		gas_rangeu(i, size)
			v_[i] -= v(i);
		return *this;
	}

	/*!
	 * @brief Self subtraction (expression version)
	 * @param exp An expression
	 * @return A reference to the current object
	 */
	template <typename left_, typename right_, typename operator_>
	inline self_t & operator-= (vector_binexp<left_, right_, size, operator_> const & exp) {
		gas_rangeu(i, size)
			v_[i] -= exp(i);
		return *this;
	}

	/*!
	 * @brief Self multiplication
	 * @param s A scalar
	 * @return A reference to the current object
	 */
	inline self_t & operator*= (double const & s) {
		gas_rangeu(i, size)
			v_[i] *= s;
		return *this;
	}

	/*!
	 * @brief Self division
	 * @param s A scalar
	 * @return A reference to the current object
	 */
	inline self_t & operator/= (double const & s) {
		gas_rangeu(i, size)
			v_[i] /= s;
		return *this;
	}

private:
	/*! @brief The real vector */
	double v_[size_];

};

/*!
 * @brief The container for binary expression
 * @param left_ The type of left operand
 * @param right_ The type of right operand
 * @param size_ The size of expression
 * @param operator_ The binary operator
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
class vector_binexp {

private:
	/*! @brief The self type */
	typedef vector_binexp<left_, right_, size_, operator_> self_t;

	/*! @brief The size of expression */
	static const unsigned int size = size_;

private:
	/*!
	 * @brief The constructor
	 * @param left The left operand
	 * @param right The right operand
	 */
	inline vector_binexp (left_ const & left, right_ const & right): l_(left), r_(right) {
	}

	/*!
	 * @brief The access operator
	 * @param i The index of the component
	 * @return The result of the evaluation of expression
	 */
	inline double operator() (unsigned int const & i) const {
		gas_assert(i < size);
		return operator_::eval(l_(i), r_(i));
	}

private:
	/*! @brief The left operand */
	left_ const & l_;

	/*! @brief The right operand */
	right_ const & r_;

	friend class vector<size_>;

	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend class vector_binexp;

	/*! @brief The sum of two vectors */
	template <unsigned int size__>
	friend vector_binexp<vector<size__>, vector<size__>, size__, gas::add>
	operator+ (vector<size__> const &, vector<size__> const &);

	/*! @brief The sum of a vector and an expression */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend vector_binexp<vector<size__>, vector_binexp<left__, right__, size__, operator__>, size__, gas::add>
	operator+ (vector<size__> const &, vector_binexp<left__, right__, size__, operator__> const &);

	/*! @brief The sum of an expression and a vector */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend vector_binexp<vector_binexp<left__, right__, size__, operator__>, vector<size__>, size__, gas::add>
	operator+ (vector_binexp<left__, right__, size__, operator__> const &, vector<size__> const &);

	/*! @brief The sum of two expressions */
	template <typename left1__, typename left2__, typename right1__, typename right2__, unsigned int size__, typename operator1__, typename operator2__>
	friend vector_binexp<vector_binexp<left1__, right1__, size__, operator1__>, vector_binexp<left2__, right2__, size__, operator2__>, size__, gas::add>
	operator+ (vector_binexp<left1__, right1__, size__, operator1__> const &, vector_binexp<left2__, right2__, size_, operator2__> const &);

	/*! @brief The subtraction of two vectors */
	template <unsigned int size__>
	friend vector_binexp<vector<size__>, vector<size__>, size__, gas::sub>
	operator- (vector<size__> const &, vector<size__> const &);

	/*! @brief The subtraction of a vector and an expression */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend vector_binexp<vector<size__>, vector_binexp<left__, right__, size__, operator__>, size__, gas::sub>
	operator- (vector<size__> const &, vector_binexp<left__, right__, size__, operator__> const &);

	/*! @brief The subtraction of an expression and a vector */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend vector_binexp<vector_binexp<left__, right__, size__, operator__>, vector<size__>, size__, gas::sub>
	operator- (vector_binexp<left__, right__, size__, operator__> const &, vector<size__> const &);

	/*! @brief The subtraction of two expressions */
	template <typename left1__, typename left2__, typename right1__, typename right2__, unsigned int size__, typename operator1__, typename operator2__>
	friend vector_binexp<vector_binexp<left1__, right1__, size__, operator1__>, vector_binexp<left2__, right2__, size__, operator2__>, size__, gas::sub>
	operator- (vector_binexp<left1__, right1__, size__, operator1__> const &, vector_binexp<left2__, right2__, size_, operator2__> const &);

	/*! @brief The dot product between a vector and an expression */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend double dot (vector<size__> const &, vector_binexp<left__, right__, size__, operator__> const &);

	/*! @brief The dot product between an expression and a vector */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend double dot (vector_binexp<left__, right__, size__, operator__> const &, vector<size__> const &);

	/*! @brief The dot product between two expressions */
	template <typename left1__, typename left2__, typename right1__, typename right2__, unsigned int size__, typename operator1__, typename operator2__>
	friend double dot (vector_binexp<left1__, right1__, size__, operator1__> const &, vector_binexp<left2__, right2__, size_, operator2__> const &);

	/*! @brief The norm of an expression */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend double norm (vector_binexp<left__, right__, size__, operator__> const &);

};

/*!
 * @brief The container for binary expression (specialization for operation with
 *        scalar)
 * @param left_ The type of left operand
 * @param size_ The size of expression
 * @param operator_ The binary operator
 */
template <typename left_, unsigned int size_, typename operator_>
class vector_binexp<left_, double, size_, operator_> {

private:
	/*! @brief The self type */
	typedef vector_binexp<left_, double, size_, operator_> self_t;

	/*! @brief The size of expression */
	static const unsigned int size = size_;

private:
	/*! @brief The constructor */
	inline vector_binexp (left_ const & left, double const & right): l_(left), r_(right) {
	}

	/*!
	 * @brief The access operator
	 * @param i The index of the component
	 * @return The result of the evaluation of expression
	 */
	inline double operator() (unsigned int const & i) const {
		gas_assert(i < size);
		return operator_::eval(l_(i), r_);
	}

private:
	/*! @brief The left operand */
	left_ const & l_;

	/*! @brief The right operand */
	double const & r_;

	friend class vector<size_>;

	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend class vector_binexp;

	/*! @brief The multiplication by a scalar */
	template <unsigned int size__>
	friend vector_binexp<vector<size__>, double, size__, gas::mul>
	operator* (vector<size__> const &, double const &);

	/*! @brief The multiplication by a scalar */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend vector_binexp<vector_binexp<left__, right__, size__, operator__>, double, size__, gas::mul>
	operator* (vector_binexp<left__, right__, size__, operator__> const &, double const &);

	/*! @brief The division by a scalar */
	template <unsigned int size__>
	friend vector_binexp<vector<size__>, double, size__, gas::div>
	operator/ (vector<size__> const &, double const &);

	/*! @brief The division by a scalar */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend vector_binexp<vector_binexp<left__, right__, size__, operator__>, double, size__, gas::div>
	operator/ (vector_binexp<left__, right__, size__, operator__> const &, double const &);

	/*! @brief The dot product between a vector and an expression */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend double dot (vector<size__> const &, vector_binexp<left__, right__, size__, operator__> const &);

	/*! @brief The dot product between an expression and a vector */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend double dot (vector_binexp<left__, right__, size__, operator__> const &, vector<size__> const &);

	/*! @brief The dot product between two expressions */
	template <typename left1__, typename left2__, typename right1__, typename right2__, unsigned int size__, typename operator1__, typename operator2__>
	friend double dot (vector_binexp<left1__, right1__, size__, operator1__> const &, vector_binexp<left2__, right2__, size_, operator2__> const &);

	/*! @brief The norm of an expression */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend double norm (vector_binexp<left__, right__, size__, operator__> const &);

};

/*!
 * @brief The sum of two vectors
 * @param a The first vector
 * @param b The second vector
 * @return An expression containing the operations
 */
template <unsigned int size_>
inline vector_binexp<vector<size_>, vector<size_>, size_, gas::add>
operator+ (vector<size_> const & a, vector<size_> const & b) {
	return vector_binexp<vector<size_>, vector<size_>, size_, gas::add>(a, b);
}

/*!
 * @brief The sum of a vector and an expression
 * @param a A vector
 * @param b An expression
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline vector_binexp<vector<size_>, vector_binexp<left_, right_, size_, operator_>, size_, gas::add>
operator+ (vector<size_> const & a, vector_binexp<left_, right_, size_, operator_> const & b) {
	return vector_binexp<vector<size_>, vector_binexp<left_, right_, size_, operator_>, size_, gas::add>(a, b);
}

/*!
 * @brief The sum of an expression and a vector
 * @param a An expression
 * @param b A vector
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline vector_binexp<vector_binexp<left_, right_, size_, operator_>, vector<size_>, size_, gas::add>
operator+ (vector_binexp<left_, right_, size_, operator_> const & a, vector<size_> const & b) {
	return vector_binexp<vector_binexp<left_, right_, size_, operator_>, vector<size_>, size_, gas::add>(a, b);
}

/*!
 * @brief The sum of two expressions
 * @param a The first expression
 * @param b The second expression
 * @return An expression containing the operations
 */
template <typename left1_, typename left2_, typename right1_, typename right2_, unsigned int size_, typename operator1_, typename operator2_>
inline vector_binexp<vector_binexp<left1_, right1_, size_, operator1_>, vector_binexp<left2_, right2_, size_, operator2_>, size_, gas::add>
operator+ (vector_binexp<left1_, right1_, size_, operator1_> const & a, vector_binexp<left2_, right2_, size_, operator2_> const & b) {
	return vector_binexp<vector_binexp<left1_, right1_, size_, operator1_>, vector_binexp<left2_, right2_, size_, operator2_>, size_, gas::add>(a, b);
}

/*!
 * @brief The subtraction of two vectors
 * @param a The first vector
 * @param b The second vector
 * @return An expression containing the operations
 */
template <unsigned int size_>
inline vector_binexp<vector<size_>, vector<size_>, size_, gas::sub>
operator- (vector<size_> const & a, vector<size_> const & b) {
	return vector_binexp<vector<size_>, vector<size_>, size_, gas::sub>(a, b);
}

/*!
 * @brief The subtraction of a vector and an expression
 * @param a A vector
 * @param b An expression
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline vector_binexp<vector<size_>, vector_binexp<left_, right_, size_, operator_>, size_, gas::sub>
operator- (vector<size_> const & a, vector_binexp<left_, right_, size_, operator_> const & b) {
	return vector_binexp<vector<size_>, vector_binexp<left_, right_, size_, operator_>, size_, gas::sub>(a, b);
}

/*!
 * @brief The subtraction of an expression and a vector
 * @param a An expression
 * @param b A vector
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline vector_binexp<vector_binexp<left_, right_, size_, operator_>, vector<size_>, size_, gas::sub>
operator- (vector_binexp<left_, right_, size_, operator_> const & a, vector<size_> const & b) {
	return vector_binexp<vector_binexp<left_, right_, size_, operator_>, vector<size_>, size_, gas::sub>(a, b);
}

/*!
 * @brief The subtraction of two expressions
 * @param a The first expression
 * @param b The second expression
 * @return An expression containing the operations
 */
template <typename left1_, typename left2_, typename right1_, typename right2_, unsigned int size_, typename operator1_, typename operator2_>
inline vector_binexp<vector_binexp<left1_, right1_, size_, operator1_>, vector_binexp<left2_, right2_, size_, operator2_>, size_, gas::sub>
operator- (vector_binexp<left1_, right1_, size_, operator1_> const & a, vector_binexp<left2_, right2_, size_, operator2_> const & b) {
	return vector_binexp<vector_binexp<left1_, right1_, size_, operator1_>, vector_binexp<left2_, right2_, size_, operator2_>, size_, gas::sub>(a, b);
}

/*!
 * @brief The multiplication by a scalar
 * @param a A vector
 * @param b A scalar
 * @return An expression containing the operations
 */
template <unsigned int size_>
inline vector_binexp<vector<size_>, double, size_, gas::mul>
operator* (vector<size_> const & a, double const & b) {
	return vector_binexp<vector<size_>, double, size_, gas::mul>(a, b);
}

/*!
 * @brief The multiplication by a scalar
 * @param a An expression
 * @param b A scalar
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline vector_binexp<vector_binexp<left_, right_, size_, operator_>, double, size_, gas::mul>
operator* (vector_binexp<left_, right_, size_, operator_> const & a, double const & b) {
	return vector_binexp<vector_binexp<left_, right_, size_, operator_>, double, size_, gas::mul>(a, b);
}

/*!
 * @brief The division by a scalar
 * @param a A vector
 * @param b A scalar
 * @return An expression containing the operations
 */
template <unsigned int size_>
inline vector_binexp<vector<size_>, double, size_, gas::div>
operator/ (vector<size_> const & a, double const & b) {
	return vector_binexp<vector<size_>, double, size_, gas::div>(a, b);
}

/*!
 * @brief The division by a scalar
 * @param a An expression
 * @param b A scalar
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline vector_binexp<vector_binexp<left_, right_, size_, operator_>, double, size_, gas::div>
operator/ (vector_binexp<left_, right_, size_, operator_> const & a, double const & b) {
	return vector_binexp<vector_binexp<left_, right_, size_, operator_>, double, size_, gas::div>(a, b);
}

} } }

#endif // _gas_numerical_tiny_vector_
