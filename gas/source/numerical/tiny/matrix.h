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
 * @file matrix.h
 * @brief Matrix of fixed size at compile time
 */

#ifndef _gas_numerical_tiny_matrix_
#define _gas_numerical_tiny_matrix_

#include "../../gas/static.h"

namespace gas { namespace numerical { namespace tiny {

template <unsigned int row_, unsigned int col_>
class matrix;

template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
class matrix_binexp;

/*!
 * @brief Matrix class with fixed size at compile time
 * @param row_ The number of rows of matrix
 * @param col_ The number of columns of matrix
 */
template <unsigned int row_, unsigned int col_>
class matrix {

public:
	/*! @brief The self type */
	typedef matrix<row_, col_> self_t;

	/*! @brief The number of rows of matrix */
	static const unsigned int row = row_;

	/*! @brief The number of columns of matrix */
	static const unsigned int col = col_;

public:
	/*!
	 *@brief The constructor by value, all components with the same value
	 *@param s The value for all components
	 */
	inline matrix (double const & s = 0.) {
		rangeu(i, row) rangeu(j, col)
			a_[i][j] = s;
	}

	/*!
	 * @brief The copy constructor
	 * @param m An other matrix
	 */
	inline matrix (self_t const & m) {
		rangeu(i, row) rangeu(j, col)
			a_[i][j] = m(i,j);
	}

	/*!
	 * @brief The expression constructor
	 * @param exp An expression
	 */
	template <typename left_, typename right_, typename operator_>
	inline matrix (matrix_binexp<left_, right_, row_, col_, operator_> const & exp) {
		rangeu(i, row) range(j, 0, col)
				a_[i][j] = exp(i,j);
	}

	/*!
	 * @brief Access to a specif component
	 * @param i The index of row
	 * @param j The index of column
	 * @return The reference to the i,j component
	 */
	inline double & operator() (unsigned int const & i, unsigned int const & j) {
		gas_pre((i < row) and (j < col));
		return a_[i][j];
	}

	/*!
	 * @brief Access to a specif component (const version)
	 * @param i The index of row
	 * @param j The index of column
	 * @return The reference to the i,j component
	 */
	inline double const & operator() (unsigned int const & i, unsigned int const & j) const {
		gas_pre((i < row) and (j < col));
		return a_[i][j];
	}

	/*!
	 * @brief The assign operator
	 * @param s The value for all components
	 * @return A reference to the current object
	 */
	inline self_t & operator= (double const & s) {
		rangeu(i, row) range(j, 0, col)
			a_[i][j] = s;
		return *this;
	}

	/*!
	 * @brief The copy operator
	 * @param m An other matrix
	 * @return A reference to the current object
	 */
	inline self_t & operator= (self_t const & m) {
		rangeu(i, row) rangeu(j, col)
			a_[i][j] = m(i,j);
		return *this;
	}

	/*!
	 * @brief The expression copy operator
	 * @param exp An expression
	 * @return A reference to the current object
	 */
	template <typename left_, typename right_, typename operator_>
	inline self_t & operator= (matrix_binexp<left_, right_, row_, col_, operator_> const & exp) {
		rangeu(i, row) rangeu(j, col)
				a_[i][j] = exp(i,j);
		return *this;
	}

	/*!
	 * @brief Self addition
	 * @param m A matrix
	 * @return A reference to the current object
	 */
	inline self_t & operator+= (self_t const & m) {
		rangeu(i, row) rangeu(j, col)
			a_[i][j] += m(i,j);
		return *this;
	}

	/*!
	 * @brief Self addition (expression version)
	 * @param exp An expression
	 * @return A reference to the current object
	 */
	template <typename left_, typename right_, typename operator_>
	inline self_t & operator+= (matrix_binexp<left_, right_, row_, col_, operator_> const & exp) {
		rangeu(i, row) rangeu(j, col)
			a_[i][j] += exp(i,j);
		return *this;
	}

	/*!
	 * @brief Self subtraction
	 * @param m A matrix
	 * @return A reference to the current object
	 */
	inline self_t & operator-= (self_t const & m) {
		rangeu(i, row) rangeu(j, col)
			a_[i][j] -= m(i,j);
		return *this;
	}

	/*!
	 * @brief Self subtraction (expression version)
	 * @param exp An expression
	 * @return A reference to the current object
	 */
	template <typename left_, typename right_, typename operator_>
	inline self_t & operator-= (matrix_binexp<left_, right_, row_, col_, operator_> const & exp) {
		rangeu(i, row) rangeu(j, col)
			a_[i][j] += exp(i,j);
		return *this;
	}

	/*!
	 * @brief Self multiplication
	 * @param s A scalar
	 * @return A reference to the current object
	 */
	inline self_t & operator*= (double const & s) {
		rangeu(i, row) rangeu(j, col)
				a_[i][j] *= s;
		return *this;
	}

	/*!
	 * @brief Self division
	 * @param s A scalar
	 * @return A reference to the current object
	 */
	inline self_t & operator/= (double const & s) {
		rangeu(i, row) rangeu(j, col)
				a_[i][j] /= s;
		return *this;
	}

private:
	/*! @brief The real matrix */
	double a_[row_][col_];

};

/*!
 * @brief The container for binary expression
 * @param left_ The type of left operand
 * @param right_ The type of right operand
 * @param row_ The number of rows of matrix
 * @param col_ The number of columns of matrix
 * @param operator_ The binary operator
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
class matrix_binexp {

private:
	/*! @brief The self type */
	typedef matrix_binexp<left_, right_, row_, col_, operator_> self_t;

	/*! @brief The number of rows of matrix */
	static const unsigned int row = row_;

	/*! @brief The number of columns of matrix */
	static const unsigned int col = col_;

private:
	/*!
	 * @brief The constructor
	 * @param left The left operand
	 * @param right The right operand
	 */
	inline matrix_binexp (left_ const & left, right_ const & right): l_(left), r_(right) {
	}

	/*!
	 * @brief The access operator
	 * @param i The index of row
	 * @param j The index of column
	 * @return The result of the evaluation of expression
	 */
	inline double const & operator() (unsigned int const & i, unsigned int const & j) const {
		gas_pre((i < row) and (j < col));
		return operator_::eval(l_(i,j), r_(i,j));
	}

private:
	/*! @brief The left operand */
	left_ const & l_;

	/*! @brief The right operand */
	right_ const & r_;

	friend class matrix<row_, col_>;

	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend class vector_binexp;

	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend class matrix_binexp;

	/*! @brief The sum of two matrixes */
	template <unsigned int row__, unsigned int col__>
	friend matrix_binexp<matrix<row__, col__>, matrix<row__, col__>, row__, col__, gas::add>
	operator+ (matrix<row__, col__> const &, matrix<row__, col__> const &);

	/*! @brief The sum of a matrix and an expression */
	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend matrix_binexp<matrix<row__, col__>, matrix_binexp<left__, right__, row__, col__, operator__>, row__, col__, gas::add>
	operator+ (matrix<row__, col__> const &, matrix_binexp<left__, right__, row__, col__, operator__> const &);

	/*! @brief The sum of an expression and a matrix */
	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend matrix_binexp<matrix_binexp<left__, right__, row__, col__, operator__>, matrix<row__, col__>, row__, col__, gas::add>
	operator+ (matrix_binexp<left__, right__, row__, col__, operator__> const &, matrix<row__, col__> const &);

	/*! @brief The sum of two expressions */
	template <typename left1__, typename left2__, typename right1__, typename right2__, unsigned int row__, unsigned int col__, typename operator1__, typename operator2__>
	friend matrix_binexp<matrix_binexp<left1__, right1__, row__, col__, operator1__>, matrix_binexp<left2__, right2__, row__, col__, operator2__>, row__, col__, gas::add>
	operator+ (matrix_binexp<left1__, right1__, row__, col__, operator1__> const &, matrix_binexp<left2__, right2__, row__, col__, operator2__> const &);

	/*! @brief The subtraction of two matrixes */
	template <unsigned int row__, unsigned int col__>
	friend matrix_binexp<matrix<row__, col__>, matrix<row__, col__>, row__, col__, gas::sub>
	operator- (matrix<row__, col__> const &, matrix<row__, col__> const &);

	/*! @brief The subtraction of a matrix and an expression */
	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend matrix_binexp<matrix<row__, col__>, matrix_binexp<left__, right__, row__, col__, operator__>, row__, col__, gas::sub>
	operator- (matrix<row__, col__> const &, matrix_binexp<left__, right__, row__, col__, operator__> const &);

	/*! @brief The subtraction of an expression and a matrix */
	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend matrix_binexp<matrix_binexp<left__, right__, row__, col__, operator__>, matrix<row__, col__>, row__, col__, gas::sub>
	operator- (matrix_binexp<left__, right__, row__, col__, operator__> const &, matrix<row__, col__> const &);

	/*! @brief The subtraction of two expressions */
	template <typename left1__, typename left2__, typename right1__, typename right2__, unsigned int row__, unsigned int col__, typename operator1__, typename operator2__>
	friend matrix_binexp<matrix_binexp<left1__, right1__, row__, col__, operator1__>, matrix_binexp<left2__, right2__, row__, col__, operator2__>, row__, col__, gas::sub>
	operator- (matrix_binexp<left1__, right1__, row__, col__, operator1__> const &, matrix_binexp<left2__, right2__, row__, col__, operator2__> const &);

	/*! @brief The determinant of a square expression */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend double det (matrix_binexp<left__, right__, size__, size__, operator__> const &);

	/*! @brief The determinant of a square expression (2x2) */
	template <typename left__, typename right__, typename operator__>
	inline double det (matrix_binexp<left__, right__, 2u, 2u, operator__> const &);

};

/*!
 * @brief The container for binary expression (specialization for operation with
 *        scalar)
 * @param left_ The type of left operand
 * @param row_ The number of rows of matrix
 * @param col_ The number of columns of matrix
 * @param operator_ The binary operator
 */
template <typename left_, unsigned int row_, unsigned int col_, typename operator_>
class matrix_binexp<left_, double, row_, col_, operator_> {

private:
	/*! @brief The self type */
	typedef matrix_binexp<left_, double, row_, col_, operator_> self_t;

	/*! @brief The number of rows of matrix */
	static const unsigned int row = row_;

	/*! @brief The number of columns of matrix */
	static const unsigned int col = col_;

private:
	/*!
	 * @brief The constructor
	 * @param left The left operand
	 * @param right The right operand
	 */
	inline matrix_binexp (left_ const & left, double const & right): l_(left), r_(right) {
	}

	/*!
	 * @brief The access operator
	 * @param i The index of row
	 * @param j The index of column
	 * @return The result of the evaluation of expression
	 */
	inline double const & operator() (unsigned int const & i, unsigned int const & j) const {
		gas_pre((i < row) and (j < col));
		return operator_::eval(l_(i,j), r_);
	}

private:
	/*! @brief The left operand */
	left_ const & l_;

	/*! @brief The right operand */
	double const & r_;

	friend class matrix<row_, col_>;

	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend class vector_binexp;

	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend class matrix_binexp;

	/*! @brief The multiplication by a scalar */
	template <unsigned int row__, unsigned int col__>
	friend matrix_binexp<matrix<row__, col__>, double, row__, col__, gas::mul>
	operator* (matrix<row__, col__> const &, double const &);

	/*! @brief The multiplication by a scalar */
	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend matrix_binexp<matrix_binexp<left__, right__, row__, col__, operator__>, double, row__, col__, gas::mul>
	operator* (matrix_binexp<left__, right__, row__, col__, operator__> const &, double const &);

	/*! @brief The division by a scalar */
	template <unsigned int row__, unsigned int col__>
	friend matrix_binexp<matrix<row__, col__>, double, row__, col__, gas::div>
	operator/ (matrix<row__, col__> const &, double const &);

	/*! @brief The division by a scalar */
	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend matrix_binexp<matrix_binexp<left__, right__, row__, col__, operator__>, double, row__, col__, gas::div>
	operator/ (matrix_binexp<left__, right__, row__, col__, operator__> const &, double const &);

	/*! @brief The determinant of a square expression */
	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend double det (matrix_binexp<left__, right__, size__, size__, operator__> const &);

	/*! @brief The determinant of a square expression (2x2) */
	template <typename left__, typename right__, typename operator__>
	inline double det (matrix_binexp<left__, right__, 2u, 2u, operator__> const &);

};


/*!
 * @brief The sum of two matrixes
 * @param a The first matrix
 * @param b The second matrix
 * @return An expression containing the operations
 */
template <unsigned int row_, unsigned int col_>
inline matrix_binexp<matrix<row_, col_>, matrix<row_, col_>, row_, col_, gas::add>
operator+ (matrix<row_, col_> const & a, matrix<row_, col_> const & b) {
	return matrix_binexp<matrix<row_, col_>, matrix<row_, col_>, row_, col_, gas::add>(a, b);
}

/*!
 * @brief The sum of a matrix and an expression
 * @param a A matrix
 * @param b An expression
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
inline matrix_binexp<matrix<row_, col_>, matrix_binexp<left_, right_, row_, col_, operator_>, row_, col_, gas::add>
operator+ (matrix<row_, col_> const & a, matrix_binexp<left_, right_, row_, col_, operator_> const & b) {
	return matrix_binexp<matrix<row_, col_>, matrix_binexp<left_, right_, row_, col_, operator_>, row_, col_, gas::add>(a, b);
}

/*!
 * @brief The sum of an expression and a matrix
 * @param a An expression
 * @param b A matrix
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
inline matrix_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, matrix<row_, col_>, row_, col_, gas::add>
operator+ (matrix_binexp<left_, right_, row_, col_, operator_> const & a, matrix<row_, col_> const & b) {
	return matrix_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, matrix<row_, col_>, row_, col_, gas::add>(a, b);
}

/*!
 * @brief The sum of two expressions
 * @param a The first expression
 * @param b The second expression
 * @return An expression containing the operations
 */
template <typename left1_, typename left2_, typename right1_, typename right2_, unsigned int row_, unsigned int col_, typename operator1_, typename operator2_>
inline matrix_binexp<matrix_binexp<left1_, right1_, row_, col_, operator1_>, matrix_binexp<left2_, right2_, row_, col_, operator2_>, row_, col_, gas::add>
operator+ (matrix_binexp<left1_, right1_, row_, col_, operator1_> const & a, matrix_binexp<left2_, right2_, row_, col_, operator2_> const & b) {
	return matrix_binexp<matrix_binexp<left1_, right1_, row_, col_, operator1_>, matrix_binexp<left2_, right2_, row_, col_, operator2_>, row_, col_, gas::add>(a, b);
}


/*!
 * @brief The subtraction of two matrixes
 * @param a The first matrix
 * @param b The second matrix
 * @return An expression containing the operations
 */
template <unsigned int row_, unsigned int col_>
inline matrix_binexp<matrix<row_, col_>, matrix<row_, col_>, row_, col_, gas::sub>
operator- (matrix<row_, col_> const & a, matrix<row_, col_> const & b) {
	return matrix_binexp<matrix<row_, col_>, matrix<row_, col_>, row_, col_, gas::sub>(a, b);
}

/*!
 * @brief The subtraction of a matrix and an expression
 * @param a A matrix
 * @param b An expression
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
inline matrix_binexp<matrix<row_, col_>, matrix_binexp<left_, right_, row_, col_, operator_>, row_, col_, gas::sub>
operator- (matrix<row_, col_> const & a, matrix_binexp<left_, right_, row_, col_, operator_> const & b) {
	return matrix_binexp<matrix<row_, col_>, matrix_binexp<left_, right_, row_, col_, operator_>, row_, col_, gas::sub>(a, b);
}

/*!
 * @brief The subtraction of an expression and a matrix
 * @param a An expression
 * @param b A matrix
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
inline matrix_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, matrix<row_, col_>, row_, col_, gas::sub>
operator- (matrix_binexp<left_, right_, row_, col_, operator_> const & a, matrix<row_, col_> const & b) {
	return matrix_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, matrix<row_, col_>, row_, col_, gas::sub>(a, b);
}

/*!
 * @brief The subtraction of two expressions
 * @param a The first expression
 * @param b The second expression
 * @return An expression containing the operations
 */
template <typename left1_, typename left2_, typename right1_, typename right2_, unsigned int row_, unsigned int col_, typename operator1_, typename operator2_>
inline matrix_binexp<matrix_binexp<left1_, right1_, row_, col_, operator1_>, matrix_binexp<left2_, right2_, row_, col_, operator2_>, row_, col_, gas::sub>
operator- (matrix_binexp<left1_, right1_, row_, col_, operator1_> const & a, matrix_binexp<left2_, right2_, row_, col_, operator2_> const & b) {
	return matrix_binexp<matrix_binexp<left1_, right1_, row_, col_, operator1_>, matrix_binexp<left2_, right2_, row_, col_, operator2_>, row_, col_, gas::sub>(a, b);
}

/*!
 * @brief The multiplication by a scalar
 * @param a A matrix
 * @param b A scalar
 * @return An expression containing the operations
 */
template <unsigned int row_, unsigned int col_>
inline matrix_binexp<matrix<row_, col_>, double, row_, col_, gas::mul>
operator* (matrix<row_, col_> const & a, double const & b) {
	return matrix_binexp<matrix<row_, col_>, double, row_, col_, gas::mul>(a, b);
}

/*!
 * @brief The multiplication by a scalar
 * @param a An expression
 * @param b A scalar
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
inline matrix_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, double, row_, col_, gas::mul>
operator* (matrix_binexp<left_, right_, row_, col_, operator_> const & a, double const & b) {
	return matrix_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, double, row_, col_, gas::mul>(a, b);
}

/*!
 * @brief The division by a scalar
 * @param a A matrix
 * @param b A scalar
 * @return An expression containing the operations
 */
template <unsigned int row_, unsigned int col_>
inline matrix_binexp<matrix<row_, col_>, double, row_, col_, gas::div>
operator/ (matrix<row_, col_> const & a, double const & b) {
	return matrix_binexp<matrix<row_, col_>, double, row_, col_, gas::div>(a, b);
}

/*!
 * @brief The division by a scalar
 * @param a An expression
 * @param b A scalar
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
inline matrix_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, double, row_, col_, gas::div>
operator/ (matrix_binexp<left_, right_, row_, col_, operator_> const & a, double const & b) {
	return matrix_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, double, row_, col_, gas::div>(a, b);
}

} } }

#endif // _gas_numerical_tiny_matrix_
