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
 * @file mul.h
 * @brief The multiplication between matrices and vectors
 */

#ifndef _gas_numerical_tiny_mul_
#define _gas_numerical_tiny_mul_

#include "matrix.h"
#include "vector.h"

namespace gas { namespace numerical { namespace tiny {

/*!
 * @brief An empty class to identify the product between matrix and vector
 */
class mul_matvec {
};

/*!
 * @brief The container for binary expression (specialization for matrix
 *        multiplication)
 * @param row_ The number of rows of matrix
 * @param size_ The size of expression
 */
template <typename left_, typename right_, unsigned int size_>
class vector_binexp<left_, right_, size_, mul_matvec> {

private:
	/*! @brief The self type */
	typedef vector_binexp<left_, right_, size_, mul_matvec> self_t;

	/*! @brief The size of expression */
	static const unsigned int size = size_;

private:
	/*! @brief The constructor */
	inline vector_binexp (left_ const & left, right_ const & right): l_(left), r_(right) {
	}

	/*!
	 * @brief The access operator
	 * @param i The index of the component
	 * @return The result of the evaluation of expression
	 */
	inline double operator() (unsigned int const & i) const {
		gas_assert(i < size);
		double r(0.);
		range(j, 0, left_::col)
			r += (l_(i,j) * r_(j));
		return r;
	}

private:
	/*! @brief The left operand */
	left_ const & l_;

	/*! @brief The right operand */
	right_ const & r_;

	friend class vector<size_>;

	template <typename left__, typename right__, unsigned int size__, typename operator__>
	friend class vector_binexp;

	/*! @brief The multiplication between matrix and vector */
	template <unsigned int row__, unsigned int col__>
	friend vector_binexp<matrix<row__, col__>, vector<col__>, row__, mul_matvec>
	operator* (matrix<row__, col__> const &, vector<col__> const &);

	/*! @brief The multiplication between matrix and expression */
	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend vector_binexp<matrix<row__, col__>, vector_binexp<left__, right__, col__, operator__>, row__, mul_matvec>
	operator* (matrix<row__, col__> const &, vector_binexp<left__, right__, col__, operator__> const &);

	/*! @brief The multiplication between expression and vector */
	template <typename left__, typename right__, unsigned int row__, unsigned int col__, typename operator__>
	friend vector_binexp<matrix_binexp<left__, right__, row__, col__, operator__>, vector<col__>, row__, mul_matvec>
	operator* (matrix_binexp<left__, right__, row__, col__, operator__> const &, vector<col__> const &);

	/*! @brief The multiplication between two expressions */
	template <typename left1__, typename left2__, typename right1__, typename right2__, unsigned int row__, unsigned int col__, typename operator1__, typename operator2__>
	friend vector_binexp<matrix_binexp<left1__, right1__, row__, col__, operator1__>, vector_binexp<left2__, right2__, col__, operator2__>, row__, mul_matvec>
	operator* (matrix_binexp<left1__, right1__, row__, col__, operator1__> const &, vector_binexp<left2__, right2__, col__, operator2__> const &);

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
 * @brief The multiplication between matrix and vector
 * @param a The matrix
 * @param b The vector
 * @return An expression containing the operations
 */
template <unsigned int row_, unsigned int col_>
inline vector_binexp<matrix<row_, col_>, vector<col_>, row_, mul_matvec>
operator* (matrix<row_, col_> const & a, vector<col_> const & b) {
	return vector_binexp<matrix<row_, col_>, vector<col_>, row_, mul_matvec>(a, b);
}

/*!
 * @brief The multiplication between matrix and expression
 * @param a The matrix
 * @param b The expression
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
inline vector_binexp<matrix<row_, col_>, vector_binexp<left_, right_, col_, operator_>, row_, mul_matvec>
operator* (matrix<row_, col_> const & a, vector_binexp<left_, right_, col_, operator_> const & b) {
	return vector_binexp<matrix<row_, col_>, vector_binexp<left_, right_, col_, operator_>, row_, mul_matvec>(a, b);
}

/*!
 * @brief The multiplication between expression and vector
 * @param a The expression
 * @param b The vector
 * @return An expression containing the operations
 */
template <typename left_, typename right_, unsigned int row_, unsigned int col_, typename operator_>
inline vector_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, vector<col_>, row_, mul_matvec>
operator* (matrix_binexp<left_, right_, row_, col_, operator_> const & a, vector<col_> const & b) {
	return vector_binexp<matrix_binexp<left_, right_, row_, col_, operator_>, vector<col_>, row_, mul_matvec>(a, b);
}

/*!
 * @brief The multiplication between two expressions
 * @param a The expression
 * @param b The expression
 * @return An expression containing the operations
 */
template <typename left1_, typename left2_, typename right1_, typename right2_, unsigned int row_, unsigned int col_, typename operator1_, typename operator2_>
inline vector_binexp<matrix_binexp<left1_, right1_, row_, col_, operator1_>, vector_binexp<left2_, right2_, col_, operator2_>, row_, mul_matvec>
operator* (matrix_binexp<left1_, right1_, row_, col_, operator1_> const & a, vector_binexp<left2_, right2_, col_, operator2_> const & b) {
	return vector_binexp<matrix_binexp<left1_, right1_, row_, col_, operator1_>, vector_binexp<left2_, right2_, col_, operator2_>, row_, mul_matvec>(a, b);
}

} } }

#endif // _gas_numerical_tiny_mul_
