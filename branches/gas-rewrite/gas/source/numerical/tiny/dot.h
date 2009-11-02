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
 * @file dot.h
 * @brief The dot product between two tiny vector (and norm)
 */

#ifndef _gas_numeric_tiny_dot_
#define _gas_numeric_tiny_dot_

#include "../../gas/macro.h"
#include "vector.h"

namespace gas { namespace numerical { namespace tiny {

/*!
 * @brief The dot product between two vectors
 * @param a The fisrt vector
 * @param b The second vector
 * @return The dot product
 */
template <unsigned int size_>
inline double dot (vector<size_> const & a, vector<size_> const & b) {
	double r(0.);
	range(i, 0, size_)
		r += (a(i) * b(i));
	return r;
}

/*!
 * @brief The dot product between a vector and an expression
 * @param a A vector
 * @param b An expression
 * @return The dot product
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline double dot (vector<size_> const & a, vector_binexp<left_, right_, size_, operator_> const & b) {
	double r(0.);
	range(i, 0, size_)
		r += (a(i) * b(i));
	return r;
}

/*!
 * @brief The dot product between an expression and a vector
 * @param a An expression
 * @param b A vector
 * @return The dot product
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline double dot (vector_binexp<left_, right_, size_, operator_> const & a, vector<size_> const & b) {
	double r(0.);
	range(i, 0, size_)
		r += (a(i) * b(i));
	return r;
}

/*!
 * @brief The dot product between two expressions
 * @param a The first expression
 * @param b The second expression
 * @return The dot product
 */
template <typename left1_, typename left2_, typename right1_, typename right2_, unsigned int size_, typename operator1_, typename operator2_>
inline double dot (vector_binexp<left1_, right1_, size_, operator1_> const & a, vector_binexp<left2_, right2_, size_, operator2_> const & b) {
	double r(0.);
	range(i, 0, size_)
		r += (a(i) * b(i));
	return r;
}

/*!
 * @brief The norm of a vector
 * @param a A vector
 * @return The norm
 */
template <unsigned int size_>
inline double norm (vector<size_> const & a) {
	double r(0.);
	range(i, 0, size_) {
		double t(a(i));
		r += (t*t);
	}
	return r;
}

/*!
 * @brief The norm of an expression
 * @param a An expression
 * @return The norm
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
inline double norm (vector_binexp<left_, right_, size_, operator_> const & a) {
	double r(0.);
	range(i, 0, size_) {
		double t(a(i));
		r += (t*t);
	}
	return r;
}

} } }

#endif // _gas_numeric_tiny_dot_
