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
 * @file det.h
 * @brief The determinant of tiny matrix
 */

#ifndef _gas_numerical_tiny_det_
#define _gas_numerical_tiny_det_

#include "matrix.h"

namespace gas { namespace numerical { namespace tiny {

/*!
 * @brief The determinant of a square matrix
 * @param A A matrix
 * @return The determinant
 */
template <unsigned int size_>
double det (matrix<size_, size_> A) {
	// TODO da implementare
	return 0.;
}

/*!
 * @brief The determinant of a square expression
 * @param exp An expression
 * @return The determinant
 */
template <typename left_, typename right_, unsigned int size_, typename operator_>
double det (matrix_binexp<left_, right_, size_, size_, operator_> const & exp) {
	// TODO da implementare
	return 0.;
}

/*!
 * @brief The determinant of a square matrix (2x2)
 * @param A A matrix
 * @return The determinant
 */
inline double det (matrix<2u, 2u> const & A) {
	return A(0,0)*A(1,1)-A(0,1)*A(1,0);
}

/*!
 * @brief The determinant of a square expression (2x2)
 * @param exp An expression
 * @return The determinant
 */
template <typename left_, typename right_, typename operator_>
inline double det (matrix_binexp<left_, right_, 2u, 2u, operator_> const & exp) {
	return exp(0,0)*exp(1,1)-exp(0,1)*exp(1,0);
}

} } }

#endif // _gas_numerical_tiny_det_
