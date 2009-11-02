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
 * @file derivative.h
 * @brief Functions for build the derivatives of function
 */

#ifndef _gas_functional_derivative_
#define _gas_functional_derivative_

namespace gas { namespace functional {

/*!
 * @brief Construct the derivative of a function
 * @param f A function
 * @return The derivative of the function
 */
template <typename function_>
typename function_::d_t d (function_ const & f) {
	typedef typename function_::d_t r_t;
	return r_t(f);
}

/*!
 * @brief Construct the derivative by first argument of a function
 * @param f A function
 * @return The derivative by first argument of the function
 */
template <typename function_>
typename function_::dx_t dx (function_ const & f) {
	typedef typename function_::dx_t rx_t;
	return rx_t(f);
}

/*!
 * @brief Construct the derivative by second argument of a function
 * @param f A function
 * @return The derivative by second argument of the function
 */
template <typename function_>
typename function_::dy_t dy (function_ const & f) {
	typedef typename function_::dy_t ry_t;
	return ry_t(f);
}

} }

#endif // _gas_functional_derivative_
