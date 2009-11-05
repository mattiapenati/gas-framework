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
 * @file static.h
 * @brief Classes with static method to evaluate some functions
 */

#ifndef _gas_static_
#define _gas_static_

namespace gas {

/*!
 * @brief The sum operation
 */
class add {

public:
	/*!
	 * @brief Evaluate the sum of two numbers
	 * @param a The first number
	 * @param b The second number
	 * @return The sum of two numbers
	 */
	inline static double eval (double const & a, double const & b) {
		return a + b;
	}

};

/*!
 * @brief The subtraction operation
 */
class sub {

public:
	/*!
	 * @brief Evaluate the difference of two numbers
	 * @param a The first number
	 * @param b The second number
	 * @return The difference of two numbers
	 */
	inline static double eval (double const & a, double const & b) {
		return a - b;
	}

};

/*!
 * @brief The multiplication operation
 */
class mul {

public:
	/*!
	 * @brief Evaluate the product of two numbers
	 * @param a The first number
	 * @param b The second number
	 * @return The product of two numbers
	 */
	inline static double eval (double const & a, double const & b) {
		return a * b;
	}

};

/*!
 * @brief The division operation
 */
class div {

public:
	/*!
	 * @brief Evaluate the quotient of two numbers
	 * @param a The first number
	 * @param b The second number
	 * @return The quotient of two numbers
	 */
	inline static double eval (double const & a, double const & b) {
		return a / b;
	}

};

}

#endif // _gas_static_
