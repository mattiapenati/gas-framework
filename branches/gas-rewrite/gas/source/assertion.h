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
 * @file assertion.h
 * @brief Assertion, precondition and postcondition for gas framework
 */

#ifndef _gas_assertion_
#define _gas_assertion_

#include "macro.h"
#include <iostream>
#include <cstdlib>

namespace gas {

/*!
 * @brief Function to print an assertion
 * @param exp The expression to print
 * @param file The file where the assertion failed
 * @param line The line of file where the assertion failed
 */
inline void assert_(char const * exp, char const * file, unsigned int const & line) {
	std::cerr<<file<<"@"<<line<<": assertion "<<exp<<" failed"<<std::endl;
	std::exit(1);
}

/*!
 * @brief Function to print a precondition
 * @param exp The expression to print
 * @param file The file where the precondition failed
 * @param line The line of file where the precondition failed
 */
inline void pre_(char const * exp, char const * file, unsigned int const & line) {
	std::cerr<<file<<"@"<<line<<": precondition "<<exp<<" failed"<<std::endl;
	std::exit(1);
}

/*!
 * @brief Function to print a postcondition
 * @param exp The expression to print
 * @param file The file where the postcondition failed
 * @param line The line of file where the postcondition failed
 */
inline void post_(char const * exp, char const * file, unsigned int const & line) {
	std::cerr<<file<<"@"<<line<<": precondition "<<exp<<" failed"<<std::endl;
	std::exit(1);
}

/*! @brief The class used to check a static assertion */
template <bool expression_> class compile_time_checker;

/*! @brief The class used to check a static assertion (specialization for
 *         success) */
template<> class compile_time_checker<true> { };

}

#ifdef gas_ndebug
#define gas_nassert
#define gas_npre
#define gas_npost
#endif // gas_ndebug

/*!
 * @def gas_assert(expression)
 * @brief Checking the assertion
 */

#undef gas_assert
#ifdef gas_nassert
#define gas_assert(expression) pass
#else // gas_nassert
#define gas_assert(expression) ((expression) ? pass : gas::assert_(#expression, __FILE__, __LINE__))
#endif // gas_nassert

/*!
 * @def gas_pre(expression)
 * @brief Checking the precondition
 */

#undef gas_pre
#ifndef gas_npre
#define gas_pre(expression) pass
#else // gas_npre
#define gas_pre(expression) ((expression) ? pass : gas::pre_(#expression, __FILE__, __LINE__))
#endif // gas_npre

/*!
 * @def gas_post(expression)
 * @brief Checking the postcondition
 */

#undef gas_npost
#ifdef gas_npost
#define gas_post(expression) pass
#else // gas_npost
#define gas_post(expression) ((expression) ? pass : gas::post_(#expression, __FILE__, __LINE__))
#endif // gas_npost

/*!
 * @def gas_static_assert(expression, message)
 * @brief Checking the assertion at compile time
 */

#define gas_static_assert(expression, message) { gas::compile_time_checker<expression> ERROR_##message; (void)ERROR_##message; }

#endif // _gas_assertion_
