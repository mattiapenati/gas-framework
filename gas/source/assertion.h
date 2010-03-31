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

#ifndef GAS_ASSERTION_H
#define GAS_ASSERTION_H

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
inline void m_assert(char const * exp, char const * file, unsigned int const & line) {
	std::cerr<<file<<"@"<<line<<": assertion "<<exp<<" failed"<<std::endl;
	std::exit(1);
}

/*!
 * @brief Function to print a precondition
 * @param exp The expression to print
 * @param file The file where the precondition failed
 * @param line The line of file where the precondition failed
 */
inline void m_pre(char const * exp, char const * file, unsigned int const & line) {
	std::cerr<<file<<"@"<<line<<": precondition "<<exp<<" failed"<<std::endl;
	std::exit(1);
}

/*!
 * @brief Function to print a postcondition
 * @param exp The expression to print
 * @param file The file where the postcondition failed
 * @param line The line of file where the postcondition failed
 */
inline void m_post(char const * exp, char const * file, unsigned int const & line) {
	std::cerr<<file<<"@"<<line<<": precondition "<<exp<<" failed"<<std::endl;
	std::exit(1);
}

/*! @brief The class used to check a static assertion */
template <bool expression_> class compile_time_checker;

/*! @brief The class used to check a static assertion (specialization for
 *         success) */
template<> class compile_time_checker<true> { };

}

#ifdef NDEBUG
#define GAS_NDEBUG
#endif

#ifdef GAS_NDEBUG
#define GAS_NASSERT
#define GAS_NPRE
#define GAS_NPOST
#endif // GAS_NDEBUG

/*!
 * @def GAS_ASSERT(expression)
 * @brief Checking the assertion
 */

#ifdef GAS_ASSERT
#undef GAS_ASSERT
#endif // GAS_ASSERT

#ifdef GAS_NASSERT
#define GAS_ASSERT(expression) GAS_PASS
#else // GAS_NASSERT
#define GAS_ASSERT(expression) ((expression) ? GAS_PASS : gas::m_assert(#expression, __FILE__, __LINE__))
#endif // GAS_NASSERT

#define gas_assert GAS_ASSERT

/*!
 * @def GAS_PRE(expression)
 * @brief Checking the precondition
 */

#ifdef GAS_PRE
#undef GAS_PRE
#endif // GAS_PRE

#ifndef GAS_NPRE
#define GAS_PRE(expression) GAS_PASS
#else // GAS_NPRE
#define GAS_PRE(expression) ((expression) ? GAS_PASS : gas::m_pre(#expression, __FILE__, __LINE__))
#endif // GAS_NPRE

#define gas_pre GAS_PRE

/*!
 * @def GAS_POST(expression)
 * @brief Checking the postcondition
 */

#ifdef GAS_POST
#undef GAS_POST
#endif // GAS_POST

#ifdef GAS_NPOST
#define GAS_POST(expression) pass
#else // GAS_NPOST
#define GAS_POST(expression) ((expression) ? pass : gas::m_post(#expression, __FILE__, __LINE__))
#endif // GAS_NPOST

#define gas_post GAS_POST

/*!
 * @def GAS_STATIC_ASSERT(expression, message)
 * @brief Checking the assertion at compile time
 */

#ifdef GAS_STATIC_ASSERT
#undef GAS_STATIC_ASSERT
#endif // GAS_STATIC_ASSERT

#define GAS_STATIC_ASSERT(expression, message) { gas::compile_time_checker<expression> ERROR_##message; (void)ERROR_##message; }

#define gas_static_assert GAS_STATIC_ASSERT

#endif // GAS_ASSERTION_H
