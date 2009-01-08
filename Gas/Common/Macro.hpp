/*
 * Copyright (c) 2008, Davide Ferrarese & Mattia Penati
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
 * 3. Neither the name of the <ORGANIZATION> nor the names of its
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

#ifndef _GAS_DEFINE_H_
#define _GAS_DEFINE_H_

/*
 * General definition
 */

#undef MIN
#define MIN(A, B) ((A < B) ? A : B)
#undef MAX
#define MAX(A, B) ((A < B) ? B : A)
#undef NULL
#define NULL 0

/*
 * Utility for Array.hpp
 */

/* The size of a small array */
#undef _GAS_SMALL_ARRAY_
#define _GAS_SMALL_ARRAY_ 1024
/* The macro to check if an array is small or big */
#undef _GAS_IS_SMALL_ARRAY_
#define _GAS_IS_SMALL_ARRAY_(N, T) ((N * sizeof(T))<=_GAS_SMALL_ARRAY_)

/*
 * Utility for Meta.Math
 */
/* These macro is used to define a mathematical function to use with meta
 * programming, it is simple, you must define the function in the namespace
 * Common::Math then use the macro. The second version is for binary operator,
 * like plus or minus. All method are inlined for better performance.
 * You can access to each function by the default method <FUNCTION_NAME>::RET(x) 
 * or <BINARY_OPERATOR_NAME>::RET(a, b).
 */
#define GAS_DEFINE_META_MATH_FUNCTION(f) \
template<typename T> \
struct f { \
	static inline T const RET(T const &x){ \
		return Common::Math::f<T>(x); \
	} \
};

#define GAS_DEFINE_META_MATH_BINARY(f) \
template<typename T> \
struct f { \
	static inline T const RET(T const &a, T const &b){ \
		return Common::Math::f<T>(a, b); \
	} \
};

#endif
