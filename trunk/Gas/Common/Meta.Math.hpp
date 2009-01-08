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

#ifndef _GAS_META_MATH_H_ /* BEGIN _GAS_META_MATH_H_ */
#define _GAS_META_MATH_H_

#include "Macro.hpp"
#include "Math.hpp"

namespace Common { namespace Meta { namespace Math {
	/* Identity */ GAS_DEFINE_META_MATH_FUNCTION(Id)

	/* Four basic operations */
	GAS_DEFINE_META_MATH_BINARY(Sum)
	GAS_DEFINE_META_MATH_BINARY(Sub)
	GAS_DEFINE_META_MATH_BINARY(Mul)
	GAS_DEFINE_META_MATH_BINARY(Div)

	/* Division mod */
	GAS_DEFINE_META_MATH_BINARY(Mod)

	/* Max  */ GAS_DEFINE_META_MATH_BINARY(Max)
	/* Min  */ GAS_DEFINE_META_MATH_BINARY(Min)

	/* Abs  */ GAS_DEFINE_META_MATH_FUNCTION(Abs)
	/* Conj */ GAS_DEFINE_META_MATH_FUNCTION(Conj)
	/* Sqrt */ GAS_DEFINE_META_MATH_FUNCTION(Sqrt)
}}}

#endif /* END _GAS_META_MATH_H_ */
