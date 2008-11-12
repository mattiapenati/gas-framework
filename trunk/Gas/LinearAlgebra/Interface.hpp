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

#ifndef _GAS_LINEARALGEBRA_INTERFACE_H_
#define _GAS_LINEARALGEBRA_INTERFACE_H_

#include <cstdlib>

namespace LinearAlgebra { namespace Interface {
	/** A vector interface for the conventional notation v(i) 
	 * 	@class Vector 
	 *  @brief An interface to define the common vector method **/
	template<size_t N, typename T>
	class Vector {
		public:
			virtual inline T &operator()(size_t const) = 0;
			virtual inline T const &operator()(size_t const) const = 0;

			inline size_t Size() { return N; };
	};

	/** A matrix interface for the conventional notation A(i,j) 
	 *  @class Matrix
	 *  @brief An interface to define the common metrix method **/
	template<size_t M, size_t N, typename T>
	class Matrix {
		public:
			virtual inline T &operator()(size_t const, size_t const) = 0;
			virtual inline T const &operator()(size_t const, size_t const) const = 0;

			inline size_t Rows() { return M; };
			inline size_t Cols() { return N; };
	};
}}

#endif
