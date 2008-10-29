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

#ifndef _GAS_VECTOR_H_
#define _GAS_VECTOR_H_

#include <cstdlib>
#include <iostream>

namespace LinearAlgebra {
	/**
	 * A simple class for Vector that use expression templates
	 */
	template<size_t N, typename T=double>
	class Vector {
		template <size_t M, typename S>
		friend std::ostream &operator<<(std::ostream const &stream, Vector<M, S> const &v);

		private:
			T data[N];
		public:
			Vector(T Value);
			T &operator[](size_t const Index);
			Vector<N, T> &operator=(Vector<N, T> const &v);

			static Vector<N, T> Factory(T *Array);
	};

	template<size_t N, typename T>
	Vector<N, T> &operator+(Vector<N, T> const &v, Vector<N, T> const &w);

	template<size_t N, typename T>
	Vector<N, T> &operator-(Vector<N, T> const &v, Vector<N, T> const &w);

	template<size_t N, typename T>
	T &operator*(Vector<N, T> const &v, Vector<N, T> const &w);

// 	template<size_t M, typename T, class V, class W>
// 	struct MetaDotProduct {
// 		static T RET = (V[M] * W[M]) + MetaDotProduct<M-1, T, V, W>::RET;
// 	};

	template<size_t N, typename T>
	Vector<N, T> &operator*(T const &a, Vector<N, T> const &w);
}

#endif
