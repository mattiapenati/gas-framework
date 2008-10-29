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

#include "Vector.hpp"
#include "Gas/Common/Common.h"

namespace LinearAlgebra{
	template <size_t M, typename S>
	std::ostream &operator<<(std::ostream const &stream, Vector<M, S> const &v) {
		stream << "[";
		range(i, 0, M) stream << " " << v.data[i];
		stream << " ]";
	}

	/**
	 * A constructor that set all values of the vector
	 */
	template<size_t N, typename T>
	Vector<N, T>::Vector(T Value) {
		range(i, 0, N) data[i] = Value;
	}
	
	/**
	 * Overloading of [] operator, to access to an element
	 */
	template<size_t N, typename T>
	T &Vector<N, T>::operator[](size_t const Index) {
		return data[Index];
	}

	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(Vector<N, T> const &v) {
		range(i, 0, N) data[i] = v[i];
		return *this;
	}

	template<size_t N, typename T>
	Vector<N, T> &operator+(Vector<N, T> const &v, Vector<N, T> const &w){
		Vector<N, T> k;
		range(i, 0, N) k[i] = v[i] + w[i];
		return k;
	}

	template<size_t N, typename T>
	Vector<N, T> &operator-(Vector<N, T> const &v, Vector<N, T> const &w){
		Vector<N, T> k;
		range(i, 0, N) k[i] = v[i] - w[i];
		return k;
	}

	template<size_t M, size_t N, typename T>
	struct MetaDot {
		static T RET(Vector<N, T> &v, Vector<N, T> &w) {
			return (v[M] * w[M]) + MetaDot<M-1, N, T>::RET(v, w);
		}
	};
	template<size_t N, typename T>
	struct MetaDot<0, N, T> {
		static T RET(Vector<N, T> &v, Vector<N, T> &w) {
			return (v[0] * w[0]);
		}
	};

	/**
	 * Scalar product between two vectors
	 */
	template<size_t N, typename T>
	T operator*(Vector<N, T> const &v, Vector<N, T> const &w){
		#ifdef _USE_META_
			return MetaDot<N, N, T>::RET(v, w);
		#else
			T x = 0;
			range(i, 0, N) x += v[i] * w[i];
			return x;
		#endif
	}

	template<size_t N, typename T>
	Vector<N, T> operator*(T const &a, Vector<N, T> const &w){
		Vector<N, T> v;
		range(i, 0, N) v[i] = a * w[i];
		return v;
	}

	template<size_t N, typename T>
	bool operator==(Vector<N, T> const &v, Vector<N, T> const &w) {
		range(i, 0, N) { if (v[i] != w[i]) return 0;}
		return 1;
	}
}
