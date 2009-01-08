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

#ifndef _GAS_ARRAY_H_ /* BEGIN _GAS_ARRAY_H_ */
#define _GAS_ARRAY_H_

/* TODO Documentation */

#include <cstring>

#include "Macro.hpp"

namespace Common {
	/* Compile-time array */
	template<typename T, size_t N=0>
	class Array {
		private:
			template<size_t M, typename S, bool C>
			struct Storage {
				S data_[M];
				Storage<M, S, C> &operator=(Storage<M, S, C> const &s) { memcpy(&data_, &(s.data_), M * sizeof(S)); return *this; }
			};
			template<size_t M, typename S>
			struct Storage<M, S, false> {
				T *data_;
				Storage() { data_ = new S[M]; }
				~Storage() { delete[] data_; }
				Storage<M, S, false> &operator=(Storage<M, S, false> const &s) { memcpy(data_, s.data_, M * sizeof(S)); return *this; }
			};

			Storage<N, T, _GAS_IS_SMALL_ARRAY_(N, T)> store_;
		public:
			inline T &operator[](size_t const i) { return store_.data_[i]; }
			inline T const &operator[](size_t const i) const { return store_.data_[i]; }
			Array<T, N> &operator=(Array<T, N> const &a) { store_ = a.store_; return *this; }
	};
	/* Run-time array */
	template<typename T>
	class Array<T, 0> {
		private:
			T *data_;
			size_t N_;
		public:
			Array(size_t const N) { N_ = N; data_ = new T[N]; }
			~Array() { delete[] data_; }
			inline T &operator[](size_t const i) { return data_[i]; }
			inline T const &operator[](size_t const i) const { return data_[i]; }
			Array<T, 0> &operator=(Array<T, 0> const &a) { memcpy(data_, a.data_, MIN(N_, a.N_) * sizeof(T)); return *this; }
	};
	/* Run-time matrix */
	template<typename T>
	class Array<Array<T, 0>, 0> {
		private:
			T *data_;
			size_t M_;
			size_t N_;
		public:
			Array(size_t const M, size_t const N) { M_ = M; N_ = N; data_ = new T[M * N]; }
			~Array() { delete[] data_; }
			inline T *operator[](size_t const i) { return &data_[i*M_]; }
			inline T const *operator[](size_t const i) const { return &data_[i*M_]; }
			Array<Array<T, 0>, 0> &operator=(Array<Array<T, 0>, 0> const &a) { size_t D_ = MIN(M_ * N_, a.M_ * a.N_); memcpy(data_, a.data_, D_ * sizeof(T)); return *this; }
	};
}

#endif /* END _GAS_ARRAY_H_ */
