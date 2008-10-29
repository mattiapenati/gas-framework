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

#include "Matrix.hpp"
#include "Gas/Common/Common.h"
#include <cblas.h>

namespace LinearAlgebra {
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T>::Matrix(T Value) {
		range(i, 0, M) range(j, 0, N) data[i][j] = Value;
	}

	template<size_t M, size_t N, typename T>
	Vector<N, T> &Matrix<M, N, T>::operator[](size_t const Index) {
		return Vector<N, T>::Factory(data[Index]);
	}

	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &operator+(Matrix<M, N, T> const &A, Matrix<M, N, T> const &B){
		Matrix<M, N, T> C;
		range(i, 0, M) range(j, 0, N) C.data[i][j] = A.data[i][j]+B.data[i][j];
		return C;
	}

	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &operator-(Matrix<M, N, T> const &A, Matrix<M, N, T> const &B){
		Matrix<M, N, T> C;
		range(i, 0, M) range(j, 0, N) C.data[i][j] = A.data[i][j]-B.data[i][j];
		return C;
	}

	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &operator*(T const &a, Matrix<M, N, T> const &A){
		Matrix<M, N, T> C;
		range(i, 0, M) range(j, 0, N) C.data[i][j] = a*A.data[i][j];
		return C;
	}

	template<size_t M, size_t N, typename T>
	Vector<M, double> &operator*(Matrix<M, N, double> const &A, Vector<N, double> const &v){
		Vector<M, double> y;
		cblas_dgemv(101, 111, M, N, 1, A, N, v, 1, 0, y, 1);
		return y;
	}

	template<size_t M, size_t N, typename T>
	Vector<M, float> &operator*(Matrix<M, N, float> const &A, Vector<N, float> const &v){
		Vector<M, float> y;
		cblas_sgemv(101, 111, M, N, 1, A, N, v, 1, 0, y, 1);
		return y;
	}

	template<size_t M, size_t N, size_t K, typename T>
	Matrix<M, N, double> &operator*(Matrix<M, K, double> const &A, Matrix<K, N, double> const &B){
		Matrix<M, N, double> C;
		cblas_dgemm(101, 111, 111, M, K, N, 1, A, K, B, N, 1, C, N);
		return C;	
	}	

	template<size_t M, size_t N, size_t K, typename T>
	Matrix<M, N, float> &operator*(Matrix<M, K, float> const &A, Matrix<K, N, float> const &B){
		Matrix<M, N, float> C;
		cblas_sgemm(101, 111, 111, M, K, N, 1, A, K, B, N, 1, C, N);
		return C;	
	}
}
