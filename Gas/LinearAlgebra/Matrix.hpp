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

#ifndef _GAS_MATRIX_H_
#define _GAS_MATRIX_H_

#include <cstdlib>
extern "C" {
#include <cblas.h>
}
#include "Vector.hpp"
#include "Gas/Common/Common.h"

namespace LinearAlgebra {
	template<size_t M, size_t N, typename T=double>
	class Matrix {
		private:
			T data[M][N];
		public:
			Matrix();
			Matrix(T const);
			T *operator[](size_t const);
			bool operator==(T const &);
			bool operator==(Matrix<M, N, T> const &);
			Matrix<M, N, T> &operator=(T const &);
			Matrix<M, N, T> &operator=(Matrix<M, N, T> const &);
			Matrix<M, N, T> &operator+=(T const &);
			Matrix<M, N, T> &operator+=(Matrix<M, N, T> const &);
			Matrix<M, N, T> &operator-=(T const &);
			Matrix<M, N, T> &operator-=(Matrix<M, N, T> const &);
			Matrix<M, N, T> &operator*=(T const &);
			Matrix<M, N, T> &operator/=(T const &);	
		
			template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator+(Matrix<P, Q, S>, Matrix<P, Q, S> const &);
	//		template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator+(Vector<P, Q, S>, S const &);		
	//		template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator+(S const &, Vector<M, S>);
			template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator-(Matrix<P, Q, S>, Matrix<P, Q, S> const &);
	//		template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator-(Vector<M, S>, S const &);
	//		template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator-(S const &, Vector<M, S>);
			template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator*(S const &, Matrix<P, Q, S>);
	//		template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator*(Vector<M, S>, S const &);
			template<size_t P, size_t Q, typename S> friend Matrix<P, Q, S> &operator/(Matrix<P, Q, S>, S const &);
			template<size_t P, size_t Q, typename S> friend Vector<P, S> operator*(Matrix<P, Q, S> const &, Vector<Q, S> const &);
			template<size_t P, size_t Q, size_t K, typename S> friend Matrix<P, Q, S> operator*(Matrix<P, K, S> const &, Matrix<K, Q, S> const &);
	//		template<size_t P, size_t Q, typename S> friend std::ostream &operator<<(std::ostream &, Vector<M, S> const &);
	};

	/** The default constructor **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T>::Matrix() {
	}

	/** The constructor to initialize the entire matrix with the same value
	 *  @param Value The vaule to use **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T>::Matrix(T const Value) {
		range(i, 0, N) range(j, 0, M) data[i][j] = Value;
	}

	/** The operator [] to access to the components of vector
	 *  @param i The index of component **/
	template<size_t M, size_t N, typename T>
	T *Matrix<M, N, T>::operator[](size_t const Index) {
		return &(data[Index]);
	}

	/** The operator == to compare all components of a matrix with a value
	 *  @param v The value to compare **/
	template<size_t M, size_t N, typename T>
	bool Matrix<M, N, T>::operator==(T const &v) {
		range(i, 0, M) range(j, 0, N){ if (data[i][j] != v) return 0; }
		return 1;
	}

	/** The operator == to compare two matrix
	 *  @param v The second vector **/
	template<size_t M, size_t N, typename T>
	bool Matrix<M, N, T>::operator==(Matrix<M, N, T> const &A) {
		if (this != &A) {
			range(i, 0, M) range(j, 0, N){ if (data[i][j] != A.data[i][j]) return 0; }
		}
		return 1;
	}

	/** The operator = to copy a value in all matrix
	 *  @param a The value to copy **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &Matrix<M, N, T>::operator=(T const &a) {
		range(i, 0, M) range(j, 0, N) data[i][j] = a;
		return *this;
	}

	/** The operator = to copy a matrix in a matrix
	 *  @param v The matrix to copy **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &Matrix<M, N, T>::operator=(Matrix<M, N, T> const &A) {
		if (this != &A) { range(i, 0, M) range(j, 0, N) data[i][j] = A.data[i][j]; }
		return *this;
	}

	/** The operator += to add a value to all components
	 *  @param a The value to add **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &Matrix<M, N, T>::operator+=(T const &a) {
		if (a != 0) { range(i, 0, M) range(j, 0, N)data[i][j] += a; }
		return *this;
	}

	/** The operator += to add a matrix
	 **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &Matrix<M, N, T>::operator+=(Matrix<M, N, T> const &A) {
		range(i, 0, M) range(j, 0, M) data[i][j] += A.data[i][j];
		return *this;
	}

	/** The operator -= to subtract a value to all components
	 *  @param a The value to add **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &Matrix<M, N, T>::operator-=(T const &a) {
		if (a != 0) { range(i, 0, M) range(j, 0, M) data[i][j] -= a; }
		return *this;
	}

	/** The operator -= to subtract a matrix
	 **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &Matrix<M, N, T>::operator-=(Matrix<M, N, T> const &A) {
		range(i, 0, M) range(j, 0, N) data[i][j] -= A.data[i][j];
		return *this;
	}

	/** The operator *= is the scalar product 
	 *  @param a The value to multiply **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &Matrix<M, N, T>::operator*=(T const &a) {
		range(i, 0, M) range(j, 0, N) data[i][j] *= a;
		return *this;
	}

	/** The operator /= to divide a matrix
	 *  @param a The value to divide **/
	template<size_t M, size_t N, typename T>
	Matrix<M, N, T> &Matrix<M, N, T>::operator/=(T const &a) {
		range(i, 0, M) range(j, 0, N) data[i][j] /= a;
		return *this;
	}
	/** The operator + to add two matrices
	 *  @param A The first matrix
	 *  @param B The second matrix **/
	template<size_t P, size_t Q, typename S>
	Matrix<P, Q, S> &operator+(Matrix<P, Q, S> A, Matrix<P, Q, S> const &B) {
		return A += B;
	}

	/** The operator - to subtract two matrices
	 *  @param A The first matrix
	 *  @param B The second matrix **/
	template<size_t P, size_t Q, typename S>
	Matrix<P, Q, S> &operator-(Matrix<P, Q, S> A, Matrix<P, Q, S> const &B) {
		return A -= B;
	}

	/** Moltiplication by a scalar
	 *  @param a The scalar
	 *  @param A The matrix **/
	template<size_t P, size_t Q, typename S> 
	Matrix<P, Q, S> &operator*(S const &a, Matrix<P, Q, S> A) {
		return A *= a;
	}

	/** Division by a scalar
	 *  @param A The matrix
	 *  @param a The scalar **/
	template<size_t P, size_t Q, typename S> 
	Matrix<P, Q, S> &operator/(Matrix<P, Q, S> A, S const &a) {
		return A /= a;
	}
	
	/** Product Matrix Vector
	 * @param A The matrix
	 * @param v The vector **/
	template<size_t P, size_t Q>
	Vector<P, double> operator*(Matrix<P, Q, double> const &A, Vector<Q, double> const &v){
		Vector<P, double> y;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, P, Q, 1, A.data, Q, v.data, 1, 0, y.data, 1);
		return y;
	}

	template<size_t P, size_t Q>
	Vector<P, float> operator*(Matrix<P, Q, float> const &A, Vector<Q, float> const &v){
		Vector<P, float> y;
		cblas_sgemv(CblasRowMajor, CblasNoTrans, P, Q, 1, A.data, Q, v.data, 1, 0, y.data, 1);
		return y;
	}

	template<size_t P, size_t Q, size_t K>
	Matrix<P, Q, double> operator*(Matrix<P, K, double> const &A, Matrix<K, Q, double> const &B){
		Matrix<P, Q, double> C;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, P, K, Q, 1, A.data, K, B.data, Q, 1, C.data, Q);
		return C;	
	}	

	template<size_t P, size_t Q, size_t K>
	Matrix<P, Q, float> operator*(Matrix<P, K, float> const &A, Matrix<K, Q, float> const &B){
		Matrix<P, Q, float> C;
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, P, K, Q, 1, A.data, K, B.data, Q, 1, C.data, Q);
		return C;	
	}
}

#endif

