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
#define _GAS_VECTOR_H_ /* BEGIN _GAS_VECTOR_H_ */

#include <cstdlib>
#include <iostream>
#include "Gas/Common/Common.h"

namespace LinearAlgebra {
	/** A simple class for Vector **/
	template<size_t N, typename T=double>
	class Vector {
		private:
			T data[N];
		public:
			Vector();
			Vector(T const);
			Vector(Vector<N, T> const &);
			T &operator[](size_t const);
			T const &operator[](size_t const) const;
			bool operator==(T const &);
			bool operator==(Vector<N, T> const &);
			Vector<N, T> &operator=(T const &);
			Vector<N, T> &operator=(Vector<N, T> const &);
			Vector<N, T> &operator+=(T const &);
			Vector<N, T> &operator+=(Vector<N, T> const &);
			Vector<N, T> &operator-=(T const &);
			Vector<N, T> &operator-=(Vector<N, T> const &);
			Vector<N, T> &operator*=(T const &);
			Vector<N, T> &operator/=(T const &);

			template<size_t M, typename S> friend Vector<M, S> &operator+(Vector<M, S>, Vector<M, S> const &);
			template<size_t M, typename S> friend Vector<M, S> &operator+(Vector<M, S>, S const &);
			template<size_t M, typename S> friend Vector<M, S> &operator+(S const &, Vector<M, S>);
			template<size_t M, typename S> friend Vector<M, S> &operator-(Vector<M, S>, Vector<M, S> const &);
			template<size_t M, typename S> friend Vector<M, S> &operator-(Vector<M, S>, S const &);
			template<size_t M, typename S> friend Vector<M, S> &operator-(S const &, Vector<M, S>);
			template<size_t M, typename S> friend Vector<M, S> &operator*(S const &, Vector<M, S>);
			template<size_t M, typename S> friend Vector<M, S> &operator*(Vector<M, S>, S const &);
			template<size_t M, typename S> friend Vector<M, S> &operator/(Vector<M, S>, S const &);
			template<size_t M, typename S> friend S &operator*(Vector<M, S> const &, Vector<M, S> const &);

			template<size_t M, typename S> friend std::ostream &operator<<(std::ostream &, Vector<M, S> const &);
	};

	/** The default constructor **/
	template<size_t N, typename T>
	Vector<N, T>::Vector() {
	}

	/** The constructor to initialize the entire vector with the same value
	 *  @param Value The vaule to use **/
	template<size_t N, typename T>
	Vector<N, T>::Vector(T const Value) {
		range(i, 0, N) data[i] = Value;
	}

	/** The constructor to copy a vector
	 *  @param v The vector to copy **/
	template<size_t N, typename T>
	Vector<N, T>::Vector(Vector<N, T> const &v) {
		range(i, 0, N) data[i] = v.data[i];
	}

	/** The operator [] to access to the components of vector
	 *  @param i The index of component **/
	template<size_t N, typename T>
	T &Vector<N, T>::operator[](size_t const i) {
		return data[i];
	}

	/** The operator [] to access to the components of vector
	 *  @param i The index of component **/
	template<size_t N, typename T>
	T const &Vector<N, T>::operator[](size_t const i) const {
		return data[i];
	}

	/** The operator == to compare all components of a vector with a value
	 *  @param v The value to compare **/
	template<size_t N, typename T>
	bool Vector<N, T>::operator==(T const &v) {
		range(i, 0, N) { if (data[i] != v) return 0; }
		return 1;
	}

	/** The operator == to compare two vectors
	 *  @param v The second vector **/
	template<size_t N, typename T>
	bool Vector<N, T>::operator==(Vector<N, T> const &v) {
		if (this != &v) {
			range(i, 0, N) { if (data[i] != v.data[i]) return 0; }
		}
		return 1;
	}

	/** The operator = to copy a value in all components of a vector
	 *  @param a The scalar value to copy **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(T const &a) {
		range(i, 0, N) data[i] = a;
		return *this;
	}

	/** The operator = to copy a vector in a vector
	 *  @param v The vector to copy **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(Vector<N, T> const &v) {
		if (this != &v) { range(i, 0, N) data[i] = v.data[i]; }
		return *this;
	}

	/** The operator += to add a value to all components
	 *  @param a The scalar value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator+=(T const &a) {
		if (a != 0) { range(i, 0, N) data[i] += a; }
		return *this;
	}

	/** The operator += to add a vector
	 *  @param v The vector to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator+=(Vector<N, T> const &v) {
		range(i, 0, N) data[i] += v.data[i];
		return *this;
	}

	/** The operator -= to subtract a value to all components
	 *  @param a The value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator-=(T const &a) {
		if (a != 0) { range(i, 0, N) data[i] -= a; }
		return *this;
	}

	/** The operator -= to subtract a vector
	 *  @param v The vector to subtract **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator-=(Vector<N, T> const &v) {
		range(i, 0, N) data[i] -= v.data[i];
		return *this;
	}

	/** The operator *= to multiply a vector by a scalar value
	 *  @param a The scalar value to multiply **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator*=(T const &a) {
		range(i, 0, N) data[i] *= a;
		return *this;
	}

	/** The operator *= to divide a vector by a scalar value
	 *  @param a The scalar value to divide **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator/=(T const &a) {
		range(i, 0, N) data[i] /= a;
		return *this;
	}

	/** The operator + to add two vectors
	 *  @param v The first vector
	 *  @param w The second vector **/
	template<size_t M, typename S>
	Vector<M, S> &operator+(Vector<M, S> v, Vector<M, S> const &w) {
		return v += w;
	}

	/** The operator + to add a vector and a scalar value
	 *  @param v The vector
	 *  @param a The scalar value **/
	template<size_t M, typename S>
	Vector<M, S> &operator+(Vector<M, S> v, S &a) {
		return v += a;
	}

	/** The operator + to add a vector and a scalar value
	 *  @param a The value
	 *  @param v The vector **/
	template<size_t M, typename S>
	Vector<M, S> &operator+(S &a, Vector<M, S> v) {
		return v += a;
	}

	/** The operator - to subtract two vectors
	 *  @param v The first vector
	 *  @param w The second vector **/
	template<size_t M, typename S>
	Vector<M, S> &operator-(Vector<M, S> v, Vector<M, S> const &w) {
		return v += w;
	}

	/** The operator - to subtract a vector and a scalar value
	 *  @param v The vector
	 *  @param a The scalar value **/
	template<size_t M, typename S>
	Vector<M, S> &operator-(Vector<M, S> v, S &a) {
		return v -= a;
	}

	/** The operator - to subtract a vector and a scalar value
	 *  @param a The scalar value
	 *  @param v The vector **/
	template<size_t M, typename S>
	Vector<M, S> &operator-(S &a, Vector<M, S> v) {
		return (-1) * (v - a);
	}

	/** Moltiplication by a scalar value
	 *  @param a The scalar value
	 *  @param v The vector **/
	template<size_t M, typename S> 
	Vector<M, S> &operator*(S const &a, Vector<M, S> v) {
		return v *= a;
	}

	/** Moltiplication by a scalar value
	 *  @param v The vector
	 *  @param a The scalar value**/
	template<size_t M, typename S> 
	Vector<M, S> &operator*(Vector<M, S> v, S const &a) {
		return v *= a;
	}

	/** Division by a scalar value
	 *  @param v The vector
	 *  @param a The scalar value**/
	template<size_t M, typename S> 
	Vector<M, S> &operator/(Vector<M, S> v, S const &a) {
		return v /= a;
	}

	/** The scalar product between two vector
	 *  @param v The first vector
	 *  @param w The second vector **/
	template<size_t M, typename S>
	S &operator*(Vector<M, S> const &v, Vector<M, S> const &w) {
		S x = 0;
		range(i, 0, M) x += (v[i] * w[i]);
		return x;
	}

	/** The operator << to print a vector in a stream
	 *  @param stream The stream
	 *  @param v The vector **/
	template<size_t M, typename S>
	std::ostream &operator<<(std::ostream &stream, Vector<M, S> const &v) {
		stream << "[";
		if (M <= 10 )
			range(i, 0, M) stream << " " << v.data[i];
		else
			stream << " " << v.data[0] << " " << v.data[1] << " ... " << v.data[M-2] << " " << v.data[M-1];
		stream << " ]";
	}

}

#endif /* END _GAS_VECTOR_H_ */
