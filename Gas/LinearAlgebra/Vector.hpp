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

#include <cstddef>
#include <cassert>
#include <iostream>
#include "Gas/Common/Common.h"
#include "Interface.hpp"

namespace LinearAlgebra {
	/** A simple class for vector
	 *  @class Vector 
	 *  @brief A simple class for vector **/
	template<size_t N, typename T=double>
	class Vector: virtual public Interface::Vector<N, T> {
		private:
			T Data_[N];

		public:
			Vector();
			Vector(T const);
			Vector(Vector<N, T> const &);
			~Vector();
			inline T &operator()(size_t const);
			inline T const &operator()(size_t const) const;
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
	 *  @param a The vaule to use **/
	template<size_t N, typename T>
	Vector<N, T>::Vector(T const a) {
		range(i, 0, N) Data_[i] = a;
	}

	/** The constructor to copy a vector
	 *  @param V The vector to copy **/
	template<size_t N, typename T>
	Vector<N, T>::Vector(Vector<N, T> const &V) {
		range(i, 0, N) Data_[i] = V.Data_[i];
	}

	/** The default destructor **/
	template<size_t N, typename T>
	Vector<N, T>::~Vector() {
	}

	/** The operator () to access to the components of vector
	 *  @param i The index of component **/
	template<size_t N, typename T>
	T &Vector<N, T>::operator()(size_t const i) {
		assert(i < N);
		return Data_[i];
	}

	/** The operator () to access to the components of vector
	 *  @param i The index of component **/
	template<size_t N, typename T>
	T const &Vector<N, T>::operator()(size_t const i) const {
		assert(i < N);
		return Data_[i];
	}

	/** The operator == to compare all components of a vector with a value
	 *  @param a The value to compare **/
	template<size_t N, typename T>
	bool Vector<N, T>::operator==(T const &a) {
		range(i, 0, N) { if (Common::Math::Abs(Data_[i] - a) > 8 * Common::Limits<T>::Epsilon) return 0; }
		return 1;
	}

	/** The operator == to compare two vectors
	 *  @param V The second vector **/
	template<size_t N, typename T>
	bool Vector<N, T>::operator==(Vector<N, T> const &V) {
		if (this != &V) {
			range(i, 0, N) { if (Common::Math::Abs(Data_[i] - V.Data_[i]) > 8 * Common::Limits<T>::Epsilon) return 0; }
		}
		return 1;
	}

	/** The operator = to copy a value in all components of a vector
	 *  @param a The scalar value to copy **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(T const &a) {
		range(i, 0, N) Data_[i] = a;
		return *this;
	}

	/** The operator = to copy a vector in a vector
	 *  @param V The vector to copy **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(Vector<N, T> const &V) {
		if (this != &V) { range(i, 0, N) Data_[i] = V.Data_[i]; }
		return *this;
	}

	/** The operator += to add a value to all components
	 *  @param a The scalar value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator+=(T const &a) {
		if (a != 0) { range(i, 0, N) Data_[i] += a; }
		return *this;
	}

	/** The operator += to add a vector
	 *  @param V The vector to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator+=(Vector<N, T> const &V) {
		range(i, 0, N) Data_[i] += V.Data_[i];
		return *this;
	}

	/** The operator -= to subtract a value to all components
	 *  @param a The value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator-=(T const &a) {
		if (a != 0) { range(i, 0, N) Data_[i] -= a; }
		return *this;
	}

	/** The operator -= to subtract a vector
	 *  @param V The vector to subtract **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator-=(Vector<N, T> const &V) {
		range(i, 0, N) Data_[i] -= V.Data_[i];
		return *this;
	}

	/** The operator *= to multiply a vector by a scalar value
	 *  @param a The scalar value to multiply **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator*=(T const &a) {
		if (a != 1) { range(i, 0, N) Data_[i] *= a; }
		return *this;
	}

	/** The operator /= to divide a vector by a scalar value
	 *  @param a The scalar value to divide **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator/=(T const &a) {
		if ((a != 1) and (a != 0)) { range(i, 0, N) Data_[i] /= a; }
		return *this;
	}

	/** The operator + to add two vectors
	 *  @param V The first vector
	 *  @param W The second vector **/
	template<size_t M, typename S>
	Vector<M, S> &operator+(Vector<M, S> V, Vector<M, S> const &W) {
		return V += W;
	}

	/** The operator + to add a vector and a scalar value
	 *  @param V The vector
	 *  @param a The scalar value **/
	template<size_t M, typename S>
	Vector<M, S> &operator+(Vector<M, S> V, S &a) {
		return V += a;
	}

	/** The operator + to add a vector and a scalar value
	 *  @param a The value
	 *  @param V The vector **/
	template<size_t M, typename S>
	Vector<M, S> &operator+(S &a, Vector<M, S> V) {
		return V += a;
	}

	/** The operator - to subtract two vectors
	 *  @param V The first vector
	 *  @param W The second vector **/
	template<size_t M, typename S>
	Vector<M, S> &operator-(Vector<M, S> V, Vector<M, S> const &W) {
		return V += W;
	}

	/** The operator - to subtract a vector and a scalar value
	 *  @param V The vector
	 *  @param a The scalar value **/
	template<size_t M, typename S>
	Vector<M, S> &operator-(Vector<M, S> V, S &a) {
		return V -= a;
	}

	/** The operator - to subtract a vector and a scalar value
	 *  @param a The scalar value
	 *  @param V The vector **/
	template<size_t M, typename S>
	Vector<M, S> &operator-(S &a, Vector<M, S> V) {
		return (-1) * (V - a);
	}

	/** Moltiplication by a scalar value
	 *  @param a The scalar value
	 *  @param V The vector **/
	template<size_t M, typename S> 
	Vector<M, S> &operator*(S const &a, Vector<M, S> V) {
		return V *= a;
	}

	/** Moltiplication by a scalar value
	 *  @param V The vector
	 *  @param a The scalar value**/
	template<size_t M, typename S> 
	Vector<M, S> &operator*(Vector<M, S> V, S const &a) {
		return V *= a;
	}

	/** Division by a scalar value
	 *  @param V The vector
	 *  @param a The scalar value**/
	template<size_t M, typename S> 
	Vector<M, S> &operator/(Vector<M, S> V, S const &a) {
		return V /= a;
	}

	/** The scalar product between two vector
	 *  @param V The first vector
	 *  @param W The second vector **/
	template<size_t M, typename S>
	S &operator*(Vector<M, S> const &V, Vector<M, S> const &W) {
		S x = 0;
		range(i, 0, M) x += (V.Data_[i] * Common::Math::Conj(W.Data_[i]));
		return x;
	}

	/** The operator << to print a vector in a stream
	 *  @param s The stream
	 *  @param V The vector **/
	template<size_t M, typename S>
	std::ostream &operator<<(std::ostream &s, Vector<M, S> const &V) {
		s << "[";
		if (M <= 10 )
			range(i, 0, M) s << " " << V.Data_[i];
		else
			s << " " << V.Data_[0] << " " << V.Data_[1] << " ... " << V.Data_[M-2] << " " << V.Data_[M-1];
		s << " ]";
	}

}

#endif /* END _GAS_VECTOR_H_ */
