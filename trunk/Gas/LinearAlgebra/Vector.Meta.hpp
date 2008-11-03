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
	namespace Meta {
		template<size_t N, typename T>
		class Expr {
			public:
				virtual T &operator[](size_t const i) = 0;
				virtual T const &operator[](size_t const i) const = 0;
		};

		template<size_t N, typename T>
		class Const: virtual public Expr<N, T> {
			private:
				T a;
			public:
				Const(T const &A) { a = A; };
				T &operator[](size_t const i) { return a; };
				T const &operator[](size_t const i) const { return a; };
		};

		template<size_t N, typename T, T F(T, T)>
		class BinExpr: virtual public Expr<N, T> {
			private:
				Expr<N, T> *l;
				Expr<N, T> *r;
			public:
				BinExpr(Expr<N, T> const &L, Expr<N, T> const &R) { l = &L; r = &R; };
				BinExpr(Expr<N, T> const *L, Expr<N, T> const *R) { l = L; r = R; };
				T &operator[](size_t const i) { return F((*l)[i], (*r)[i]); };
				T const operator[](size_t const i) const { return F((*l)[i], (*r)[i]); };
		};

		template<size_t M, size_t N, typename T>
		struct Dot {
			static T &RET(Expr<N, T> &l, Expr<N, T> &r) {
				return (l[M-1] * r[M-1]) + Dot<M-1, N, T>::RET(l, r);
			}
		};
		template<size_t N, typename T>
		struct Dot<1, N, T> {
			static T &RET(Expr<N, T> &l, Expr<N, T> &r) {
				return l[0] * r[0];
			}
		};
	}

	/** A simple class for Vector **/
	template<size_t N, typename T=double>
	class Vector: virtual public Meta::Expr<N, T> {
		private:
			T data[N];
		public:
			Vector();
			Vector(T const);
			Vector(Meta::Expr<N, T> const &);
			T &operator[](size_t const);
			T const &operator[](size_t const) const;
			bool operator==(T const &);
			bool operator==(Meta::Expr<N, T>const &);
			Vector<N, T> &operator=(T const &);
			Vector<N, T> &operator=(Meta::Expr<N, T> const &);
			Vector<N, T> &operator+=(T const &);
			Vector<N, T> &operator+=(Meta::Expr<N, T> const &);
			Vector<N, T> &operator-=(T const &);
			Vector<N, T> &operator-=(Meta::Expr<N, T> const &);
			Vector<N, T> &operator*=(T const &);
			Vector<N, T> &operator/=(T const &);

			template<size_t M, typename S> friend bool operator==(Meta::Expr<M, S> const &, Meta::Expr<M, S> const &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Sum> &operator+(Meta::Expr<M, S> const &, Meta::Expr<M, S> const &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Sum> &operator+(Meta::Expr<M, S> const &, S const &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Sum> &operator+(S const &, Meta::Expr<M, S> const &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Sub> &operator-(Meta::Expr<M, S> const &, Meta::Expr<M, S> const &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Sub> &operator-(Meta::Expr<M, S> const &, S const &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Sub> &operator-(S const &, Meta::Expr<M, S> const &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Mul> &operator*(S const &, Meta::Expr<M, S> &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Mul> &operator*(Meta::Expr<M, S> &, S const &);
			template<size_t M, typename S> friend Meta::BinExpr<M, S, Common::Meta::Div> &operator/(Meta::Expr<M, S> &, S const &);
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
	Vector<N, T>::Vector(Meta::Expr<N, T> const &e) {
		range(i, 0, N) data[i] = e[i];
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
	bool Vector<N, T>::operator==(Meta::Expr<N, T> const &e) {
		if (this != &e) {
			range(i, 0, N) { if (data[i] != e[i]) return 0; }
		}
		return 1;
	}

	/** The operator = to copy a value in all vector
	 *  @param a The value to copy **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(T const &a) {
		range(i, 0, N) data[i] = a;
		return *this;
	}

	/** The operator = to copy a vector in all vector
	 *  @param v The vector to copy **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(Meta::Expr<N, T> const &e) {
		if (this != &e) { range(i, 0, N) data[i] = e[i]; }
		return *this;
	}

	/** The operator += to add a value to all components
	 *  @param a The value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator+=(T const &a) {
		if (a != 0) { range(i, 0, N) data[i] += a; }
		return *this;
	}

	/** The operator += to add a vector
	 *  @param v The value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator+=(Meta::Expr<N, T> const &e) {
		range(i, 0, N) data[i] += e[i];
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
	 *  @param v The value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator-=(Meta::Expr<N, T> const &e) {
		range(i, 0, N) data[i] -= e[i];
		return *this;
	}

	/** The operator *= to multiply a vector
	 *  @param a The value to multiply **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator*=(T const &a) {
		range(i, 0, N) data[i] *= a;
		return *this;
	}

	/** The operator *= to divide a vector
	 *  @param a The value to divide **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator/=(T const &a) {
		range(i, 0, N) data[i] /= a;
		return *this;
	}

	/** The operator == to compare two vectors
	 *  @param v The first vector
	 *  @param w The second vector **/
	template<size_t M, typename S>
	bool operator==(Meta::Expr<M, S> const &v, Meta::Expr<M, S> const &w) {
		if (&v != &w ) {
			range(i, 0, M) { if (v[i] != w[i]) return 0; }
		}
		return 1;
	}

	/** The operator + to add two vectors
	 *  @param l The first vector
	 *  @param r The second vector **/
	template<size_t M, typename S>
	Meta::BinExpr<M, S, Common::Meta::Sum> &operator+(Meta::Expr<M, S> const &l, Meta::Expr<M, S> const &r) {
		return Meta::BinExpr<M, S, Common::Meta::Sum>(&l, &r);
	}

	/** The operator + to add a vector and a value
	 *  @param v The vector
	 *  @param a The value **/
	template<size_t M, typename S>
	Meta::BinExpr<M, S, Common::Meta::Sum> &operator+(Meta::Expr<M, S> const &l, S const &r) {
		return Meta::BinExpr<M, S, Common::Meta::Sum>(&l, Meta::Const<M, S>(&r));
	}

	/** The operator + to add a vector and a value
	 *  @param a The value
	 *  @param v The vector **/
	template<size_t M, typename S>
	Meta::BinExpr<M, S, Common::Meta::Sum> &operator+(S const &l, Meta::Expr<M, S> const &r) {
		return Meta::BinExpr<M, S, Common::Meta::Sum>(Meta::Const<M, S>(&l), &r);
	}

	/** The operator - to subtract two vectors
	 *  @param v The first vector
	 *  @param w The second vector **/
	template<size_t M, typename S>
	Meta::BinExpr<M, S, Common::Meta::Sub> &operator-(Meta::Expr<M, S> const &l, Meta::Expr<M, S> const &r) {
		return Meta::BinExpr<M, S, Common::Meta::Sub>(l, r);
	}

	/** The operator - to subtract a vector and a value
	 *  @param v The vector
	 *  @param a The value **/
	template<size_t M, typename S>
	Meta::BinExpr<M, S, Common::Meta::Sub> &operator-(Meta::Expr<M, S> const &l, S const &r) {
		return Meta::BinExpr<M, S, Common::Meta::Sub>(l, Meta::Const<M, S>(r));
	}

	/** The operator - to subtract a vector and a value
	 *  @param a The value
	 *  @param v The vector **/
	template<size_t M, typename S>
	Meta::BinExpr<M, S, Common::Meta::Sub> &operator-(S const &l, Meta::Expr<M, S> const &r) {
		return Meta::BinExpr<M, S, Common::Meta::Sub>(Meta::Const<M, S>(l), r);
	}

	/** Moltiplication by a scalar
	 *  @param a The scalar
	 *  @param v The vector **/
	template<size_t M, typename S> 
	Meta::BinExpr<M, S, Common::Meta::Mul> &operator*(S const &a, Meta::Expr<M, S> &v) {
		return Meta::BinExpr<M, S, Common::Meta::Mul>(Meta::Const<M ,S>(a), v);
	}

	/** Moltiplication by a scalar
	 *  @param v The vector
	 *  @param a The scalar **/
	template<size_t M, typename S> 
	Meta::BinExpr<M, S, Common::Meta::Mul> &operator*(Meta::Expr<M, S> &v, S const &a) {
		return Meta::BinExpr<M, S, Common::Meta::Mul>(v, Meta::Const<M, S>(a));
	}

	/** Division by a scalar
	 *  @param v The vector
	 *  @param a The scalar **/
	template<size_t M, typename S>
	Meta::BinExpr<M, S, Common::Meta::Mul> &operator/(Meta::Expr<M, S> &v, S const &a) {
		return Meta::BinExpr<M, S, Common::Meta::Div>(v, Meta::Const<M, S>(a));
	}

	/** The scalar product between two vector
	 *  @param v The first vector
	 *  @param w The second vector **/
	template<size_t M, typename S>
	S &operator*(Vector<M, S> const &v, Vector<M, S> const &w) {
		return Meta::Dot<M, M, S>(v, w);
	}

	/** The operator << to prin a vector in a stream
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
