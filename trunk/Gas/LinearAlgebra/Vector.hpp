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

#ifndef _GAS_LINEARALGEBRA_VECTOR_H_ /* BEGIN _GAS_LINEARALGEBRA_VECTOR_H_ */
#define _GAS_LINEARALGEBRA_VECTOR_H_

#include "Gas/Gas.h"

/* Definitions of classes methods */
namespace LinearAlgebra {

	/* A simple class for vector that use the expressions templates for evaluating operations */
	template<size_t N, typename T, typename C>
	class Vector {
		public:
			/* Constructor */
			Vector();
			Vector(T const);
			Vector(Vector<N, T, C> const &);
			template<typename D>
			Vector(Vector<N, T, D> const &);
			template<typename E, typename F>
			explicit Vector(Meta::Expression<N, T, E, F> const &);

			/* Desctructor */
			~Vector();

			/* Question */
			inline size_t Size();
			inline T &operator()(size_t const);
			inline T const &operator()(size_t const) const;

			/* Check */
			bool operator==(T const &);
			template<typename D>
			bool operator==(Vector<N, T, D> const &);
			template<typename E, typename F>
			bool operator==(Meta::Expression<N, T, E, F> const &);

			/* Copy */
			Vector<N, T, C> &operator=(T const &);
			template<typename D>
			Vector<N, T, C> &operator=(Vector<N, T, D> const &);
			template<typename E, typename F>
			Vector<N, T, C> &operator=(Meta::Expression<N, T, E, F> const &);

			/* Sum */
			Vector<N, T, C> &operator+=(T const &);
			template<typename D>
			Vector<N, T, C> &operator+=(Vector<N, T, D> const &);
			template<typename E, typename F>
			Vector<N, T, C> &operator+=(Meta::Expression<N, T, E, F> const &);

			/* Subtract */
			Vector<N, T, C> &operator-=(T const &);
			template<typename D>
			Vector<N, T, C> &operator-=(Vector<N, T, D> const &);
			template<typename E, typename F>
			Vector<N, T, C> &operator-=(Meta::Expression<N, T, E, F> const &);

			/* Multiply */
			Vector<N, T, C> &operator*=(T const &);

			/* Divide */
			Vector<N, T, C> &operator/=(T const &);

		private:
			C v_;
	};

	/* Vector<N, T, C> + Vector<N, T, D> */
	template<size_t N, typename T, typename C, typename D>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T, C> &, Vector<N, T, D> &);
	/* T + Vector<N, T, C> */
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Sum<T> > >
	operator+(T const &, Vector<N, T, C> &);
	/* Vector<N, T> + T */
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T, C> &, T const &);
	/* Expression<N, T, E, F> + Vector<N, T, C> */
	template<size_t N, typename T, typename E, typename F, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, F>, Vector<N, T, C>, Common::Meta::Math::Sum<T> > >
	operator+(Meta::Expression<N, T, E, F> const &, Vector<N, T, C> &);
	/* Vector<N, T, C> + Expression<N, T, E, F> */
	template<size_t N, typename T, typename C, typename E, typename F>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, Meta::Expression<N, T, E, F>, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T, C> &, Meta::Expression<N, T, E, F> const &);

	/* Vector<N, T, C> - Vector<N, T, D> */
	template<size_t N, typename T, typename C, typename D>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T, C> &, Vector<N, T, D> &);
	/* T - Vector<N, T, C> */
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Sub<T> > >
	operator-(T const &, Vector<N, T, C> &);
	/* Vector<N, T, C> - T */
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T, C> &, T const &);
	/* Expression<N, T, E, F> - Vector<N, T, C> */
	template<size_t N, typename T, typename E, typename F, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, F>, Vector<N, T, C>, Common::Meta::Math::Sub<T> > >
	operator-(Meta::Expression<N, T, E, F> const &, Vector<N, T, C> &);
	/* Vector<N, T, C> - Expression<N, T, E, F> */
	template<size_t N, typename T, typename C, typename E, typename F>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, Meta::Expression<N, T, E, F>, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T, C> &, Meta::Expression<N, T, E, F> const &);

	/* T * Vector<N, T, C> */
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Mul<T> > >
	operator*(T const &, Vector<N, T, C> &);
	/* Vector<N, T, C> * T */
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Mul<T> > >
	operator*(Vector<N, T, C> &, T const &);

	/* Vector<N, T, C> / T */
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Div<T> > >
	operator/(Vector<N, T, C> &, T const &);

	/* Vector<N, T, C> * Vector<N, T, D> */
	template<size_t N, typename T, typename C, typename D>
	T operator*(Vector<N, T, C> const &, Vector<N, T, D> const &);
	/* Vector<N, T, C> * Expression<N, T, E, F> */
	template<size_t N, typename T, typename C, typename E, typename F>
	T operator*(Vector<N, T, C> const &, Meta::Expression<N, T, E, F> const &);
	/* Expression<N, T, E, G> * Vector<N, T, C> */
	template<size_t N, typename T, typename E, typename F, typename C> 
	T operator*(Meta::Expression<N, T, E, F> const &, Vector<N, T, C> const &);
	
	namespace Meta {
		/* A container for expression */
		template<size_t N, typename T, typename F, typename C>
		class Expression<N, T, Vector<N, T, C>, F> {
			public:
				explicit Expression(Vector<N, T, C> const &);
				~Expression();
				inline T const operator()(size_t const &) const;
			private:
				Vector<N, T, C> const &e_;
		};
		/* A container for binary expression: F(Vector<N, T, C>, Vector<N, T, D>) */
		template<size_t N, typename T, typename C, typename D, typename F>
		class BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, F> {
			public:
				explicit BinaryExpression(Vector<N, T, C> const &, Vector<N, T, D> const &);
				~BinaryExpression();
				inline T const operator()(size_t const &) const;
			private:
				Vector<N, T, C> const &l_;
				Vector<N, T, D> const &r_;
		};
		/* A container for binary expression: F(Vector<N, T, C>, R) */
		template<size_t N, typename T, typename C, typename R, typename F>
		class BinaryExpression<N, T, Vector<N, T, C>, R, F> {
			public:
				explicit BinaryExpression(Vector<N, T, C> const &, R const &);
				~BinaryExpression();
				inline T const operator()(size_t const &) const;
			private:
				Vector<N, T, C> const &l_;
				R r_;
		};
		/* A container for binary expression: F(L, Vector<N, T, C>) */
		template<size_t N, typename T, typename L, typename C, typename F>
		class BinaryExpression<N, T, L, Vector<N, T, C>, F> {
			public:
				explicit BinaryExpression(L const &, Vector<N, T, C> const &);
				~BinaryExpression();
				inline T const operator()(size_t const &) const;
			private:
				L l_;
				Vector<N, T, C> const &r_;
		};
		/* A container for binary expression: F(Vector<N, T, C>, T) */
		template<size_t N, typename T, typename C, typename F>
		class BinaryExpression<N, T, Vector<N, T, C>, T, F> {
			public:
				explicit BinaryExpression(Vector<N, T, C> const &, T const &);
				~BinaryExpression();
				inline T const operator()(size_t const &) const;
			private:
				Vector<N, T, C> const &l_;
				T r_;
		};
		/* A container for binary expression: F(T, Vector<N, T, C>) */
		template<size_t N, typename T, typename C, typename F>
		class BinaryExpression<N, T, T, Vector<N, T, C>, F> {
			public:
				explicit BinaryExpression(T const &, Vector<N, T, C> const &);
				~BinaryExpression();
				inline T const operator()(size_t const &) const;
			private:
				T l_;
				Vector<N, T, C> const &r_;
		};
	}
}

#include "Vector.cpp"

#endif /* END _GAS_LINEARALGEBRA_VECTOR_H_ */
