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

#ifndef _GAS_LINEARALGEBRA_META_H_ /* BEGIN _GAS_LINEARALGEBRA_META_H_ */
#define _GAS_LINEARALGEBRA_META_H_

#include "Gas/Gas.h"

namespace LinearAlgebra {
	namespace Meta {
		/* A container for expression */
		template<size_t N, typename T, typename E, typename F>
		class Expression {
			public:
				explicit Expression(E const &);
				~Expression();
				inline T const operator()(size_t const &) const;
			private:
				E e_;
		};
		
		/* A container for binary expression: F(L, R) */
		template<size_t N, typename T, typename L, typename R, typename F>
		class BinaryExpression {
			public:
				explicit BinaryExpression(L const &, R const &);
				~BinaryExpression();
				inline T const operator()(size_t const &) const;
			private:
				L l_;
				R r_;
		};
		/* A container for binary expression: F(L, T) */
		template<size_t N, typename T, typename L, typename F>
		class BinaryExpression<N, T, L, T, F> {
			public:
				explicit BinaryExpression(L const &, T const &);
				~BinaryExpression();
				inline T const operator()(size_t const &) const;
			private:
				L l_;
				T r_;
		};
		/* A container for binary expression: F(T, R) */
		template<size_t N, typename T, typename R, typename F>
		class BinaryExpression<N, T, T, R, F> {
			public:
				explicit BinaryExpression(T const &, R const &);
				~BinaryExpression();
				inline T const operator()(size_t const &) const;
			private:
				T l_;
				R r_;
		};
		
		/* Expression<N, T, E, F> + Expression<N, T, G, H> */
		template<size_t N, typename T, typename E, typename F, typename G, typename H>
		inline Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, Expression<N, T, G, H>, Common::Meta::Math::Sum<T> > >
		operator+(Expression<N, T, E, F> const &, Expression<N, T, G, H> const &);
		/* Expression<N, T, E, F> + T */
		template<size_t N, typename T, typename E, typename F>
		inline Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Sum<T> > >
		operator+(Expression<N, T, E, F> const &, T const &);
		/* T + Expression<N, T, E, F> */
		template<size_t N, typename T, typename E, typename F>
		inline Expression<N, T, BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Sum<T> > >
		operator+(T const &, Expression<N, T, E, F> const &);
		
		/* Expression<N, T, E, F> - Expression<N, T, G, H> */
		template<size_t N, typename T, typename E, typename F, typename G, typename H>
		inline Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, Expression<N, T, G, H>, Common::Meta::Math::Sub<T> > >
		operator-(Expression<N, T, E, F> const &, Expression<N, T, G, H> const &);
		/* Expression<N, T, E, F> - T */
		template<size_t N, typename T, typename E, typename F>
		inline Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Sub<T> > >
		operator-(Expression<N, T, E, F> const &, T const &);
		/* T - Expression<N, T, E, F> */
		template<size_t N, typename T, typename E, typename F>
		inline Expression<N, T, BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Sub<T> > >
		operator-(T const &, Expression<N, T, E, F> const &);
		
		/* Expression<N, T, E, F> * T */
		template<size_t N, typename T, typename E, typename F>
		inline Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Mul<T> > >
		operator*(Expression<N, T, E, F> const &, T const &);
		/* T * Expression<N, T, E, F> */
		template<size_t N, typename T, typename E, typename F>
		inline Expression<N, T, BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Mul<T> > >
		operator*(T const &, Expression<N, T, E, F> const &);
		
		/* Expression<N, T, E, F> / T */
		template<size_t N, typename T, typename E, typename F>
		inline Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Div<T> > >
		operator/(Expression<N, T, E, F> const &, T const &);
		
		/* Expression<N, T, E, F> * Expression<N, T, G, H> */
		template<size_t N, typename T, typename E, typename F, typename G, typename H> 
		T operator*(Expression<N, T, E, F> const &, Expression<N, T, G, H> const &);
	}
}

#include "Meta.cpp"

#endif /* END _GAS_LINEARALGEBRA_META_H_ */
