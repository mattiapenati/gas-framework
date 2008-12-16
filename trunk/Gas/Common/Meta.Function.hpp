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

#ifndef _GAS_META_FUNCTION_H_ /* BEGIN _GAS_META_FUNCTION_H_ */
#define _GAS_META_FUNCTION_H_

#include "Meta.Math.hpp"

namespace Common { namespace Meta { namespace Function {
	
	/* Classes */
	template<typename T> class Function;
	template<typename T> class Variable;
	template<typename T> class VirtualExpression;
	template<typename T, typename E, typename F=Common::Meta::Math::Id<T> > class Expression;
	template<typename T, typename L, typename R, typename F> class BinaryExpression;

	template<typename T>
	class Function {
		private:
			VirtualExpression<T> e_;
		public:
			template<typename E, typename F> Function(Expression<T, E, F> e): e_(e) {}
			inline T const operator()(T const &x) { return e_(x); }
	};

	/** @class Variable
	 *  @brief An incubator for function **/
	template<typename T>
	class Variable {
		public:
			inline T const operator()(T const &x) { return x; }
	};

	template<typename T>
	class VirtualExpression {
		public:
			virtual T const operator()(T const &) = 0;
	};

	/** @class Expression
	 *  @brief The main structure for evalutation tree **/
	template<typename T, typename E, typename F>
	class Expression: public VirtualExpression<T> {
		private:
			E e_;
		public:
			Expression(E const &e): e_(e) {}
			inline T const operator()(T const &x) { return F::RET(e_(x));	 } 
	};
	template<typename T, typename F>
	class Expression<T, Variable<T>, F>: public VirtualExpression<T> {
		private:
			Variable<T> &e_;
		public:
			Expression(Variable<T> const &e): e_(e) {}
			inline T const operator()(T const &x) { return F::RET(e_(x));	 } 
	};

	/** @class Binary expression
	 *  @brief The structure for binary expression **/
	template<typename T, typename L, typename R, typename F>
	class BinaryExpression {
		private:
			L l_;
			R r_;
		public:
			BinaryExpression(L const &l, R const &r): l_(l), r_(r) {}
			inline T const operator()(T const &x) { return F::RET(l_(x), r_(x)); }
	};
	template<typename T, typename L, typename F>
	class BinaryExpression<T, L, T, F> {
		private:
			L l_;
			T r_;
		public:
			BinaryExpression(L const &l, T const &r): l_(l), r_(r) {}
			inline T const operator()(T const &x) { return F::RET(l_(x), r_); }
	};
	template<typename T, typename R, typename F>
	class BinaryExpression<T, T, R, F> {
		private:
			T l_;
			R r_;
		public:
			BinaryExpression(T const &l, R const &r): l_(l), r_(r) {}
			inline T const operator()(T const &x) { return F::RET(l_, r_(x)); }
	};
	template<typename T, typename R, typename F>
	class BinaryExpression<T, Variable<T>, R, F> {
		private:
			Variable<T> &l_;
			R r_;
		public:
			BinaryExpression(Variable<T> &l, R const &r): l_(l), r_(r) {}
			inline T const operator()(T const &x) { return F::RET(l_(x), r_(x)); }
	};
	template<typename T, typename L, typename F>
	class BinaryExpression<T, L, Variable<T>, F> {
		private:
			L l_;
			Variable<T> &r_;
		public:
			BinaryExpression(L const &l, Variable<T> &r): l_(l), r_(r) {}
			inline T const operator()(T const &x) { return F::RET(l_(x), r_(x)); }
	};
	template<typename T, typename F>
	class BinaryExpression<T, Variable<T>, T, F> {
		private:
			Variable<T> &l_;
			T r_;
		public:
			BinaryExpression(Variable<T> &l, T const &r): l_(l), r_(r) {}
			inline T const operator()(T const &x) { return F::RET(l_(x), r_); }
	};
	template<typename T, typename F>
	class BinaryExpression<T, T, Variable<T>, F> {
		private:
			T l_;
			Variable<T> &r_;
		public:
			BinaryExpression(T const &l, Variable<T> &r): l_(l), r_(r) {}
			inline T const operator()(T const &x) { return F::RET(l_, r_(x)); }
	};
	template<typename T, typename F>
	class BinaryExpression<T, Variable<T>, Variable<T>, F> {
		private:
			Variable<T> &l_;
			Variable<T> &r_;
		public:
			BinaryExpression(Variable<T> &l, Variable<T> &r): l_(l), r_(r) {}
			inline T const operator()(T const &x) { return F::RET(l_(x), r_(x)); }
	};

	/**
	 * Defining the four basic operations
	 */
	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Sum<T> > >
	operator+(Variable<T> &a, Variable<T> &b) {
		typedef BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Sum<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Sum<T> > >
	operator+(T const &a, Variable<T> &b) {
		typedef BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Sum<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Sum<T> > >
	operator+(Variable<T> &a, T const &b) {
		typedef BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Sum<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Sum<T> > >
	operator+(Expression<T, E, G> const &a, Variable<T> &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Sum<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Sum<T> > >
	operator+(Variable<T> &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Sum<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename F, typename G, typename H>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Sum<T> > >
	operator+(Expression<T, E, G> const &a, Expression<T, F, H> const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Sum<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Sum<T> > >
	operator+(Expression<T, E, G> const &a, T const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Sum<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Sum<T> > >
	operator+(T const &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Sum<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}

	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Sub<T> > >
	operator-(Variable<T> &a, Variable<T> &b) {
		typedef BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Sub<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Sub<T> > >
	operator-(T const &a, Variable<T> &b) {
		typedef BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Sub<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Sub<T> > >
	operator-(Variable<T> &a, T const &b) {
		typedef BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Sub<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Sub<T> > >
	operator-(Expression<T, E, G> const &a, Variable<T> &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Sub<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Sub<T> > >
	operator-(Variable<T> &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Sub<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename F, typename G, typename H>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Sub<T> > >
	operator-(Expression<T, E, G> const &a, Expression<T, F, H> const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Sub<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Sub<T> > >
	operator-(Expression<T, E, G> const &a, T const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Sub<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Sub<T> > >
	operator-(T const &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Sub<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}

	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Mul<T> > >
	operator*(Variable<T> &a, Variable<T> &b) {
		typedef BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Mul<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Mul<T> > >
	operator*(T const &a, Variable<T> &b) {
		typedef BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Mul<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Mul<T> > >
	operator*(Variable<T> &a, T const &b) {
		typedef BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Mul<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Mul<T> > >
	operator*(Expression<T, E, G> const &a, Variable<T> &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Mul<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Mul<T> > >
	operator*(Variable<T> &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Mul<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename F, typename G, typename H>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Mul<T> > >
	operator*(Expression<T, E, G> const &a, Expression<T, F, H> const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Mul<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Mul<T> > >
	operator*(Expression<T, E, G> const &a, T const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Mul<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Mul<T> > >
	operator*(T const &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Mul<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}

	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Div<T> > >
	operator/(Variable<T> &a, Variable<T> &b) {
		typedef BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Div<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Div<T> > >
	operator/(T const &a, Variable<T> &b) {
		typedef BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Div<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Div<T> > >
	operator/(Variable<T> &a, T const &b) {
		typedef BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Div<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Div<T> > >
	operator/(Expression<T, E, G> const &a, Variable<T> &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Div<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Div<T> > >
	operator/(Variable<T> &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Div<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename F, typename G, typename H>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Div<T> > >
	operator/(Expression<T, E, G> const &a, Expression<T, F, H> const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Div<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Div<T> > >
	operator/(Expression<T, E, G> const &a, T const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Div<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Div<T> > >
	operator/(T const &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Div<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}

	/* Division mod */
	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Mod<T> > >
	operator%(Variable<T> &a, Variable<T> &b) {
		typedef BinaryExpression<T, Variable<T>, Variable<T>, Common::Meta::Math::Mod<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Mod<T> > >
	operator%(T const &a, Variable<T> &b) {
		typedef BinaryExpression<T, T, Variable<T>, Common::Meta::Math::Mod<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T>
	Expression<T, BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Mod<T> > >
	operator%(Variable<T> &a, T const &b) {
		typedef BinaryExpression<T, Variable<T>, T, Common::Meta::Math::Mod<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Mod<T> > >
	operator%(Expression<T, E, G> const &a, Variable<T> &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Variable<T>, Common::Meta::Math::Mod<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Mod<T> > >
	operator%(Variable<T> &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, Variable<T>, Expression<T, E, G>, Common::Meta::Math::Mod<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename F, typename G, typename H>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Mod<T> > >
	operator%(Expression<T, E, G> const &a, Expression<T, F, H> const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, Expression<T, F, H>, Common::Meta::Math::Mod<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Mod<T> > >
	operator%(Expression<T, E, G> const &a, T const &b) {
		typedef BinaryExpression<T, Expression<T, E, G>, T, Common::Meta::Math::Mod<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}
	template<typename T, typename E, typename G>
	Expression<T, BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Mod<T> > >
	operator%(T const &a, Expression<T, E, G> const &b) {
		typedef BinaryExpression<T, T, Expression<T, E, G>, Common::Meta::Math::Mod<T> > TempType;
		return Expression<T, TempType>(TempType(a, b));
	}

	/* Absolute value */
	template<typename T>
	Expression<T, Variable<T>, Common::Meta::Math::Abs<T> >
	Abs(Variable<T> &x) {
		return Expression<T, Variable<T>, Common::Meta::Math::Abs<T> >(x);
	}
	template<typename T, typename E>
	Expression<T, Expression<T, E>, Common::Meta::Math::Abs<T> >
	Abs(Expression<T, E> const &x) {
		return Expression<T, Expression<T, E>, Common::Meta::Math::Abs<T> >(x);
	}

	/* Complex conjugate */
	template<typename T>
	Expression<T, Variable<T>, Common::Meta::Math::Conj<T> >
	Conj(Variable<T> &x) {
		return Expression<T, Variable<T>, Common::Meta::Math::Conj<T> >(x);
	}
	template<typename T, typename E>
	Expression<T, Expression<T, E>, Common::Meta::Math::Conj<T> >
	Conj(Expression<T, E> const &x) {
		return Expression<T, Expression<T, E>, Common::Meta::Math::Conj<T> >(x);
	}

	/* Square root */
	template<typename T>
	Expression<T, Variable<T>, Common::Meta::Math::Sqrt<T> >
	Sqrt(Variable<T> &x) {
		return Expression<T, Variable<T>, Common::Meta::Math::Sqrt<T> >(x);
	}
	template<typename T, typename E>
	Expression<T, Expression<T, E>, Common::Meta::Math::Sqrt<T> >
	Sqrt(Expression<T, E> const &x) {
		return Expression<T, Expression<T, E>, Common::Meta::Math::Sqrt<T> >(x);
	}
}}}

#endif /* END _GAS_META_FUNCTION_H_ */
