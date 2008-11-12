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

#ifndef _GAS_META_H_
#define _GAS_META_H_

namespace Common { namespace Meta {
	template<typename T>
	class Expr {
		public:
			virtual T &operator()(T const &x) = 0 ;
	};

	template<typename T>
	class Const: virtual public Expr<T> {
		private:
			T &a;
		public:
			Const(T const &a): a(a) {};
			T &operator()(T const &x) { return a; };
	};

	template<typename T, T F(T, T)>
	class BinExpr: virtual public Expr<T> {
		private:
			Expr<T> &l;
			Expr<T> &r;
		public:
			BinExpr(Expr<T> const &l, Expr<T> const &r): l(l), r(r) {};
			T &operator()(T const &x) { return F(l(x), r(x)); }
	};

	template<typename T>
	class Identity: virtual public Expr<T> {
		public:
			T &operator()(T const &v) { return v; }

			template<typename S> friend BinExpr<S, Sum> &operator+(Expr<S> const &l, Expr<S> const &r) { return BinExpr<S, Sum>(l, r); }
			template<typename S> friend BinExpr<S, Sum> &operator+(S const &l, Expr<S> const &r) { return BinExpr<S, Sum>(Const<S>(l), r); }
			template<typename S> friend BinExpr<S, Sum> &operator+(Expr<S> const &l, S const &r) { return BinExpr<S, Sum>(l, Const<S>(r)); }

			template<typename S> friend BinExpr<S, Sub> &operator-(Expr<S> const &l, Expr<S> const &r) { return BinExpr<S, Sub>(l, r); }
			template<typename S> friend BinExpr<S, Sub> &operator-(S const &l, Expr<S> const &r) { return BinExpr<S, Sub>(Const<S>(l), r); }
			template<typename S> friend BinExpr<S, Sub> &operator-(Expr<S> const &l, S const &r) { return BinExpr<S, Sub>(l, Const<S>(r)); }

			template<typename S> friend BinExpr<S, Mul> &operator*(Expr<S> const &l, Expr<S> const &r) { return BinExpr<S, Mul>(l, r); }
			template<typename S> friend BinExpr<S, Mul> &operator*(S const &l, Expr<S> const &r) { return BinExpr<S, Mul>(Const<S>(l), r); }
			template<typename S> friend BinExpr<S, Mul> &operator*(Expr<S> const &l, S const &r) { return BinExpr<S, Mul>(l, Const<S>(r)); }

			template<typename S> friend BinExpr<S, Div> &operator/(Expr<S> const &l, Expr<S> const &r) { return BinExpr<S, Div>(l, r); }
			template<typename S> friend BinExpr<S, Div> &operator/(S const &l, Expr<S> const &r) { return BinExpr<S, Div>(Const<S>(l), r); }
			template<typename S> friend BinExpr<S, Div> &operator/(Expr<S> const &l, S const &r) { return BinExpr<S, Div>(l, Const<S>(r)); }
	};
} }

#endif
