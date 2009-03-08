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

namespace LinearAlgebra {
	namespace Meta {
		/* Expression<N, T, E, F> */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, E, F>::Expression(E const &e): e_(e) {
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::Expression@"<<this<<": ===Create an Expression from Object@"<<&e<<"==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, E, F>::~Expression() {
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::Expression@"<<this<<": ===Destroyed==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename E, typename F>
		T const Expression<N, T, E, F>::operator()(size_t const &i) const {
			#if _GAS_VERBOSITY_ >= 5
			std::cerr<<"LinearAlgebra::Meta::Expression@"<<this<<": ===Evaluate at "<<i<<"==="<<std::endl;
			#endif
			return F::RET(e_(i));
		}
		
		/* BinaryExpression<N, T, L, R, F> */
		template<size_t N, typename T, typename L, typename R, typename F>
		BinaryExpression<N, T, L, R, F>::BinaryExpression(L const &l, R const &r): l_(l), r_(r){
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Create a BinaryExpression from Object@"<<&l<<" and Object@"<<&r<<"==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename L, typename R, typename F>
		BinaryExpression<N, T, L, R, F>::~BinaryExpression() {
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Destroyed==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename L, typename R, typename F>
		T const BinaryExpression<N, T, L, R, F>::operator()(size_t const &i) const {
			#if _GAS_VERBOSITY_ >= 5
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Evaluate at "<<i<<"==="<<std::endl;
			#endif
			return F::RET(l_(i), r_(i));
		}
		
		/* BinaryExpression<N, T, L, T, F> */
		template<size_t N, typename T, typename L, typename F>
		BinaryExpression<N, T, L, T, F>::BinaryExpression(L const &l, T const &r): l_(l), r_(r){
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Create a BinaryExpression from Object@"<<&l<<" and "<<r<<"==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename L, typename F>
		BinaryExpression<N, T, L, T, F>::~BinaryExpression() {
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Destroyed==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename L, typename F>
		T const BinaryExpression<N, T, L, T, F>::operator()(size_t const &i) const {
			#if _GAS_VERBOSITY_ >= 5
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Evaluate at "<<i<<"==="<<std::endl;
			#endif
			return F::RET(l_(i), r_);
		}
		
		/* BinaryExpression<N, T, T, R, F> */
		template<size_t N, typename T, typename R, typename F>
		BinaryExpression<N, T, T, R, F>::BinaryExpression(T const &l, R const &r): l_(l), r_(r){
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Create a BinaryExpression from "<<l<<" and Object@"<<&r<<"==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename R, typename F>
		BinaryExpression<N, T, T, R, F>::~BinaryExpression() {
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Destroyed==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename R, typename F>
		T const BinaryExpression<N, T, T, R, F>::operator()(size_t const &i) const {
			#if _GAS_VERBOSITY_ >= 5
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Evaluate at "<<i<<"==="<<std::endl;
			#endif
			return F::RET(l_, r_(i));
		}
		
		/* Expression<N, T, E, F> + Expression<N, T, G, H> */
		template<size_t N, typename T, typename E, typename F, typename G, typename H>
		Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, Expression<N, T, G, H>, Common::Meta::Math::Sum<T> > >
		operator+(Expression<N, T, E, F> const &l, Expression<N, T, G, H> const &r) {
			#if _GAS_VERBOSITY_ >= 4
			std::cerr<<"LinearAlgebra::Meta: ===Expression@"<<&l <<" + Expression@"<<&r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, Expression<N, T, E, F>, Expression<N, T, G, H>, Common::Meta::Math::Sum<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		/* Expression<N, T, E, F> + T */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Sum<T> > >
		operator+(Expression<N, T, E, F> const &l, T const &r) {
			#if _GAS_VERBOSITY_ >= 4
			std::cerr<<"LinearAlgebra::Meta: ===Expression@"<<&l <<" + "<<r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Sum<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		/* T + Expression<N, T, E, F> */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Sum<T> > >
		operator+(T const &l, Expression<N, T, E, F> const &r) {
			#if _GAS_VERBOSITY_ >= 4
			std::cerr<<"LinearAlgebra::Meta: ==="<<l <<" + Expression@"<<&r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Sum<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		
		/* Expression<N, T, E, F> - Expression<N, T, G, H> */
		template<size_t N, typename T, typename E, typename F, typename G, typename H>
		Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, Expression<N, T, G, H>, Common::Meta::Math::Sub<T> > >
		operator-(Expression<N, T, E, F> const &l, Expression<N, T, G, H> const &r) {
			#if _GAS_VERBOSITY_ >= 4
			std::cerr<<"LinearAlgebra::Meta: ===Expression@"<<&l <<" - Expression@"<<&r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, Expression<N, T, E, F>, Expression<N, T, G, H>, Common::Meta::Math::Sub<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		/* Expression<N, T, E, F> - T */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Sub<T> > >
		operator-(Expression<N, T, E, F> const &l, T const &r) {
			#if _GAS_VERBOSITY_ >=4
			std::cerr<<"LinearAlgebra::Meta: ===Expression@"<<&l <<" - "<<r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Sub<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		/* T - Expression<N, T, E, F> */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Sub<T> > >
		operator-(T const &l, Expression<N, T, E, F> const &r) {
			#if _GAS_VERBOSITY_ >=4
			std::cerr<<"LinearAlgebra::Meta: ==="<<l <<" - Expression@"<<&r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Sub<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		
		/* Expression<N, T, E, F> * T */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Mul<T> > >
		operator*(Expression<N, T, E, F> const &l, T const &r) {
			#if _GAS_VERBOSITY_ >=4
			std::cerr<<"LinearAlgebra::Meta: ===Expression@"<<&l<<" * "<<r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Mul<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		/* T * Expression<N, T, E, F> */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Mul<T> > >
		operator*(T const &l, Expression<N, T, E, F> const &r) {
			#if _GAS_VERBOSITY_ >=4
			std::cerr<<"LinearAlgebra::Meta: ==="<<l <<" * Expression@"<<&r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, T, Expression<N, T, E, F>, Common::Meta::Math::Mul<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		
		/* Expression<N, T, E, F> / T */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Div<T> > >
		operator/(Expression<N, T, E, F> const &l, T const &r) {
			#if _GAS_VERBOSITY_ >=4
			std::cerr<<"LinearAlgebra::Meta: ===Expression@"<<&l<<" / "<<r<<"==="<<std::endl;
			#endif
			typedef BinaryExpression<N, T, Expression<N, T, E, F>, T, Common::Meta::Math::Div<T> > TempType;
			return Expression<N, T, TempType>(TempType(l, r));
		}
		
		/* Expression<N, T, E, F> * Expression<N, T, G, H> */
		template<size_t N, typename T, typename E, typename F, typename G, typename H> 
		T operator*(Expression<N, T, E, F> const &l, Expression<N, T, G, H> const &r) {
			#if _GAS_VERBOSITY_ >= 4
			std::cerr<<"LinearAlgebra: ===Expression@"<<&l<<" * Expression@"<<&r<<"==="<<std::endl;
			#endif
			T s = T();
			for(int i=0; i<N; i+=4) {
				s += (l(i) * Common::Math::Conj(r(i)));
				s += (l(i+1) * Common::Math::Conj(r(i+1)));
				s += (l(i+2) * Common::Math::Conj(r(i+2)));
				s += (l(i+3) * Common::Math::Conj(r(i+3)));
			}
			switch (N % 4) {
				case 3: s += (l(N-3) * Common::Math::Conj(r(N-3)));
				case 2: s += (l(N-2) * Common::Math::Conj(r(N-2)));
				case 1: s += (l(N-1) * Common::Math::Conj(r(N-1)));
			}
			return s;
		}
		
	}
}
