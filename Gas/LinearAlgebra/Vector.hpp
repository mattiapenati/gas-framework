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

/* Thanks to Todd Veldhuizen for Expression Templates */

/* FLAGS
 *  - _GAS_CHECK_INDEX_: check if is a valid index when you use v(i)
 *  - _GAS_BE_VERBOSE_: first level of verbosity
 *  - _GAS_BE_VERY_VERBOSE_: second level of verbosity
 */

#ifndef _GAS_VECTOR_H_ /* BEGIN _GAS_VECTOR_H_ */
#define _GAS_VECTOR_H_

#ifdef _GAS_CHECK_INDEX_
#include <cassert>
#endif

	/* only for size_t, problably I move it in Gas/Macro (TODO) */
#include <cstddef>
	/* only for cout and ostream */
#include <iostream>

#include "Gas/Common/Common.h"

/* Definition of classe names */
namespace LinearAlgebra {
	template<size_t N, typename T=double> class Vector;
	namespace Meta {
		template<size_t N, typename T, typename E, typename F=Common::Meta::Math::Id<T> > class Expression;
		template<size_t N, typename T, typename L, typename R, typename F> class BinaryExpression;
	}
}

/* Definition of classe methods */
namespace LinearAlgebra {

	/** A simple class for vector
	 *  @class Vector 
	 *  @brief A simple class for vector that use the expressions templates for operations **/
	template<size_t N, typename T>
	class Vector {
		public:
			/* Constructor */
			Vector();
			Vector(T const);
			Vector(Vector<N, T> const &);
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
			bool operator==(Vector<N, T> const &);
			template<typename E, typename F>
			bool operator==(Meta::Expression<N, T, E, F> const &);

			/* Copy */
			Vector<N, T> &operator=(T const &);
			Vector<N, T> &operator=(Vector<N, T> const &);
			template<typename E, typename F>
			Vector<N, T> &operator=(Meta::Expression<N, T, E, F> const &);

			/* Sum */
			Vector<N, T> &operator+=(T const &);
			Vector<N, T> &operator+=(Vector<N, T> const &);
			template<typename E, typename F>
			Vector<N, T> &operator+=(Meta::Expression<N, T, E, F> const &);

			/* Subtract */
			Vector<N, T> &operator-=(T const &);
			Vector<N, T> &operator-=(Vector<N, T> const &);
			template<typename E, typename F>
			Vector<N, T> &operator-=(Meta::Expression<N, T, E, F> const &);

			/* Multiply */
			Vector<N, T> &operator*=(T const &);

			/* Divide */
			Vector<N, T> &operator/=(T const &);


		private:
			Common::Array<T, N> v_; /* If you want to use without GasFramewrok is simple, change in something with [] access */
			
			static T const Zero = T();
	};

	/* Vector<N, T> + Vector<N, T> */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T> &, Vector<N, T> &);
	/* T + Vector<N, T> */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Sum<T> > >
	operator+(T const &, Vector<N, T> &);
	/* Vector<N, T> + T */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T> &, T const &);
	/* Expression<N, T,...> + Vector<N, T> */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Vector<N, T>, Common::Meta::Math::Sum<T> > >
	operator+(Meta::Expression<N, T, E, G> const &, Vector<N, T> &);
	/* Vector<N, T> + Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T> &, Meta::Expression<N, T, E, G> const &);
	/* Expression<N, T,...> + Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename F, typename G, typename H>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Meta::Expression<N, T, F, H>, Common::Meta::Math::Sum<T> > >
	operator+(Meta::Expression<N, T, E, G> const &, Meta::Expression<N, T, F, H> const &);
	/* Expression<N, T,...> + T */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Sum<T> > >
	operator+(Meta::Expression<N, T, E, G> const &, T const &);
	/* T + Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sum<T> > >
	operator+(T const &, Meta::Expression<N, T, E, G> const &);

	/* Vector<N, T> - Vector<N, T> */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T> &, Vector<N, T> &);
	/* T - Vector<N ,T> */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Sub<T> > >
	operator-(T const &, Vector<N, T> &);
	/* Vector<N ,T> - T */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T> &, T const &);
	/* Expression<N, T,...> - Vector<N, T> */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Vector<N, T>, Common::Meta::Math::Sub<T> > >
	operator-(Meta::Expression<N, T, E, G> const &, Vector<N, T> &);
	/* Vector<N, T> - Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T> &, Meta::Expression<N, T, E, G> const &);
	/* Expression<N, T,....> - Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename F, typename G, typename H>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Meta::Expression<N, T, F, H>, Common::Meta::Math::Sub<T> > >
	operator-(Meta::Expression<N, T, E, G> const &, Meta::Expression<N, T, F, H> const &);
	/* Expression<N, T,...> - T */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Sub<T> > >
	operator-(Meta::Expression<N, T, E, G> const &, T const &);
	/* T - Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sub<T> > >
	operator-(T const &, Meta::Expression<N, T, E, G> const &);

	/* T * Vector<N, T> */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Mul<T> > >
	operator*(T const &, Vector<N, T> &);
	/* Vector<N, T> * T */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Mul<T> > >
	operator*(Vector<N, T> &, T const &);
	/* Expression<N, T,...> * T */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Mul<T> > >
	operator*(Meta::Expression<N, T, E, G> const &, T const &);
	/* T * Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Mul<T> > >
	operator*(T const &, Meta::Expression<N, T, E, G> const &);

	/* Vector<N, T> / T */
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Div<T> > >
	operator/(Vector<N, T> &, T const &);
	/* Expression<N, T,...> / T */
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Div<T> > >
	operator/(Meta::Expression<N, T, E, G> const &, T const &);

	/* Vector<N, T> * Vector<N, T> */
	template<size_t N, typename T>
	T operator*(Vector<N, T> const &, Vector<N, T> const &);
	/* Vector<N, T> * Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename F>
	T operator*(Vector<N, T> const &, Meta::Expression<N, T, E, F> const &);
	/* Expression<N, T,...> * Vector<N, T> */
	template<size_t N, typename T, typename E, typename F> 
	T operator*(Meta::Expression<N, T, E, F> const &, Vector<N, T> const &);
	/* Expression<N, T,...> * Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename F, typename G, typename H> 
	T operator*(Meta::Expression<N, T, E, G> const &, Meta::Expression<N, T, F, H> const &);

	/* ostream << Vector<N, T> */
	template<size_t N, typename T>
	std::ostream &operator<<(std::ostream &, Vector<N, T> const &);
	/* ostream << Expression<N, T,...> */
	template<size_t N, typename T, typename E, typename G>
	std::ostream &operator<<(std::ostream &,Meta::Expression<N, T, E, G> const &);

	namespace Meta {
		/** A container for expression
		 *  @class Expression
		 *  @brief A container for expression, usefull for unary expression **/
		template<size_t N, typename T, typename E, typename F>
		class Expression {
			public:
				explicit Expression(E const &);
				inline T const operator()(size_t const &) const;
			private:
				E e_;
		};
		template<size_t N, typename T, typename F>
		class Expression<N, T, Vector<N, T>, F> {
			public:
				explicit Expression(Vector<N, T> &);
				inline T const operator()(size_t const &) const;
			private:
				Vector<N, T> &e_;
		};

		/** A container for binary expression
		 *  @class BinaryExpression
		 *  @brief A container for binary expression, it isn't used directly, but in an other expression **/
		template<size_t N, typename T, typename L, typename R, typename F>
		class BinaryExpression {
			public:
				explicit BinaryExpression(L const &, R const &);
				inline T const operator()(size_t const &) const;
			private:
				L l_;
				R r_;
		};
		template<size_t N, typename T, typename L, typename F>
		class BinaryExpression<N, T, L, T, F> {
			public:
				explicit BinaryExpression(L const &, T const &);
				inline T const operator()(size_t const &) const;
			private:
				L l_;
				T r_;
		};
		template<size_t N, typename T, typename R, typename F>
		class BinaryExpression<N, T, T, R, F> {
			public:
				explicit BinaryExpression(T const &, R const &);
				inline T const operator()(size_t const &) const;
			private:
				T l_;
				R r_;
		};
		template<size_t N, typename T, typename R, typename F>
		class BinaryExpression<N, T, Vector<N, T>, R, F> {
			public:
				explicit BinaryExpression(Vector<N, T> &, R const &);
				inline T const operator()(size_t const &) const;
			private:
				Vector<N, T> &l_;
				R r_;
		};
		template<size_t N, typename T, typename L, typename F>
		class BinaryExpression<N, T, L, Vector<N, T>, F> {
			public:
				explicit BinaryExpression(L const &, Vector<N, T> &);
				inline T const operator()(size_t const &) const;
			private:
				L l_;
				Vector<N, T> &r_;
		};
		template<size_t N, typename T, typename F>
		class BinaryExpression<N, T, Vector<N, T>, T, F> {
			public:
				explicit BinaryExpression(Vector<N, T> &, T const &);
				inline T const operator()(size_t const &) const;
			private:
				Vector<N, T> &l_;
				T r_;
		};
		template<size_t N, typename T, typename F>
		class BinaryExpression<N, T, T, Vector<N, T>, F> {
			public:
				explicit BinaryExpression(T const &, Vector<N, T> &);
				inline T const operator()(size_t const &) const;
			private:
				T l_;
				Vector<N, T> &r_;
		};
		template<size_t N, typename T, typename F>
		class BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, F> {
			public:
				explicit BinaryExpression(Vector<N, T> &, Vector<N, T> &);
				inline T const operator()(size_t const &) const;
			private:
				Vector<N, T> &l_;
				Vector<N, T> &r_;
		};
	}
}

/*	 _                 _                           _        _   _
 *	(_)               | |                         | |      | | (_)
 *	 _ _ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_ _  ___  _ __
 *	| | '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *	| | | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 *	|_|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *	            | |
 *	            |_|
 */

namespace LinearAlgebra {

	/** The default constructor **/
	template<size_t N, typename T>
	Vector<N, T>::Vector() {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Create an empty Vector===" << std::endl;
		#endif
	}

	/** The constructor to initialize the entire vector with the same value
	 *  @param a The vaule to use **/
	template<size_t N, typename T>
	Vector<N, T>::Vector(T const a) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Create a Vector from value " << a << "===" << std::endl;
		#endif
		for(int i=0; i<N; i+=4) {
			v_[i] = a;
			v_[i+1] = a;
			v_[i+2] = a;
			v_[i+3] = a;
		}
		switch (N % 4) {
			case 3: v_[N-3] = a;
			case 2: v_[N-2] = a;
			case 1: v_[N-1] = a;
		}
	}

	/** The constructor to copy a vector
	 *  @param v The vector to copy **/
	template<size_t N, typename T>
	Vector<N, T>::Vector(Vector<N, T> const &v) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this <<": ===Create a Vector from Vector@" << &v << "===" << std::endl;
		#endif
		v_ = v.v_;
	}

	/** The constructor to copy a vector expression
	 *  @param e The expression to copy **/
	template<size_t N, typename T>
	template<typename E, typename F>
	Vector<N, T>::Vector(Meta::Expression<N, T, E, F> const &e) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Create a Vector from Expression@" << &e << "===" << std::endl;
		#endif
		for(int i=0; i<N; i+=4) {
			v_[i] = e(i);
			v_[i+1] = e(i+1);
			v_[i+2] = e(i+2);
			v_[i+3] = e(i+3);
		}
		switch (N % 4) {
			case 3: v_[N-3] = e(N-3);
			case 2: v_[N-2] = e(N-2);
			case 1: v_[N-1] = e(N-1);
		}
	}

	/** The default destructor **/
	template<size_t N, typename T>
	Vector<N, T>::~Vector() {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Destroyed===" << std::endl;
		#endif
	}

	/** The vector size **/
	template<size_t N, typename T>
	size_t Vector<N, T>::Size() {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Reading the size===" << std::endl;
		#endif
		return N;
	}

	/** The operator () to access to the components of vector
	 *  @param i The index of component **/
	template<size_t N, typename T>
	T &Vector<N, T>::operator()(size_t const i) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Accessing to the element at position " << i << "===" << std::endl;
		#endif
		#ifdef _GAS_CHECK_INDEX_
		assert(i < N);
		#endif
		return v_[i];
	}

	/** The operator () to access to the components of vector
	 *  @param i The index of component **/
	template<size_t N, typename T>
	T const &Vector<N, T>::operator()(size_t const i) const {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Accessing to the element at position " << i << "===" << std::endl;
		#endif
		#ifdef _GAS_CHECK_INDEX_
		assert(i < N);
		#endif
		return v_[i];
	}

	/** The operator == to compare all components of a vector with a value
	 *  @param a The value to compare **/
	template<size_t N, typename T>
	bool Vector<N, T>::operator==(T const &a) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Compare with the value " << a << "===" << std::endl;
		#endif
		range(i, 0, N) { if (!Common::Limits<T>::Equal(v_[i], a)) return false; }
		return true;
	}

	/** The operator == to compare two vectors
	 *  @param v The second vector **/
	template<size_t N, typename T>
	bool Vector<N, T>::operator==(Vector<N, T> const &v) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Compare with the Vector@" << &v << "===" << std::endl;
		#endif
		if (this != &v) {
			range(i, 0, N) { if (!Common::Limits<T>::Equal(v_[i], v.v_[i])) return false; }
		}
		return true;
	}

	/** The operator == to compare a vector with an expression
	 *  @param e An expression **/
	template<size_t N, typename T>
	template<typename E, typename F>
	bool Vector<N, T>::operator==(Meta::Expression<N, T, E, F> const &e) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Compare with the Expression@" << &e << "===" << std::endl;
		#endif
		range(i, 0, N) { if (!Common::Limits<T>::Equal(v_[i], e(i))) return false; }
		return true;
	}

	/** The operator = to copy a value in all components of a vector
	 *  @param a The scalar value to copy **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(T const &a) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Copy the value " << a << "===" << std::endl;
		#endif
		for(int i=0; i<N; i+=4) {
			v_[i] = a;
			v_[i+1] = a;
			v_[i+2] = a;
			v_[i+3] = a;
		}
		switch (N % 4) {
			case 3: v_[N-3] = a;
			case 2: v_[N-2] = a;
			case 1: v_[N-1] = a;
		}
		return *this;
	}

	/** The operator = to copy a vector in a vector
	 *  @param v The vector to copy **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator=(Vector<N, T> const &v) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Copy the Vector@" << &v << "===" << std::endl;
		#endif
		if (this != &v) v_ = v.v_;
		return *this;
	}

	/** The operator = to copy an expression in a vector
	 *  @param e An expression **/
	template<size_t N, typename T>
	template<typename E, typename F>
	Vector<N, T> &Vector<N, T>::operator=(Meta::Expression<N, T, E, F> const &e) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Copy the Expression@" << &e << "===" << std::endl;
		#endif
		for(int i=0; i<N; i+=4) {
			v_[i] = e(i);
			v_[i+1] = e(i+1);
			v_[i+2] = e(i+2);
			v_[i+3] = e(i+3);
		}
		switch (N % 4) {
			case 3: v_[N-3] = e(N-3);
			case 2: v_[N-2] = e(N-2);
			case 1: v_[N-1] = e(N-1);
		}
		return *this;
	}

	/** The operator += to add a value to all components
	 *  @param a The scalar value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator+=(T const &a) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Add the value " << a << "===" << std::endl;
		#endif
		if (a != 0.) {
			for(int i=0; i<N; i+=4) {
				v_[i] += a;
				v_[i+1] += a;
				v_[i+2] += a;
				v_[i+3] += a;
			}
			switch (N % 4) {
				case 3: v_[N-3] += a;
				case 2: v_[N-2] += a;
				case 1: v_[N-1] += a;
			}
		}
		return *this;
	}

	/** The operator += to add a vector
	 *  @param v The vector to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator+=(Vector<N, T> const &v) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Add the Vector@" << &v << "===" << std::endl;
		#endif
		for(int i=0; i<N; i+=4) {
			v_[i] += v.v_[i];
			v_[i+1] += v.v_[i+1];
			v_[i+2] += v.v_[i+2];
			v_[i+3] += v.v_[i+3];
		}
		switch (N % 4) {
			case 3: v_[N-3] += v.v_[N-3];
			case 2: v_[N-2] += v.v_[N-2];
			case 1: v_[N-1] += v.v_[N-1];
		}
		return *this;
	}

	/** The operator += to add an expression
	 *  @param e An expression **/
	template<size_t N, typename T>
	template<typename E, typename F>
	Vector<N, T> &Vector<N, T>::operator+=(Meta::Expression<N, T, E, F> const &e) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Add the Expression@" << &e << "===" << std::endl;
		#endif
		for(int i=0; i<N; i+=4) {
			v_[i] += e(i);
			v_[i+1] += e(i+1);
			v_[i+2] += e(i+2);
			v_[i+3] += e(i+3);
		}
		switch (N % 4) {
			case 3: v_[N-3] += e(N-3);
			case 2: v_[N-2] += e(N-2);
			case 1: v_[N-1] += e(N-1);
		}
		return *this;
	}

	/** The operator -= to subtract a value to all components
	 *  @param a The value to add **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator-=(T const &a) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Subtract the value " << a << "===" << std::endl;
		#endif
		if (a != 0.) {
			for(int i=0; i<N; i+=4) {
				v_[i] -= a;
				v_[i+1] -= a;
				v_[i+2] -= a;
				v_[i+3] -= a;
			}
			switch (N % 4) {
				case 3: v_[N-3] -= a;
				case 2: v_[N-2] -= a;
				case 1: v_[N-1] -= a;
			}
		}
		return *this;
	}

	/** The operator -= to subtract a vector
	 *  @param V The vector to subtract **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator-=(Vector<N, T> const &v) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Subtract the Vector@" << &v << "===" << std::endl;
		#endif
		for(int i=0; i<N; i+=4) {
			v_[i] -= v.v_[i];
			v_[i+1] -= v.v_[i+1];
			v_[i+2] -= v.v_[i+2];
			v_[i+3] -= v.v_[i+3];
		}
		switch (N % 4) {
			case 3: v_[N-3] -= v.v_[N-3];
			case 2: v_[N-2] -= v.v_[N-2];
			case 1: v_[N-1] -= v.v_[N-1];
		}
		return *this;
	}

	/** The operator -= to subtract an expression
	 *  @param e An expression **/
	template<size_t N, typename T>
	template<typename E, typename F>
	Vector<N, T> &Vector<N, T>::operator-=(Meta::Expression<N, T, E, F> const &e) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Subtract the Expression@" << &e << "===" << std::endl;
		#endif
		for(int i=0; i<N; i+=4) {
			v_[i] -= e(i);
			v_[i+1] -= e(i+1);
			v_[i+2] -= e(i+2);
			v_[i+3] -= e(i+3);
		}
		switch (N % 4) {
			case 3: v_[N-3] -= e(N-3);
			case 2: v_[N-2] -= e(N-2);
			case 1: v_[N-1] -= e(N-1);
		}
		return *this;
	}

	/** The operator *= to multiply a vector by a scalar value
	 *  @param a The scalar value to multiply **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator*=(T const &a) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Multiply by value" << a << "===" << std::endl;
		#endif
		if (a != 1.) {
			for(int i=0; i<N; i+=4) {
				v_[i] *= a;
				v_[i+1] *= a;
				v_[i+2] *= a;
				v_[i+3] *= a;
			}
			switch (N % 4) {
				case 3: v_[N-3] *= a;
				case 2: v_[N-2] *= a;
				case 1: v_[N-1] *= a;
			}
		}
		return *this;
	}

	/** The operator /= to divide a vector by a scalar value
	 *  @param a The scalar value to divide **/
	template<size_t N, typename T>
	Vector<N, T> &Vector<N, T>::operator/=(T const &a) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra::Vector@" << this << ": ===Divide by value" << a << "===" << std::endl;
		#endif
		if ((a != 1.) and (a != 0)) {
			for(int i=0; i<N; i+=4) {
				v_[i] /= a;
				v_[i+1] /= a;
				v_[i+2] /= a;
				v_[i+3] /= a;
			}
			switch (N % 4) {
				case 3: v_[N-3] /= a;
				case 2: v_[N-2] /= a;
				case 1: v_[N-1] /= a;
			}
		}
		return *this;
	}

	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T> &l, Vector<N, T> &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Vector@" << &l <<" + Vector@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Sum<T> > >
	operator+(T const &l, Vector<N, T> &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===" << l <<" + Vector@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T> &l, T const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Vector@" << &l <<" + " << r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Vector<N, T>, Common::Meta::Math::Sum<T> > >
	operator+(Meta::Expression<N, T, E, G> const &l, Vector<N, T> &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Expression@" << &l <<" + Vector@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Vector<N, T>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T> &l, Meta::Expression<N, T, E, G> const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Vector@" << &l <<" + Expression@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T>, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename F, typename G, typename H>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Meta::Expression<N, T, F, H>, Common::Meta::Math::Sum<T> > >
	operator+(Meta::Expression<N, T, E, G> const &l, Meta::Expression<N, T, F, H> const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Expression@" << &l <<" + Expression@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Meta::Expression<N, T, F, H>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Sum<T> > >
	operator+(Meta::Expression<N, T, E, G> const &l, T const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Expression@" << &l <<" + " << r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sum<T> > >
	operator+(T const &l, Meta::Expression<N, T, E, G> const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===" << l <<" + Expression@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T> &l, Vector<N, T> &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Vector@" << &l <<" - Vector@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Sub<T> > >
	operator-(T const &l, Vector<N, T> &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===" << l <<" - Vector@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T> &l, T const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Vector@" << &l <<" - " << r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Vector<N, T>, Common::Meta::Math::Sub<T> > >
	operator-(Meta::Expression<N, T, E, G> const &l, Vector<N, T> &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Expression@" << &l <<" - Vector@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Vector<N, T>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T> &l, Meta::Expression<N, T, E, G> const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Vector@" << &l <<" - Expression@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T>, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename F, typename G, typename H>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Meta::Expression<N, T, F, H>, Common::Meta::Math::Sub<T> > >
	operator-(Meta::Expression<N, T, E, G> const &l, Meta::Expression<N, T, F, H> const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Expression@" << &l <<" - Expression@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, Meta::Expression<N, T, F, H>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Sub<T> > >
	operator-(Meta::Expression<N, T, E, G> const &l, T const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Expression@" << &l <<" - " << r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sub<T> > >
	operator-(T const &l, Meta::Expression<N, T, E, G> const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===" << l <<" - Expression@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}

	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Mul<T> > >
	operator*(T const &l, Vector<N, T> &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===" << l <<" * Vector@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Vector<N, T>, Common::Meta::Math::Mul<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Mul<T> > >
	operator*(Vector<N, T> &l, T const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Vector@" << &l << " * " << r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Mul<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Mul<T> > >
	operator*(Meta::Expression<N, T, E, G> const &l, T const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Expression@" << &l << " * " << r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Mul<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Mul<T> > >
	operator*(T const &l, Meta::Expression<N, T, E, G> const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===" << l <<" * Expression@" << &r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Meta::Expression<N, T, E, G>, Common::Meta::Math::Mul<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}

	template<size_t N, typename T>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Div<T> > >
	operator/(Vector<N, T> &l, T const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Vector@" << &l << " / " << r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T>, T, Common::Meta::Math::Div<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename G>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Div<T> > >
	operator/(Meta::Expression<N, T, E, G> const &l, T const &r) {
		#ifdef _GAS_BE_VERY_VERBOSE_
		std::cerr << "LinearAlgebra::Meta: ===Expression@" << &l << " / " << r << "===" << std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, G>, T, Common::Meta::Math::Div<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}

	/** The scalar product between two vector
	 *  @param v The first vector
	 *  @param w The second vector **/
	template<size_t N, typename T>
	T operator*(Vector<N, T> const &v, Vector<N, T> const &w) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra: ===Vector@" << &v << " * Vector@" << &w << "===" << std::endl;
		#endif
		T r = T();
		for(int i=0; i<N; i+=4) {
			r += (v(i) * Common::Math::Conj(w(i)));
			r += (v(i+1) * Common::Math::Conj(w(i+1)));
			r += (v(i+2) * Common::Math::Conj(w(i+2)));
			r += (v(i+3) * Common::Math::Conj(w(i+3)));
		}
		switch (N % 4) {
			case 3: r += (v(N-3) * Common::Math::Conj(w(N-3)));
			case 2: r += (v(N-2) * Common::Math::Conj(w(N-2)));
			case 1: r += (v(N-1) * Common::Math::Conj(w(N-1)));
		}
		return r;
	}
	template<size_t N, typename T, typename E, typename F>
	T operator*(Vector<N, T> const &v, Meta::Expression<N, T, E, F> const &w) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra: ===Vector@" << &v << " * Expression@" << &w << "===" << std::endl;
		#endif
		T r = T();
		for(int i=0; i<N; i+=4) {
			r += (v(i) * Common::Math::Conj(w(i)));
			r += (v(i+1) * Common::Math::Conj(w(i+1)));
			r += (v(i+2) * Common::Math::Conj(w(i+2)));
			r += (v(i+3) * Common::Math::Conj(w(i+3)));
		}
		switch (N % 4) {
			case 3: r += (v(N-3) * Common::Math::Conj(w(N-3)));
			case 2: r += (v(N-2) * Common::Math::Conj(w(N-2)));
			case 1: r += (v(N-1) * Common::Math::Conj(w(N-1)));
		}
		return r;
	}
	template<size_t N, typename T, typename E, typename F> 
	T operator*(Meta::Expression<N, T, E, F> const &v, Vector<N, T> const &w) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra: ===Expression@" << &v << " * Vector@" << &w << "===" << std::endl;
		#endif
		T r = T();
		for(int i=0; i<N; i+=4) {
			r += (v(i) * Common::Math::Conj(w(i)));
			r += (v(i+1) * Common::Math::Conj(w(i+1)));
			r += (v(i+2) * Common::Math::Conj(w(i+2)));
			r += (v(i+3) * Common::Math::Conj(w(i+3)));
		}
		switch (N % 4) {
			case 3: r += (v(N-3) * Common::Math::Conj(w(N-3)));
			case 2: r += (v(N-2) * Common::Math::Conj(w(N-2)));
			case 1: r += (v(N-1) * Common::Math::Conj(w(N-1)));
		}
		return r;
	}
	template<size_t N, typename T, typename E, typename F, typename G, typename H> 
	T operator*(Meta::Expression<N, T, E, G> const &v, Meta::Expression<N, T, F, H> const &w) {
		#ifdef _GAS_BE_VERBOSE_
		std::cerr << "LinearAlgebra: ===Expression@" << &v << " * Expression@" << &w << "===" << std::endl;
		#endif
		T r = T();
		for(int i=0; i<N; i+=4) {
			r += (v(i) * Common::Math::Conj(w(i)));
			r += (v(i+1) * Common::Math::Conj(w(i+1)));
			r += (v(i+2) * Common::Math::Conj(w(i+2)));
			r += (v(i+3) * Common::Math::Conj(w(i+3)));
		}
		switch (N % 4) {
			case 3: r += (v(N-3) * Common::Math::Conj(w(N-3)));
			case 2: r += (v(N-2) * Common::Math::Conj(w(N-2)));
			case 1: r += (v(N-1) * Common::Math::Conj(w(N-1)));
		}
		return r;
	}

	/** The operator << to print a vector in a stream
	 *  @param s The stream
	 *  @param v The vector **/
	template<size_t N, typename T>
	std::ostream &operator<<(std::ostream &s, Vector<N, T> const &v) {
		s << "[";
		if (N <= 5)
			range(i, 0, N) s << " " << v(i);
		else
			s << " " << v(0) << " " << v(1) << " ... " << v(N-2) << " " << v(N-1);
		s << " ]";
		return s;
	}
	/** The operator << to print a vector expression in a stream
	 *  @param s The stream
	 *  @param e The expression **/
	template<size_t N, typename T, typename E, typename G>
	std::ostream &operator<<(std::ostream &s,Meta::Expression<N, T, E, G> const &e) {
		s << "[";
		if (N <= 5)
			range(i, 0, N) s << " " << e(i);
		else
			s << " " << e(0) << " " << e(1) << " ... " << e(N-2) << " " << e(N-1);
		s << " ]";
		return s;
	}

	namespace Meta {
		/* Expression<N, T, E, F> */
		template<size_t N, typename T, typename E, typename F>
		Expression<N, T, E, F>::Expression(E const &e): e_(e) {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta:Expression@" << this << ": ===Create an Expression from Object@" << &e << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename E, typename F>
		T const Expression<N, T, E, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::Expression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(e_(i));
		}

		/* Expression<N, T, Vector<N, T>, F> */
		template<size_t N, typename T, typename F>
		Expression<N, T, Vector<N, T>, F>::Expression(Vector<N, T> &e): e_(e) {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::Expression@" << this << ": ===Create an Expression from Vector@" << &e << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename F>
		T const Expression<N, T, Vector<N, T>, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::Expression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(e_(i));
		}

		/* BinaryExpression<N, T, L, R, F> */
		template<size_t N, typename T, typename L, typename R, typename F>
		BinaryExpression<N, T, L, R, F>::BinaryExpression(L const &l, R const &r): l_(l), r_(r){
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Create a BinaryExpression from Object@" << &l << \
				" and Object@" << &r << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename L, typename R, typename F>
		T const BinaryExpression<N, T, L, R, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(l_(i), r_(i));
		}

		/* BinaryExpression<N, T, L, T, F> */
		template<size_t N, typename T, typename L, typename F>
		BinaryExpression<N, T, L, T, F>::BinaryExpression(L const &l, T const &r): l_(l), r_(r){
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Create a BinaryExpression from Object@" << &l << \
				" and " << r << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename L, typename F>
		T const BinaryExpression<N, T, L, T, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(l_(i), r_);
		}

		/* BinaryExpression<N, T, T, R, F> */
		template<size_t N, typename T, typename R, typename F>
		BinaryExpression<N, T, T, R, F>::BinaryExpression(T const &l, R const &r): l_(l), r_(r){
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Create a BinaryExpression from " << l << \
				" and Object@" << &r << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename R, typename F>
		T const BinaryExpression<N, T, T, R, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(l_, r_(i));
		}

		/* BinaryExpression<N, T, Vector<N, T>, R, F> */
		template<size_t N, typename T, typename R, typename F>
		BinaryExpression<N, T, Vector<N, T>, R, F>::BinaryExpression(Vector<N, T> &l, R const &r): l_(l), r_(r){
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Create a BinaryExpression from Vector@" << &l << \
				" and Object@" << &r << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename R, typename F>
		T const BinaryExpression<N, T, Vector<N, T>, R, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(l_(i), r_(i));
		}

		/* BinaryExpression<N, T, L, Vector<N, T>, F> */
		template<size_t N, typename T, typename L, typename F>
		BinaryExpression<N, T, L, Vector<N, T>, F>::BinaryExpression(L const &l, Vector<N, T> &r): l_(l), r_(r){
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Create a BinaryExpression from Object@" << &l << \
				" and Vector@" << &r << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename L, typename F>
		T const BinaryExpression<N, T, L, Vector<N, T>, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(l_(i), r_(i));
		}

		/* BinaryExpression<N, T, Vector<N, T>, T, F> */
		template<size_t N, typename T, typename F>
		BinaryExpression<N, T, Vector<N, T>, T, F>::BinaryExpression(Vector<N, T> &l, T const &r): l_(l), r_(r){
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Create a BinaryExpression from Vector@" << &l << \
				" and " << r << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename F>
		T const BinaryExpression<N, T, Vector<N, T>, T, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(l_(i), r_);
		}

		/* BinaryExpression<N, T, T, Vector<N, T>, F> */
		template<size_t N, typename T, typename F>
		BinaryExpression<N, T, T, Vector<N, T>, F>::BinaryExpression(T const &l, Vector<N, T> &r): l_(l), r_(r){
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Create a BinaryExpression from " << l << \
				" and Vector@" << &r << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename F>
		T const BinaryExpression<N, T, T, Vector<N, T>, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(l_, r_(i));
		}

		/* BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, F> */
		template<size_t N, typename T, typename F>
		BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, F>::BinaryExpression(Vector<N, T> &l, Vector<N, T> &r): l_(l), r_(r){
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Create a BinaryExpression from Vector@" << &l << \
				" and Vector@" << &r << "===" << std::endl;
			#endif
		}

		template<size_t N, typename T, typename F>
		T const BinaryExpression<N, T, Vector<N, T>, Vector<N, T>, F>::operator()(size_t const &i) const {
			#ifdef _GAS_BE_VERY_VERBOSE_
			std::cerr << "LinearAlgebra::Meta::BinaryExpression@" << this << ": ===Evaluate at " << i << "===" << std::endl;
			#endif
			return F::RET(l_(i), r_(i));
		}
	}
}

#endif /* END _GAS_VECTOR_H_ */
