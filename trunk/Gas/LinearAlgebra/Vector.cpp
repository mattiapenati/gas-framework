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

	/* The default constructor */
	template<size_t N, typename T, typename C>
	Vector<N, T, C>::Vector() {
		#if _GAS_VERBOSITY_ >= 1
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Create an empty Vector==="<<std::endl;
		#endif
	}

	/* The constructor to initialize the entire vector with the same value */
	template<size_t N, typename T, typename C>
	Vector<N, T, C>::Vector(T const a) {
		#if _GAS_VERBOSITY_ >= 1
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Create a Vector from value "<<a<<"==="<<std::endl;
		#endif
		if (N >= 4) {
			for(size_t i = 0; i < N; i += 4) {
				v_[i] = a;
				v_[i+1] = a;
				v_[i+2] = a;
				v_[i+3] = a;
			}
		}
		switch (N % 4) {
			case 3: v_[N-3] = a;
			case 2: v_[N-2] = a;
			case 1: v_[N-1] = a;
		}
	}

	/* The constructor to copy a vector */
	template<size_t N, typename T, typename C>
	Vector<N, T, C>::Vector(Vector<N, T, C> const &v) {
		#if _GAS_VERBOSITY_ >= 1
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Create a Vector from Vector@"<<&v<<"==="<<std::endl;
		#endif
		v_ = v.v_;
	}

	/* The constructor to copy a vector */
	template<size_t N, typename T, typename C>
	template<typename D>
	Vector<N, T, C>::Vector(Vector<N, T, D> const &v) {
		#if _GAS_VERBOSITY_ >= 1
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Create a Vector from Vector@"<<&v<<"==="<<std::endl;
		#endif
		v_ = v.v_;
	}

	/* The constructor to copy a vector expression */
	template<size_t N, typename T, typename C>
	template<typename E, typename F>
	Vector<N, T, C>::Vector(Meta::Expression<N, T, E, F> const &e) {
		#if _GAS_VERBOSITY_ >= 1
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Create a Vector from Expression@"<<&e<<"==="<<std::endl;
		#endif
		if (N >= 4) {
			for(size_t i = 0; i < N; i += 4) {
				v_[i] = e(i);
				v_[i+1] = e(i+1);
				v_[i+2] = e(i+2);
				v_[i+3] = e(i+3);
			}
		}
		switch (N % 4) {
			case 3: v_[N-3] = e(N-3);
			case 2: v_[N-2] = e(N-2);
			case 1: v_[N-1] = e(N-1);
		}
	}

	/* The default destructor */
	template<size_t N, typename T, typename C>
	Vector<N, T, C>::~Vector() {
		#if _GAS_VERBOSITY_ >= 1
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Destroyed==="<<std::endl;
		#endif
	}

	/* The vector size */
	template<size_t N, typename T, typename C>
	size_t Vector<N, T, C>::Size() {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Reading the size==="<<std::endl;
		#endif
		return N;
	}

	/* The operator () to access to the components of vector */
	template<size_t N, typename T, typename C>
	T &Vector<N, T, C>::operator()(size_t const i) {
		#if _GAS_VERBOSITY_ >= 5
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Accessing to the element at position "<<i<<"==="<<std::endl;
		#endif
		#ifdef _GAS_CHECK_INDEX_
		assert(i < N);
		#endif
		return v_[i];
	}

	/* The operator () to access to the components of vector */
	template<size_t N, typename T, typename C>
	T const &Vector<N, T, C>::operator()(size_t const i) const {
		#if _GAS_VERBOSITY_ >= 5
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Accessing to the element at position "<<i<<"==="<<std::endl;
		#endif
		#ifdef _GAS_CHECK_INDEX_
		assert(i < N);
		#endif
		return v_[i];
	}

	/* The operator == to compare all components of a vector with a value */
	template<size_t N, typename T, typename C>
	bool Vector<N, T, C>::operator==(T const &a) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Compare with the value "<<a<<"==="<<std::endl;
		#endif
		range(i, 0, N) { if (!Common::Limits<T>::Equal(v_[i], a)) return false; }
		return true;
	}

	/* The operator == to compare two vectors */
	template<size_t N, typename T, typename C>
	template<typename D>
	bool Vector<N, T, C>::operator==(Vector<N, T, D> const &v) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Compare with the Vector@"<<&v<<"==="<<std::endl;
		#endif
		if (this != &v) {
			range(i, 0, N) { if (!Common::Limits<T>::Equal(v_[i], v.v_[i])) return false; }
		}
		return true;
	}

	/* The operator == to compare a vector with an expression */
	template<size_t N, typename T, typename C>
	template<typename E, typename F>
	bool Vector<N, T, C>::operator==(Meta::Expression<N, T, E, F> const &e) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Compare with the Expression@"<<&e<<"==="<<std::endl;
		#endif
		range(i, 0, N) { if (!Common::Limits<T>::Equal(v_[i], e(i))) return false; }
		return true;
	}

	/* The operator = to copy a value in all components of a vector */
	template<size_t N, typename T, typename C>
	Vector<N, T, C> &Vector<N, T, C>::operator=(T const &a) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Copy the value "<<a<<"==="<<std::endl;
		#endif
		if (N >= 4) {
			for(int i=0; i<N; i+=4) {
				v_[i] = a;
				v_[i+1] = a;
				v_[i+2] = a;
				v_[i+3] = a;
			}
		}
		switch (N % 4) {
			case 3: v_[N-3] = a;
			case 2: v_[N-2] = a;
			case 1: v_[N-1] = a;
		}
		return *this;
	}

	/* The operator = to copy a vector in a vector */
	template<size_t N, typename T, typename C>
	template<typename D>
	Vector<N, T, C> &Vector<N, T, C>::operator=(Vector<N, T, D> const &v) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Copy the Vector@"<<&v<<"==="<<std::endl;
		#endif
		if (this != &v) v_ = v.v_;
		return *this;
	}

	/* The operator = to copy an expression in a vector */
	template<size_t N, typename T, typename C>
	template<typename E, typename F>
	Vector<N, T, C> &Vector<N, T, C>::operator=(Meta::Expression<N, T, E, F> const &e) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Copy the Expression@"<<&e<<"==="<<std::endl;
		#endif
		if (N >= 4) {
			for(int i=0; i<N; i+=4) {
				v_[i] = e(i);
				v_[i+1] = e(i+1);
				v_[i+2] = e(i+2);
				v_[i+3] = e(i+3);
			}
		}
		switch (N % 4) {
			case 3: v_[N-3] = e(N-3);
			case 2: v_[N-2] = e(N-2);
			case 1: v_[N-1] = e(N-1);
		}
		return *this;
	}

	/* The operator += to add a value to all components */
	template<size_t N, typename T, typename C>
	Vector<N, T, C> &Vector<N, T, C>::operator+=(T const &a) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Add the value "<<a<<"==="<<std::endl;
		#endif
		if (a != 0.) {
			if (N >= 4) {
				for(int i=0; i<N; i+=4) {
					v_[i] += a;
					v_[i+1] += a;
					v_[i+2] += a;
					v_[i+3] += a;
				}
			}
			switch (N % 4) {
				case 3: v_[N-3] += a;
				case 2: v_[N-2] += a;
				case 1: v_[N-1] += a;
			}
		}
		return *this;
	}

	/* The operator += to add a vector */
	template<size_t N, typename T, typename C>
	template<typename D>
	Vector<N, T, C> &Vector<N, T, C>::operator+=(Vector<N, T, D> const &v) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Add the Vector@"<<&v<<"==="<<std::endl;
		#endif
		if (N >= 4) {
			for(int i=0; i<N; i+=4) {
				v_[i] += v.v_[i];
				v_[i+1] += v.v_[i+1];
				v_[i+2] += v.v_[i+2];
				v_[i+3] += v.v_[i+3];
			}
		}
		switch (N % 4) {
			case 3: v_[N-3] += v.v_[N-3];
			case 2: v_[N-2] += v.v_[N-2];
			case 1: v_[N-1] += v.v_[N-1];
		}
		return *this;
	}

	/* The operator += to add an expression */
	template<size_t N, typename T, typename C>
	template<typename E, typename F>
	Vector<N, T, C> &Vector<N, T, C>::operator+=(Meta::Expression<N, T, E, F> const &e) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Add the Expression@"<<&e<<"==="<<std::endl;
		#endif
		if (N >= 4) {
			for(int i=0; i<N; i+=4) {
				v_[i] += e(i);
				v_[i+1] += e(i+1);
				v_[i+2] += e(i+2);
				v_[i+3] += e(i+3);
			}
		}
		switch (N % 4) {
			case 3: v_[N-3] += e(N-3);
			case 2: v_[N-2] += e(N-2);
			case 1: v_[N-1] += e(N-1);
		}
		return *this;
	}

	/* The operator -= to subtract a value to all components */
	template<size_t N, typename T, typename C>
	Vector<N, T, C> &Vector<N, T, C>::operator-=(T const &a) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Subtract the value "<<a<<"==="<<std::endl;
		#endif
		if (a != 0.) {
			if (N >= 4) {
				for(int i=0; i<N; i+=4) {
					v_[i] -= a;
					v_[i+1] -= a;
					v_[i+2] -= a;
					v_[i+3] -= a;
				}
			}
			switch (N % 4) {
				case 3: v_[N-3] -= a;
				case 2: v_[N-2] -= a;
				case 1: v_[N-1] -= a;
			}
		}
		return *this;
	}

	/* The operator -= to subtract a vector */
	template<size_t N, typename T, typename C>
	template<typename D>
	Vector<N, T, C> &Vector<N, T, C>::operator-=(Vector<N, T, D> const &v) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Subtract the Vector@"<<&v<<"==="<<std::endl;
		#endif
		if (N >= 4) {
			for(int i = 0; i < N; i += 4) {
				v_[i] -= v.v_[i];
				v_[i+1] -= v.v_[i+1];
				v_[i+2] -= v.v_[i+2];
				v_[i+3] -= v.v_[i+3];
			}
		}
		switch (N % 4) {
			case 3: v_[N-3] -= v.v_[N-3];
			case 2: v_[N-2] -= v.v_[N-2];
			case 1: v_[N-1] -= v.v_[N-1];
		}
		return *this;
	}

	/* The operator -= to subtract an expression */
	template<size_t N, typename T, typename C>
	template<typename E, typename F>
	Vector<N, T, C> &Vector<N, T, C>::operator-=(Meta::Expression<N, T, E, F> const &e) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Subtract the Expression@"<<&e<<"==="<<std::endl;
		#endif
		if (N >= 4) {
			for(int i=0; i<N; i+=4) {
				v_[i] -= e(i);
				v_[i+1] -= e(i+1);
				v_[i+2] -= e(i+2);
				v_[i+3] -= e(i+3);
			}
		}
		switch (N % 4) {
			case 3: v_[N-3] -= e(N-3);
			case 2: v_[N-2] -= e(N-2);
			case 1: v_[N-1] -= e(N-1);
		}
		return *this;
	}

	/* The operator *= to multiply a vector by a scalar value */
	template<size_t N, typename T, typename C>
	Vector<N, T, C> &Vector<N, T, C>::operator*=(T const &a) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Multiply by value"<<a<<"==="<<std::endl;
		#endif
		if (a != 1.) {
			if (N >= 4) {
				for(int i=0; i<N; i+=4) {
					v_[i] *= a;
					v_[i+1] *= a;
					v_[i+2] *= a;
					v_[i+3] *= a;
				}
			}
			switch (N % 4) {
				case 3: v_[N-3] *= a;
				case 2: v_[N-2] *= a;
				case 1: v_[N-1] *= a;
			}
		}
		return *this;
	}

	/* The operator /= to divide a vector by a scalar value */
	template<size_t N, typename T, typename C>
	Vector<N, T, C> &Vector<N, T, C>::operator/=(T const &a) {
		#if _GAS_VERBOSITY_ >= 2
		std::cerr<<"LinearAlgebra::Vector@"<<this<<": ===Divide by value"<<a<<"==="<<std::endl;
		#endif
		if ((a != 1.) and (a != 0)) {
			if (N >= 4) {
				for(int i=0; i<N; i+=4) {
					v_[i] /= a;
					v_[i+1] /= a;
					v_[i+2] /= a;
					v_[i+3] /= a;
				}
			}
			switch (N % 4) {
				case 3: v_[N-3] /= a;
				case 2: v_[N-2] /= a;
				case 1: v_[N-1] /= a;
			}
		}
		return *this;
	}

	template<size_t N, typename T, typename C, typename D>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T, C> &l, Vector<N, T, D> &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Vector@"<<&l<<" + Vector@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Sum<T> > >
	operator+(T const &l, Vector<N, T, C> &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ==="<<l<<" + Vector@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T, C> &l, T const &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Vector@"<<&l<<" + "<<r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename F, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, F>, Vector<N, T, C>, Common::Meta::Math::Sum<T> > >
	operator+(Meta::Expression<N, T, E, F> const &l, Vector<N, T, C> &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Expression@"<<&l<<" + Vector@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, F>, Vector<N, T, C>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename C, typename E, typename F>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, Meta::Expression<N, T, E, F>, Common::Meta::Math::Sum<T> > >
	operator+(Vector<N, T, C> &l, Meta::Expression<N, T, E, F> const &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Vector@"<<&l<<" + Expression@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T, C>, Meta::Expression<N, T, E, F>, Common::Meta::Math::Sum<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	
	template<size_t N, typename T, typename C, typename D>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T, C> &l, Vector<N, T, D> &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Vector@"<<&l<<" - Vector@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Sub<T> > >
	operator-(T const &l, Vector<N, T, C> &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ==="<<l<<" - Vector@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T, C> &l, T const &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Vector@"<<&l<<" - "<<r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename E, typename F, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, F>, Vector<N, T, C>, Common::Meta::Math::Sub<T> > >
	operator-(Meta::Expression<N, T, E, F> const &l, Vector<N, T, C> &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Expression@"<<&l<<" - Vector@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Meta::Expression<N, T, E, F>, Vector<N, T, C>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename C, typename E, typename F>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, Meta::Expression<N, T, E, F>, Common::Meta::Math::Sub<T> > >
	operator-(Vector<N, T, C> &l, Meta::Expression<N, T, E, F> const &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Vector@"<<&l<<" - Expression@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T, C>, Meta::Expression<N, T, E, F>, Common::Meta::Math::Sub<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}

	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Mul<T> > >
	operator*(T const &l, Vector<N, T, C> &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ==="<<l<<" * Vector@"<<&r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, T, Vector<N, T, C>, Common::Meta::Math::Mul<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}
	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Mul<T> > >
	operator*(Vector<N, T, C> &l, T const &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Vector@"<<&l<<" * "<<r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Mul<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}

	template<size_t N, typename T, typename C>
	Meta::Expression<N, T, Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Div<T> > >
	operator/(Vector<N, T, C> &l, T const &r) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra::Meta: ===Vector@"<<&l<<" / "<<r<<"==="<<std::endl;
		#endif
		typedef Meta::BinaryExpression<N, T, Vector<N, T, C>, T, Common::Meta::Math::Div<T> > TempType;
		return Meta::Expression<N, T, TempType>(TempType(l, r));
	}

	/* The scalar product between two vector */
	template<size_t N, typename T, typename C, typename D>
	T operator*(Vector<N, T, C> const &v, Vector<N, T, D> const &w) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra: ===Vector@"<<&v<<" * Vector@"<<&w<<"==="<<std::endl;
		#endif
		T r = T();
		if (N >= 4) {
			for(int i=0; i<N; i+=4) {
				r += (v(i) * Common::Math::Conj(w(i)));
				r += (v(i+1) * Common::Math::Conj(w(i+1)));
				r += (v(i+2) * Common::Math::Conj(w(i+2)));
				r += (v(i+3) * Common::Math::Conj(w(i+3)));
			}
		}
		switch (N % 4) {
			case 3: r += (v(N-3) * Common::Math::Conj(w(N-3)));
			case 2: r += (v(N-2) * Common::Math::Conj(w(N-2)));
			case 1: r += (v(N-1) * Common::Math::Conj(w(N-1)));
		}
		return r;
	}
	template<size_t N, typename T, typename C, typename E, typename F>
	T operator*(Vector<N, T, C> const &v, Meta::Expression<N, T, E, F> const &w) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra: ===Vector@"<<&v<<" * Expression@"<<&w<<"==="<<std::endl;
		#endif
		T r = T();
		if (N >= 4) {
			for(int i=0; i<N; i+=4) {
				r += (v(i) * Common::Math::Conj(w(i)));
				r += (v(i+1) * Common::Math::Conj(w(i+1)));
				r += (v(i+2) * Common::Math::Conj(w(i+2)));
				r += (v(i+3) * Common::Math::Conj(w(i+3)));
			}
		}
		switch (N % 4) {
			case 3: r += (v(N-3) * Common::Math::Conj(w(N-3)));
			case 2: r += (v(N-2) * Common::Math::Conj(w(N-2)));
			case 1: r += (v(N-1) * Common::Math::Conj(w(N-1)));
		}
		return r;
	}
	template<size_t N, typename T, typename E, typename F, typename C> 
	T operator*(Meta::Expression<N, T, E, F> const &v, Vector<N, T, C> const &w) {
		#if _GAS_VERBOSITY_ >= 3
		std::cerr<<"LinearAlgebra: ===Expression@"<<&v<<" * Vector@"<<&w<<"==="<<std::endl;
		#endif
		T r = T();
		if (N >= 4) {
			for(int i=0; i<N; i+=4) {
				r += (v(i) * Common::Math::Conj(w(i)));
				r += (v(i+1) * Common::Math::Conj(w(i+1)));
				r += (v(i+2) * Common::Math::Conj(w(i+2)));
				r += (v(i+3) * Common::Math::Conj(w(i+3)));
			}
		}
		switch (N % 4) {
			case 3: r += (v(N-3) * Common::Math::Conj(w(N-3)));
			case 2: r += (v(N-2) * Common::Math::Conj(w(N-2)));
			case 1: r += (v(N-1) * Common::Math::Conj(w(N-1)));
		}
		return r;
	}
	
	namespace Meta {
		/* Expression<N, T, Vector<N, T, C>, F> */
		template<size_t N, typename T, typename F, typename C>
		Expression<N, T, Vector<N, T, C>, F>::Expression(Vector<N, T, C> const &e): e_(e) {
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::Expression@"<<this<<": ===Create an Expression from Vector@"<<&e<<"==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename F, typename C>
		Expression<N, T, Vector<N, T, C>, F>::~Expression() {
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::Expression@"<<this<<": ===Destroyed==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename F, typename C>
		T const Expression<N, T, Vector<N, T, C>, F>::operator()(size_t const &i) const {
			#if _GAS_VERBOSITY_ >= 5
			std::cerr<<"LinearAlgebra::Meta::Expression@"<<this<<": ===Evaluate at "<<i<<"==="<<std::endl;
			#endif
			return F::RET(e_(i));
		}
		/* BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, F> */
		template<size_t N, typename T, typename C, typename D, typename F>
		BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, F>::BinaryExpression(Vector<N, T, C> const &l, Vector<N, T, D> const &r): l_(l), r_(r){
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Create a BinaryExpression from Vector@"<<&l<<" and Vector@"<<&r<<"==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename C, typename D, typename F>
		BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, F>::~BinaryExpression() {
			#if _GAS_VERBOSITY_ >= 3
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Destroyed==="<<std::endl;
			#endif
		}
		
		template<size_t N, typename T, typename C, typename D, typename F>
		T const BinaryExpression<N, T, Vector<N, T, C>, Vector<N, T, D>, F>::operator()(size_t const &i) const {
			#if _GAS_VERBOSITY_ >= 5
			std::cerr<<"LinearAlgebra::Meta::BinaryExpression@"<<this<<": ===Evaluate at "<<i<<"==="<<std::endl;
			#endif
			return F::RET(l_(i), r_(i));
		}
		// TODO Completare il travaso di codice
	}
}
