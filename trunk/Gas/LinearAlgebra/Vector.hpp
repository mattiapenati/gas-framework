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

#ifndef _GAS_LINEARALGEBRA_VECTOR_
#define _GAS_LINEARALGEBRA_VECTOR_

namespace LinearAlgebra {

template <typename Type, unsigned int Dimension, typename Operand, typename Function=id<Type> >
struct VectorExpression {
	Operand const & o_;
	
	inline VectorExpression (Operand const & o): o_(o) {}
	inline Type const operator() (unsigned int const & i) const { return Function::RET( o_(i) ); }
};

template <typename Type ,unsigned int Dimension ,typename LeftOperand, typename RightOperand, typename Function>
struct VectorBinaryExpression {
	LeftOperand const & l_;
	RightOperand const & r_;
	
	inline VectorBinaryExpression (LeftOperand const & l, RightOperand const & r): l_(l) , r_(r) {}
	inline Type const operator() (unsigned int const & i) const { return Function::RET(l_(i) , r_(i)); }
};

/* specialization for scalar types */
template <typename Type,  unsigned int Dimension, typename RightOperand, typename Function>
struct VectorBinaryExpression<Type, Dimension, Type, RightOperand, Function> {
	Type const & l_;
	RightOperand const & r_;
	
	inline VectorBinaryExpression (Type const & l, RightOperand const & r): l_(l) , r_(r) {}
	inline Type const operator() (unsigned int const & i) const { return Function::RET(l_ , r_(i)); }
};
template <typename Type ,unsigned int Dimension, typename LeftOperand, typename Function>
struct VectorBinaryExpression<Type, Dimension, LeftOperand, Type, Function> {
	LeftOperand const & l_;
	Type const & r_;
	
	inline VectorBinaryExpression (LeftOperand const & l, Type const & r) : l_(l) , r_(r) {}
	inline Type const operator() (unsigned int const & i) const { return Function::RET(l_(i) , r_); }
};
/* specialization to write less code */
template <typename Type, unsigned int Dimension, typename LeftOperand, typename RightOperand, typename ExpressionFunction, typename BinaryFunction>
struct VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, LeftOperand, RightOperand, BinaryFunction>, ExpressionFunction> {
	VectorBinaryExpression<Type, Dimension, LeftOperand, RightOperand, BinaryFunction> const o_;
	
	inline VectorExpression (LeftOperand const & l, RightOperand const & r): o_(l , r) {}
	inline Type const operator() (unsigned int const & i) const { return ExpressionFunction::RET(o_(i)); }
};

/* main class */
template <typename Type, unsigned int Dimension>
class Vector {

	public:
		/* constructor */
		Vector ();
		Vector (Type const &);
		Vector (Vector<Type, Dimension> const &);
		template <typename Operand, typename Function> Vector (VectorExpression<Type, Dimension, Operand, Function> const &);
		
		
		/* destructor */
		~Vector ();
		
		/* assign */
		inline Vector<Type, Dimension> & operator= (Vector<Type, Dimension> const &);
		template <typename Operand, typename Function> inline Vector<Type, Dimension> & operator= (VectorExpression<Type, Dimension, Operand, Function> const &);
		
		/* access */
		inline Type & operator() (unsigned int const &);
		inline Type const & operator() (unsigned int const &) const;
	
	private:
		Type value_[Dimension];

};

/* default constructor */
template <typename Type ,unsigned int Dimension>
Vector<Type, Dimension>::Vector () {
}

/* constructor by value */
template <typename Type ,unsigned int Dimension>
Vector<Type, Dimension>::Vector (Type const & scalar) {
	for (unsigned int i = 0u; i < Dimension; ++i)
		value_[i] = scalar;
}

/* copy constructor */
template <typename Type, unsigned int Dimension>
Vector<Type, Dimension>::Vector (Vector<Type, Dimension> const & vector) {
	for (unsigned int i = 0u; i < Dimension; ++i)
		value_[i] = vector(i);
}

/* exrpession constructor */
template <typename Type, unsigned int Dimension>
template <typename Operand, typename Function>
Vector<Type, Dimension>::Vector (VectorExpression<Type, Dimension, Operand, Function> const & expression) {
	for (unsigned int i = 0u; i < Dimension; ++i)
		value_[i] = expression(i);
}

/* destructor */
template <typename Type, unsigned int Dimension>
Vector<Type, Dimension>::~Vector () {
}

/* assign from vector */
template <typename Type, unsigned int Dimension>
Vector<Type, Dimension> & Vector<Type, Dimension>::operator= (Vector<Type, Dimension> const & vector) {
	for (unsigned int i = 0u; i < Dimension; ++i)
		value_[i] = vector(i);
	return *this;
}

/* assign from exrpession */
template <typename Type, unsigned int Dimension>
template <typename Operand, typename Function>
Vector<Type, Dimension> & Vector<Type, Dimension>::operator= (VectorExpression<Type, Dimension, Operand, Function> const & expression) {
	for (unsigned int i = 0u; i < Dimension; ++i)
		value_[i] = expression(i);
	return *this;
}

/* access */
template <typename Type, unsigned int Dimension>
Type & Vector<Type, Dimension>::operator() (unsigned int const & i) {
	return value_[i];
}
template <typename Type, unsigned int Dimension>
Type const & Vector<Type, Dimension>::operator() (unsigned int const & i) const {
	return value_[i];
}

/* add */
template <typename Type, unsigned int Dimension>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, Vector<Type, Dimension>, add<Type> >, id<Type> >
operator+ (Vector<Type, Dimension> const & x, Vector<Type, Dimension> const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, Vector<Type, Dimension>, add<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, Operand, Function>, Vector<Type, Dimension>, add<Type> >, id<Type> >
operator+ (VectorExpression<Type, Dimension, Operand, Function> const & x, Vector<Type, Dimension> const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, Operand, Function>, Vector<Type, Dimension>, add<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, VectorExpression<Type, Dimension, Operand, Function>, add<Type> >, id<Type> >
operator+ (Vector<Type, Dimension> const & x, VectorExpression<Type, Dimension, Operand, Function> const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, VectorExpression<Type, Dimension, Operand, Function>, add<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Dimension, typename LeftOperand, typename LeftFunction, typename RightOperand, typename RightFunction>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, LeftOperand, LeftFunction>, VectorExpression<Type, Dimension, RightOperand, RightFunction>, add<Type> >, id<Type> >
operator+ (VectorExpression<Type, Dimension, LeftOperand, LeftFunction> const & x, VectorExpression<Type, Dimension, RightOperand, RightFunction> const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, LeftOperand, LeftFunction>, VectorExpression<Type, Dimension, RightOperand, RightFunction>, add<Type> >, id<Type> >(x, y);
}

/* sub */
template <typename Type, unsigned int Dimension>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, Vector<Type, Dimension>, sub<Type> >, id<Type> >
operator- (Vector<Type, Dimension> const & x, Vector<Type, Dimension> const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, Vector<Type, Dimension>, sub<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, Operand, Function>, Vector<Type, Dimension>, sub<Type> >, id<Type> >
operator- (VectorExpression<Type, Dimension, Operand, Function> const & x, Vector<Type, Dimension> const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, Operand, Function>, Vector<Type, Dimension>, sub<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, VectorExpression<Type, Dimension, Operand, Function>, sub<Type> >, id<Type> >
operator- (Vector<Type, Dimension> const & x, VectorExpression<Type, Dimension, Operand, Function> const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, VectorExpression<Type, Dimension, Operand, Function>, sub<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Dimension, typename LeftOperand, typename LeftFunction, typename RightOperand, typename RightFunction>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, LeftOperand, LeftFunction>, VectorExpression<Type, Dimension, RightOperand, RightFunction>, sub<Type> >, id<Type> >
operator- (VectorExpression<Type, Dimension, LeftOperand, LeftFunction> const & x, VectorExpression<Type, Dimension, RightOperand, RightFunction> const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, LeftOperand, LeftFunction>, VectorExpression<Type, Dimension, RightOperand, RightFunction>, sub<Type> >, id<Type> >(x, y);
}

/* mul */
template <typename Type, unsigned int Dimension>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, Type, mul<Type> >, id<Type> >
operator* (Vector<Type, Dimension> const & x, Type const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, Type, mul<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, Operand, Function>, Type, mul<Type> >, id<Type> >
operator* (VectorExpression<Type, Dimension, Operand, Function> const & x, Type const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, Operand, Function>, Type, mul<Type> >, id<Type> >(x, y);
}

/* div */
template <typename Type, unsigned int Dimension>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, Type, div<Type> >, id<Type> >
operator/ (Vector<Type, Dimension> const & x, Type const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, Vector<Type, Dimension>, Type, div<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, Operand, Function>, Type, div<Type> >, id<Type> >
operator/ (VectorExpression<Type, Dimension, Operand, Function> const & x, Type const & y) {
	return VectorExpression<Type, Dimension, VectorBinaryExpression<Type, Dimension, VectorExpression<Type, Dimension, Operand, Function>, Type, div<Type> >, id<Type> >(x, y);
}

/* dot product */
template <typename Type, unsigned int Dimension>
inline Type dot (Vector<Type, Dimension> const & v, Vector<Type, Dimension> const & w) {
	Type r = (v(0) * w(0));
	for (unsigned int i = 1u; i < Dimension; ++i)
		r += (v(i) * w(i));
	return r;
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline Type dot (VectorExpression<Type, Dimension, Operand, Function> const & v, Vector<Type, Dimension> const & w) {
	Type r = (v(0) * w(0));
	for (unsigned int i = 1u; i < Dimension; ++i)
		r += (v(i) * w(i));
	return r;
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline Type dot (Vector<Type, Dimension> const & v, VectorExpression<Type, Dimension, Operand, Function> const & w) {
	Type r = (v(0) * w(0));
	for (unsigned int i = 1u; i < Dimension; ++i)
		r += (v(i) * w(i));
	return r;
}
template <typename Type, unsigned int Dimension, typename LeftOperand, typename LeftFunction, typename RightOperand, typename RightFunction>
inline Type dot (VectorExpression<Type, Dimension, LeftOperand, LeftFunction> const & v, VectorExpression<Type, Dimension, RightOperand, RightFunction> const & w) {
	Type r = (v(0) * w(0));
	for (unsigned int i = 1u; i < Dimension; ++i)
		r += (v(i) * w(i));
	return r;
}

/* dot product */
template <typename Type, unsigned int Dimension>
inline Type operator* (Vector<Type, Dimension> const & v, Vector<Type, Dimension> const & w) {
	Type r = (v(0) * w(0));
	for (unsigned int i = 1u; i < Dimension; ++i)
		r += (v(i) * w(i));
	return r;
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline Type operator* (VectorExpression<Type, Dimension, Operand, Function> const & v, Vector<Type, Dimension> const & w) {
	Type r = (v(0) * w(0));
	for (unsigned int i = 1u; i < Dimension; ++i)
		r += (v(i) * w(i));
	return r;
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline Type operator* (Vector<Type, Dimension> const & v, VectorExpression<Type, Dimension, Operand, Function> const & w) {
	Type r = (v(0) * w(0));
	for (unsigned int i = 1u; i < Dimension; ++i)
		r += (v(i) * w(i));
	return r;
}
template <typename Type, unsigned int Dimension, typename LeftOperand, typename LeftFunction, typename RightOperand, typename RightFunction>
inline Type operator* (VectorExpression<Type, Dimension, LeftOperand, LeftFunction> const & v, VectorExpression<Type, Dimension, RightOperand, RightFunction> const & w) {
	Type r = (v(0) * w(0));
	for (unsigned int i = 1u; i < Dimension; ++i)
		r += (v(i) * w(i));
	return r;
}

/* norm */
template <typename Type, unsigned int Dimension>
inline Type norm (Vector<Type, Dimension> const & v) {
	Type r = dot(v, v); 
	return std::sqrt(r);
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline Type norm (VectorExpression<Type, Dimension, Operand, Function> const & v) {
	Type r = dot(v, v);
	return std::sqrt(r);
}

}

#endif // _GAS_LINEARALGEBRA_VECTOR_
