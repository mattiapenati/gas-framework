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

#ifndef _GAS_LINEARALGEBRA_MATRIX_
#define _GAS_LINEARALGEBRA_MATRIX_

namespace LinearAlgebra {

template <typename Type, unsigned int Row, unsigned int Column, typename Operand, typename Function=id<Type> >
struct MatrixExpression {
	Operand const & o_;
	
	inline MatrixExpression (Operand const & o) : o_(o) {}
	inline Type const operator() (unsigned int const & i, unsigned int const & j) const { return Function::RET(o_(i,j)); }
};

template <typename Type, unsigned int Row, unsigned int Column, typename LeftOperand, typename RightOperand, typename Function>
struct MatrixBinaryExpression {
	LeftOperand const & l_;
	RightOperand const & r_;
	
	inline MatrixBinaryExpression (LeftOperand const & l, RightOperand const & r) : l_(l), r_(r) {}
	inline Type const operator() (unsigned int const & i, unsigned int const & j) const { return Function::RET(l_(i,j), r_(i,j)); }
};

/* specialization for scalar types */
template <typename Type, unsigned int Row, unsigned int Column, typename RightOperand, typename Function>
struct MatrixBinaryExpression<Type, Row, Column, Type, RightOperand, Function> {
	Type const & l_;
	RightOperand const & r_;
	
	inline MatrixBinaryExpression (Type const & l, RightOperand const & r) : l_(l), r_(r) {}
	inline Type const operator() (unsigned int const & i, unsigned int const & j) const { return Function::RET(l_, r_(i,j)); }
};
template <typename Type, unsigned int Row, unsigned int Column, typename LeftOperand, typename Function>
struct MatrixBinaryExpression<Type, Row, Column, LeftOperand, Type, Function> {
	LeftOperand const & l_;
	Type const & r_;
	
	inline MatrixBinaryExpression (LeftOperand const & l, Type const & r) : l_(l), r_(r) {}
	inline Type const operator() (unsigned int const & i, unsigned int const & j) const { return Function::RET(l_(i,j), r_); }
};
/* specialization to write less code */
template <typename Type, unsigned int Row, unsigned int Column, typename LeftOperand, typename RightOperand, typename ExpressionFunction, typename BinaryFunction>
struct MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, LeftOperand, RightOperand, BinaryFunction>, ExpressionFunction> {
	MatrixBinaryExpression<Type, Row, Column, LeftOperand, RightOperand, BinaryFunction> const o_;
	
	inline MatrixExpression (LeftOperand const & l, RightOperand const & r) : o_(l, r) {}
	inline Type const operator() (unsigned int const & i, unsigned int const & j) const { return ExpressionFunction::RET(o_(i,j)); }
};

/* main class */
template <typename Type, unsigned int Row, unsigned int Column>
class Matrix {

	public:
		/* constructor */
		Matrix ();
		Matrix (Type const &);
		Matrix (Matrix <Type, Row, Column> const &);
		
		/* destructor */
		~Matrix ();
		
		/* assign */
		inline Matrix<Type, Row, Column> & operator= (Matrix<Type, Row, Column> const &);
		template <typename Operand, typename Function> inline Matrix<Type, Row, Column> & operator= (MatrixExpression<Type, Row, Column, Operand, Function> const &);
		
		/* access */
		inline Type & operator() (unsigned int const &, unsigned int const &);
		inline Type const & operator() (unsigned int const &, unsigned int const &) const ;
	
	private:
		Type value_[Row][Column];

};

/* default constructor */
template <typename Type, unsigned int Row, unsigned int Column>
Matrix<Type, Row, Column>::Matrix () {
}

/* constructor by value */
template <typename Type, unsigned int Row, unsigned int Column>
Matrix<Type, Row, Column>::Matrix (Type const & scalar) {
	for (unsigned int i = 0u; i <Row; ++i) {
		for (unsigned int j = 0u; j <Column; ++j) {
			value_[i][j] = scalar;
		}
	}
}
/* copy constructor */
template <typename Type, unsigned int Row, unsigned int Column>
Matrix<Type, Row, Column>::Matrix (Matrix <Type, Row, Column> const & matrix) {
	for (unsigned int i = 0u; i <Row; ++i) {
		for (unsigned int j = 0u; j <Column; ++j) {
			value_[i][j] = matrix(i,j);
		}
	}
}

/* assign from matrix */
template <typename Type, unsigned int Row, unsigned int Column>
Matrix<Type, Row, Column> & Matrix<Type, Row, Column>::operator= (Matrix<Type, Row, Column> const & matrix) {
	for (unsigned int i = 0u; i <Row; ++i) {
		for (unsigned int j = 0u; j <Column; ++j) {
			value_[i][j] = matrix(i,j);
		}
	}
	return *this;
}

/* assign from exrpession */
template <typename Type, unsigned int Row, unsigned int Column>
template <typename Operand, typename Function>
Matrix<Type, Row, Column> & Matrix<Type, Row, Column>::operator= (MatrixExpression<Type, Row, Column, Operand, Function> const & expression) {
	for (unsigned int i = 0u; i <Row; ++i) {
		for (unsigned int j = 0u; j <Column; ++j) {
			value_[i][j] = expression(i,j);
		}
	}
	return *this;
}

/* destructor */
template <typename Type, unsigned int Row, unsigned int Column>
Matrix<Type, Row, Column>::~Matrix () {
}

/* access */
template <typename Type, unsigned int Row, unsigned int Column>
Type & Matrix<Type, Row, Column>::operator() (unsigned int const & i, unsigned int const & j) {
	return value_[i][j];
}
template <typename Type, unsigned int Row, unsigned int Column>
Type const & Matrix<Type, Row, Column>::operator() (unsigned int const & i, unsigned int const & j) const {
	return value_[i][j];
}

/* add */
template <typename Type, unsigned int Row, unsigned int Column>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Matrix<Type, Row, Column>, add<Type> >, id<Type> >
operator+ (Matrix<Type, Row, Column> const & x, Matrix<Type, Row, Column> const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Matrix<Type, Row, Column>, add<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Row, unsigned int Column, typename Operand, typename Function>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, Operand, Function>, Matrix<Type, Row, Column>, add<Type> >, id<Type> >
operator+ (MatrixExpression<Type, Row, Column, Operand, Function> const & x, Matrix<Type, Row, Column> const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, Operand, Function>, Matrix<Type, Row, Column>, add<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Row, unsigned int Column, typename Operand, typename Function>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, MatrixExpression<Type, Row, Column, Operand, Function>, add<Type> >, id<Type> >
operator+ (Matrix<Type, Row, Column> const & x, MatrixExpression<Type, Row, Column, Operand, Function> const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, MatrixExpression<Type, Row, Column, Operand, Function>, add<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Row, unsigned int Column, typename LeftOperand, typename LeftFunction, typename RightOperand, typename RightFunction>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, LeftOperand, LeftFunction>, MatrixExpression<Type, Row, Column, RightOperand, RightFunction>, add<Type> >, id<Type> >
operator+ (MatrixExpression<Type, Row, Column, LeftOperand, LeftFunction> const & x, MatrixExpression<Type, Row, Column, RightOperand, RightFunction> const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, LeftOperand, LeftFunction>, MatrixExpression<Type, Row, Column, RightOperand, RightFunction>, add<Type> >, id<Type> >(x, y);
}

/* sub */
template <typename Type, unsigned int Row, unsigned int Column>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Matrix<Type, Row, Column>, sub<Type> >, id<Type> >
operator- (Matrix<Type, Row, Column> const & x, Matrix<Type, Row, Column> const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Matrix<Type, Row, Column>, sub<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Row, unsigned int Column, typename Operand, typename Function>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, Operand, Function>, Matrix<Type, Row, Column>, sub<Type> >, id<Type> >
operator- (MatrixExpression<Type, Row, Column, Operand, Function> const & x, Matrix<Type, Row, Column> const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, Operand, Function>, Matrix<Type, Row, Column>, sub<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Row, unsigned int Column, typename Operand, typename Function>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, MatrixExpression<Type, Row, Column, Operand, Function>, sub<Type> >, id<Type> >
operator- (Matrix<Type, Row, Column> const & x, MatrixExpression<Type, Row, Column, Operand, Function> const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, MatrixExpression<Type, Row, Column, Operand, Function>, sub<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Row, unsigned int Column, typename LeftOperand, typename LeftFunction, typename RightOperand, typename RightFunction>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, LeftOperand, LeftFunction>, MatrixExpression<Type, Row, Column, RightOperand, RightFunction>, sub<Type> >, id<Type> >
operator- (MatrixExpression<Type, Row, Column, LeftOperand, LeftFunction> const & x, MatrixExpression<Type, Row, Column, RightOperand, RightFunction> const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, LeftOperand, LeftFunction>, MatrixExpression<Type, Row, Column, RightOperand, RightFunction>, sub<Type> >, id<Type> >(x, y);
}

/* mul */
template <typename Type, unsigned int Row, unsigned int Column>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Type, mul<Type> >, id<Type> >
operator* (Matrix<Type, Row, Column> const & x, Type const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Type, mul<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Row, unsigned int Column, typename Operand, typename Function>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, Operand, Function>, Type, mul<Type> >, id<Type> >
operator* (MatrixExpression<Type, Row, Column, Operand, Function> const & x, Type const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, Operand, Function>, Type, mul<Type> >, id<Type> >(x, y);
}

/* div */
template <typename Type, unsigned int Row, unsigned int Column>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Type, div<Type> >, id<Type> >
operator/ (Matrix<Type, Row, Column> const & x, Type const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Type, div<Type> >, id<Type> >(x, y);
}
template <typename Type, unsigned int Row, unsigned int Column, typename Operand, typename Function>
inline MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, Operand, Function>, Type, div<Type> >, id<Type> >
operator/ (MatrixExpression<Type, Row, Column, Operand, Function> const & x, Type const & y) {
	return MatrixExpression<Type, Row, Column, MatrixBinaryExpression<Type, Row, Column, MatrixExpression<Type, Row, Column, Operand, Function>, Type, div<Type> >, id<Type> >(x, y);
}

/* det */
template <typename Type, unsigned int Dimension>
inline Type
det (Matrix<Type, Dimension, Dimension> const & matrix) {
	Matrix<Type, Dimension, Dimension> F;
	F = matrix;
	/* factorization */
	for (unsigned int k = 0; k < Dimension; ++k) {
		for (unsigned int j = k+1; j < Dimension; ++j)
			F(j,k) /= F(k,k);
		for (unsigned int j = k+1; j < Dimension; ++j) {
			for (unsigned int i = k+1; i < Dimension; ++i)
				F(i,j) -= F(i,k) * F(k,j);
		}
	}
	/* det */
	Type r = F(0,0);
	for (unsigned int i = 1; i < Dimension; ++i)
		r *= F(i,i);
	return r; 
}
template <typename Type, unsigned int Dimension, typename Operand, typename Function>
inline Type
det (MatrixExpression<Type, Dimension, Dimension, Operand, Function> const & matrix) {
	Matrix<Type, Dimension, Dimension> F;
	F = matrix;
	/* factorization */
	for (unsigned int k = 0; k < Dimension; ++k) {
		for (unsigned int j = k+1; j < Dimension; ++j)
			F(j,k) /= F(k,k);
		for (unsigned int j = k+1; j < Dimension; ++j) {
			for (unsigned int i = k+1; i < Dimension; ++i)
				F(i,j) -= F(i,k) * F(k,j);
		}
	}
	/* det */
	Type r = F(0,0);
	for (unsigned int i = 1; i < Dimension; ++i)
		r *= F(i,i);
	return r; 
}

/* det 2x2 */
template <typename Type >
inline Type
det (Matrix<Type, 2, 2> const & matrix) {
	return matrix(0,0)*matrix(1,1)-matrix(1,0)*matrix(0,1);
}
template <typename Type, typename Operand, typename Function>
inline Type
det (MatrixExpression<Type, 2, 2, Operand, Function> const & matrix) {
	return matrix(0,0)*matrix(1,1)-matrix(1,0)*matrix(0,1);
}

}

#endif // _GAS_LINEARALGEBRA_MATRIX_
