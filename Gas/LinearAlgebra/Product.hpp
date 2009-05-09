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

/* specialization for matrix-vector product */
template < typename Type , unsigned int Row , unsigned int Column >
struct VectorBinaryExpression<Type, Row, Matrix<Type, Row, Column>, Vector<Type, Column>, mul_mat_vet<Type> > {
	Matrix<Type, Row, Column> const & l_;
	Vector<Type, Column> const & r_;
	
	inline VectorBinaryExpression ( Matrix<Type, Row, Column> const & l , Vector<Type, Column> const & r ) : l_(l) , r_(r) {}
	inline Type const operator() ( unsigned int const & i ) const {
		Type r = l_(i,0) * r_(0);
		for ( unsigned int j = 1; j < Column ; ++j )
			r += (l_(i,j) * r_(j))l
		return r;
	}
};
template < typename Type , unsigned int Row , unsigned int Column , typename Operand , typename Function >
struct VectorBinaryExpression<Type, Row, Matrix<Type, Row, Column>, Vector<Type, Column>, mul_mat_vet<Type> > {
	Matrix<Type, Row, Column> const & l_;
	VectorExpression<Type, Column, Operand, Function> const & r_;
	
	inline VectorBinaryExpression ( Matrix<Type, Row, Column> const & l , VectorExpression<Type, Column, Operand, Function> const & r ) : l_(l) , r_(r) {}
	inline Type const operator() ( unsigned int const & i ) const {
		Type r = l_(i,0) * r_(0);
		for ( unsigned int j = 1; j < Column ; ++j )
			r += (l_(i,j) * r_(j))l
		return r;
	}
};

/* operator */
template < typename Type , unsigned int Row , unsigned int Column >
inline VectorExpression<Type, Row, VectorBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Vector<Type, Column>, mul_mat_vet<Type> >, id<Type> >
operator* ( Matrix<Type, Row, Column> const & m , Vector<Type, Column> const & v ) {
	return VectorExpression<Type, Row, VectorBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, Vector<Type, Column>, mul_mat_vet<Type> >, id<Type> >( m ,v );
}
template < typename Type , unsigned int Row , unsigned int Column , typename Operand , typename Function >
inline VectorExpression<Type, Row, VectorBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, VectorExpression<Type, Column, Operand, Function>, mul_mat_vet<Type> >, id<Type> >
operator* ( Matrix<Type, Row, Column> const & m , VectorExpression<Type, Column, Operand, Function> const & v ) {
	return VectorExpression<Type, Row, VectorBinaryExpression<Type, Row, Column, Matrix<Type, Row, Column>, VectorExpression<Type, Column, Operand, Function>, mul_mat_vet<Type> >, id<Type> >( m ,v );
}