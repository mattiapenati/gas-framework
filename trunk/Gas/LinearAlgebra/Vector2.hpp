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

/* main class */
template < typename Type >
class Vector<Type, 2> {

	public:
		/* constructor */
		Vector ( );
		Vector ( Type const & );
		Vector ( Vector<Type, 2> const & );
		Vector ( Type const & , Type const & );
		template < typename Operand , typename Function > Vector ( VectorExpression<Type, 2, Operand, Function> const & );
		
		/* destructor */
		~Vector ( );
		
		/* assign */
		inline Vector<Type, 2> & operator= ( Vector<Type, 2> const & );
		template < typename Operand , typename Function > inline Vector<Type, 2> & operator= ( VectorExpression<Type, 2, Operand, Function> const & );
		
		/* access */
		inline Type & operator() ( unsigned int const & );
		inline Type const & operator() ( unsigned int const & ) const;
	
	private:
		Type value_[2];

};

/* default constructor */
template < typename Type >
Vector < Type , 2 >::Vector ( ) {
}

/* constructor by value */
template < typename Type >
Vector < Type , 2 >::Vector ( Type const & scalar ) {
	value_[0] = scalar;
	value_[1] = scalar;
}

/* copy constructor */
template < typename Type >
Vector < Type , 2 >::Vector ( Vector<Type, 2> const & vector ) {
	value_[0] = vector(0);
	value_[1] = vector(1);
}
/* exrpession constructor */
template < typename Type >
template < typename Operand , typename Function >
Vector < Type , 2 >::Vector ( VectorExpression<Type, 2, Operand, Function> const & expression ) {
	value_[0] = expression(0);
	value_[1] = expression(1);
}

/* constructor by two values */
template < typename Type >
Vector < Type , 2 >::Vector ( Type const & v0 , Type const & v1 ) {
	value_[0] =  v0;
	value_[1] =  v1;
}

/* destructor */
template < typename Type >
Vector < Type , 2 >::~Vector ( ) {
}

/* assign from vector */
template < typename Type >
Vector<Type, 2> & Vector < Type , 2 >::operator= ( Vector<Type, 2> const & vector ) {
	value_[0] = vector(0);
	value_[1] = vector(1);
	return *this;
}

/* assign from exrpession */
template < typename Type >
template < typename Operand , typename Function >
Vector<Type, 2> & Vector < Type , 2 >::operator= ( VectorExpression<Type, 2, Operand, Function> const & expression ) {
	value_[0] = expression(0);
	value_[1] = expression(1);
	return *this;
}

/* access */
template < typename Type >
Type & Vector < Type , 2 >::operator() ( unsigned int const & i ) {
	return value_[i];
}
template < typename Type >
Type const & Vector < Type , 2 >::operator() ( unsigned int const & i ) const {
	return value_[i];
}