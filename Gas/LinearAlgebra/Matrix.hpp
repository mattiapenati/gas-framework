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

template < typename Type , unsigned int Row , unsigned int Column >
class Matrix {

	public:
		/* constructor */
		Matrix ( );
		Matrix ( Type const & );
		Matrix ( Matrix < Type , Row , Column > const & );
		
		/* destructor */
		~Matrix ( );
		
		/* access */
		inline Type & operator() ( unsigned int const & , unsigned int const & );
		inline Type const & operator() ( unsigned int const & , unsigned int const & ) const ;
	
	private:
		Type value_[Row][Column];

};

/* default constructor */
template < typename Type , unsigned int Row , unsigned int Column >
Matrix < Type , Row , Column >::Matrix ( ) {
}

/* constructor by value */
template < typename Type , unsigned int Row , unsigned int Column >
Matrix < Type , Row , Column >::Matrix ( Type const & scalar ) {
	for ( unsigned int i = 0u ; i < Row ; ++i ) {
		for ( unsigned int j = 0u ; j < Column ; ++j ) {
			value_[i][j] = scalar;
		}
	}
}
/* copy constructor */
template < typename Type , unsigned int Row , unsigned int Column >
Matrix < Type , Row , Column >::Matrix ( Matrix < Type , Row , Column > const & matrix ) {
	for ( unsigned int i = 0u ; i < Row ; ++i ) {
		for ( unsigned int j = 0u ; j < Column ; ++j ) {
			value_[i][j] = matrix(i,j);
		}
	}
}

/* destructor */
template < typename Type , unsigned int Row , unsigned int Column >
Matrix < Type , Row , Column >::~Matrix ( ) {
}

/* access */
template < typename Type , unsigned int Row , unsigned int Column >
Type & Matrix < Type , Row , Column >::operator() ( unsigned int const & i , unsigned int const & j ) {
	return value_[i][j];
}
template < typename Type , unsigned int Row , unsigned int Column >
Type const & Matrix < Type , Row , Column >::operator() ( unsigned int const & i , unsigned int const & j ) const {
	return value_[i][j];
}