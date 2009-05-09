/*                                                          
 * Copyright (c) 2008, Alfonso Fascì, Davide Ferrarese, Mattia Penati     
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

/*!
 * Questa classe contiene le informazioni geometriche di una geometria 
 * triangolare e permette di trasformare le coordinate dei punti dal triangolo
 * di riferimento di vertici \f$(0,0)\f$, \f$(1,0)\f$, \f$(0,1)\f$ a quello
 * attuale.
 * La trasformazione usata è una trasformazione affine della forma 
 * \f$\mathbf{x} = J \hat{\mathbf{x}} + \mathbf{b}\f$.
 */

class Triangle {

	public:
		
		/*!
		 * Costruttore di default, crea una geometria vuota.
		 */
		inline Triangle ( );
		/*!
		 * Unico costruttore, partendo da una faccia della triangolazione
		 * calcola i paramentri delle trasformazioni. L'utilizzo del template
		 * permette di specializzare la classe per un tipo di faccia
		 * personalizzata.
		 * @param f La faccia da cui costruire la geometria.
		 */
		template < typename Face > inline Triangle ( Face const & f );
		
		/*!
		 * Trasforma le coordinate dal triangolo di riferimento, partendo da
		 * quelle sul triangolo attuale.
		 * @param x L'ascissa nel riferimento.
		 * @param y L'ordinata nel riferimento.
		 * @return L'ascissa nel triangolo attuale.
		 */
		inline double xTransform ( double const & x , double const & y ) const;
		
		/*!
		 * Trasforma le coordinate dal triangolo di riferimento, partendo da
		 * quelle sul triangolo attuale.
		 * @param x L'ascissa nel riferimento.
		 * @param y L'ordinata nel riferimento.
		 * @return L'ordinata nel triangolo attuale.
		 */
		inline double yTransform ( double const & x , double const & y ) const;
		
		/*!
		 * Operatore di copia.
		 * @param t Un altro triangolo.
		 * @return Ritorna il riferimento al lvalue.
		 */
		inline Triangle & operator= ( Triangle const & t );
		
		/*!
		 * Il determinante della trasformazione, il suo valore diviso per 2
		 * permette di ottenere l'area del triangolo.
		 * @return Il determinante della trasformazione.
		 */
		inline double det ( ) const;
	
	private:
		/* parametri della trasformazione */
		LinearAlgebra::Matrix<double, 2, 2> J;
		LinearAlgebra::Vector<double, 2> b;
		/* determinante della trasformazione */
		double detJ;

};

/* costruttore vuoto */
Triangle::Triangle ( ) {
}

/* costruttore di copia dalla faccia */
template < typename Face >
Triangle::Triangle ( Face const & f ) {
	
	/* coordinate dei vertici del triangolo */
	LinearAlgebra::Vector<double, 2> p0( f.vertex(0)->point().x() , f.vertex(0)->point().y() );
	LinearAlgebra::Vector<double, 2> p1( f.vertex(1)->point().x() , f.vertex(1)->point().y() );
	LinearAlgebra::Vector<double, 2> p2( f.vertex(2)->point().x() , f.vertex(2)->point().y() );
	
	/* temporanei */
	LinearAlgebra::Vector<double, 2> v1(p1-p0);
	LinearAlgebra::Vector<double, 2> v2(p2-p0);
	
	/* calcolo di J e b */
	J(0,0) = v1(0);
	J(1,0) = v1(1);
	J(0,1) = v2(0);
	J(1,1) = v2(1);
	b = p0;
	
	/* det J */
	detJ = LinearAlgebra::det(J);
}
		
/* trasforma le coordinate dal riferimento al triangolo attuale */
double Triangle::xTransform ( double const & x , double const & y ) const {
	return J(0,0) * x + J(0,1) * y + b(0);
}
double Triangle::yTransform ( double const & x , double const & y ) const {
	return J(1,0) * x + J(1,1) * y + b(1);
}

/* copia della geometria */
Triangle & Triangle::operator= ( Triangle const & t ) {
	J = t.J;
	b = t.b;
	detJ = t.detJ;
	
	return *this;
};

/* determinante di J */
double Triangle::det ( ) const {
	return detJ; 
}
