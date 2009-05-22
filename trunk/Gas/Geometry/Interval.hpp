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

#ifndef _GAS_GEOMETRY_INTERVAL_
#define _GAS_GEOMETRY_INTERVAL_

namespace Geometry {

/*!
 * Questa classe contiene le informazioni geometriche di un intervallo e 
 * permette di trasformare le coordinate dei punti dall'intervallo di 
 * riferimento \f$[-1,1]\f$ a quello attuale \f$[a,b]\f$.
 * La trasformazione usata è una trasformazione affine della forma \f$x=
 * \alpha\hat{x}+\beta\f$.
 */ 

class Interval {
	
	public:
		
		/*!
		 * Costruttore di default, crea una geometria vuota.
		 */
		inline Interval ();
		/*!
		 * Partendo dagli estremi dell'intervallo viene costruita la geometria.
		 * @param a Estremo inferiore
		 * @param b Estremo superiore
		 * @pre a < b
		 */
		inline Interval (double const & a , double const & b);
		
		/*!
		 * Trasforma la coordinata dall'estremo di riferimento a quello attuale
		 * @param x L'ascissa di riferimento
		 * @return L'ascissa nel triangolo attuale.
		 */
		inline double xTransform (double const & x) const;
		
		/*!
		 * Operatore di copia.
		 * @param i Un altro intervallo.
		 * @return Ritorna il riferimento al lvalue.
		 */
		inline Interval & operator= (Interval const & i);
		
		/*!
		 * Il determinante della trasformazione locale, equivalente alla 
		 * metà della lunghezza dell'intervallo.
		 * @param x L'ascissa di riferimento dove valutare il determinante 
		 *          della trasformazione.
		 * @return Il determinante della trasformazione.
		 */
		inline double det (double const & x = 0.) const;
		
	private:
		double a_;
		double b_;
	
};

/* costruttore vuoto */
Interval::Interval () {
}

/* costruttore dagli estremi */
Interval::Interval (double const & a ,double const & b) {
	a_ = a;
	b_ = b;
}

/* trasforma le coordinate dei punti */
double Interval::xTransform (double const & x) const {
	return ( ( b_ - a_ ) * x + ( b_ + a_ ) ) / 2.;
}

/* operatore di copia */
Interval & Interval::operator= (Interval const & i) {
	a_ = i.a_;
	b_ = i.b_;
	return *this;
}

/* determinante della trasformazione */
double Interval::det (double const & x) const {
	return ( b_ - a_ ) / 2.;
}

}

#endif // _GAS_GEOMETRY_INTERVAL_
