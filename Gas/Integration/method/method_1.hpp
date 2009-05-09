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

template < typename Method_ , typename Geometry_ >
class method_1 {
	
	public:
		typedef Geometry_ Geometry;
	
	private:
		/* punti trasformati dal riferimento */
		struct TransformPolicy {
			Geometry_ const *gp;
			inline TransformPolicy ( Geometry_ const & g ) : gp( &g ) { };
			inline double x ( double const & x_ ) { return gp->xTransform( x_ ); }
		};
		/* punti senza trasformazione */
		struct NoTransformPolicy {
			Geometry_ const *gp;
			inline NoTransformPolicy ( Geometry_ const & g ) : gp( &g ) { };
			inline double x ( double const & x_ ) { return x_; }
		};
		
	public:
		typedef TransformPolicy Transform;
		typedef NoTransformPolicy NoTransform;
	
	public:
		
		/* costruttore */
		method_1 ( );
		
		/* assegnamento geometria */
		void domain ( Geometry_ const & );
		
		/* integratori */
		template < typename TransformationPolicy , typename FunctionType >
		double integrate ( FunctionType const & ) const ;
		
		template < typename TransformationPolicy1 , typename TransformationPolicy2 , typename FunctionType1 , typename FunctionType2 >
		double integrateMul ( FunctionType1 const & , FunctionType2 const & ) const ;
	
	private:
		Geometry_ g_;
	
};

/* costruttore di default */
template < typename Method_ , typename Geometry_ >
method_1 < Method_ , Geometry_ >::method_1 ( ) {
}

/* assegnamento della geometria */
template < typename Method_ , typename Geometry_ >
void method_1 < Method_ , Geometry_ >::domain ( Geometry_ const & g ) {
	g_ = g;
}

/* integratore semplice */
template < typename Method_ , typename Geometry_ >
template < typename TransformationPolicy , typename FunctionType >
double method_1 < Method_ , Geometry_ >::integrate ( FunctionType const & f ) const {
	
	TransformationPolicy t(g_);
	double r = 0.;
	unsigned int n  = Method_::nPoints;
	
	for ( unsigned int i = 0 ; i < n ; ++i )
		r += std::abs(g_.det(Method_::x[i])) * Method_::w[i] * f(t.x(Method_::x[i]));
	
	return r;
}

/* integratore di moltiplicazione */
template < typename Method_ , typename Geometry_ >
template < typename TransformationPolicy1 , typename TransformationPolicy2 , typename FunctionType1 , typename FunctionType2 >
double method_1 < Method_ , Geometry_ >::integrateMul ( FunctionType1 const & f , FunctionType2 const & g ) const {
	
	TransformationPolicy1 t1(g_);
	TransformationPolicy2 t2(g_);
	double r = 0.;
	unsigned int n  = Method_::nPoints;
	
	for ( unsigned int i = 0 ; i < n ; ++i )
		r += std::abs(g_.det(Method_::x[i])) * Method_::w[i] * f(t1.x(Method_::x[i])) * g(t2.x(Method_::x[i]));
	
	return r;
}