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

/* un esempio di utilizzo
 * - costruisco l'integratore con il metodo richiesto, Newton-Cotes
 *   su un triangolo (basato sulla faccia CGAL) e di ordine 2
 *
 *    typedef Integrator< Method::NewtonCotes_2<Geometry::Triangle<Fb>, 1> > Int
 *    Int it;
 *
 * - definisco il dominio su cui viene applicato partendo dalla faccia
 *
 *    it.domain(Int::Geometry(face));
 *
 * - applico il metodo ad una funzione definita sul triangolo di 
 *   partenza
 *
 *    it.integrate<Int::Transform>(f);
 *
 * - oppure se lo vogliamo applicare per calcolare il termine noto
 *   di un problema ellittico, sapendo la base sul riferimento
 *
 *    it.integrate<Int::Transform, Int::NoTransform>(f, phi0);
 *
 * Si puo sintetizzare la chiamata in
 *
 *   it.domain(Int::Geometry(face)).integrate<Int::Transform>(f);
 */

template<typename IntegrationMethod>
class Integrator {
	public:
		typedef Integrator<IntegrationMethod> self;
	
		typedef typename IntegrationMethod::Geometry Geometry;
		
		/* queste sono le policies che definiscono se la funzione
		 * da integrare è definita sulla geometria oppure sulla
		 * geometria di riferimento */
		typedef typename IntegrationMethod::Transform Transform;
		typedef typename IntegrationMethod::NoTransform NoTransform;
		
	private:
		IntegrationMethod m_;
		
	public:
		/* constructor */
		Integrator () : m_ () {
		}
		
		/* impostazione del dominio di integrazione, 
		 * la geometria del problema viene definita 
		 * dal metodo scelto */
		inline self & domain (Geometry const & g) {
			m_.domain(g);
			return *this;
		}
		
		/* applicazione del metodo di integrazione su 
		 * una singola funzione */
		template<typename TransformationPolicy, typename FunctionType>
		inline double integrate (FunctionType const & f) {
			return m_.template integrate<TransformationPolicy>(f);
		}
		
		/* questo metodo viene usato per la moltiplicazione
		 * tra due funzioni: utile per il termine noto e il
		 * termine di reazione */
		template<typename TransformationPolicy1, typename TransformationPolicy2, typename FunctionType1, typename FunctionType2>
		inline double integrateMul (FunctionType1 const & f, FunctionType2 const & g) {
			return m_.template integrateMul<TransformationPolicy1, TransformationPolicy2>(f, g);
		}
};
