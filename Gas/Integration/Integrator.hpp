#ifndef GAS_INTEGRATOR_HPP
#define GAS_INTEGRATOR_HPP

/* un esempio di utilizzo
 * - costruisco l'integratore con il metodo richiesto, Newton-Cotes
 *   su un triangolo (basato sulla faccia CGAL) e di ordine 2
 *
 *    typedef Integrator< NewtonCotes_2<Geometry::Triangle<Fb>, 1> > Int
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
		typedef typename IntegrationMethod::Function Function;
		
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
		self & domain (Geometry const & g) {
			m_.setGeometry(g);
			return *this;
		}
		
		/* applicazione del metodo di integrazione su 
		 * una singola funzione */
		template<typename TransformationPolicy>
		double integrate (Function const & f) {
			return m_.template apply<TransformationPolicy>(f);
		}
		
		/* questo metodo viene usato per la moltiplicazione
		 * tra due funzioni: utile per il termine noto e il
		 * termine di reazione */
		template<typename TransformationPolicy1, typename TransformationPolicy2>
		double integrateMul (Function const & f, Function const & g) {
			return m_.template applyMul<TransformationPolicy1, TransformationPolicy2>(f, g);
		}
};

#endif // GAS_INTEGRATOR_HPP
