#ifndef GAS_NEWTONC2_HPP
#define GAS_NEWTONC2_HPP

template<typename Geometry, unsigned int Order> class NewtonCotes_2;

/* Questa classe definisce il metodo di Newton-Cotes su un triangolo
 * di ordine 1. In modo del tutto identico si definisce il metodo di
 * ordine successivo */
template<>
template<typename Fb>
class NewtonCotes_2<Geometry::Triangle<Fb>, 1> {

	public:
		typedef Geometry::Triangle<Fb> Geometry;
		typedef double (*Function)(double const &, double const &);

	private:
		/* nodi */
		double x[3];
		double y[3];
		double w[3];
		/* geometria */
		Geometry g_;
	
	private:
		struct TransformPolicy {
			inline TransformPolicy () {}
			inline double x (double const & x) { return g_.xTransform(x); }
			inline double y (double const & y) { return g_.yTransform(y); }
		};
		struct NoTransformPolicy {
			inline NoTransformPolicy () {}
			inline double x (double const & x) { return x; }
			inline double y (double const & y) { return y; }
		};
		
	public:
		typedef TransformPolicy Transform;
		typedef NoTransformPolicy NoTransform;

	public:
		NewtonCotes_2 () : g_() {
			/* init nodes and weight */
			x[0] = 0.; x[1] = 1.; x[2] = 0.;
			y[0] = 0.; y[1] = 0.; y[2] = 1.;
			w[0] = w[1] = w[2] = 1./6.;
		}
		
		void setGeometry(Geometry const & g) {
			g_ = g;
		}
		
		/* la definizione di questo template permette
		 * di dichiararare se bisogna, oppure no, 
		 * trasformare le coordinate dei nodi nel
		 * triangolo definito come dominio */ 
		template<typename TransformationPolicy>
		double apply (Function const & f) {
			/* chiamando t.x() si ottiene la coordinata
			 * nel triangolo definito, invece che nel
			 * riferimento */
			TransformationPolicy t;
			double r = 0.;
			r += w[0] * f( t.x(x[0]) , t.y(y[0]) );
			r += w[1] * f( t.x(x[1]) , t.y(y[1]) );
			r += w[2] * f( t.x(x[2]) , t.y(y[2]) );
			
			return r;
		}
		template<typename TransformationPolicy1, typename TransformationPolicy2>
		double applyMul (Function const & f, Function const & g) {
			TransformationPolicy1 t1;
			TransformationPolicy2 t2;
			double r = 0.;
			
			r += w[0] * f( t1.x(x[0]) , t1.y(y[0]) ) * g( t2.x(x[0]) , t2.y(y[0]) );
			r += w[1] * f( t1.x(x[1]) , t1.y(y[1]) ) * g( t2.x(x[1]) , t2.y(y[1]) );
			r += w[2] * f( t1.x(x[2]) , t1.y(y[2]) ) * g( t2.x(x[2]) , t2.y(y[2]) );
			
			r *= g_.area();
			
			return r;
		}
		
};

#endif // GAS_NEWTONC2_HPP