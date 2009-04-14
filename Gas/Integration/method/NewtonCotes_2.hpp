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
		/* pesi */
		double w[3];
		/* geometria */
		Geometry g_;
	
	private:
		struct TransformPolicy {
			inline TransformPolicy () {}
			inline double x (double const & x) { return g_.xTransform(x, y); }
			inline double y (double const & y) { return g_.yTransform(x, y); }
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

template<typename Fb>
class NewtonCotes_2<Geometry::Triangle<Fb>, 2> {

	public:
		typedef Geometry::Triangle<Fb> Geometry;
		typedef double (*Function)(double const &, double const &);

	private:
		/* nodi */
		double x[6];
		double y[6];
		/* pesi */
		double w[6];
		/* geometria */
		Geometry g_;
	
	private:
		struct TransformPolicy {
			inline TransformPolicy () {}
			inline double x (double const & x) { return g_.xTransform(x, y); }
			inline double y (double const & y) { return g_.yTransform(x, y); }
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
			x[0] = 0.; x[1] = 1.; x[2] = 0.; x[3] = 0.; x[4] = 0.5; x[5] = 0.5;
			y[0] = 0.; y[1] = 0.; y[2] = 1.; y[3] = 0.5; y[4] = 0.; y[5] = 0.5;
			w[0] = w[1] = w[2] = 0.2199034873106;
			w[3] = w[4] = w[5] = 0.446763179356;
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
			r += w[3] * f( t.x(x[3]) , t.y(y[3]) );
			r += w[4] * f( t.x(x[4]) , t.y(y[4]) );
			r += w[5] * f( t.x(x[5]) , t.y(y[5]) );
			
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
			r += w[3] * f( t1.x(x[3]) , t1.y(y[3]) ) * g( t2.x(x[3]) , t2.y(y[3]) );
			r += w[4] * f( t1.x(x[4]) , t1.y(y[4]) ) * g( t2.x(x[4]) , t2.y(y[4]) );
			r += w[5] * f( t1.x(x[5]) , t1.y(y[5]) ) * g( t2.x(x[5]) , t2.y(y[5]) );
			
			r *= g_.area();
			
			return r;
		}
		
};
#endif // GAS_NEWTONC2_HPP
