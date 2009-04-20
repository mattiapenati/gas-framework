#ifndef GAS_NEWTONC2_HPP
#define GAS_NEWTONC2_HPP

template<typename Geometry, unsigned int Order> class NewtonCotes_2;

// TODO Generalizzare costruendo una classe method che dichiara il
//      funzionamento base, poi i derivati devono solo definire il
//      numero dei nodi, le loro coordinate e i pesi.

/* Questa classe definisce il metodo di Newton-Cotes su un triangolo
 * di ordine 1. In modo del tutto identico si definisce il metodo di
 * ordine successivo */
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
	
	// Transformation Policies
	private:
		struct TransformPolicy {
			Geometry const *local_g;
			inline TransformPolicy (Geometry const & g) : local_g(&g) {};
			inline double x (double const & x_, double const & y_) { return local_g->xTransform(x_, y_); }
			inline double y (double const & x_, double const & y_) { return local_g->yTransform(x_, y_); }
		};
		struct NoTransformPolicy {
			Geometry const *local_g;
			inline NoTransformPolicy (Geometry const & g) : local_g(&g) {};
			inline double x (double const & x_, double const & y_) { return x_; }
			inline double y (double const & x_, double const & y_) { return y_; }
		};
		
	public:
		typedef TransformPolicy Transform;
		typedef NoTransformPolicy NoTransform;

	public:
		NewtonCotes_2 () : g_() {
			/* init nodes and weight */
			x[0] = 0.1666666666667; x[1] = 0.6666666666667; x[2] = 0.1666666666667;
			y[0] = 0.6666666666667; y[1] = 0.1666666666667; y[2] = 0.1666666666667;
			w[0] = w[1] = w[2] = 0.6666666666667;
		}
		
		void setGeometry(Geometry const & g) {
			g_ = g;
		}
		
		/* la definizione di questo template permette di dichiararare se bisogna, oppure no, 
		 * trasformare le coordinate dei nodi nel triangolo definito come dominio */ 
		template<typename TransformationPolicy>
		double apply (Function const & f) {
			/* chiamando t.x() si ottiene la coordinata nel triangolo definito, 
			 * invece che nel riferimento */
			TransformationPolicy t(g_);
			double r = 0.;
			
			r += w[0] * f( t.x(x[0], y[0]) , t.y(x[0], y[0]) );
			r += w[1] * f( t.x(x[1], y[1]) , t.y(x[1], y[1]) );
			r += w[2] * f( t.x(x[2], y[2]) , t.y(x[2], y[2]) );
			
			r *= g_.det();
			
			return r;
		}
		template<typename TransformationPolicy1, typename TransformationPolicy2>
		double applyMul (Function const & f, Function const & g) {
			TransformationPolicy1 t1(g_);
			TransformationPolicy2 t2(g_);
			double r = 0.;
			
			r += w[0] * f( t1.x(x[0], y[0]) , t1.y(x[0], y[0]) ) * g( t2.x(x[0], y[0]) , t2.y(x[0], y[0]) );
			r += w[1] * f( t1.x(x[1], y[1]) , t1.y(x[1], y[1]) ) * g( t2.x(x[1], y[1]) , t2.y(x[1], y[1]) );
			r += w[2] * f( t1.x(x[2], y[2]) , t1.y(x[2], y[2]) ) * g( t2.x(x[2], y[2]) , t2.y(x[2], y[2]) );
			
			r *= g_.det();
			
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
	
	// Transformation Policies
	private:
		struct TransformPolicy {
			Geometry const *local_g;
			inline TransformPolicy (Geometry const & g) : local_g(&g) {};
			inline double x (double const & x_, double const & y_) { return local_g->xTransform(x_, y_); }
			inline double y (double const & x_, double const & y_) { return local_g->yTransform(x_, y_); }
		};
		struct NoTransformPolicy {
			Geometry const *local_g;
			inline NoTransformPolicy (Geometry const & g) : local_g(&g) {};
			inline double x (double const & x_, double const & y_) { return x_; }
			inline double y (double const & x_, double const & y_) { return y_; }
		};
		
	public:
		typedef TransformPolicy Transform;
		typedef NoTransformPolicy NoTransform;

	public:
		NewtonCotes_2 () : g_() {
			/* init nodes and weight */
			x[0] = 0.0915762135098; x[1] = 0.8168475729805; x[2] = 0.0915762135098; x[3] = 0.1081030181681; x[4] = 0.4459484909160; x[5] = 0.4459484909160;
			y[0] = 0.0915762135098; y[1] = 0.0915762135098; y[2] = 0.8168475729805; y[3] = 0.4459484909160; y[4] = 0.1081030181681; y[5] = 0.4459484909160;
			w[0] = w[1] = w[2] = 0.2199034873106;
			w[3] = w[4] = w[5] = 0.446763179356;
		}
		
		void setGeometry(Geometry const & g) {
			g_ = g;
		}
		
		template<typename TransformationPolicy>
		double apply (Function const & f) {
			TransformationPolicy t(g_);
			double r = 0.;
			
			r += w[0] * f( t.x(x[0], y[0]) , t.y(x[0], y[0]) );
			r += w[1] * f( t.x(x[1], y[1]) , t.y(x[1], y[1]) );
			r += w[2] * f( t.x(x[2], y[2]) , t.y(x[2], y[2]) );
			r += w[3] * f( t.x(x[3], y[3]) , t.y(x[3], y[3]) );
			r += w[4] * f( t.x(x[4], y[4]) , t.y(x[4], y[4]) );
			r += w[5] * f( t.x(x[5], y[5]) , t.y(x[5], y[5]) );
			
			r *= g_.det();
			
			return r;
		}
		template<typename TransformationPolicy1, typename TransformationPolicy2>
		double applyMul (Function const & f, Function const & g) {
			TransformationPolicy1 t1(g_);
			TransformationPolicy2 t2(g_);
			double r = 0.;
			
			r += w[0] * f( t1.x(x[0], y[0]) , t1.y(x[0], y[0]) ) * g( t2.x(x[0], y[0]) , t2.y(x[0], y[0]) );
			r += w[1] * f( t1.x(x[1], y[1]) , t1.y(x[1], y[1]) ) * g( t2.x(x[1], y[1]) , t2.y(x[1], y[1]) );
			r += w[2] * f( t1.x(x[2], y[2]) , t1.y(x[2], y[2]) ) * g( t2.x(x[2], y[2]) , t2.y(x[2], y[2]) );
			r += w[3] * f( t1.x(x[3], y[3]) , t1.y(x[3], y[3]) ) * g( t2.x(x[3], y[3]) , t2.y(x[3], y[3]) );
			r += w[4] * f( t1.x(x[4], y[4]) , t1.y(x[4], y[4]) ) * g( t2.x(x[4], y[4]) , t2.y(x[4], y[4]) );
			r += w[5] * f( t1.x(x[5], y[5]) , t1.y(x[5], y[5]) ) * g( t2.x(x[5], y[5]) , t2.y(x[5], y[5]) );
			
			r *= g_.det();
			
			return r;
		}
		
};
#endif // GAS_NEWTONC2_HPP
