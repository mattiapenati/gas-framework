#ifndef GAS_TRIANGLE_HPP
#define GAS_TRIANGLE_HPP

/* questo template definisce il tipo che viene usato per memorizzare una faccia */
template<typename Fb>
class Triangle {
	public:
	
	private:
		/* parametri della trasformazione */
		double J[2][2];
		double b[2];
		double detJ;
	public:
		Triangle () {
		}
		
		Triangle (Fb const & face) {
			/* coordinate dei punti: non sono necessari per i metodi numerici */
			double x[3];
			double y[3];
				
			x[0] = face.vertex(0)->point().x(); y[0] = face.vertex(0)->point().y();
			x[1] = face.vertex(1)->point().x(); y[1] = face.vertex(1)->point().y();
			x[2] = face.vertex(2)->point().x(); y[2] = face.vertex(2)->point().y();
			
			/* calcolare J e b: questi sono necessari */
			J[0][0] = x[1] - x[0]; 
			J[0][1] = x[2] - x[0];
			J[1][0] = y[1] - y[0];
			J[1][1] = y[2] - y[0];
			b[0] = x[0];
			b[1] = y[0];
		}
		
		/* copia */
		Triangle<Fb> & operator= (Triangle<Fb> const & t) {
			/* trasformazione */
			J[0][0] = t.J[0][0]; J[0][1] = t.J[0][1];
			J[1][0] = t.J[1][0]; J[1][1] = t.J[1][1];
			/* traslazione */
			b[0] = t.b[0]; b[1] = t.b[1];
			
			return *this;
		};
		
		/* trasforma le coordinate dal riferimento 
		 * a questo triangolo */
		double xTransform (double const & x, double const & y) const {
			return J[0][0] * x + J[0][1] * y + b[0];
		}
		double yTransform (double const & x, double const & y) const {
			return J[1][0] * x + J[1][1] * y + b[1];
		}
		/* il determinante di J */
		double area () const {
			return J[0][0]*J[1][1] - J[1][0]*J[0][1];
		}
};

#endif // GAS_TRIANGLE_HPP
