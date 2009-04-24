/* Classe per l'utilizzo della geometria triangolare, con trasformazione
 * affine x = Jy + b */
class Triangle {

	public:
		
		/* costruttori */
		Triangle ( );
		template < typename Face > Triangle ( Face const & );
		
		/* trasformazione delle coordinate */
		double xTransform ( double const & , double const & ) const;
		double yTransform ( double const & , double const & ) const;
		
		/* copia della geometria */
		Triangle & operator= ( Triangle const & );
		
		/* determinante della trasformazione */
		
		double det ( ) const;
	
	private:
		/* parametri della trasformazione */
		double J[2][2];
		double b[2];
		/* determinante della trasformazione */
		double detJ;

};

/* costruttore di default */
Triangle::Triangle ( ) {
}

/* costruttore di copia dalla faccia */
template < typename Face >
Triangle::Triangle ( Face const & f ) {
	double x[3];
	double y[3];
	
	x[0] = f.vertex(0)->point().x(); y[0] = f.vertex(0)->point().y();
	x[1] = f.vertex(1)->point().x(); y[1] = f.vertex(1)->point().y();
	x[2] = f.vertex(2)->point().x(); y[2] = f.vertex(2)->point().y();
	
	/* calcolo di J e b */
	J[0][0] = x[1] - x[0]; 
	J[0][1] = x[2] - x[0];
	J[1][0] = y[1] - y[0];
	J[1][1] = y[2] - y[0];
	b[0] = x[0];
	b[1] = y[0];
	
	/* det J */
	detJ = J[0][0]*J[1][1] - J[1][0]*J[0][1];
}
		
/* trasforma le coordinate dal riferimento a questo triangolo */
double Triangle::xTransform ( double const & x , double const & y ) const {
	return J[0][0] * x + J[0][1] * y + b[0];
}
double Triangle::yTransform ( double const & x , double const & y ) const {
	return J[1][0] * x + J[1][1] * y + b[1];
}

/* copia della geometria */
Triangle & Triangle::operator= ( Triangle const & t ) {
	/* trasformazione */
	J[0][0] = t.J[0][0]; J[0][1] = t.J[0][1];
	J[1][0] = t.J[1][0]; J[1][1] = t.J[1][1];
	/* traslazione */
	b[0] = t.b[0]; b[1] = t.b[1];
	/* determinante */
	detJ = t.detJ;
	
	return *this;
};

/* determinante di J */
double Triangle::det ( ) const {
	return detJ; 
}
