template < typename Method_ , typename Geometry_ >
class Method_2 {
	
	public:
		typedef Geometry_ Geometry;
		typedef double (*Function) ( double const & , double const & );
	
	private:
		/* punti trasformati dal riferimento */
		struct TransformPolicy {
			Geometry_ const *gp;
			inline TransformPolicy ( Geometry_ const & g ) : gp( &g ) { };
			inline double x ( double const & x_ , double const & y_ ) { return gp->xTransform( x_ , y_ ); }
			inline double y ( double const & x_ , double const & y_ ) { return gp->yTransform( x_ , y_ ); }
		};
		/* punti senza trasformazione */
		struct NoTransformPolicy {
			Geometry_ const *gp;
			inline NoTransformPolicy ( Geometry_ const & g ) : gp( &g ) { };
			inline double x ( double const & x_ , double const & y_ ) { return x_; }
			inline double y ( double const & x_ , double const & y_ ) { return y_; }
		};
		
	public:
		typedef TransformPolicy Transform;
		typedef NoTransformPolicy NoTransform;
	
	public:
		
		/* costruttore */
		Method_2 ( );
		
		/* assegnamento geometria */
		void domain ( Geometry_ const & );
		
		/* integratori */
		template < typename TransformationPolicy >
		double integrate ( Function const & ) const ;
		
		template < typename TransformationPolicy1 , typename TransformationPolicy2 >
		double integrateMul ( Function const & , Function const & ) const ;
	
	private:
		Geometry_ g_;
	
};

/* costruttore di default */
template < typename Method_ , typename Geometry_ >
Method_2 < Method_ , Geometry_ >::Method_2 ( ) {
}

/* assegnamento della geometria */
template < typename Method_ , typename Geometry_ >
void Method_2 < Method_ , Geometry_ >::domain ( Geometry_ const & g ) {
	g_ = g;
}

/* integratore semplice */
template < typename Method_ , typename Geometry_ >
template < typename TransformationPolicy >
double Method_2 < Method_ , Geometry_ >::integrate ( Function const & f ) const {
	
	TransformationPolicy t(g_);
	double r = 0.;
	unsigned int n  = Method_::nPoints;
	
	for ( unsigned int i = 0 ; i < n ; ++i ) {
		r += Method_::w[i] * f( t.x(Method_::x[i], Method_::y[i]) , t.y(Method_::x[i], Method_::y[i]) );
	}
	
	r *= std::abs(g_.det());
	
	return r;
}

/* integratore di moltiplicazione */
template < typename Method_ , typename Geometry_ >
template < typename TransformationPolicy1 , typename TransformationPolicy2 >
double Method_2 < Method_ , Geometry_ >::integrateMul ( Function const & f , Function const & g ) const {
	
	TransformationPolicy1 t1(g_);
	TransformationPolicy2 t2(g_);
	double r = 0.;
	unsigned int n  = Method_::nPoints;
	
	for ( unsigned int i = 0 ; i < n ; ++i ) {
		r += Method_::w[i] * 
			f( t1.x(Method_::x[i], Method_::y[i]) , t1.y(Method_::x[i], Method_::y[i]) ) * 
			g( t2.x(Method_::x[i], Method_::y[i]) , t2.y(Method_::x[i], Method_::y[i]) );
	}
	
	r *= std::abs(g_.det());
	
	return r;
}