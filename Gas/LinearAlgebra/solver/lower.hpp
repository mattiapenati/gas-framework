struct Lower : public solver < Lower > {

	template < typename Type , unsigned int Dimension >
	static void apply ( Matrix < Type , Dimension , Dimension > const & L , Vector < Type , Dimension > & x , Vector < Type , Dimension > const & b ) {
	
		x(0) = b(0) / L(0,0) ;
		
		for ( unsigned int i = 1u ; i < Dimension ; ++i ) {
			x(i) = b(i);
			for ( unsigned int j = 0u ; j < i ; ++j ) {
				x(i) -= L(i,j) * x(j);
			}
			x(i) /= L(i,i);
		}
	
	}
	
};