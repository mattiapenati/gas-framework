struct Diagonal : public solver < Diagonal > {

	template < typename Type , unsigned int Dimension >
	static void apply ( Matrix < Type , Dimension , Dimension > const & D , Vector < Type , Dimension > & x , Vector < Type , Dimension > const & b ) {
	
		for ( unsigned int i = 0u ; i < Dimension ; ++i )
			x(i) = b(i) / D(i,i);
	
	}
	
};