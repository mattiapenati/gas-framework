struct Upper : public solver < Upper > {

	template < typename Type , unsigned int Dimension >
	static void apply ( Matrix < Type , Dimension , Dimension > const & U , Vector < Type , Dimension > & x , Vector < Type , Dimension > const & b ) {
	
		unsigned int N = Dimension - 1u;
		
		x(N) = b(N) /  U(N,N);
		
		for ( int i = N-1 ; i >= 0 ; --i ) {
			x(i) = b(i);
			for ( unsigned int j = i+1 ; j < Dimension ; ++j ) {
				x(i) -= U(i,j) * x(j);
			}
			x(i) /= U(i,i);
		}
	
	}
	
};