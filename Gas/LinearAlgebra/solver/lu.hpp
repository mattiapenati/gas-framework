struct LU : public solver < LU > {

	template < typename Type , unsigned int Dimension >
	static void apply ( Matrix < Type , Dimension , Dimension > const & A , Vector < Type , Dimension > & x , Vector < Type , Dimension > const & b ) {
		Matrix < Type , Dimension , Dimension > F(A);
		Vector < Type , Dimension > y;
		
		/* factorization */
		for ( unsigned int k = 0 ; k < Dimension ; ++k ) {
			for ( unsigned int j = k+1 ; j < Dimension ; ++j )
				F(j,k) /= F(k,k);
			for ( unsigned int j = k+1; j < Dimension; ++j ) {
				for ( unsigned int i = k+1 ; i < Dimension ; ++i )
					F(i,j) -= F(i,k) * F(k,j);
			}
		}
		
		Lower1::solve( F , y , b );
		Upper::solve( F , x , y );
	}

};