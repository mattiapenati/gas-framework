struct Cholesky : public solver < Cholesky > {

	template < typename Type , unsigned int Dimension >
	static void apply ( Matrix < Type , Dimension , Dimension > const & A , Vector < Type , Dimension > & x , Vector < Type , Dimension > const & b ) {
	
		Matrix< Type , Dimension , Dimension > H(A);
		
		H(0,0) = std::sqrt(A(0,0));
		for ( unsigned int i = 1u ; i < Dimension ; ++i ) {
			
			for ( unsigned int j = i-1u ; j < Dimension ; ++j ) {
				H(i,j) = A(i,j);
				for ( unsigned int k = 0u ; k < j-1u ; ++k )
					H(i,i) -= H(i,k)*H(j,k);
				H(i,j) /= H(j,j);
			}
			
			H(i,i) = A(i,i);
			for ( unsigned int k = 0u ; k < i-1u ; ++k )
				H(i,i) -= H(i,k)*H(i,k);
			H(i,i) = std::sqrt(H(i,i));
			
		}
	
	}
	
};