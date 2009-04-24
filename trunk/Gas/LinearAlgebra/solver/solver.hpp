template < typename Method >
struct solver {

	template < typename Type , unsigned int Dimension >
	static void solve ( Matrix < Type , Dimension , Dimension > const & A , Vector < Type , Dimension > & x , Vector < Type , Dimension > const & b ) {
		Method::template apply< Type , Dimension >(A , x , b);
	}

};