#include "Gas/Gas.h"
#include <sstream>

double u ( double const & x , double const & y ) {
	return std::sin( M_PI * x ) * std::sin( M_PI * y );
}

double v ( double const & x , double const & y ) {
	return std::exp( 10 * x );
} 

double f(double const & x, double const & y) {
	return (M_PI * M_PI - 100) * u( x , y ) * v( x , y ) - 20 * M_PI * std::cos( M_PI * x ) * std::sin( M_PI * y ) * std::exp( 10 * x );
}

int main(int argc, char **argv) {
	Poisson::PointList bound;
	bound.push_back(Poisson::Point(+1, +1));
	bound.push_back(Poisson::Point(-1, +1));
	bound.push_back(Poisson::Point(-1, -1));
	bound.push_back(Poisson::Point(+1, -1));
	Poisson p(bound, f, 0.2);
	p.saveToSVG("grid.svg", true);
	unsigned int inseriti;
	unsigned int iterata = 1;
	do {
		inseriti = p.refine(6.e-3);
		std::cout << "Iterata " << iterata << ": " << inseriti << " punti inseriti/eliminati" << std::endl;
		std::ostringstream s1;
		s1 << "grid" << iterata << ".svg";		
		p.saveToSVG(s1.str().data(), true);
		++iterata;
	} while(inseriti);
	return 0;
}
