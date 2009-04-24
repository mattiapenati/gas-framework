#include "Gas/Gas.h"

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
	Poisson p(bound, f, 0.15);
	p.saveToSVG("grid.svg", true);
	return 0;
}
