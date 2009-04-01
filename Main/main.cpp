#include "Gas/Poisson/Poisson.hpp"

#define NP 25

#include <cmath>


double f(double const & x, double const & y) {
	return 4.;
}

int main(int argc, char **argv) {
	Poisson::PointList bound;
	double th = 0.;
	double x, y;
	for (int i = 0; i < NP; ++i) {
		th = 2. * i * M_PI / NP; 
		x = std::cos(th);
		y = std::sin(th);
		bound.push_back(Poisson::Point(x, y));
	}
	Poisson p(bound, f);
	p.saveToSVG("grid.svg");
	return 0;
}
