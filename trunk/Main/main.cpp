#include "Gas/Poisson/Poisson.hpp"

#define NP 25

#include <cmath>


double f(double const & x, double const & y) {
	return std::sin(2 * x)*std::sin(2 * y);
}

int main(int argc, char **argv) {
	Poisson::PointList bound;
	bound.push_back(Poisson::Point(+M_PI, +M_PI));
	bound.push_back(Poisson::Point(-M_PI, +M_PI));
	bound.push_back(Poisson::Point(-M_PI, -M_PI));
	bound.push_back(Poisson::Point(+M_PI, -M_PI));
	Poisson p(bound, f, 0.2);
	p.saveToSVG("grid.svg", false);
	return 0;
}
