#define CGAL_NO_PRECONDITIONS

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

double sol(double const & x , double const & y) {
	return u(x,y) * v(x,y);
	//return u(x,y) / (M_PI * M_PI);
}

double sol_x(double const & x, double const & y) {
	return M_PI * std::cos(M_PI*x) * std::sin(M_PI*y) * v(x,y) + u(x,y) * 10 * v(x,y);
	//return std::cos( M_PI * x ) * std::sin( M_PI * y ) / M_PI;
}

double sol_y(double const & x, double const & y) {
	return M_PI * std::sin(M_PI*x) * std::cos(M_PI*y) * v(x,y);
	//return std::sin( M_PI * x ) * std::cos( M_PI * y ) / M_PI;
}

double zero(double const & x, double const & y) {
	return 0.;
}

int main(int argc, char **argv) {
	Poisson::PointList bound;
	bound.push_back(Poisson::Point(+1, +1));
	bound.push_back(Poisson::Point(-1, +1));
	bound.push_back(Poisson::Point(-1, -1));
	bound.push_back(Poisson::Point(+1, -1));

	Poisson p(bound, f, 1.);

	// Errori
	double errL2 = p.errL2(sol);
	double norL2 = p.normL2();
	double errH1 = p.errH1(sol, sol_x, sol_y);
	double norH1 = p.normH1();

	// Output
	std::cout << "====== Iterata " << 0 << " ======" << std::endl;
	std::cout << "Errore L2: " << errL2 << std::endl;
	std::cout << "Errore relativo L2: " << errL2/norL2 << std::endl;
	std::cout << "Errore H1: " << errH1 << std::endl;
	std::cout << "Errore relativo H1: " << errH1/norH1 << std::endl;

	// Grafico
	p.saveToSVG("grid0.svg", true);

	unsigned int inseriti;
	unsigned int iterata = 1;

	do {
		// Raffinamento
		inseriti = p.refine<PosteriorH1<Poisson::CDT> >(1.e-2);

		// Errori
		errL2 = p.errL2(sol);
		norL2 = p.normL2();
		errH1 = p.errH1(sol, sol_x, sol_y);
		norH1 = p.normH1();

		// Output
		std::cout << "====== Iterata " << iterata << " ======" << std::endl;
		std::cout << "Raffinamento: " << inseriti << " punti inseriti/eliminati" << std::endl;
		std::cout << "Errore L2: " << errL2 << std::endl;
		std::cout << "Errore relativo L2: " << errL2/norL2 << std::endl;
		std::cout << "Errore H1: " << errH1 << std::endl;
		std::cout << "Errore relativo H1: " << errH1/norH1 << std::endl;

		// Grafico
		std::ostringstream s1;
		s1 << "grid" << iterata << ".svg";
		p.saveToSVG(s1.str().data(), true);

		// Avanzamento
		++iterata;
	} while(inseriti);

	// p.saveToSVG("grid.svg", true);

	return 0;
}
