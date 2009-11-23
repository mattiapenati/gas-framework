/*
 * Copyright (c) 2009, Politecnico di Milano
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the Politecnico di Milano nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <gas>

inline double u1 (double const & x, double const & y) { return std::sin(M_PI * x) * std::sin(M_PI * y); }
inline double u2 (double const & x, double const & y) { return std::cos(M_PI * x) * std::sin(M_PI * y); }

inline double v (double const & x, double const & y) { return std::exp(10 * x); }

class forzante: public gas::functional::function<2u, forzante> {

public:
	inline double operator() (double const & x, double const & y) const {
		return ((M_PI * M_PI - 100) * u1(x, y) - 20 * M_PI * u2(x, y)) * v(x, y);
	}

};

class costante: public gas::functional::function<2u, costante> {

public:
	inline double operator() (double const & x, double const & y) const {
		return 1.;
	}

};

typedef costante forzante_t;

#include "poisson.h"

int main (int argc, char * argv[]) {

	typedef poisson::triangulation::point_t point_t;

	gas::chrono timer;
	gas::chrono local;

	timer.start();

	/* definizione del bordo */
	std::vector<point_t> boundary;
	boundary.push_back(point_t(+1., +1.));
	boundary.push_back(point_t(-1., +1.));
	boundary.push_back(point_t(-1., -1.));
	boundary.push_back(point_t(+1., -1.));

	/* costruzione della triangolazione */
	std::cerr << "-- Triangolazione (";
	local.start();
	poisson::triangulation mesh(boundary.begin(), boundary.end(), 0.1);
	local.stop();
	std::cerr << ") " << local.elapsed() << std::endl;

	/* costruzione del problema */
	std::cerr << "-- Problema (";
	local.start();
	poisson::problem problem(mesh);
	local.stop();
	std::cerr << ") " << local.elapsed() << std::endl;

	/* risoluzione del problema */
	std::cerr << "-- Risoluzione (";
	local.start();
	problem.solve();
	local.stop();
	std::cerr << ") " << local.elapsed() << std::endl;

	timer.stop();

	/* raffinamento */

	/* stampa della soluzione */
	std::cerr << "-- File di output";
	local.start();

	poisson::svg svg(problem);
	poisson::ps ps(problem);
	poisson::vtk vtk(problem);

	std::ofstream out_svg("solution.svg");
	out_svg << svg;

	std::ofstream out_ps("solution.eps");
	out_ps << ps;

	std::ofstream out_vtk("solution.vtk");
	out_vtk << vtk;

	local.stop();
	std::cerr << " " << local.elapsed() << std::endl;

	/* stampo informazioni */
	std::cerr << std::endl;
	std::cerr << "-- Tempo impiegato: " << timer.elapsed() << std::endl;
	std::cerr << "-- Numero di nodi: " << mesh.nodes() << std::endl;
	std::cerr << "-- Numero delle facce: " << mesh.faces() << std::endl;
	std::cerr << "-- Elementi non nulli: " << problem.no_zeros() << std::endl;

	return 0;
}
