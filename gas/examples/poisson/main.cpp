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

#include "poisson.h"

int main (int argc, char * argv[]) {

	typedef poisson::triangulation::point_t point_t;

	gas::chrono timer;

	timer.start();

	/* definizione del bordo */
	std::vector<point_t> boundary;
	boundary.push_back(point_t(+1., +1.));
	boundary.push_back(point_t(-1., +1.));
	boundary.push_back(point_t(-1., -1.));
	boundary.push_back(point_t(+1., -1.));

	/* costruzione della triangolazione */
	poisson::triangulation mesh(boundary.begin(), boundary.end(), 0.05);

	/* costruzione del problema */
	poisson::problem problem(mesh);

	/* risoluzione del problema */
	problem.solve();

	/* raffinamento */
	unsigned int n(0);
	unsigned int i(0);

	unsigned int const max(200);
	double const eps(1.e-2);

	poisson::posteriorH1 criteria(problem, eps);

	do {
		n = mesh.refine(criteria);
		problem.solve();
		++i;
	} while((n > 0) && (i < max));

	timer.stop();

	/* stampa della soluzione */
	poisson::svg svg(problem);
	poisson::ps ps(problem);
	poisson::vtk vtk(problem);

	std::ofstream out_svg("solution.svg");
	out_svg << svg;

	std::ofstream out_ps("solution.eps");
	out_ps << ps;

	std::ofstream out_vtk("solution.vtk");
	out_vtk << vtk;

	/* stampo informazioni */
	std::cout << "-- Tempo impiegato: " << timer.elapsed() << std::endl;
	std::cout << "-- Numero di nodi: " << mesh.nodes() << std::endl;
	std::cout << "-- Numero delle facce: " << mesh.faces() << std::endl;
	std::cout << "-- Elementi non nulli: " << problem.no_zeros() << std::endl;
	std::cout << "-- Norma L2: " << problem.normL2() << std::endl;
	std::cout << "-- Norma H1: " << problem.normH1() << std::endl;
	std::cout << "-- Errore L2: " << problem.errL2(poisson::soluzione()) << std::endl;

	return 0;
}
