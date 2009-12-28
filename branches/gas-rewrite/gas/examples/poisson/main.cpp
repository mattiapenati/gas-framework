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

	unsigned const N(400);
	for (unsigned i(0); i < N; ++i)
		boundary.push_back(point_t(
				std::cos(2 * M_PI * i / N),
				std::sin(2 * M_PI * i / N)
				));

	/* costruzione della triangolazione */
	poisson::triangulation mesh(boundary.begin(), boundary.end(), 0.025);

	/* costruzione del problema */
	poisson::problem problem(mesh);

	/* risoluzione del problema */
	problem.solve();

	{
		poisson::ps ps(problem);
		std::ofstream out_ps("iniziale.ps");
		out_ps << ps;
	}

	/* raffinamento */
	std::pair<unsigned, unsigned> n = std::make_pair(0u, 0u);
	std::pair<unsigned, unsigned> nold;
	unsigned int i(0);

	unsigned int const max(100);
	double const eps(0.5);

	do {

		nold.first = n.second;
		nold.second = n.first;

		poisson::posteriorH1 stimator(problem, eps);
		n = mesh.refine(stimator);
		problem.solve();

		std::stringstream ss;
		ss << "soluzione" << i << ".eps";

		poisson::ps ps(problem);
		std::ofstream out_ps(ss.str().c_str());
		out_ps << ps;

		std::cout << i << ", " << n.first << ", " << n.second << std::endl;

		++i;
	} while(
			((n.first + n.second) > 0) &&
			(i < max) &&
			(nold != n)
		);

	timer.stop();

	/* stampa della soluzione */
	{
		poisson::svg svg(problem);
		poisson::ps ps(problem);
		poisson::vtk vtk(problem);

		std::ofstream out_svg("solution.svg");
		out_svg << svg;

		std::ofstream out_ps("solution.eps");
		out_ps << ps;

		std::ofstream out_vtk("solution.vtk");
		out_vtk << vtk;
	}

	/* stampo informazioni */
	std::cerr << "-- Tempo impiegato: " << timer.elapsed() << std::endl;
	std::cerr << "-- Numero di nodi: " << mesh.nodes() << std::endl;
	std::cerr << "-- Numero delle facce: " << mesh.faces() << std::endl;
	std::cerr << "-- Elementi non nulli: " << problem.no_zeros() << std::endl;
	std::cerr << "-- Norma L2: " << problem.normL2() << std::endl;
	std::cerr << "-- Norma H1: " << problem.normH1() << std::endl;
	std::cerr << "-- Errore L2: " << problem.errL2(poisson::soluzione()) << std::endl;

	return 0;
}
