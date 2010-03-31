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

	int const N(400);
	boundary.reserve(N);
	for (int i(0); i < N; ++i)
		boundary.push_back(point_t(
				std::cos(2 * M_PI * i / N),
				std::sin(2 * M_PI * i / N)
				));

	/* costruzione della triangolazione */
	poisson::triangulation mesh(boundary.begin(), boundary.end(), 0.25);

	/* costruzione del problema */
	poisson::problem problem(mesh);

	/* risoluzione del problema */
	problem.solve();

	{	/* stampa della soluzione */
		poisson::svg out(problem);
		std::ofstream out_file("iniziale.svg");
		out_file << out;
	}

	/* raffinamento */
	std::pair<int, int> n = std::make_pair(0, 0);
	std::pair<int, int> nold;
	int i(0);

	int const max(100);
	double const eps(0.5);

	do {

		nold.first = n.second;
		nold.second = n.first;

		poisson::posteriorH1 stimator(problem, eps);
		n = mesh.refine(stimator);
		problem.solve();

		std::stringstream ss;
		ss << "soluzione" << i << ".svg";

		poisson::svg out(problem);
		std::ofstream out_file(ss.str().c_str());
		out_file << out;

		std::cout << i << ", " << n.first << ", " << n.second << std::endl;

		++i;
	} while(
			((n.first + n.second) > 0) &&
			(i < max) &&
			(nold != n)
		);

	timer.stop();

	{	/* stampa della soluzione */
		poisson::svg out(problem);
		std::ofstream out_file("solution.svg");
		out_file << out;
	}

	/* stampo informazioni */
	std::cerr << "-- Tempo impiegato: " << timer << std::endl;
	std::cerr << "-- Numero di nodi: " << mesh.nodes() << std::endl;
	std::cerr << "-- Numero delle facce: " << mesh.faces() << std::endl;
	std::cerr << "-- Elementi non nulli: " << problem.no_zeros() << std::endl;
	std::cerr << "-- Norma L2: " << problem.normL2() << std::endl;
	std::cerr << "-- Norma H1: " << problem.normH1() << std::endl;
	std::cerr << "-- Errore L2: " << problem.errL2(poisson::soluzione()) << std::endl;

	return 0;
}
