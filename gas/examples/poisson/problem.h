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

#ifndef _poisson_problem_
#define _poisson_problem_

#include <gas>

#include "triangulation.h"
#include "printer.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace poisson {

class problem {

private:
	typedef Eigen::SparseMatrix<double> matrix_t;
	typedef Eigen::VectorXd vector_t;

public:
	inline problem (triangulation const & cdt): cdt_(cdt), no_zeros_(0u) {
		std::cerr << "costruttore";
	}

	void solve();

	inline unsigned int no_zeros() const {
		return no_zeros_;
	}

private:
	/*! @brief Triangolazione */
	triangulation const & cdt_;

	/*! @brief Soluzione */
	vector_t x;

	/*! @brief Elementi non nulli */
	unsigned int no_zeros_;

	friend class svg;
	friend class ps;
	friend class vtk;

};

void problem::solve() {

	/* dimensione del problema */
	unsigned int const n(cdt_.nodes());

	/* strutture algebriche */
	matrix_t A(n,n);
	vector_t b(n);

	/* riempimento delle strutture */
	std::cerr << "strutture ";
	{

		/* tipi */
		typedef triangulation::face_iterator_t iterator_t;

		typedef gas::geometry::unit::triangle triangle_t;
		typedef gas::numerical::quadrature::newton_cotes<triangle_t, 6u> method_t;
		typedef gas::numerical::quadrature::formula<method_t> integrator_t;

		typedef gas::functional::base::P1<triangle_t> base_t;
		typedef gas::functional::element<base_t> element_t;

		using gas::functional::dx;
		using gas::functional::dy;

		/* A(u,v) = F(v) \forall v */
		#define A(u,v) dx(u) * dx(v) + dy(u) * dy(v)
		#define F(v) f * v

		Eigen::DynamicSparseMatrix<double> Atmp(n,n);
		integrator_t s;

		forzante_t f;

		for (iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) {
			s(*it);
			element_t e(*it);
			gas_rangeu(i, 3u) {
				gas_rangeu(j, 3u) {
					/* forma bilineare */
					double const r(s.integrate(A(e.b(j), e.b(i))));
					Atmp.coeffRef(it->i(i), it->i(j)) += r;
				}
				/* termine noto */
				double const r(s.integrate(F(e.b(i))));
				b.coeffRef(it->i(i)) += r;
			}
		}

		A = Atmp;
	}

	no_zeros_ = A.nonZeros();

	/* condizioni al bordo */
	std::cerr << "condizioni ";
	{
		typedef triangulation::boundary_circulator_t circulator_t;

		circulator_t circ = cdt_.boundary();
		circulator_t begin = circ;
		do {
			unsigned int const i(circ->i());
			A.coeffRef(i, i) = 1.e20;
			b.coeffRef(i) = 0.;
			++circ;
		} while (circ != begin);
	}

	/* risoluzione del sistema */
	std::cerr << "risoluzione";
	{
		x.resize(n);
#ifdef EIGEN_UMFPACK_SUPPORT
		Eigen::SparseLU<matrix_t,Eigen::UmfPack> lu(A);
#else
		Eigen::SparseLU<matrix_t> lu(A);
#endif
		lu.solve(b, &x);
	}

}

}

#endif // _poisson_problem_
