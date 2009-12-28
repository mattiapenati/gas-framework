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

namespace poisson {

class soluzione: public gas::functional::function<2u, soluzione> {

public:
	inline double operator() (double const & x, double const & y) const {
		double const k(3.);
		double const n(10.);

		double const r(std::sqrt(x * x + y * y));

		return std::exp(k * r) * std::pow(r, n) * (1 - r);
	}

};

class problem {

private:
	class forzante: public gas::functional::function<2u, forzante> {

	public:
		inline double operator() (double const & x, double const & y) const {
			double const k(3.);
			double const n(10.);

			double const r(std::sqrt(x * x + y * y));

			double const a(n * n);
			double const b(2 * k * n + k - n * n - 2 * n - 1);
			double const c(k * k - 2 * k * n - 3 * k);
			double const d(-k * k);
			return -std::exp(k * r) * std::pow(r, n - 2) * (a + r * (b + r * (c + d * r)));
		}

	};

public:
	typedef forzante forzante_t;

private:
	typedef Eigen::SparseMatrix<double> matrix_t;
	typedef Eigen::VectorXd vector_t;

private:
	typedef gas::geometry::unit::triangle triangle_t;
	typedef gas::numerical::quadrature::newton_cotes<triangle_t, 6u> method_t;

public:
	typedef gas::numerical::quadrature::formula<method_t> integrator_t;

public:
	typedef gas::functional::base::P1<triangle_t> base_t;
	typedef gas::functional::element<base_t> element_t;

public:
	inline problem (triangulation const & cdt): cdt_(cdt), no_zeros_(0u) {
	}

	void solve ();

	double normL2 () const;

	double normH1 () const;

	template <typename function_>
	double errL2 (function_ const & f) const;

	template <typename function_>
	double errH1 (function_ const & f) const;

	inline unsigned int no_zeros () const {
		return no_zeros_;
	}

	inline unsigned int faces () const {
		return cdt_.faces();
	}

	inline double operator() (unsigned int const & i) const {
		return x_(i);
	}

	inline triangulation const & mesh() const {
		return cdt_;
	}

private:
	/*! @brief Soluzione */
	vector_t x_;

	/*! @brief Elementi non nulli */
	unsigned int no_zeros_;

	/*! @brief Pointer to the mesh */
	triangulation const & cdt_;

};

void problem::solve () {

	/* tipi */

	typedef triangulation::face_iterator_t iterator_t;
	typedef triangulation::boundary_circulator_t circulator_t;

	using gas::functional::dx;
	using gas::functional::dy;

	/* dimensione del problema */

	unsigned int const n(cdt_.nodes());
	unsigned int const N(element_t::base_t::n);

	/* riempimento delle strutture */

	#define bilinear_form(u,v) (dx(u) * dx(v)) + (dy(u) * dy(v))
	//#define bilinear_form(u,v) u * v
	#define functional(v) f * v

	Eigen::DynamicSparseMatrix<double> Atmp(n,n);
	vector_t b;
	b.setZero(n);

	integrator_t s;
	forzante_t f;

	for (iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) {
		s(*it);
		element_t e(*it);

		/* forma bilineare */
		gas_rangeu(i, N) {
			gas_rangeu(j, N) {
				unsigned int const ii(it->i(i));
				unsigned int const jj(it->i(j));

				double const r(s.integrate(bilinear_form(e.b(j), e.b(i))));

				Atmp.coeffRef(ii, jj) += r;
			}
		}

		/* termine noto */
		gas_rangeu(i, N) {
			unsigned int const ii(it->i(i));

			double const r(s.integrate(functional(e.b(i))));

			b.coeffRef(ii) += r;
		}
	}

	/* condizioni al bordo */

	circulator_t circ = cdt_.boundary();
	circulator_t begin = circ;
	do {
		unsigned int const i(circ->i());
		Atmp.coeffRef(i, i) = 1.e30;
		b(i) = 0.;
		++circ;
	} while (circ != begin);


	/* strutture algebriche */

	matrix_t A(Atmp);
	no_zeros_ = A.nonZeros();
	x_.resize(n);

	/* risoluzione del sistema con LU */
	Eigen::SparseLU<matrix_t, Eigen::UmfPack> lu(A);
	gas_assert(lu.succeeded());
	lu.solve(b, &x_);

}

double problem::normL2() const {
	double n(0.);
	{
		typedef triangulation::face_iterator_t iterator_t;

		integrator_t s;

		unsigned int const N(element_t::base_t::n);

		for (iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) {
			/* soluzione */
			Eigen::Vector3d const u(
					x_(it->i(0)),
					x_(it->i(1)),
					x_(it->i(2))
					);

			/* matrice locale */
			s(*it);
			element_t e(*it);

			Eigen::Matrix3d S;
			gas_rangeu(i, N) { gas_rangeu(j, N) {
				S(i,j) = s.integrate(e.b(i) * e.b(j));
			} }

			/* norma locale */
			n += u.dot(S * u);
		}
	}
	return std::sqrt(n);
}

double problem::normH1() const {
	double n(0.);

	typedef triangulation::face_iterator_t iterator_t;

	integrator_t s;

	unsigned int const N(element_t::base_t::n);

	for (iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) {
		/* soluzione */
		Eigen::Vector3d const u(
				x_(it->i(0)),
				x_(it->i(1)),
				x_(it->i(2))
				);

		/* matrice locale */
		s(*it);
		element_t e(*it);

		Eigen::Matrix3d S;
		gas_rangeu(i, N) { gas_rangeu(j, N) {
			using gas::functional::dx;
			using gas::functional::dy;

			S(i,j) = s.integrate(
					dx(e.b(i)) * dx(e.b(j)) +
					dy(e.b(i)) * dy(e.b(j)) +
					e.b(i) * e.b(j)
					);
		} }

		/* norma locale */
		n += u.dot(S * u);
	}

	return std::sqrt(n);
}

template <typename function_>
double problem::errL2 (function_ const & f) const {
	double n(0.);

	typedef triangulation::face_iterator_t iterator_t;

	integrator_t s;

	unsigned int const N(element_t::base_t::n);

	for (iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) {
		/* soluzione */
		Eigen::Vector3d const u(
				x_(it->i(0)),
				x_(it->i(1)),
				x_(it->i(2))
				);

		s(*it);
		element_t e(*it);

		/* matrice locale */

		Eigen::Matrix3d S;
		gas_rangeu(i, N) { gas_rangeu(j, N) {
			S(i,j) = s.integrate(e.b(i) * e.b(j));
		} }

		/* vettore locale */
		Eigen::Vector3d b;
		gas_rangeu(i, N) {
			b(i) = s.integrate(f * e.b(i));
		}

		double const f2(s.integrate(f * f));
		double const fu(u.dot(b));
		double const u2(u.dot(S * u));

		/* norma locale */
		n += std::abs(u2 + f2 - 2. * fu);
	}

	return std::sqrt(n);
}

}

#endif // _poisson_problem_
