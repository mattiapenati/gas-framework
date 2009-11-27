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

#ifndef _poisson_posterior_
#define _poisson_posterior_

namespace poisson {

class posteriorH1 {

public:
	inline posteriorH1 (problem const & p, double const & eps)
			: p_(p),
			  eps_(eps),
			  max_res((3. * eps) / (2. * std::sqrt(p_.faces()))),
			  min_res(max_res / 3.) {
		/* calcolo norma H1 */
		// TODO
	}

	template <typename face_>
	double residue (face_ const & face) const;

	inline bool insert (double const & res) const {
		return (res > max_res);
	}

	inline bool remove (double const & res) const {
		return (res < min_res);
	}

private:
	problem const & p_;

	double const & eps_;

	double const max_res;
	double const min_res;

};

template <typename face_>
double posteriorH1::residue (face_ const & face) const {

	/* coordinate dei nodi */
	Eigen::Vector2d const p0(face.x(0), face.y(0));
	Eigen::Vector2d const p1(face.x(1), face.y(1));
	Eigen::Vector2d const p2(face.x(2), face.y(2));

	/* punto medio */
	Eigen::Vector2d const m(
			(p0(0) + p1(0) + p2(0)) / 3.,
			(p0(1) + p1(1) + p2(1)) / 3.
			);

	/* area */
	Eigen::Matrix3d tri;
	tri << p0(0), p0(1), 1., p1(0), p1(1), 1., p2(0), p2(1), 1.;
	double const area(0.5 * std::abs(tri.determinant()));

	/* distanze */
	double const d01((p0-p1).norm());
	double const d12((p1-p2).norm());
	double const d02((p0-p2).norm());
	double const h(std::max(d01, std::max(d12, d02)));

	/* soluzione */
	Eigen::Vector3d const u(
			p_(face.i(0)),
			p_(face.i(1)),
			p_(face.i(2))
			);

	/* matrice del gradiente */
	problem::element_t e(face);
	Eigen::Matrix<double, 2, 3> GG;
	gas_rangeu(i, 3)
		GG.col(i) << dx(e.b(i))(m(0), m(1)), dy(e.b(i))(m(0), m(1));

	/* gradiente locale */
	Eigen::Vector2d const g(GG * u);

	/* gradiente Clement */
	Eigen::Vector2d const gC0(0., 0.); // TODO
	Eigen::Vector2d const gC1(0., 0.); // TODO
	Eigen::Vector2d const gC2(0., 0.); // TODO

	/* laplaciano */
	Eigen::Vector3d const g0(gC0(0), gC1(0), gC2(0));
	Eigen::Vector3d const g1(gC0(1), gC1(1), gC2(1));
	Eigen::Vector2d const l0(GG * g0);
	Eigen::Vector2d const l1(GG * g1);

	double const lap(l0(0) + l1(1));

	/* residuo sulla faccia */
	problem::integrator_t s;
	problem::forzante_t f;
	s(face);

	double const res_faccia = h * std::sqrt(
			s.integrate(f * f)  +
			2 * lap * s.integrate(f) +
			lap * lap * area
			);

	/* residuo sui lati */
	// TODO
	double res_lato(0.);
	res_lato = std::sqrt(h * res_lato) / 2.;

	/* residuo totale */
	double const res_totale(res_faccia + res_lato);

	// TODO
	return (max_res+min_res)*0.5;

}

}

#endif // _poisson_posterior_
