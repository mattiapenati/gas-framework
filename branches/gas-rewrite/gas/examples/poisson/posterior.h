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

#ifndef POISSON_POSTERIOR_H
#define POISSON_POSTERIOR_H

namespace poisson {

/*!
 * @brief Posterior estimator in H1 norm
 */
class posteriorH1 {

public:
	/*!
	 * @brief The constructor
	 * @param p The problem
	 * @param eps The cutoff to the error
	 */
	inline posteriorH1 (problem const & p, double const eps)
			: p_(p),
			  eps_(eps),
			  max_res_((3. * eps) / (2. * std::sqrt(p_.faces()))),
			  min_res_(max_res_ / 3.),
			  h1_norm_(p.normH1()) {
	}

	/*!
	 * @brief Evaluate the residue on a face
	 * @param face The face
	 * @return The residue
	 */
	template <typename face_>
	double residue (face_ const & face) const;

	/*!
	 * A test to check the residue is small enough to insert a new point
	 * @param res The residue
	 * @return The result of test
	 */
	inline bool insert (double const & res) const {
		return (res > max_res_);
	}

	/*!
	 * A test to check the residue is small enough to remove a point
	 * @param res The residue
	 * @return The result of test
	 */
	inline bool remove (double const & res) const {
		return (res < min_res_);
	}

private:
	/*
	 * Compute the gradient vector of a function on a triangle, from the
	 * coordinate of vertices and the values of the function in these vertices
	 */
	Eigen::Vector2d grad(
			Eigen::Vector2d const & p0, double const & u0,
			Eigen::Vector2d const & p1, double const & u1,
			Eigen::Vector2d const & p2, double const & u2) const;

	template < typename VertexPointer >
	Eigen::Vector2d gradClement(VertexPointer const & v) const;

private:
	problem const & p_;

	double const eps_;

	double const max_res_;
	double const min_res_;

	double const h1_norm_;

private:
	/*
	 * Wrap a circulator to be used with integrator
	 */
	class faceCirculatorWrapper {
	private:
		typedef triangulation::cdt_t::Face_circulator circulator;
	public:
		inline faceCirculatorWrapper(circulator const & c): c_(c) {}
		inline double x (int const i) const {
			return c_->vertex(i)->point().x();
		}
		inline double y (int const i) const {
			return c_->vertex(i)->point().y();
		}
		inline int i (int const i) const {
			return c_->vertex(i)->info().index;
		}
	private:
		circulator const & c_;
	};

private:
	struct baseClement: public gas::functional::function<2u, baseClement> {
		inline baseClement (int const & i): i(i) { }
		inline double operator() (double const & x, double const & y) const {
			switch(i) {
			case 0: return 1.;
			case 1: return x;
			case 2: return y;
			}
			return 0.;
		}
	private:
		int const i;
	};

};


/* calcolo del gradiente */
Eigen::Vector2d posteriorH1::grad(
		Eigen::Vector2d const & p0, double const & u0,
		Eigen::Vector2d const & p1, double const & u1,
		Eigen::Vector2d const & p2, double const & u2) const {

	/* matrice */
	Eigen::Matrix2d A;
	A.row(0) = (p0 - p2);
	A.row(1) = (p1 - p2);

	/* termine noto */
	Eigen::Vector2d const b(
			(u0 - u2),
			(u1 - u2));

	/* risolvo */
	return A.inverse() * b;
}


/* calcolo del gradiente di clement */
template < typename VertexPointer >
Eigen::Vector2d posteriorH1::gradClement(VertexPointer const & v) const {

	/* coordinate del punto */
	double const x(v->point().x());
	double const y(v->point().y());

	/* matrice di stiffness */
	Eigen::Matrix3d M;

	/* termine noto */
	Eigen::Vector3d bx;
	Eigen::Vector3d by;

	/* soluzione */
	Eigen::Vector3d ax, ay;

	/* integratore */
	problem::integrator_t s;

	/* circolatore */
	typedef triangulation::cdt_t::Face_circulator circulator_t;

	circulator_t circ(p_.mesh().cdt_.incident_faces(v));
	circulator_t end(circ);

	/* inizializzazione delle strutture */
	M.setZero();
	bx.setZero();
	by.setZero();

	do {
		if (!p_.mesh().cdt_.is_infinite(circ)) {
			/* wrapper */
			faceCirculatorWrapper w(circ);

			/* integratore */
			s(w);

			/* matrice M */
			gas_rangeu(i, 3) {
				gas_rangeu(j, 3) {
					M(i,j) += s.integrate(baseClement(j) * baseClement(i));
				}
			}

			/* coordinate dei nodi */
			Eigen::Vector2d const p0(w.x(0), w.y(0));
			Eigen::Vector2d const p1(w.x(1), w.y(1));
			Eigen::Vector2d const p2(w.x(2), w.y(2));

			/* soluzione */
			Eigen::Vector3d const u(p_(w.i(0)), p_(w.i(1)), p_(w.i(2)));

			/* gradiente locale */
			Eigen::Vector2d const g(grad(p0, u(0), p1, u(1), p2, u(2)));

			/* termini noti */
			gas_rangeu(i, 3) {
				double const tmp(s.integrate(baseClement(i)));
				bx(i) += g(0) * tmp;
				by(i) += g(1) * tmp;
			}
		}
		++circ;
	} while(circ != end);

	/* soluzione del sistema */
	M = M.inverse();
	ax = M * bx;
	ay = M * by;

	/* gradiente */
	Eigen::Vector2d gC(0., 0.);
	gas_rangeu(i, 3) {
		gC(0) += ax(i) * baseClement(i)(x, y);
		gC(1) += ay(i) * baseClement(i)(x, y);
	}

	return gC;
}

/* calcolo del residuo */
template <typename face_>
double posteriorH1::residue (face_ const & face) const {

	/* iteratore CGAL */
	typedef triangulation::face::cgal_face_iterator_t cgal_face_iterator_t;
	cgal_face_iterator_t const cgal_it(face.it_);

	/* coordinate dei nodi della faccia */
	Eigen::Vector2d const p0(face.x(0), face.y(0));
	Eigen::Vector2d const p1(face.x(1), face.y(1));
	Eigen::Vector2d const p2(face.x(2), face.y(2));

	/* area del triangolo */
	Eigen::Matrix3d tri;
	tri << p0(0), p0(1), 1., p1(0), p1(1), 1., p2(0), p2(1), 1.;
	double const area(0.5 * std::abs(tri.determinant()));

	/* lunghezze dei lati */
	double const d01((p0-p1).norm());
	double const d12((p1-p2).norm());
	double const d02((p0-p2).norm());
	double const h(std::max(d01, std::max(d12, d02)));

	/* soluzione sulla faccia */
	Eigen::Vector3d const u(
			p_(face.i(0)),
			p_(face.i(1)),
			p_(face.i(2)));

	/* gradiente locale */
	Eigen::Vector2d const g(grad(p0, u(0), p1, u(1), p2, u(2)));

	/* gradiente Clement su ogni vertice */
	Eigen::Vector2d const gC0(gradClement(cgal_it->vertex(0)));
	Eigen::Vector2d const gC1(gradClement(cgal_it->vertex(1)));
	Eigen::Vector2d const gC2(gradClement(cgal_it->vertex(2)));

	/* laplaciano */
	Eigen::Vector2d const l0(grad(p0, gC0(0), p1, gC1(0), p2, gC2(0)));
	Eigen::Vector2d const l1(grad(p0, gC0(1), p1, gC1(1), p2, gC2(1)));

	double const lap(l0(0) + l1(1));

	/* residuo sulla faccia */
	problem::integrator_t s;
	problem::forzante_t f;
	s(face);

	/* -lap(u) = f --> r = f + lap(u) --> r^2 */
	double const res_faccia = h * std::sqrt(
			s.integrate(f * f)  +
			2 * lap * s.integrate(f) +
			lap * lap * area
			);

	/* residuo sui lati */
	double res_lato(0.);

	/* giro sui lati */
	gas_rangeu(i, 3) {
		/* controllo se il lato è di bordo */
		if (!p_.mesh().cdt_.is_infinite(cgal_it->neighbor(i))) {
			/* coordinate dei vertici del triangolo vicino */
			Eigen::Vector2d const pn0(
					cgal_it->neighbor(i)->vertex(0)->point().x(),
					cgal_it->neighbor(i)->vertex(0)->point().y());
			Eigen::Vector2d const pn1(
					cgal_it->neighbor(i)->vertex(1)->point().x(),
					cgal_it->neighbor(i)->vertex(1)->point().y());
			Eigen::Vector2d const pn2(
					cgal_it->neighbor(i)->vertex(2)->point().x(),
					cgal_it->neighbor(i)->vertex(2)->point().y());

			/* soluzione del vicino */
			Eigen::Vector3d const un(
					p_(cgal_it->neighbor(i)->vertex(0)->info().index),
					p_(cgal_it->neighbor(i)->vertex(1)->info().index),
					p_(cgal_it->neighbor(i)->vertex(2)->info().index)
					);

			/* gradiente del vicino */
			Eigen::Vector2d const gn(grad(pn0, un(0), pn1, un(1), pn2, un(2)));

			/* normale */
			switch(i) {
			case 0:
			{
				Eigen::Vector2d const n(
						p1(1) - p2(1),
						p2(0) - p1(0));
				double const t(n.dot(g - gn));
				res_lato += (t * t) / d12;
			}
				break;
			case 1:
			{
				Eigen::Vector2d const n(
						p2(1) - p0(1),
						p0(0) - p2(0));
				double const t(n.dot(g - gn));
				res_lato += (t * t) / d02;
			}
				break;
			case 2:
			{
				Eigen::Vector2d const n(
						p0(1) - p1(1),
						p1(0) - p0(0));
				double const t(n.dot(g - gn));
				res_lato += (t * t) / d01;
			}
				break;
			}
		}
	}
	res_lato = std::sqrt(h * res_lato) / 2.;

	/* residuo totale */
	double const res_totale(res_faccia + res_lato);
	return res_totale / h1_norm_;

}

}

#endif // POISSON_POSTERIOR_H
