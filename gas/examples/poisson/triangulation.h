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

#ifndef _poisson_triangulation_
#define _poisson_triangulation_

namespace poisson {

class triangulation {

private:
	/*! @brief The information to save onto the node of triangulation */
	class point_info {

	public:
		/*! @brief The constructor */
		inline point_info (): index(0u), removable(true) {}

		/*! @brief The index of node */
		unsigned int index;

		/*! @brief A flag to check if the node is removable */
		bool removable;

	};

private:
	/*! @brief The information to save onto the face of triangulation */
	class face_info {

	public:
		/*! @brief The constructor */
		inline face_info (): index(0u) {}

		/*! @brief The index of face */
		unsigned int index;

	};

private:
	/* Kernel */
	struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

	/* Triangolazione (vertice) */
	typedef CGAL::Triangulation_vertex_base_with_info_2<point_info, K> vb_t;

	/* Triangolazione (faccia) */
	typedef CGAL::Delaunay_mesh_face_base_2<K,
	            CGAL::Constrained_triangulation_face_base_2<K,
	                 CGAL::Triangulation_face_base_with_info_2<face_info , K>
	            >
	        > fb_t;

	/* Struttura dati */
	typedef CGAL::Triangulation_data_structure_2<vb_t, fb_t> tds_t;

	/* Triangolazione */
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, tds_t> cdt_t;

private:
	class face;
	class face_iterator;

	class boundary;
	class boundary_circulator;

private:
	/*! @brief A wrapper to CGAL triangulation's face */
	class face {

	private:
		typedef cdt_t::Finite_faces_iterator cgal_face_iterator_t;

		/*! @brief The internal iterator */
		cgal_face_iterator_t it_;

		/*! @brief The constructor */
		inline face (): it_() {}
		inline face (cgal_face_iterator_t const & it): it_(it) {}

	public:
		/*! @brief The first coordinate of nodes */
		inline double x (unsigned int const & i) const {
			gas_assert(i < 3);
			return it_->vertex(i)->point().x();
		}

		/*! @brief The second coordinate of nodes */
		inline double y (unsigned int const & i) const {
			gas_assert(i < 3);
			return it_->vertex(i)->point().y();
		}

		inline unsigned int i (unsigned int const & i) const {
			gas_assert(i < 3);
			return it_->vertex(i)->info().index;
		}

		inline unsigned int i () const {
			return it_->info().index;
		}

		friend class face_iterator;
		friend class posteriorH1;

	};

	class boundary {

	private:
		typedef cdt_t::Vertex_circulator cgal_vertex_circulator_t;

		/*! @brief The internal circulator */
		cgal_vertex_circulator_t circ_;

		inline boundary (): circ_() {}
		inline boundary (cgal_vertex_circulator_t const & circ): circ_(circ) {}

	public:
		inline unsigned int i () const {
			return circ_->info().index;
		}

		friend class boundary_circulator;

	};

public:
	typedef face face_t;
	typedef boundary boundary_t;

private:
	class face_iterator {

	private:
		typedef face_iterator self_t;

	public:
		typedef face_t const * pointer_t;
		typedef face_t const & reference_t;

	public:
		inline face_iterator (): f_() {}
		inline face_iterator (self_t const & i): f_(i.f_.it_) {}

		inline self_t & operator= (self_t const & i) {
			f_.it_ = i.f_.it_;
			return *this;
		}

		inline bool operator== (self_t const & i) const { return (f_.it_ == i.f_.it_); }
		inline bool operator!= (self_t const & i) const { return (f_.it_ != i.f_.it_); }

		inline self_t & operator++ () {
			++(f_.it_);
			return *this;
		}
		inline self_t operator++ (int) {
			self_t _it(*this);
			++(f_.it_);
			return _it;
		}

		inline reference_t operator* () const { return f_; }
		inline pointer_t operator-> () const { return &f_; }

	private:
		/*! @brief The internal iterator */
		face f_;

	private:
		typedef cdt_t::Finite_faces_iterator cgal_face_iterator_t;

		inline face_iterator (cgal_face_iterator_t const & it): f_(it) {}

		friend class triangulation;
		friend class posteriorH1;

	};

	class boundary_circulator {

	private:
		typedef boundary_circulator self_t;

	public:
		typedef boundary_t const * pointer_t;
		typedef boundary_t const & reference_t;

	public:
		inline boundary_circulator (): b_() {}
		inline boundary_circulator (self_t const & i): b_(i.b_.circ_) {}

		inline self_t & operator= (self_t const & i) {
			b_.circ_ = i.b_.circ_;
			return *this;
		}

		inline bool operator== (self_t const & i) const { return (b_.circ_ == i.b_.circ_); }
		inline bool operator!= (self_t const & i) const { return (b_.circ_ != i.b_.circ_); }

		inline self_t & operator++ () {
			++(b_.circ_);
			return *this;
		}
		inline self_t operator++ (int) {
			self_t _it(*this);
			++(b_.circ_);
			return _it;
		}

		inline reference_t operator* () const { return b_; }
		inline pointer_t operator-> () const { return &b_; }

	private:
		/*! @brief The internal iterator */
		boundary b_;

	private:
		typedef cdt_t::Vertex_circulator cgal_vertex_circulator_t;

		inline boundary_circulator (cgal_vertex_circulator_t const & it): b_(it) {}

		friend class triangulation;

	};


public:
	/*! @brief The point's type */
	typedef cdt_t::Point point_t;

	/*! @brief The iterator on faces */
	typedef face_iterator face_iterator_t;
	typedef boundary_circulator boundary_circulator_t;

public:
	/*!
	 * @brief The constructor from the boundary nodes of domains
	 * @param begin The iterator the point to the begin of list of nodes
	 * @param end The iterator the point to the end of list of nodes
	 * @param h
	 */
	template <typename iterator_>
	triangulation (iterator_ begin, iterator_ end, double h = 0.);

	template <typename stimator_>
	std::pair<unsigned, unsigned> refine (stimator_ const & stimator);

	/*!
	 * @brief Number of nodes
	 * @return
	 */
	inline unsigned int nodes () const {
		return n_nodes_;
	}

	/*!
	 * @brief Number of faces
	 * @return
	 */
	inline unsigned int faces () const {
		return n_faces_;
	}

	/*!
	 * @brief The begin of const iterator on face
	 * @return
	 */
	inline face_iterator_t face_begin() const {
		return face_iterator_t(cdt_.finite_faces_begin());
	}

	/*!
	 * @brief The end of const iterator on face
	 * @return
	 */
	inline face_iterator_t face_end() const {
		return face_iterator_t(cdt_.finite_faces_end());
	}

	boundary_circulator_t boundary() const {
		return boundary_circulator_t(cdt_.incident_vertices(cdt_.infinite_vertex()));
	}

private:
	/*! @brief The internal structure for the triangulation */
	cdt_t cdt_;

	/*! @brief The number of nodes */
	unsigned int n_nodes_;

	/*! @brief The number of faces */
	unsigned int n_faces_;

	friend class posteriorH1;
};

template <typename iterator_>
triangulation::triangulation (iterator_ begin, iterator_ end, double criteria) {

	/* inserimento dei nodi e dei lati */
	{
		typedef std::vector<cdt_t::Vertex_handle> vertex_list_t;
		typedef vertex_list_t::const_iterator iterator_t;

		vertex_list_t _vertex_list;

		/* inserimento dei nodi */
		for (iterator_ _i = begin; _i != end; ++_i)
			_vertex_list.push_back(cdt_.insert(*_i));

		/* lati constrained */
		for (iterator_t _i = _vertex_list.begin() + 1; _i != _vertex_list.end(); ++_i)
			cdt_.insert_constraint(*(_i-1), *_i);
		cdt_.insert_constraint(_vertex_list.front(), _vertex_list.back());
	}

	/* calcolo lunghezza minima */
	{
		double _minimo(0.);

		for (iterator_ i = begin + 1; i != end; ++i) {
			double const _dx((i->x() - (i-1)->x()));
			double const _dy((i->y() - (i-1)->y()));
			double const _dist(_dx * _dx + _dy * _dy);
			_minimo = std::min(_dist, _minimo);
		}

		double const _dx((begin->x() - (end-1)->x()));
		double const _dy((begin->y() - (end-1)->y()));
		double const _dist(_dx * _dx + _dy * _dy);
		_minimo = std::min(_dist, _minimo);

		if (criteria <= 0. ) criteria = _minimo;
	}

	/* selezione i nodi non eliminabili */
	{
		typedef cdt_t::Vertex_circulator circulator_t;

		circulator_t _c(cdt_.incident_vertices (cdt_.infinite_vertex()));
		circulator_t _begin(_c);

		do {
			_c->info().removable = false;
			++_c;
		} while(_c != _begin);
	}

	/* raffinamento */
	{
		typedef CGAL::Delaunay_mesh_size_criteria_2<cdt_t> Criteria;

		CGAL::refine_Delaunay_mesh_2(cdt_, Criteria(0.125, criteria));
	}

	/* numerazione dei nodi */
	{
		typedef cdt_t::Finite_vertices_iterator iterator_t;

		unsigned int n(0);

		for (iterator_t _i(cdt_.finite_vertices_begin()); _i != cdt_.finite_vertices_end(); ++_i)
			_i->info().index = n++;

		n_nodes_ = n;
	}

	/* numerazione delle facce */
	{
		typedef cdt_t::Finite_faces_iterator iterator_t;

		unsigned int n(0);

		for (iterator_t i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
			i->info().index = n++;

		n_faces_ = n;
	}

}

template <typename stimator_>
std::pair<unsigned, unsigned> triangulation::refine (stimator_ const & stimator) {

	/* calcolo dei residui */
	Eigen::VectorXd res(faces());

	for (face_iterator_t i(face_begin()); i != face_end(); ++i)
		res(i->i()) = stimator.residue(*i);

	std::ofstream resFile("residui.dat", std::ios_base::app);
	resFile << res << std::endl << "=============================" << std::endl;
	resFile.close();

	/* elenco dei nodi da inserire e rimuovere */
	typedef cdt_t::Vertex_handle vertex_t;

	std::list<point_t> da_inserire;
	std::list<vertex_t> da_rimuovere;

	/* da inserire */
	for (face_iterator_t i(face_begin()); i != face_end(); ++i) {
		double const r(res(i->i()));

		if (stimator.insert(r)) {
			Eigen::Vector3d x(i->x(0), i->x(1), i->x(2));
			Eigen::Vector3d y(i->y(0), i->y(1), i->y(2));

			double const px(x.sum() / 3.);
			double const py(y.sum() / 3.);

			da_inserire.push_back(point_t(px, py));
		}
	}

	/* da rimuovere */
	for (cdt_t::Finite_vertices_iterator i(cdt_.finite_vertices_begin()); i != cdt_.finite_vertices_end(); ++i) {
		if (i->info().removable) {
			/* calcolo residuo medio sulla patch */
			double r(0.);
			unsigned int n(0);

			cdt_t::Face_circulator face_c = cdt_.incident_faces(i);
			cdt_t::Face_circulator end_c = face_c;

			do {
				if (!cdt_.is_infinite(face_c)) {
					r += res(face_c->info().index);
					++n;
				}
				++face_c;
			} while(face_c != end_c);

			r /= n;

			if (stimator.remove(r))
				da_rimuovere.push_back(i);
		}
	}

	/* rimozione */
	unsigned int const rimossi(da_rimuovere.size());
	if (rimossi) {
		for (std::list<vertex_t>::iterator it(da_rimuovere.begin()); it != da_rimuovere.end(); ++it)
			cdt_.remove(*it);
	}

	/* inserimento */
	unsigned int inseriti(da_inserire.size());
	if (inseriti) {
		cdt_.insert(da_inserire.begin(), da_inserire.end());
	}

	/* numerazione dei nodi */
	{
		typedef cdt_t::Finite_vertices_iterator iterator_t;

		unsigned int n(0);

		for (iterator_t _i(cdt_.finite_vertices_begin()); _i != cdt_.finite_vertices_end(); ++_i)
			_i->info().index = n++;

		n_nodes_ = n;
	}

	/* numerazione delle facce */
	{
		typedef cdt_t::Finite_faces_iterator iterator_t;

		unsigned int n(0);

		for (iterator_t i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
			i->info().index = n++;

		n_faces_ = n;
	}

	/* ritorno */
	return std::make_pair(rimossi, inseriti);

}

}

#endif // _poisson_triangulation_
