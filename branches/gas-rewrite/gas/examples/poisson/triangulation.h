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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <iterator>
#include <algorithm>
#include <vector>

#include "printer.h"

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
	class face_info {
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

		inline unsigned int index (unsigned int const & i) const {
			gas_assert(i < 3);
			return it_->vertex(i)->info().index;
		}

		friend class face_iterator;

	};

public:
	typedef face face_t;

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

	};


public:
	/*! @brief The point's type */
	typedef cdt_t::Point point_t;

	/*! @brief The iterator on faces */
	typedef face_iterator face_iterator_t;

public:
	/*!
	 * @brief The constructor from the boundary nodes of domains
	 * @param begin The iterator the point to the begin of list of nodes
	 * @param end The iterator the point to the end of list of nodes
	 * @param criteria
	 */
	template <typename iterator_>
	triangulation (iterator_ begin, iterator_ end, double criteria = 0.);

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
	face_iterator_t face_end() const {
		return face_iterator_t(cdt_.finite_faces_end());
	}

private:
	/*! @brief The internal structure for the triangulation */
	cdt_t cdt_;

	friend class svg;

};

template <typename iterator_>
triangulation::triangulation (iterator_ begin, iterator_ end, double criteria) {
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

		unsigned int _n(0);

		for (iterator_t _i(cdt_.finite_vertices_begin()); _i != cdt_.finite_vertices_end(); ++_i)
			_i->info().index = _n++;
	}
}

}

#endif // _poisson_triangulation_
