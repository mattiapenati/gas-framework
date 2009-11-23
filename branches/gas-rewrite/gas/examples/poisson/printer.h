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

#ifndef _poisson_printer_
#define _poisson_printer_

#include <gas>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "problem.h"
#include "triangulation.h"

namespace poisson {

class svg {

public:
	svg (problem const & p): p_(p), x_(1024u), y_(768u) {
	}

private:
	problem const & p_;

	unsigned int x_, y_;

	std::string to_string ();

	friend std::ostream & operator<< (std::ostream & out, svg & text);

};

std::ostream & operator<< (std::ostream & out, svg & text) {
	out << text.to_string();
	return out;
}

std::string svg::to_string () {

	/* file di stile */
	gas::css style;

	style.new_id("grid");
		style.property("stroke", "white");
		style.property("stroke-width", 2);
		style.property("fill", "none");

	style.new_id("solution");
		style.property("fill", "red");

	gas_rangeu(i,256) {
		std::ostringstream name;
		std::ostringstream color;

		name << "gray" << i;
		color << "rgb(" << i << ", " << i << ", " << i << ")";

		style.new_class(name.str().c_str());
			style.property("fill", color.str().c_str());
	}

	/* dimensioni quadro */
	double _xmin, _xmax, _ymin, _ymax, _rmin, _rmax;
	{
		typedef triangulation::cdt_t cdt_t;
		typedef triangulation::cdt_t::Finite_vertices_iterator iterator_t;
		typedef problem::vector_t vector_t;

		cdt_t const & _cdt(p_.cdt_.cdt_);
		vector_t const & x(p_.x);

		iterator_t _i(_cdt.finite_vertices_begin());
		_xmin = _xmax = _i->point().x();
		_ymin = _ymax = _i->point().y();
		_rmin = _rmax = x.coeff(_i->info().index);

		for (; _i != _cdt.finite_vertices_end(); ++_i) {
			double const _x(_i->point().x());
			_xmin = std::min(_xmin, _x);
			_xmax = std::max(_xmax, _x);

			double const _y(_i->point().y());
			_ymin = std::min(_ymin, _y);
			_ymax = std::max(_ymax, _y);

			double const _r(x.coeff(_i->info().index));
			_rmin = std::min(_rmin, _r);
			_rmax = std::max(_rmax, _r);
		}
	}

	/* rapporto e risoluzione */
	unsigned int _x, _bx;
	unsigned int _y, _by;
	{
		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);

		double const _ratio(_dx / _dy);

		if (_ratio < 4./3.) {
			_x = 762u * _ratio;
			_y = 762u;
		} else {
			_x = 1016u;
			_y = 1016u / _ratio;
		}

		_bx = 3u * _ratio;
		_by = 3u;
	}

	/* file svg */
	gas::svg svg(_x + (2u * _bx), _y + (2u * _by));

	svg.title("Soluzione");
	svg.description("Soluzione del problema di Poisson");

	svg.style(style);

	/* soluzione */
	{
		typedef triangulation::face_iterator_t iterator_t;
		typedef problem::vector_t vector_t;

		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);
		double const _dr(_rmax - _rmin);

		vector_t const & x(p_.x);

		svg.open_group("solution");

		for (iterator_t _i(p_.cdt_.face_begin()); _i != p_.cdt_.face_end(); ++_i) {

			unsigned int const x0(_bx + _x * (_i->x(0) - _xmin) / _dx);
			unsigned int const y0(_by + _y * (_xmax - _i->y(0)) / _dy);

			unsigned int const x1(_bx + _x * (_i->x(1) - _xmin) / _dx);
			unsigned int const y1(_by + _y * (_xmax - _i->y(1)) / _dy);

			unsigned int const x2(_bx + _x * (_i->x(2) - _xmin) / _dx);
			unsigned int const y2(_by + _y * (_xmax - _i->y(2)) / _dy);

			double const r0(x(_i->i(0)));
			double const r1(x(_i->i(0)));
			double const r2(x(_i->i(0)));

			double const r((r0 + r1 + r2) / 3.);

			unsigned int c(255 * ((r-_rmin)/_dr));

			std::ostringstream color;
			color << "gray" << c;

			svg.triangle(x0, y0, x1 ,y1, x2, y2, color.str().c_str());
		}

		svg.close_group();
	}

	/* griglia */
	{
		typedef triangulation::face_iterator_t iterator_t;

		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);

		svg.open_group("grid");

		for (iterator_t _i(p_.cdt_.face_begin()); _i != p_.cdt_.face_end(); ++_i) {

			unsigned int const x0(_bx + _x * (_i->x(0) - _xmin) / _dx);
			unsigned int const y0(_by + _y * (_ymax - _i->y(0)) / _dy);

			unsigned int const x1(_bx + _x * (_i->x(1) - _xmin) / _dx);
			unsigned int const y1(_by + _y * (_ymax - _i->y(1)) / _dy);

			unsigned int const x2(_bx + _x * (_i->x(2) - _xmin) / _dx);
			unsigned int const y2(_by + _y * (_ymax - _i->y(2)) / _dy);

			svg.triangle(x0, y0, x1 ,y1, x2, y2);
		}

		svg.close_group();
	}

	/* ritorno */
	std::ostringstream out;
	out << svg;
	return out.str();

}


class ps {

public:
	ps (problem const & p): p_(p) {
	}

private:
	problem const & p_;

	std::string to_string ();

	friend std::ostream & operator<< (std::ostream & out, ps & text);

};

std::ostream & operator<< (std::ostream & out, ps & text) {
	out << text.to_string();
	return out;
}

std::string ps::to_string () {
	std::ostringstream out;

	out << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;

	/* dimensioni quadro */
	double _xmin, _xmax, _ymin, _ymax, _rmin, _rmax;
	{
		typedef triangulation::cdt_t cdt_t;
		typedef triangulation::cdt_t::Finite_vertices_iterator iterator_t;
		typedef problem::vector_t vector_t;

		cdt_t const & _cdt(p_.cdt_.cdt_);
		vector_t const & x(p_.x);

		iterator_t _i(_cdt.finite_vertices_begin());
		_xmin = _xmax = _i->point().x();
		_ymin = _ymax = _i->point().y();
		_rmin = _rmax = x.coeff(_i->info().index);

		for (; _i != _cdt.finite_vertices_end(); ++_i) {
			double const _x(_i->point().x());
			_xmin = std::min(_xmin, _x);
			_xmax = std::max(_xmax, _x);

			double const _y(_i->point().y());
			_ymin = std::min(_ymin, _y);
			_ymax = std::max(_ymax, _y);

			double const _r(x.coeff(_i->info().index));
			_rmin = std::min(_rmin, _r);
			_rmax = std::max(_rmax, _r);
		}
	}

	/* rapporto e risoluzione */
	unsigned int _x, _bx;
	unsigned int _y, _by;
	{
		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);

		double const _ratio(_dx / _dy);

		if (_ratio < 4./3.) {
			_x = 762u * _ratio;
			_y = 762u;
		} else {
			_x = 1016u;
			_y = 1016u / _ratio;
		}

		_bx = 3u * _ratio;
		_by = 3u;
	}

	out << "%%BoundingBox: 0 0 " << _x + (2u * _bx) << " " << _y + (2u * _by) << std::endl;

	/* soluzione */
	{
		typedef triangulation::face_iterator_t iterator_t;
		typedef problem::vector_t vector_t;

		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);
		double const _dr(_rmax - _rmin);

		vector_t const & x(p_.x);

		for (iterator_t _i(p_.cdt_.face_begin()); _i != p_.cdt_.face_end(); ++_i) {

			unsigned int const x0(_bx + _x * (_i->x(0) - _xmin) / _dx);
			unsigned int const y0(_by + _y * (_i->y(0) - _ymin) / _dy);

			unsigned int const x1(_bx + _x * (_i->x(1) - _xmin) / _dx);
			unsigned int const y1(_by + _y * (_i->y(1) - _ymin) / _dy);

			unsigned int const x2(_bx + _x * (_i->x(2) - _xmin) / _dx);
			unsigned int const y2(_by + _y * (_i->y(2) - _ymin) / _dy);

			double const r0(x(_i->i(0)));
			double const r1(x(_i->i(0)));
			double const r2(x(_i->i(0)));

			double const r((r0 + r1 + r2) / 3.);

			float const c(((r - _rmin) / _dr));

			out << "newpath" << std::endl;
			out << x0 << " " << y0 << " moveto" << std::endl;
			out << x1 << " " << y1 << " lineto" << std::endl;
			out << x2 << " " << y2 << " lineto" << std::endl;
			out << "closepath" << std::endl;
			out << "gsave" << std::endl;
			out << c << " " << c << " " << c << " " << "setrgbcolor" << std::endl;
			out << "fill" << std::endl;
		}
	}

	return out.str();

}


class vtk {

public:
	vtk (problem const & p): p_(p) {
	}

private:
	problem const & p_;

	std::string to_string ();

	friend std::ostream & operator<< (std::ostream & out, vtk & text);

};

std::ostream & operator<< (std::ostream & out, vtk & text) {
	out << text.to_string();
	return out;
}

std::string vtk::to_string () {
	std::ostringstream out;

	out << "# vtk DataFile Version 2.0" << std::endl;
	out << "Soluzione" << std::endl;
	out << "ASCII" << std::endl;
	out << "DATASET UNSTRUCTURED_GRID" << std::endl;

	/* soluzione */
	{
		typedef triangulation::cdt_t cdt_t;
		typedef cdt_t::Finite_vertices_iterator iterator_t;

		cdt_t const & cdt(p_.cdt_.cdt_);

		out << "POINTS " << p_.cdt_.nodes() << " float" << std::endl;
		for (iterator_t it(cdt.finite_vertices_begin()); it != cdt.finite_vertices_end(); ++it)
			out << it->point().x() << " " << it->point().y() << " " << p_.x.coeff(it->info().index) << std::endl;
	}
	{
		typedef triangulation::face_iterator_t iterator_t;

		triangulation const & cdt(p_.cdt_);
		unsigned int const N(p_.cdt_.faces());

		out << std::endl  << "CELLS " << N << " " << 4 * N << std::endl;

		for (iterator_t it(cdt.face_begin()); it != cdt.face_end(); ++it)
			out << "3 " << it->i(0) << " " << it->i(1) << " " << it->i(2) << std::endl;

		out << std::endl << "CELL_TYPES " << N << std::endl;
		gas_rangeu(i, N)
			out << "5" << std::endl;
	}

	return out.str();

}

}

#endif // _poisson_printer_
