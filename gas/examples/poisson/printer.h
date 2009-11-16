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

	friend std::ostream & operator<< (std::ostream & out, svg & style);

};

std::ostream & operator<< (std::ostream & out, svg & style) {
	out << style.to_string();
	return out;
}

std::string svg::to_string () {
	/* file di stile */
	gas::css style;

	style.new_id("grid");
		style.property("stroke", "white");
		style.property("stroke-width", 2);

	style.new_class("solution");

	/* dimensioni quadro */
	double _xmin, _xmax, _ymin, _ymax;
	{
		typedef triangulation::cdt_t cdt_t;
		typedef triangulation::cdt_t::Finite_vertices_iterator iterator_t;

		cdt_t const & _cdt(p_.cdt_.cdt_);

		iterator_t _i(_cdt.finite_vertices_begin());
		_xmin = _xmax = _i->point().x();
		_ymin = _ymax = _i->point().y();

		for (; _i != _cdt.finite_vertices_end(); ++_i) {
			double const _x(_i->point().x());
			_xmin = std::min(_xmin, _x);
			_xmax = std::max(_xmax, _x);

			double const _y(_i->point().y());
			_ymin = std::min(_ymin, _y);
			_ymax = std::max(_ymax, _y);
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

	/* griglia */
	{
		svg.open_group("grid");

		typedef triangulation::face_iterator_t iterator_t;
		typedef gas::numerical::tiny::vector<6u> vector_t;

		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);

		for (iterator_t _i(p_.cdt_.face_begin()); _i != p_.cdt_.face_end(); ++_i) {
			vector_t _v;
			_v(0) = _i->x(0);
			_v(1) = _i->y(0);
			_v(2) = _i->x(1);
			_v(3) = _i->y(1);
			_v(4) = _i->x(2);
			_v(5) = _i->y(2);

			unsigned int const x0(_bx + _x * (_v(0) - _xmin) / _dx);
			unsigned int const y0(_by + _y * (_xmax - _v(1)) / _dy);

			unsigned int const x1(_bx + _x * (_v(2) - _xmin) / _dx);
			unsigned int const y1(_by + _y * (_xmax - _v(3)) / _dy);

			unsigned int const x2(_bx + _x * (_v(4) - _xmin) / _dx);
			unsigned int const y2(_by + _y * (_xmax - _v(5)) / _dy);

			svg.triangle(x0, y0, x1 ,y1, x2, y2);
		}
	}

	/* ritorno */
	std::ostringstream out;
	out << svg;
	return out.str();
}

}

#endif // _poisson_printer_
