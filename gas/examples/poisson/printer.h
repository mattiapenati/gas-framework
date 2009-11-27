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
		typedef triangulation::face_iterator_t iterator_t;

		triangulation const & cdt(p_.mesh());

		iterator_t i(cdt.face_begin());

		Eigen::Vector3d const x(i->x(0), i->x(1), i->x(2));
		Eigen::Vector3d const y(i->y(0), i->y(1), i->y(2));
		Eigen::Vector3d const r(p_(i->i(0)), p_(i->i(1)), p_(i->i(2)));

		_xmin = x.minCoeff();
		_xmax = x.maxCoeff();

		_ymin = y.minCoeff();
		_ymax = y.maxCoeff();

		_rmin = r.minCoeff();
		_rmax = r.maxCoeff();

		++i;

		for (; i != cdt.face_end(); ++i) {
			Eigen::Vector3d const x(i->x(0), i->x(1), i->x(2));
			Eigen::Vector3d const y(i->y(0), i->y(1), i->y(2));
			Eigen::Vector3d const r(p_(i->i(0)), p_(i->i(1)), p_(i->i(2)));

			_xmin = std::min(_xmin, x.minCoeff());
			_xmax = std::max(_xmax, x.maxCoeff());

			_ymin = std::min(_ymin, y.minCoeff());
			_ymax = std::max(_ymax, y.maxCoeff());

			_rmin = std::min(_rmin, r.minCoeff());
			_rmax = std::max(_rmax, r.maxCoeff());
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

		triangulation const & cdt(p_.mesh());

		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);
		double const _dr(_rmax - _rmin);

		svg.open_group("solution");

		for (iterator_t i(cdt.face_begin()); i != cdt.face_end(); ++i) {
			unsigned int x[3];
			unsigned int y[3];

			x[0] = _bx + _x * (i->x(0) - _xmin) / _dx;
			x[1] = _bx + _x * (i->x(1) - _xmin) / _dx;
			x[2] = _bx + _x * (i->x(2) - _xmin) / _dx;

			y[0] = _by + _y * (_ymax - i->y(0)) / _dy;
			y[1] = _by + _y * (_ymax - i->y(1)) / _dy;
			y[2] = _by + _y * (_ymax - i->y(2)) / _dy;

			double const r((p_(i->i(0)) + p_(i->i(1)) + p_(i->i(2))) / 3.);
			unsigned int c(255 * (r - _rmin) / _dr);

			std::ostringstream color;
			color << "gray" << c;

			svg.triangle(x[0], y[0], x[1] ,y[1], x[2], y[2], color.str().c_str());
		}

		svg.close_group();
	}

	/* griglia */
	{
		typedef triangulation::face_iterator_t iterator_t;

		triangulation const & cdt(p_.mesh());

		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);

		svg.open_group("grid");

		for (iterator_t i(cdt.face_begin()); i != cdt.face_end(); ++i) {
			unsigned int x[3];
			unsigned int y[3];

			x[0] = _bx + _x * (i->x(0) - _xmin) / _dx;
			x[1] = _bx + _x * (i->x(1) - _xmin) / _dx;
			x[2] = _bx + _x * (i->x(2) - _xmin) / _dx;

			y[0] = _by + _y * (_ymax - i->y(0)) / _dy;
			y[1] = _by + _y * (_ymax - i->y(1)) / _dy;
			y[2] = _by + _y * (_ymax - i->y(2)) / _dy;

			svg.triangle(x[0], y[0], x[1] ,y[1], x[2], y[2]);
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

	double _xmin, _xmax, _ymin, _ymax, _rmin, _rmax;
	{
		typedef triangulation::face_iterator_t iterator_t;

		triangulation const & cdt(p_.mesh());

		iterator_t i(cdt.face_begin());

		Eigen::Vector3d const x(i->x(0), i->x(1), i->x(2));
		Eigen::Vector3d const y(i->y(0), i->y(1), i->y(2));
		Eigen::Vector3d const r(p_(i->i(0)), p_(i->i(1)), p_(i->i(2)));

		_xmin = x.minCoeff();
		_xmax = x.maxCoeff();

		_ymin = y.minCoeff();
		_ymax = y.maxCoeff();

		_rmin = r.minCoeff();
		_rmax = r.maxCoeff();

		++i;

		for (; i != cdt.face_end(); ++i) {
			Eigen::Vector3d const x(i->x(0), i->x(1), i->x(2));
			Eigen::Vector3d const y(i->y(0), i->y(1), i->y(2));
			Eigen::Vector3d const r(p_(i->i(0)), p_(i->i(1)), p_(i->i(2)));

			_xmin = std::min(_xmin, x.minCoeff());
			_xmax = std::max(_xmax, x.maxCoeff());

			_ymin = std::min(_ymin, y.minCoeff());
			_ymax = std::max(_ymax, y.maxCoeff());

			_rmin = std::min(_rmin, r.minCoeff());
			_rmax = std::max(_rmax, r.maxCoeff());
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

		triangulation const & cdt(p_.mesh());

		double const _dx(_xmax - _xmin);
		double const _dy(_ymax - _ymin);
		double const _dr(_rmax - _rmin);

		for (iterator_t i(cdt.face_begin()); i != cdt.face_end(); ++i) {
			unsigned int x[3];
			unsigned int y[3];

			x[0] = _bx + _x * (i->x(0) - _xmin) / _dx;
			x[1] = _bx + _x * (i->x(1) - _xmin) / _dx;
			x[2] = _bx + _x * (i->x(2) - _xmin) / _dx;

			y[0] = _by + _y * (i->y(0) - _ymin) / _dy;
			y[1] = _by + _y * (i->y(1) - _ymin) / _dy;
			y[2] = _by + _y * (i->y(2) - _ymin) / _dy;

			double const r((p_(i->i(0)) + p_(i->i(1)) + p_(i->i(2))) / 3.);
			float const c(((r - _rmin) / _dr));

			out << "newpath" << std::endl;
			out << x[0] << " " << y[0] << " moveto" << std::endl;
			out << x[1] << " " << y[1] << " lineto" << std::endl;
			out << x[2] << " " << y[2] << " lineto" << std::endl;
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
		typedef triangulation::face_iterator_t iterator_t;

		triangulation const & cdt(p_.mesh());

		unsigned int const M(cdt.nodes());
		unsigned int const N(cdt.faces());

		Eigen::VectorXd x(M);
		Eigen::VectorXd y(M);

		for (iterator_t i(cdt.face_begin()); i != cdt.face_end(); ++i) {
			gas_rangeu(j, 3) {
				x(i->i(j)) = i->x(j);
				y(i->i(j)) = i->y(j);
			}
		}

		out << "POINTS " << M << " float" << std::endl;
		gas_rangeu(i, M)
			out << x(i) << " " << y(i) << " " <<  p_(i) << std::endl;

		out << std::endl;

		out << "CELLS " << N << " " << 4 * N << std::endl;
		for (iterator_t it(cdt.face_begin()); it != cdt.face_end(); ++it)
			out << "3 " << it->i(0) << " " << it->i(1) << " " << it->i(2) << std::endl;

		out << std::endl;

		out << "CELL_TYPES " << N << std::endl;
		gas_rangeu(i, N)
			out << "5" << std::endl;
	}

	return out.str();

}

}

#endif // _poisson_printer_
