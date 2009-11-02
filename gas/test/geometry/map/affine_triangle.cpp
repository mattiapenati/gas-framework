/*
 * Copyright (c) 2008, Politecnico di Milano
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

#undef gas_ndebug
#include "gas.h"

struct fake_triangle {
	inline fake_triangle () { }
	inline double const x (unsigned int const & i) const {
		switch (i) {
		case 0: return 3.;
		case 1: return 9.;
		case 2: return 6.;
		}
	}
	inline double const y (unsigned int const & i) const {
		switch (i) {
		case 0: return 3.;
		case 1: return 1.;
		case 2: return 8.;
		}
	}
};

#define TEST gas_geometry_map_affine_triangle

class TEST {
	public:
		TEST ();
		void execute ();
		void check ();
	private:
		gas::geometry::map::affine<gas::geometry::unit::triangle> g;
		double x, y, X, Y, det, DET, dxdX, dxdY, dydX, dydY, dXdx, dXdy, dYdx, dYdy;
};

TEST::TEST (): g(fake_triangle()) {
}

void TEST::execute () {
	x = g.x(1./3., 1./3.);
	y = g.y(1./3., 1./3.);
	X = g.X(6., 4.);
	Y = g.Y(6., 4.);

	det = g.det(1./3., 1./3.);
	DET = g.DET(6., 4.);

	dxdX = g.dxdX(1./3., 1./3.);
	dxdY = g.dxdY(1./3., 1./3.);
	dydX = g.dydX(1./3., 1./3.);
	dydY = g.dydY(1./3., 1./3.);

	dXdx = g.dXdx(6., 4.);
	dXdy = g.dXdy(6., 4.);
	dYdx = g.dYdx(6., 4.);
	dYdy = g.dYdy(6., 4.);
}

void TEST::check () {
	gas_assert(x == 6.);
	gas_assert(y == 4.);
	gas_assert(X == 1./3.);
	gas_assert(Y == 1./3.);

	gas_assert(det == 36.);
	gas_assert(DET == 1./36.);

	gas_assert(dxdX == 6.);
	gas_assert(dxdY == 3.);
	gas_assert(dydX == -2.);
	gas_assert(dydY == 5.);

	gas_assert(dXdx == 5./36.);
	gas_assert(dXdy == -1./12.);
	gas_assert(dYdx == 1./18.);
	gas_assert(dYdy == 1./6.);
}

gas_unit(TEST)
