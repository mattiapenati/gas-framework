/* Copyright (c) 2009, Politecnico di Milano
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
#include <gas>
#include <cmath>

#define TEST gas_geometry_map_affine_tetra

struct fake_tetra {
	inline fake_tetra () {
	}
	inline double const x (unsigned int const & i) const {
		switch (i) {
		case 0:
			return 3.;
		case 1:
			return 9.;
		case 2:
			return 6.;
		}
		return 6.;
	}
	inline double const y (unsigned int const & i) const {
		switch (i) {
		case 0:
			return 3.;
		case 1:
			return 1.;
		case 2:
			return 8.;
		}
		return 4.;
	}
	inline double const z (unsigned int const & i) const {
		switch (i) {
		case 0:
			return 2.;
		case 1:
			return 4.;
		case 2:
			return 1.;
		}
		return 6.;
	}
};

class TEST {
public:
	TEST ();
	void execute ();
	void check ();
private:
	typedef gas::geometry::map::affine<gas::geometry::unit::tetra> map;

	map g;

	gas::numerical::tiny::matrix<3u, 3u> A;
	gas::numerical::tiny::matrix<3u, 3u> invA;
	double x, y, z, X, Y, Z, det, DET;
};

TEST::TEST () :
	g(fake_tetra()) {
}

void TEST::execute () {
	x = g.x(0.25, 0.25, 0.25);
	y = g.y(0.25, 0.25, 0.25);
	z = g.z(0.25, 0.25, 0.25);

	X = g.X(6., 4., 3.25);
	Y = g.Y(6., 4., 3.25);
	Z = g.Z(6., 4., 3.25);

	A(0, 0) = g.dxdX(0.25, 0.25, 0.25);
	A(0, 1) = g.dxdY(0.25, 0.25, 0.25);
	A(0, 2) = g.dxdZ(0.25, 0.25, 0.25);
	A(1, 0) = g.dydX(0.25, 0.25, 0.25);
	A(1, 1) = g.dydY(0.25, 0.25, 0.25);
	A(1, 2) = g.dydZ(0.25, 0.25, 0.25);
	A(2, 0) = g.dzdX(0.25, 0.25, 0.25);
	A(2, 1) = g.dzdY(0.25, 0.25, 0.25);
	A(2, 2) = g.dzdZ(0.25, 0.25, 0.25);

	invA(0, 0) = g.dXdx(6., 4., 3.25);
	invA(0, 1) = g.dXdy(6., 4., 3.25);
	invA(0, 2) = g.dXdz(6., 4., 3.25);
	invA(1, 0) = g.dYdx(6., 4., 3.25);
	invA(1, 1) = g.dYdy(6., 4., 3.25);
	invA(1, 2) = g.dYdz(6., 4., 3.25);
	invA(2, 0) = g.dZdx(6., 4., 3.25);
	invA(2, 1) = g.dZdy(6., 4., 3.25);
	invA(2, 2) = g.dZdz(6., 4., 3.25);

	det = g.det(0.25, 0.25, 0.25);
	DET = g.DET(6., 4., 3.25);
}

void TEST::check () {
	typedef gas::geometry::map::info<map> info;
	using gas::geometry::unit::tetra;

	gas_assert(info::d == 3u);
	gas_assert((gas::same_type<info::unit_t, tetra>::value));

	double const eps(1.e-14);

	gas_assert(std::abs(x - 6.) < eps);
	gas_assert(std::abs(y - 4.) < eps);
	gas_assert(std::abs(z - 3.25) < eps);

	gas_assert(std::abs(X - 0.25) < eps);
	gas_assert(std::abs(Y - 0.25) < eps);
	gas_assert(std::abs(Z - 0.25) < eps);

	gas_assert(std::abs(A(0,0) - 6.) < eps);
	gas_assert(std::abs(A(0,1) - 3.) < eps);
	gas_assert(std::abs(A(0,2) - 3.) < eps);
	gas_assert(std::abs(A(1,0) + 2.) < eps);
	gas_assert(std::abs(A(1,1) - 5.) < eps);
	gas_assert(std::abs(A(1,2) - 1.) < eps);
	gas_assert(std::abs(A(2,0) - 2.) < eps);
	gas_assert(std::abs(A(2,1) + 1.) < eps);
	gas_assert(std::abs(A(2,2) - 4.) < eps);

	gas_assert(std::abs(invA(0,0) - 7./44.) < eps);
	gas_assert(std::abs(invA(0,1) + 5./44.) < eps);
	gas_assert(std::abs(invA(0,2) + 1./11.) < eps);
	gas_assert(std::abs(invA(1,0) - 5./66.) < eps);
	gas_assert(std::abs(invA(1,1) - 3./22.) < eps);
	gas_assert(std::abs(invA(1,2) + 1./11.) < eps);
	gas_assert(std::abs(invA(2,0) + 2./33.) < eps);
	gas_assert(std::abs(invA(2,1) - 1./11.) < eps);

	gas_assert(std::abs(invA(2,2) - 3./11.) < eps);
	gas_assert(std::abs(det - 132.) < eps);
	gas_assert(std::abs(DET - 1./132.) < eps);
}
gas_unit(TEST)

