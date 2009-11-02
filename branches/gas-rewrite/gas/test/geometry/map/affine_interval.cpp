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

struct fake_interval {
	inline fake_interval () { }
	inline double const a() const { return 4.; }
	inline double const b() const { return 6.; }
};

#define TEST gas_geometry_map_affine_interval

class TEST {
	public:
		TEST ();
		void execute ();
		void check ();
	private:
		gas::geometry::map::affine<gas::geometry::unit::interval> i;

		double x, X, det, DET, dxdX, dXdx;
};

TEST::TEST (): i(fake_interval()) {
}

void TEST::execute () {
	X = i.X(5.);
	x = i.x(0.);
	dXdx = i.dXdx(5.);
	dxdX = i.dxdX(0.);
	DET = i.DET(5.);
	det = i.det(0.);
}

void TEST::check () {
	gas_assert(X == 0.);
	gas_assert(x == 5.);
	gas_assert(dXdx == 1.);
	gas_assert(dxdX == 1.);
	gas_assert(DET == 1.);
	gas_assert(det == 1.);
}

gas_unit(TEST)
