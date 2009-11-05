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

#undef gas_ndebug
#include "gas.h"

#define TEST gas_numerical_quadrature_newtoncotes_test

struct fake_triangle {
	inline fake_triangle () { }
	inline double const x (unsigned int const & i) const {
		switch (i) {
		case 0: return 3.;
		case 1: return 9.;
		}
		return 6.;
	}
	inline double const y (unsigned int const & i) const {
		switch (i) {
		case 0: return 3.;
		case 1: return 1.;
		}
		return 8.;
	}
};

inline double f1(double const & x, double const & y) {
	return x*y;
}

inline double f2(double const & x, double const & y) {
	return f1(x,y)*f1(x,y);
}

class TEST {
	public:
		TEST ();
		void execute ();
		void check ();
	private:
		typedef gas::geometry::unit::triangle triangle;

		typedef gas::numerical::quadrature::newton_cotes<triangle, 3u> method3;
		typedef gas::numerical::quadrature::newton_cotes<triangle, 6u> method6;

		typedef gas::numerical::quadrature::formula<method3> formula3;
		typedef gas::numerical::quadrature::formula<method6> formula6;

		formula3 q3;
		formula6 q6;

		double r[2];
};

TEST::TEST () {
}

void TEST::execute () {
	fake_triangle i;

	r[0] = q3(i).integrate(f1);
	r[1] = q6(i).integrate(f2);
}

void TEST::check () {
	double e[2];

	e[0] = std::abs(r[0] - 423.)/423.;
	e[1] = std::abs(r[1] - 11394.)/11394.;

	double const eps(1.e-12);

	gas_assert(e[0] < eps);
	gas_assert(e[1] < eps);
}

gas_unit(TEST)
