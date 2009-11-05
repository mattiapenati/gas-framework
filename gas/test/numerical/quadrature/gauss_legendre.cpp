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

#include <cmath>

#undef gas_ndebug
#include "gas.h"

template <int n>
inline double f (double const & x) {
	return x * f<n-1>(x);
}

template <>
inline double f<0> (double const & x) {
	return 1.;
}

#define TEST gas_numerical_quadrature_gauss_legendre_test

class TEST {
	public:
		TEST ();
		void execute ();
		void check ();
	private:
		typedef gas::geometry::unit::interval interval;
		typedef gas::numerical::quadrature::gauss_legendre<interval, 2u> method2;
		typedef gas::numerical::quadrature::gauss_legendre<interval, 3u> method3;
		typedef gas::numerical::quadrature::formula<method2> formula2;
		typedef gas::numerical::quadrature::formula<method3> formula3;

		formula2 q2;
		formula3 q3;

		double r[2];
};

TEST::TEST () {
	r[0] = q2.integrate(f<2>);
	r[1] = q3.integrate(f<4>);
}

void TEST::execute () {
}

void TEST::check () {
}

gas_unit(TEST)
