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

#include <cmath>

#undef gas_ndebug
#include "gas.h"

#define TEST gas_numerical_quadrature_gauss_legendre_test

#define A 0.
#define B 2.

struct fake_interval {
	inline fake_interval () {
	}
	inline double const a () const {
		return A;
	}
	inline double const b () const {
		return B;
	}
};

template<int n> inline double f (double const & x) {
	return x * f<n - 1> (x);
}
template<> inline double f<0> (double const & x) {
	return 1.;
}
template<int n> inline double exact () {
	return (f<n + 1> (B) - f<n + 1> (A)) / (n + 1);
}

class TEST {
public:
	TEST ();
	void execute ();
	void check ();
private:
	typedef gas::geometry::unit::interval interval;

	typedef gas::numerical::quadrature::gauss_legendre<interval, 2u> method2;
	typedef gas::numerical::quadrature::gauss_legendre<interval, 3u> method3;
	typedef gas::numerical::quadrature::gauss_legendre<interval, 4u> method4;
	typedef gas::numerical::quadrature::gauss_legendre<interval, 5u> method5;
	typedef gas::numerical::quadrature::gauss_legendre<interval, 6u> method6;
	typedef gas::numerical::quadrature::gauss_legendre<interval, 7u> method7;

	typedef gas::numerical::quadrature::formula<method2> formula2;
	typedef gas::numerical::quadrature::formula<method3> formula3;
	typedef gas::numerical::quadrature::formula<method4> formula4;
	typedef gas::numerical::quadrature::formula<method5> formula5;
	typedef gas::numerical::quadrature::formula<method6> formula6;
	typedef gas::numerical::quadrature::formula<method7> formula7;

	formula2 q2;
	formula3 q3;
	formula4 q4;
	formula5 q5;
	formula6 q6;
	formula7 q7;

	double r[6];
};

TEST::TEST () {
}

void TEST::execute () {
	fake_interval i;

	r[0] = q2(i).integrate(f<3> );
	r[1] = q3(i).integrate(f<5> );
	r[2] = q4(i).integrate(f<7> );
	r[3] = q5(i).integrate(f<9> );
	r[4] = q6(i).integrate(f<11> );
	r[5] = q7(i).integrate(f<13> );
}

void TEST::check () {
	double e[6];

	e[0] = std::abs((r[0] - exact<3>()) / (exact<3>()));
	e[1] = std::abs((r[1] - exact<5>()) / (exact<5>()));
	e[2] = std::abs((r[2] - exact<7>()) / (exact<7>()));
	e[3] = std::abs((r[3] - exact<9>()) / (exact<9>()));
	e[4] = std::abs((r[4] - exact<11>()) / (exact<11>()));
	e[5] = std::abs((r[5] - exact<13>()) / (exact<13>()));

	double const eps(1.e-12);

	gas_assert(e[0]< eps);
	gas_assert(e[1]< eps);
	gas_assert(e[2]< eps);
	gas_assert(e[3]< eps);
	gas_assert(e[4]< eps);
	gas_assert(e[5]< eps);
}

gas_unit(TEST)
