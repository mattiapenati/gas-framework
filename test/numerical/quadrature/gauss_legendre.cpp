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
#include <gas>
#include <cmath>

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

	template <unsigned int nodes_>
	class method: public gas::numerical::quadrature::gauss_legendre<interval, nodes_> {};

	typedef gas::numerical::quadrature::formula< method<2u> > formula2;
	typedef gas::numerical::quadrature::formula< method<3u> > formula3;
	typedef gas::numerical::quadrature::formula< method<4u> > formula4;
	typedef gas::numerical::quadrature::formula< method<5u> > formula5;
	typedef gas::numerical::quadrature::formula< method<6u> > formula6;
	typedef gas::numerical::quadrature::formula< method<7u> > formula7;

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
	using namespace gas::numerical::quadrature;

	gas_assert((gas::same_type<info<formula2>::unit_t, interval>::value));
	gas_assert((gas::same_type<info<formula3>::unit_t, interval>::value));
	gas_assert((gas::same_type<info<formula4>::unit_t, interval>::value));
	gas_assert((gas::same_type<info<formula5>::unit_t, interval>::value));
	gas_assert((gas::same_type<info<formula6>::unit_t, interval>::value));
	gas_assert((gas::same_type<info<formula7>::unit_t, interval>::value));

	gas_assert(info<formula2>::d == 1u);
	gas_assert(info<formula3>::d == 1u);
	gas_assert(info<formula4>::d == 1u);
	gas_assert(info<formula5>::d == 1u);
	gas_assert(info<formula6>::d == 1u);
	gas_assert(info<formula7>::d == 1u);

	gas_assert((gas::same_type<info<formula2>::method_t, method<2u> >::value));
	gas_assert((gas::same_type<info<formula3>::method_t, method<3u> >::value));
	gas_assert((gas::same_type<info<formula4>::method_t, method<4u> >::value));
	gas_assert((gas::same_type<info<formula5>::method_t, method<5u> >::value));
	gas_assert((gas::same_type<info<formula6>::method_t, method<6u> >::value));
	gas_assert((gas::same_type<info<formula7>::method_t, method<7u> >::value));

	gas_assert(info<formula2>::n == 2u);
	gas_assert(info<formula3>::n == 3u);
	gas_assert(info<formula4>::n == 4u);
	gas_assert(info<formula5>::n == 5u);
	gas_assert(info<formula6>::n == 6u);
	gas_assert(info<formula7>::n == 7u);

	gas_assert(info<formula2>::degree == 3u);
	gas_assert(info<formula3>::degree == 5u);
	gas_assert(info<formula4>::degree == 7u);
	gas_assert(info<formula5>::degree == 9u);
	gas_assert(info<formula6>::degree == 11u);
	gas_assert(info<formula7>::degree == 13u);

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
