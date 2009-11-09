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
#include <cmath>

#define TEST gas_functional_base_P1_tetra

class TEST {
public:
	TEST ();
	void execute ();
	void check ();

	typedef gas::functional::base::P1<gas::geometry::unit::tetra> base;

	gas::numerical::tiny::vector<4u> b;
	gas::numerical::tiny::vector<4u> dbdX;
	gas::numerical::tiny::vector<4u> dbdY;
	gas::numerical::tiny::vector<4u> dbdZ;
};

TEST::TEST () {
}

void TEST::execute () {

	rangeu(i, 4)
		b(i) = base::b(i, 1./4., 1./4., 1./4.);

	rangeu(i, 4)
		dbdX(i) = base::dbdX(i, 1./4., 1./4., 1./4.);

	rangeu(i, 4)
		dbdY(i) = base::dbdY(i, 1./4., 1./4., 1./4.);

	rangeu(i, 4)
		dbdZ(i) = base::dbdZ(i, 1./4., 1./4., 1./4.);
}

void TEST::check () {
	typedef gas::functional::base::info<base> info;
	using gas::geometry::unit::tetra;

	gas_assert(info::d == 3u);
	gas_assert(info::n == 4u);
	gas_assert((gas::same_type<info::unit_t, tetra>::value));

	double const eps(1.e-14);

	gas_assert(std::abs(b(0) - 1./4.) < eps);
	gas_assert(std::abs(b(1) - 1./4.) < eps);
	gas_assert(std::abs(b(2) - 1./4.) < eps);
	gas_assert(std::abs(b(3) - 1./4.) < eps);

	gas_assert(std::abs(dbdX(0) + 1.0) < eps);
	gas_assert(std::abs(dbdX(1) - 1.0) < eps);
	gas_assert(std::abs(dbdX(2) - 0.0) < eps);
	gas_assert(std::abs(dbdX(3) - 0.0) < eps);

	gas_assert(std::abs(dbdY(0) + 1.0) < eps);
	gas_assert(std::abs(dbdY(1) - 0.0) < eps);
	gas_assert(std::abs(dbdY(2) - 1.0) < eps);
	gas_assert(std::abs(dbdY(3) - 0.0) < eps);

	gas_assert(std::abs(dbdZ(0) + 1.0) < eps);
	gas_assert(std::abs(dbdZ(1) - 0.0) < eps);
	gas_assert(std::abs(dbdZ(2) - 0.0) < eps);
	gas_assert(std::abs(dbdZ(3) - 1.0) < eps);
}

gas_unit(TEST)
