/*
 * Copyright (c) 2008, Davide Ferrarese & Mattia Penati
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

#ifndef _GAS_LINEARALGEBRA_
#define _GAS_LINEARALGEBRA_

#include <cmath>

namespace LinearAlgebra {

template < typename Type > struct id { static inline Type RET ( Type const & x ) { return x; } }; 

template < typename Type > struct add { static inline Type RET ( Type const & x , Type const & y ) { return x + y; } };
template < typename Type > struct sub { static inline Type RET ( Type const & x , Type const & y ) { return x - y; } };
template < typename Type > struct mul { static inline Type RET ( Type const & x , Type const & y ) { return x * y; } };
template < typename Type > struct div { static inline Type RET ( Type const & x , Type const & y ) { return x / y; } };

template < typename Type > struct mul_mat_vet {};

}

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Product.hpp"

#include "Vector2.hpp"

namespace LinearAlgebra { namespace Solver {

#include "solver/solver.hpp"
#include "solver/upper.hpp"
#include "solver/lower.hpp"
#include "solver/lower1.hpp"
#include "solver/lu.hpp"
#include "solver/diagonal.hpp"
#include "solver/cholesky.hpp"

} }

#endif // _GAS_LINEARALGEBRA_