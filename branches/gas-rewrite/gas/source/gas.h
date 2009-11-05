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

/*!
 * @file gas.h
 * @brief The main header, includes all other files
 */

#ifndef _gas_
#define _gas_

/*!
 * @namespace gas
 * @brief The main namespace
 *
 * @namespace gas::functional
 * @brief Classes and functions to manage functional elements
 *
 * @namespace gas::geometry
 * @brief Classes and functions to manage geometric elements
 *
 * @namespace gas::geometry::map
 * @brief The maps to change the coordinates
 *
 * @namespace gas::geometry::unit
 * @brief The basic shapes on which you can define base function and quadrature
 *        formulae
 *
 * @namespace gas::numerical
 * @brief Classes and function for numerical methods
 *
 * @namespace gas::numerical::quadrature
 * @brief The quadrature formulae to integrate the function
 *
 * @namespace gas::numerical::tiny
 * @brief Linear algebra structure with fixed size at compile time
 */

#include "functional/derivative.h"

#include "gas/assertion.h"
#include "gas/chrono.h"
#include "gas/macro.h"
#include "gas/static.h"
#include "gas/test.h"
#include "gas/type.h"

#include "geometry/map/affine.h"
#include "geometry/unit/interval.h"
#include "geometry/unit/square.h"
#include "geometry/unit/triangle.h"

#include "numerical/quadrature/formula.h"
#include "numerical/quadrature/gauss_legendre.h"
#include "numerical/quadrature/method.h"
#include "numerical/quadrature/newton_cotes.h"
#include "numerical/tiny/det.h"
#include "numerical/tiny/dot.h"
#include "numerical/tiny/matrix.h"
#include "numerical/tiny/mul.h"
#include "numerical/tiny/utility.h"
#include "numerical/tiny/vector.h"

#endif // _gas_
