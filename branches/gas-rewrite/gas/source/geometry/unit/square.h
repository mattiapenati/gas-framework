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
 * @file square.h
 * @brief The square obtained by cartesian product of two intervals,
 *        \f$(-1,1)\times(-1,1)\f$
 */

#ifndef _gas_geometry_unit_square_
#define _gas_geometry_unit_square_

#include "interval.h"

namespace gas { namespace geometry { namespace unit {

/*!
 * @brief The square obtained by cartesian product of two intervals,
 *        \f$(-1,1)\times(-1,1)\f$
 */
class square {

public:
	/*! @brief The self type */
	typedef square self_t;

public:
	/*! @brief Dimension of geometry */
	static unsigned int const d = 2u;

	/*!
	 * @brief Check the membership of a point by its coordinates
	 * @param x The first coordinate of point
	 * @param y The second coordinate of point
	 * @return True if the points is locate in the square
	 */
	static inline bool in (double const & x, double const & y) {
		return (interval::in(x) and interval::in(y));
	}

};

} } }

#endif // _gas_geometry_unit_square_
