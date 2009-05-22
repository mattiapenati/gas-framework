/*                                                          
 * Copyright (c) 2008, Alfonso Fascì, Davide Ferrarese, Mattia Penati     
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

template < typename Geometry , unsigned int Degree > struct NewtonCotes;

/* integrazione sull'intervallo */
template < unsigned int Degree >
struct NewtonCotes<Geometry::Interval, Degree>: public method_1<NewtonCotes<Geometry::Interval, Degree>, Geometry::Interval> {
	/* nodi e pesi */
	static const double x[];
	static const double w[];
	/* numbero di punti */
	static const unsigned int nPoints;
};

/* Ordine 1 su intervallo */
template < > const unsigned int NewtonCotes < Geometry::Interval , 1 >::nPoints = 2;
template < > const double NewtonCotes < Geometry::Interval , 1 >::x[2] = {
	-1.0000000000000000, +1.0000000000000000
};
template < > const double NewtonCotes < Geometry::Interval , 1 >::w[2] = {
	1.0000000000000000, 1.0000000000000000
};
/* Ordine 2 su intervallo */
template < > const unsigned int NewtonCotes < Geometry::Interval , 2 >::nPoints = 3;
template < > const double NewtonCotes < Geometry::Interval , 2 >::x[3] = {
	-1.0000000000000000, +0.0000000000000000, +1.0000000000000000
};
template < > const double NewtonCotes < Geometry::Interval , 2 >::w[3] = {
	0.3333333333333333, 1.3333333333333333, 0.3333333333333333
};
/* Ordine 3 su intervallo */
template < > const unsigned int NewtonCotes < Geometry::Interval , 3 >::nPoints = 4;
template < > const double NewtonCotes < Geometry::Interval , 3 >::x[4] = {
	-1.0000000000000000, -0.3333333333333333, +0.3333333333333333, +1.0000000000000000
};
template < > const double NewtonCotes < Geometry::Interval , 3 >::w[4] = {
	0.2500000000000000, 0.7500000000000000, 0.7500000000000000, 0.2500000000000000
};

/* For more details on Newton-Cotes quadrature formulae on triangule see:
 * @article{0501496v2,
 *	eprint = {math/0501496v2},
 *	author = {Mark A. Taylor, Beth A. Wingate, Len P. Bos},
 *	title = {Several new quadrature formulas for polynomial integration in the triangle},
 *	year = {2007}
 * }
*/

/* Integrazione sul triangolo */
template < unsigned int Degree >
struct NewtonCotes<Geometry::Triangle, Degree>: public method_2<NewtonCotes<Geometry::Triangle, Degree>, Geometry::Triangle> {
	/* nodi e pesi */
	static const double x[];
	static const double y[];
	static const double w[];
	/* numbero di punti */
	static const unsigned int nPoints;
};

/* Ordine 1 su triangolo : grado di esattezza 2 */
template < > const unsigned int NewtonCotes < Geometry::Triangle , 1 >::nPoints = 3;
template < > const double NewtonCotes < Geometry::Triangle , 1 >::x[3] = {
	0.1666666666667, 0.6666666666667, 0.1666666666667
};
template < > const double NewtonCotes < Geometry::Triangle , 1 >::y[3] = {
	0.6666666666667, 0.1666666666667, 0.1666666666667
};
template < > const double NewtonCotes < Geometry::Triangle , 1 >::w[3] = {
	0.1666666666667, 0.1666666666667, 0.1666666666667
};

/* Ordine 2 su triangolo : grado di esattezza 4 */
template < > const unsigned int NewtonCotes < Geometry::Triangle , 2 >::nPoints = 6;
template < > const double NewtonCotes < Geometry::Triangle , 2 >::x[6] = {
	0.0915762135098, 0.8168475729805, 0.0915762135098, 0.1081030181681, 0.4459484909160, 0.4459484909160
};
template < > const double NewtonCotes < Geometry::Triangle , 2 >::y[6] = {
	0.0915762135098, 0.0915762135098, 0.8168475729805, 0.4459484909160, 0.1081030181681, 0.4459484909160
};
template < > const double NewtonCotes < Geometry::Triangle , 2 >::w[6] = {
	0.0549758718277, 0.0549758718277, 0.0549758718277, 0.1116907948390, 0.1116907948390, 0.1116907948390
};

/* Ordine 3 su triangolo : grado di esattezza 5 */
template < > const unsigned int NewtonCotes < Geometry::Triangle , 3 >::nPoints = 10;
template < > const double NewtonCotes < Geometry::Triangle , 3 >::x[10] = {
	0.0000000000000, 1.0000000000000, 0.0000000000000, 0.2673273531185, 0.6728175529461, 0.0649236350054,
	0.6716498539042, 0.0654032456800, 0.2693767069140, 0.3386738503896
};
template < > const double NewtonCotes < Geometry::Triangle , 3 >::y[10] = {
	1.0000000000000, 0.0000000000000, 0.0000000000000, 0.6728199218710, 0.2673288599482, 0.6716530111494, 
	0.0649251690029, 0.2693789366453, 0.0654054874919, 0.3386799893027
};
template < > const double NewtonCotes < Geometry::Triangle , 3 >::w[10] = {
	0.0065678024876, 0.0065679153017, 0.0068540986900, 0.0587095966456, 0.0587103059567, 0.0620062948279,
	0.0620076230630, 0.0629651151382, 0.0629665133414, 0.1126447345479
};
