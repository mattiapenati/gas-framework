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

template < typename Geometry , unsigned int Degree > struct Gauss;

/* integrazione su un intervallo */
template < unsigned int Degree >
struct Gauss<Geometry::Interval, Degree>: public method_1<Gauss<Geometry::Interval, Degree>, Geometry::Interval> {
	/* nodi e pesi */
	static const double x[];
	static const double w[];
	/* numbero di punti */
	static const unsigned int nPoints;
};

/* ordine 1 su intervallo */
template < > const unsigned int Gauss<Geometry::Interval, 1>::nPoints = 2;
template < > const double Gauss<Geometry::Interval, 1>::x[2] = {
	-0.5773502691896258, +0.5773502691896258
};
template < > const double Gauss<Geometry::Interval, 1>::w[2] = {
	1.0000000000000000, 1.0000000000000000
};

/* ordine 2 su intervallo */
template < > const unsigned int Gauss<Geometry::Interval, 2>::nPoints = 3;
template < > const double Gauss<Geometry::Interval, 2>::x[3] = {
	-0.7745966692414834, 0.0000000000000000, +0.7745966692414834
};
template < > const double Gauss<Geometry::Interval, 2>::w[3] = {
	0.5555555555555556, 0.8888888888888889, 0.5555555555555556
};

/* ordine 3 su intervallo */
template < > const unsigned int Gauss<Geometry::Interval, 3>::nPoints = 4;
template < > const double Gauss<Geometry::Interval, 3>::x[4] = {
	-0.8611363115940526, -0.3399810435848563, +0.3399810435848563, +0.8611363115940526
};
template < > const double Gauss<Geometry::Interval, 3>::w[4] = {
	0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539
};