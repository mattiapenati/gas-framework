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

#ifndef _GAS_ELEMENT_P1_TRIANGLE_
#define _GAS_ELEMENT_P1_TRIANGLE_

#include "../LinearAlgebra/LinearAlgebra.h"
#include "../Geometry/Geometry.h"
#include "../Integration/Integration.h"

namespace Element {

template <>
class P1<Geometry::Triangle> {
	private:
		typedef P1<Geometry::Triangle> self;

	public:
		enum {
			Dimension = 2,
			DegreeOfFreedom = 3
		};

	public:
		static inline double PHI0 (double const & x, double const & y) { return 1-x-y; }
		static inline double PHI1 (double const & x, double const & y) { return x; }
		static inline double PHI2 (double const & x, double const & y) { return y; }

		static inline LinearAlgebra::Vector<double, 2> GRAD0 (double const & x, double const & y) { return LinearAlgebra::Vector<double, 2>(-1.,-1.); }
		static inline LinearAlgebra::Vector<double, 2> GRAD1 (double const & x, double const & y) { return LinearAlgebra::Vector<double, 2>(+1.,+0.); }
		static inline LinearAlgebra::Vector<double, 2> GRAD2 (double const & x, double const & y) { return LinearAlgebra::Vector<double, 2>(+0.,+1.); }

	private:
		struct base_function {
			unsigned int index_;
			Geometry::Triangle const *g_;
			inline base_function (unsigned int const & index, Geometry::Triangle const & g): index_(index), g_(&g) {}
			inline double operator() (double const & x, double const & y) const {
				switch (index_) {
					case 0: return 1 - g_->invx(x,y) - g_->invy(x,y);
					case 1: return g_->invx(x,y);
					case 2: return g_->invy(x,y);
				}
			} 
		};
		struct grad_function {
			unsigned int index_;
			Geometry::Triangle const *g_;
			inline grad_function (unsigned int const & index, Geometry::Triangle const & g): index_(index), g_(&g) {}
			inline LinearAlgebra::Vector<double, 2> operator() (double const & x, double const & y) const {
				switch (index_) {
					case 0: return LinearAlgebra::Vector<double, 2>(
						(g_->y(+1.,-1.) - g_->y(0.,0.)) / g_->det(0.,0.),
						(g_->x(-1.,+1.) - g_->x(0.,0.)) / g_->det(0.,0.)
					);
					case 1: return LinearAlgebra::Vector<double, 2>(
						(g_->y(0.,+1.) - g_->y(0.,0.)) / g_->det(0.,0.),
						(g_->x(0.,-1.) - g_->x(0.,0.)) / g_->det(0.,0.)
					);
					case 2: return LinearAlgebra::Vector<double, 2>(
						(g_->y(-1.,0.) - g_->y(0.,0.)) / g_->det(0.,0.),
						(g_->x(+1.,0.) - g_->x(0.,0.)) / g_->det(0.,0.)
					);
				}
			}
		};

	public:
		P1(): g_(), u_(1.) {
		}

		P1(Geometry::Triangle const & g): g_(g), u_(1.) {
		}

		P1(Geometry::Triangle const & g, LinearAlgebra::Vector<double, 3> const & u): g_(g), u_(u) {
		}

		P1(Geometry::Triangle const & g, double const & u0, double const & u1, double const & u2): g_(g), u_() {
			u_(0) = u0;
			u_(1) = u1;
			u_(2) = u2;
		}

		inline double & operator() (unsigned int const & i) { return u_(i); }
		inline double const & operator() (unsigned int const & i) const { return u_(i); }

		inline double operator() (double const & x, double const & y) {
			return u_(0) + (u_(1) - u_(0)) * g_.invx(x,y) + (u_(2) - u_(0)) * g_.invy(x,y);
		}
		
		inline base_function base(unsigned int const & i) {
			return base_function(i, g_);
		} 

		inline LinearAlgebra::Vector<double, 2> grad (double const & x, double const & y) {
			LinearAlgebra::Vector<double, 2> v0(
						(g_.y(+1.,-1.) - g_.y(0.,0.)) / g_.det(0.,0.),
						(g_.x(-1.,+1.) - g_.x(0.,0.)) / g_.det(0.,0.)
			);
			LinearAlgebra::Vector<double, 2> v1(
						(g_.y(0.,+1.) - g_.y(0.,0.)) / g_.det(0.,0.),
						(g_.x(0.,-1.) - g_.x(0.,0.)) / g_.det(0.,0.)
			);
			LinearAlgebra::Vector<double, 2> v2(
						(g_.y(-1.,0.) - g_.y(0.,0.)) / g_.det(0.,0.),
						(g_.x(+1.,0.) - g_.x(0.,0.)) / g_.det(0.,0.)
			);
			LinearAlgebra::Vector<double, 2> r = (v0 * u_(0)) + (v1 * u_(1)) + (v2 * u_(2));
			return r;
		}

		inline grad_function grad(unsigned int const & i) {
			return grad_function(i, g_);
		} 

		LinearAlgebra::Matrix<double, 3, 3> stiffness () {
			/* integrator */
			typedef Integrator<Method::NewtonCotes<Geometry::Triangle, 1> > IntegratorType;
			typedef IntegratorType::Transform Policy;
			IntegratorType integ(g_);
			/* matrice */
			LinearAlgebra::Matrix<double, 3, 3> M;
			/* calcolo */
			for (unsigned int i = 0; i < 3; ++i)
				for (unsigned int j = 0; j < 3; ++j)
					M(i,j) = integ.integrateMul<Policy, Policy>(base(i), base(j));
			/* return */
			return M;
		}

		double normL2 () {
			return std::sqrt(u_*(stiffness()*u_));
		}

		double normH1() {
			LinearAlgebra::Vector<double, 2> v = grad(0.,0.);
			double gu = g_.area() * (v*v);
			return std::sqrt(u_*(stiffness()*u_) + gu);
		}

		template <typename Function>
		double normL2 (Function const & f) {
			/* integrator */
			typedef Integrator<Method::NewtonCotes<Geometry::Triangle, 3> > IntegratorType;
			typedef IntegratorType::Transform Policy;
			IntegratorType integ(g_);
			/* double product */
			LinearAlgebra::Vector<double, 3> b;
			for (unsigned int i = 0; i < 3; ++i)
				b(i) = integ.integrateMul<Policy, Policy>(base(i), f);
			/* norm f */
			double fl2 = integ.integrateMul<Policy, Policy>(f, f);
			/* norm u */
			double ul2 = u_ * (stiffness()*u_);
			/* return */
			return std::sqrt(std::abs(fl2 + ul2 - 2. * (b*u_)));
		}

		template <typename Function1, typename Function2, typename Function3>
		double normH1 (Function1 const & f, Function2 const & fx, Function3 const & fy) {
			/* integrator */
			typedef Integrator<Method::NewtonCotes<Geometry::Triangle, 3> > IntegratorType;
			IntegratorType integ(g_);
			/* L2 norm */
			double nl2 = normL2(f); nl2 *= nl2;
			/* L2 norm of grad f */
			double fxl2 = integ.integrateMul<IntegratorType::Transform, IntegratorType::Transform>(fx,fx);
			double fyl2 = integ.integrateMul<IntegratorType::Transform, IntegratorType::Transform>(fy,fy);
			/* integral of grad f */
			LinearAlgebra::Vector<double, 2> b(
				integ.integrate<IntegratorType::Transform>(fx),
				integ.integrate<IntegratorType::Transform>(fy)
			);
			/* grad of u */
			LinearAlgebra::Vector<double, 2> v = grad(0.,0.);
			/* return */
			return std::sqrt(std::abs(nl2 + fxl2 + fyl2 + (g_.area()*(v*v)) - 2. * (b*v)));
		}
	
	private:
		Geometry::Triangle g_;
		LinearAlgebra::Vector<double, 3> u_;

};

}

#endif // _GAS_ELEMENT_P1_TRIANGLE_