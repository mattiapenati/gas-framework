/*
 * Copyright (c) 2013, Politecnico di Milano
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
 * @file posterior.h
 * @brief Includes the implementation of the class to compute error estimate and mark elements for refinement
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_POSTERIOR_H
#define _IMG_SEG_POSTERIOR_H


namespace img_seg {
  
  
/*!@brief Class to compute error estimation and to mark triangles for refinement */
class posterior {
  
  public:
	  /*! @brief Constructor
	  *   @param p The problem
	  *   @param data_r Struct to store refinement data
	  */
	  inline posterior (problem & p, const struct ref_data & data_r): p_(p), data_r_(data_r) {};

	  /*! @brief Destructor */
	  ~posterior() {};
	  
	  /*! @brief Compute the initial error on a given face
	  *  @param _it Iterator to the face we are treating
	  *  eta_0(S) = || u_0 -Iu_0 || in L^2(Omega) norm
	  */
	  double initial_error (triangulation::face_iterator_t & _it);
	  
	  /*! @brief Compute the residual error on a given face
	  *  @param _it Iterator to the face we are treating
	  *  eta_S(Phi^m)^2 = A + B + C
	  *  A = C_0^2 h_S^4 \Bigg\| \frac{\Phi_h^m(\mathbf{x}) - I_h^m\Phi_h^{m-1}(\mathbf{x})}{\tau \delta_{\rho}(\Phi_h^m(\mathbf{x}))} - G(u_0(\mathbf{x})) \Bigg\|_{L^2(S)}^2
	  *  B = C_1^2 h_S^3 \sum_{e \subset \partial S \cap \Omega}{\Bigg\| \Bigg\llbracket \frac{\nabla \Phi_h^m(\mathbf{x})}{Q_{\epsilon}(\Phi_h^m(\mathbf{x}))} \Bigg\rrbracket \Bigg\|_{L^2(e)}^2}
	  *  C = C_1^2 h_S^3 \sum_{e \subset \partial S \cap \partial \Omega}{\Bigg\| \frac{\nabla \Phi_h^m(\mathbf{x})}{Q_{\epsilon}(\Phi_h^m(\mathbf{x}))} \cdot \mathbf{n} \Bigg\|_{L^2(e)}^2}
	  */
	  double loop_error (triangulation::face_iterator_t & _it);
	  
	  /*! @brief GERS marking algorithm
	  *  @param loc_err Vector that contains previously computed local error
	  *  @param glob_err Global error on the mesh
	  *  @param marked Vector of booleans to store if an element has to be refined
	  */
	  void mark (Eigen::VectorXd & loc_err, double const & glob_err, std::vector<bool> & marked);

	  /*! @brief Return reference to the problem */
	  inline problem & pb_return() const { return p_; }
	  
	  /*! @brief Return the struct that contains refinement data */
	  inline struct ref_data data() const { return data_r_; }	
	  
	  /*! @brief Return tolerance for initial adaptivity process */
	  inline double init_tol() const { return data_r_.initial_tolerance; }
	  
	  /*! @brief Return tolerance for loop adaptivity process */
	  inline double loop_tol() const { return data_r_.loop_tolerance; }	
	  
	  /*! @brief Return GERS parameter nu */
	  inline double nu() const { return data_r_.nu; }
	  
	  /*! @brief Return GERS parameter theta_star */
	  inline double theta_star() const { return data_r_.theta_star; }		

  private:
	  /*! @brief Reference to the problem */
	  problem & p_;

	  /*! @brief Reference to the struct that contains refinement data */
	  const struct ref_data & data_r_;
	
};


double posterior::initial_error (triangulation::face_iterator_t & _it)
{
    // Get initial image information
    initial_image u0(p_.init_image());
    
    // Get the restriction of the initial image to the mesh
    Eigen::VectorXd u0_vec;
    u0_vec.setZero(p_.mesh().nodes());
    p_.get_init_img_vec(u0_vec);
    integral_argument u0h(u0_vec, _it);    

    problem::integrator_t s;
    s(*_it);

    // Compute initial error on the given element
    double quad_err(s.integrate((u0.fun()-u0h.fun())*(u0.fun()-u0h.fun())));
    
    /*
    // Export element where quad_err is nan
    if(quad_err != quad_err)
    {
      std::ofstream file("result/__ERROR__.txt",std::ios::app);
      file << _it->i() << "\t" << quad_err << std::endl;
      file << "\t" << _it->i(0) << "\t" << _it->x(0) << "\t" << _it->y(0) << std::endl;
      file << "\t" << _it->i(1) << "\t" << _it->x(1) << "\t" << _it->y(1) << std::endl;
      file << "\t" << _it->i(2) << "\t" << _it->x(2) << "\t" << _it->y(2) << std::endl;	  
      file.close();	
    }
    */

    return std::sqrt(quad_err);
}


double posterior::loop_error (triangulation::face_iterator_t & _it)
{  
    typedef triangulation::face::cgal_face_iterator_t cgal_face_iterator_t;
    typedef triangulation::face_iterator_t iterator_f;
    
    triangulation::face face_(*_it);
    cgal_face_iterator_t const cgal_it(face_.it_);
    
    // Weights for error estimation
    double C0(0.01);
    double C1(0.495);

    // Nodes coordinates
    Eigen::Vector2d const P0(_it->x(0), _it->y(0));
    Eigen::Vector2d const P1(_it->x(1), _it->y(1));
    Eigen::Vector2d const P2(_it->x(2), _it->y(2));

    // Edges information
    double const E0d((P1-P2).norm());
    double const E1d((P2-P0).norm());
    double const E2d((P0-P1).norm());
    Eigen::Vector3d edge_dim(E0d,E1d,E2d);
    double const hS(std::max(E0d, std::max(E1d, E2d)));

    // Get solution at current time
    Eigen::VectorXd sol;
    sol.setZero(p_.mesh().nodes());
    p_.get_solution(sol);

    // Get solution at previous time
    Eigen::VectorXd old_sol;
    old_sol.setZero(p_.mesh().nodes());
    p_.get_old_solution(old_sol);		
    
    // Get the restriction of the initial image to the mesh
    Eigen::VectorXd u0_vec;
    u0_vec.setZero(p_.mesh().nodes());
    p_.get_init_img_vec(u0_vec);
    
    // Get the initial image
    initial_image img = p_.init_image();

    // Get tho
    double rho(p_.rho());
	    
    // Create force term
    ForceFactory* myFac = new ForceFactory(p_.data().model, p_.mesh(), sol, img, u0_vec, rho);
    Force* ff = myFac->create();
    
    // Compute coefficients from previous time step and build linear approximating functions 
    Eigen::VectorXd onesn;
    onesn.setZero(p_.mesh().nodes());
    for(int k=0; k<p_.mesh().nodes();++k) { onesn(k)=1.0; }
    
    // Gradient module, delta_rho and gradient
    Eigen::VectorXd grad_mod, q_e, d_r;
    Eigen::DynamicSparseMatrix<double> grad(p_.mesh().faces(),2);
    
    grad_mod.setZero(p_.mesh().faces());
    q_e.setZero(p_.mesh().faces());	
    d_r.setZero(p_.mesh().nodes());
    
    gradient_module(p_.mesh(), sol, grad_mod);	
    radq((q_eps2(p_.eps(),grad_mod)),q_e);    
    d_r=delta_rho(p_.rho(),old_sol);	
    gradient(p_.mesh(), sol, grad);		
    
    // Coefficients for the equation 
    Eigen::VectorXd v1, v2, v3;
    Eigen::DynamicSparseMatrix<double> v4(p_.mesh().faces(),2);

    v1.setZero(p_.mesh().nodes());
    v2.setZero(p_.mesh().nodes());
    v3.setZero(p_.mesh().nodes());
    
    v1=(sol.cwise()/(p_.tau()*delta_rho(p_.rho(),sol)));
    v2=(old_sol.cwise()/(p_.tau()*delta_rho(p_.rho(),sol)));	

    if (p_.data().model == 0)
    {
      v3=p_.lambda()*ff->vector(); // ChanVeseForce	  
    }
    else 
    {
      v3=ff->vector(); // RegionGaussianForce and RegionNonParametricForce
    }

    for (iterator_f k(p_.mesh().face_begin()); k != p_.mesh().face_end(); ++k) 
    {
	v4.coeffRef(k->i(),0)=(grad.coeffRef(k->i(),0)/q_e(k->i()));
	v4.coeffRef(k->i(),1)=(grad.coeffRef(k->i(),1)/q_e(k->i()));	   
    }
    
    integral_argument f1(v1, _it);
    integral_argument f2(v2, _it);
    integral_argument f3(v3, _it);	
    
    problem::integrator_t s;
    s(*_it);

    // Compute residue on the face
    double face_guesstimator = std::pow(C0*hS*hS,2)*(s.integrate((f1.fun()-f2.fun()-f3.fun())*(f1.fun()-f2.fun()-f3.fun())));

    // Compute residue on the edges
    double edge_guesstimator(0.0);

    // Loop over the edges
    gas_rangeu(i, 3) 
    {
	// If the edge doesn't belong to the boundary of Omega
	if (!p_.mesh().cdt_.is_infinite(cgal_it->neighbor(i))) 
	{

	  // Compute the jump of the gradient along the edge
	  Eigen::Vector2d const g(v4.coeffRef(_it->i(),0)-v4.coeffRef(_it->neighbor(i),0),v4.coeffRef(_it->i(),1)-v4.coeffRef(_it->neighbor(i),1));
	  Eigen::Vector2d n;
	  
	  // Extern normal vector
	  switch(i)
	  {
	      case 0:
	      {
		  double theta = std::atan((P2(0)-P1(0))/(P1(1)-P2(1)));
		  n(std::cos(theta),std::sin(theta));				
	      }
	      break;
	      case 1:
	      {
		  double theta = std::atan((P2(0)-P0(0))/(P0(1)-P2(1)));
		  n(std::cos(theta),std::sin(theta));							  
	      }
	      break;
	      case 2:
	      {
		  double theta = std::atan((P0(0)-P1(0))/(P1(1)-P0(1)));
		  n(std::cos(theta),std::sin(theta));						  
	      }
	      break;
	  }
	  
	  double const t(n.dot(g));		    
	  edge_guesstimator += std::pow(hS,3)*std::pow(C1*t,2)*edge_dim(i)/2.;		    

	}
	else // If the edge is a boundary edge
	{
	  // Compute gradient
	  Eigen::Vector2d const g(v4.coeffRef(_it->i(),0),v4.coeffRef(_it->i(),1));
	  Eigen::Vector2d n;
	  
	  // Extern normal vector		  
	  switch(i)
	  {
	      case 0:
	      {
		  double theta = std::atan((P2(0)-P1(0))/(P1(1)-P2(1)));
		  n(std::cos(theta),std::sin(theta));			
	      }
	      break;
	      case 1:
	      {
		  double theta = std::atan((P2(0)-P0(0))/(P0(1)-P2(1)));
		  n(std::cos(theta),std::sin(theta));							  			  
	      }
	      break;
	      case 2:
	      {
		  double theta = std::atan((P0(0)-P1(0))/(P1(1)-P0(1)));
		  n(std::cos(theta),std::sin(theta));						  			  
	      }
	      break;
	  }
	  
	  double const t(n.dot(g));
	  edge_guesstimator += std::pow(hS,3)*std::pow(C1*t,2)*edge_dim(i);		    		    
	} 
    }
    
    
    // Total residue
    double guesstimator(face_guesstimator+edge_guesstimator);	
    
    // Clean up memory
    delete ff;
    ff = NULL;
    delete myFac;
    myFac = NULL;
    
    return std::sqrt(guesstimator);
}


void posterior::mark (Eigen::VectorXd & loc_err, double const & glob_err, std::vector<bool> & marked)
{
      double sum(0.0);
      double t(1.0);	
      double max_err(loc_err.maxCoeff());

      // Loop over the elements to identify optimal subset of triangles to be refined
      do
      {
	t -= this->nu();

	// Loop over all the elements
	for (int i(0); i != loc_err.size(); ++i)
	{
	    // If the element isn't already marked
	    if(!marked[i])
	    {
		// Checking condition to choose elements to be refined
		// See report for details about Guaranteed Error Reduction Strategy
		if(loc_err(i)>std::abs(t)*max_err)
		{
		  marked[i]=true;
		  sum += std::pow(loc_err(i),2);
		}
	    }
	}
	
      } while ( (sum < std::pow((1-this->theta_star())*glob_err,2)) && (t > 0) );	
}


} //namespace


#endif 