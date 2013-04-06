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
 * @file adapt_mesh.h
 * @brief Includes the implementation of the class to manage the problem formulation
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_PROBLEM_H
#define _IMG_SEG_PROBLEM_H


namespace img_seg {
  
/*! @brief Function that computes the gradient
 *  @param cdt_ Triangulation where the function is computed
 *  @param z_vec_ Function to be derived
 *  @param grad Gradient of z_vec_ computed on the triangulation cdt_
 */  
void gradient(triangulation const & cdt_, const Eigen::VectorXd & z_vec_, Eigen::DynamicSparseMatrix<double> & grad)
{
	typedef triangulation::face_iterator_t iterator_f;
	
	for (iterator_f _it(cdt_.face_begin()); _it != cdt_.face_end(); ++_it)
	{
	  
	  double xA(_it->x(0));
	  double yA(_it->y(0));
	  double zA(z_vec_(_it->i(0)));
	  double xB(_it->x(1));
	  double yB(_it->y(1));
	  double zB(z_vec_(_it->i(1)));		
	  double xC(_it->x(2));
	  double yC(_it->y(2));
	  double zC(z_vec_(_it->i(2)));
		
	  Eigen::Vector3d row0(xA,yA,1);
	  Eigen::Vector3d row1(xB,yB,1);
	  Eigen::Vector3d row2(xC,yC,1);
	
	  Eigen::Matrix3d A;
	  Eigen::Vector3d coeff;
	  Eigen::Vector3d zz(zA,zB,zC);

	  A.row(0) = row0;
	  A.row(1) = row1;
	  A.row(2) = row2;
	
	  coeff = A.inverse()*zz;
		
	  grad.coeffRef(_it->i(),0)=coeff(0);			
	  grad.coeffRef(_it->i(),1)=coeff(1);					
	}
	
}


/*! @brief Function that computes the module of the gradient 
 *  @param cdt_ Triangulation where the function is computed
 *  @param z_vec_ Function to be derived
 *  @param grad_mod Module of the gradient of z_vec_ computed on the triangulation cdt_
 */
void gradient_module(triangulation const & cdt_, const Eigen::VectorXd & z_vec_, Eigen::VectorXd & grad_mod)
{
	typedef triangulation::face_iterator_t iterator_f;
	
	for (iterator_f _it(cdt_.face_begin()); _it != cdt_.face_end(); ++_it)
	{  
	    double xA(_it->x(0));
	    double yA(_it->y(0));
	    double zA(z_vec_(_it->i(0)));
	    double xB(_it->x(1));
	    double yB(_it->y(1));
	    double zB(z_vec_(_it->i(1)));		
	    double xC(_it->x(2));
	    double yC(_it->y(2));
	    double zC(z_vec_(_it->i(2)));
		
	    Eigen::Vector3d row0(xA,yA,1);
	    Eigen::Vector3d row1(xB,yB,1);
	    Eigen::Vector3d row2(xC,yC,1);
	    
	    Eigen::Matrix3d A;
	    Eigen::Vector3d coeff;
	    Eigen::Vector3d zz(zA,zB,zC);

	    A.row(0) = row0;
	    A.row(1) = row1;
	    A.row(2) = row2;
	    
	    coeff = A.inverse()*zz;
	
	    grad_mod(_it->i())=std::sqrt((coeff(0)*coeff(0) + coeff(1)*coeff(1)));		
	}

}


/*! @brief Class that describes the problem */
class problem
{
  private:
	typedef gas::geometry::unit::triangle triangle_t;
	typedef gas::numerical::quadrature::newton_cotes<triangle_t, 6u> method_t;
  
  public:
	  // Defining element and basis function
	  typedef gas::functional::base::P1<triangle_t> base_t;
	  typedef gas::functional::element<base_t> element_t;  
	  typedef gas::numerical::quadrature::formula<method_t> integrator_t;
	  
	  /*! @brief Constructor 
	  *  @param cdt Triangulation
	  *  @param phi Initial level-set curve
	  *  @param u0 Initial image (Image to be processed)
	  *  @param data Struct thta contains problem data
	  */
	  problem (triangulation const & cdt, initial_curve &phi, initial_image &u0, const struct pb_data & data);
	  
	  /*! @brief Destructor */
	  ~problem() {};
	  
	  /*! @brief Return the triangulation */ 
	  inline triangulation const & mesh() const { return cdt_; }
	  
	  /*! @brief Return the struct containing the problem data */
	  inline struct pb_data data() const { return data_; }

	  /*! @brief Return the number of faces in the triangulation */
	  inline int faces () const { return cdt_.faces(); }	
	  
	  /*! @brief Return the number of non-zeros elements */
	  inline int no_zeros () const { return no_zeros_; }

	  /*! @brief Return lambda */
	  inline double lambda() const { return data_.lambda; }
	  
	  /*! @brief Return epsilon */
	  inline double eps() const { return data_.epsilon; }

	  /*! @brief Return rho */
	  inline double rho() const { return data_.rho; }
	  
	  /*! @brief Return mu */
	  inline double mu() const { return data_.mu; }
	  
	  /*! @brief Return tau */
	  inline double tau() const { return data_.tau; }
	  
	  /*! @brief Return true if the problem has already been solved */
	  inline bool already_solved() const { return solved; }	
	  
	  /*! @brief Return true if the original image has been modified 
	  *  If false the printer class plots the initial image interpolated over the triangulation
	  */
	  inline bool not_original() const { return modified; }
	  
	  /*! @brief Set the modified bool to true
	  *  Called to allow printer class to plot the solution
	  */
	  inline void setNotOriginal() { modified=true; }		
	  
	  /*! @brief Return true if the current step is a refinement step*/
	  inline bool is_refinement() const { return refinement; }	
	  
	  /*! @brief Set the refinement value equal to the given boolean 
	  *  This script allows to plot the mesh during refinement steps 
	  *  and to plot the solution during solver steps 
	  */
	  inline void update_refinement(bool value) { refinement = value; }		
	  
	  /*! @brief Return the solution at current time 
	  *  @param vec Vector to store the current solution in
	  */
	  void get_solution(Eigen::VectorXd & vec) const { vec = x_; }
	  
	  /*! @brief Return the solution at current time at node i */
	  inline double operator() (int const & i) const { return x_(i); }	
	  
	  /*! @brief Return the solution at previous time step
	  *  @param vec Vector to store the previous solution in
	  */
	  void get_old_solution(Eigen::VectorXd & vec) const { vec = x_old_; }	

	  /*! @brief Return the solution at previous time step at node i */
	  inline double old_solution(int const & i) const { return x_old_(i); }
	  
	  /*! @brief Return the restriction of initial image over the triangulation
	  *  @param vec Vector to store the image interpolation in
	  */
	  void get_init_img_vec(Eigen::VectorXd & vec) const { vec = u_0_; }
	  
	  /*! @brief Return initial image at node i*/
	  inline double init_img_vec(int const & i) const { return u_0_(i); }		
	  
	  /*! @brief Return initial image*/
	  inline initial_image init_image() const { return img_in_; }		

	  /*! @brief Solver*/
	  void solve();
	  
	  /*! @brief Compute the evolutionary error to evaluate stopping criterion*/
	  double stopping_res();
	  
	  /*! @brief Update data after initial refinement */
	  void init_update();
	  
	  /*! @brief Update data after loop refinement or solver*/
	  void loop_update();

	  /*! @brief Update the solution with newly interpolated values
	  *  @param vec Vector of pairs containing the information about the coordinates (x,y) 
	  *	       and the solution at previous and current time step (x_old_,x_)
	  */
	  void updateInterpolator(std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & vec);	
	  
	  /*! @brief Resize data structure to fit new dimensions */
	  void InterpolationResize();
	  
	  /*! @brief Interpolate the solution in a new point (x,y)
	  *  @param it Iterator to the face we are working on
	  *  @param x x value for the new point
	  *  @param y y value for the new point
	  *  @param vec Vector that stores the coordinates of the new point (x,y) and the values of 
	  *             the solution at previous and current time step (x_old_,x_)
	  *  @param new_point True if the point is newly inserted and false if it's the vertex of an 
	  * 		     existing triangle
	  */
	  void interpolate(triangulation::face_iterator_t & it, double const & x, double const & y, std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & vec, bool new_point) const ;
	  
  private:
	  /*! @brief Mesh */
	  triangulation const & cdt_;	  
    
	  /*! @brief Solution at current time step */
	  Eigen::VectorXd x_;
	  
	  /*! @brief Solution at previous time step */
	  Eigen::VectorXd x_old_;

	  /*! @brief Initial curve */
	  initial_curve curve_;		
	  
	  /*! @brief Initial image */
	  initial_image img_in_;
	  
	  /*! @brief Restriction of the initial curve to the triangulation */
	  Eigen::VectorXd u_0_;
	  
	  /*! @brief Struct that contains the problem data */
	  const struct pb_data & data_;
	  
	  /*! @brief Elementi non nulli */
	  int no_zeros_;

	  /*! @brief  Boolean variable "Problem already solved"*/	
	  bool solved;
	  
	  /*! @brief Boolean variable "Problem is not at first time step" */	
	  bool modified;
	  
	  /*! @brief Boolean variable "Refinement step" */
	  bool refinement;
	
};


problem::problem (triangulation const & cdt, initial_curve &phi, initial_image &u0, const struct pb_data & data): 
cdt_(cdt), curve_(phi), img_in_(u0), data_(data), no_zeros_(0u), solved(false), modified(false), refinement(false)
{
    u_0_.setZero(cdt_.nodes());
    img_in_.Trestrict(cdt_,u_0_);  

    x_old_.setZero(cdt_.nodes());
    phi.vectorize(x_old_);
}


void problem::solve()
{
    typedef triangulation::face_iterator_t iterator_f;
    typedef triangulation::nodes_iterator_t iterator_n;	

    // Using derivatives from gas
    using gas::functional::dx;
    using gas::functional::dy;

    // Functions
    #define delta_rho(r,z) ((r*onesn).cwise()/(M_PI*(r*r + (z*z).cwise())))
    #define q_eps2(e,g) (e*e + (g*g).cwise())        
    // Generalized mass matrix
    #define generalized_mass(u,v) ((u*v)*coef_f1.fun())
    // Generalized stiffness matrix
    #define generalized_stiffness(u,v) ((dx(u)*dx(v) + dy(u)*dy(v))*coef_f2.fun())
    // Generalized source
    #define source(v) (coef_f4.fun()*v)
    // Generalized mass matrix for previous time step
    #define old_mass(v) (coef_f3.fun()*v)    

    // Create force term
    ForceFactory* myFactory = new ForceFactory(this->data().model, cdt_, x_old_, img_in_, u_0_, this->rho());
    Force* force_t = myFactory->create();     
   
    // Dimensions of the problem
    int const N(element_t::base_t::n);

    // Temporary vectors
    Eigen::VectorXd ones, onesn;

    ones.setZero(cdt_.faces());
    onesn.setZero(cdt_.nodes());
    for(int k=0; k<cdt_.faces();++k) { ones(k)=1.0; }
    for(int k=0; k<cdt_.nodes();++k) { onesn(k)=1.0; }
    
    // Coefficients for the equation
    Eigen::VectorXd gradient_mod, qe, coef_v1, coef_v2, coef_v3, coef_v4;
    
    gradient_mod.setZero(cdt_.faces());
    qe.setZero(cdt_.faces());
    coef_v1.setZero(cdt_.nodes());
    coef_v2.setZero(cdt_.faces());
    coef_v3.setZero(cdt_.nodes());
    coef_v4.setZero(cdt_.nodes());

    gradient_module(cdt_, x_old_, gradient_mod);	    
    radq((q_eps2(this->eps(),gradient_mod)),qe);
    coef_v1=(onesn.cwise()/(delta_rho(this->rho(),x_old_)));
    coef_v2=(ones.cwise()/qe);    
    coef_v3=x_old_.cwise()*coef_v1;
    coef_v4=force_t->vector();	    
    
    // Algebraic structures to solve the problem
    Eigen::DynamicSparseMatrix<double> Mtmp(cdt_.nodes(),cdt_.nodes());
    Eigen::DynamicSparseMatrix<double> Atmp(cdt_.nodes(),cdt_.nodes());
    Eigen::VectorXd b;

    b.setZero(cdt_.nodes());
    
    integrator_t s;    
    
    // Loop over the faces
    for (iterator_f it(cdt_.face_begin()); it != cdt_.face_end(); ++it)
    {
	s(*it);
	element_t e(*it);
	
	// Approximating the coefficients with linear functions to evaluate them in Gauss points
	integral_argument coef_f1(coef_v1, it);
	int_arg_mod coef_f2(coef_v2, it);		
	integral_argument coef_f3(coef_v3, it);
	integral_argument coef_f4(coef_v4, it);
	
	// Generalized mass matrix
	gas_rangeu(i, N)
	{
	    gas_rangeu(j, N)
	    {
		int const ii(it->i(i));
		int const jj(it->i(j));
		
		double const r(s.integrate(generalized_mass(e.b(j),e.b(i))));

		Mtmp.coeffRef(ii, jj) += r;
	    }
	}

	// Generalized stiffness matrix
	gas_rangeu(i, N)
	{
	    gas_rangeu(j, N)
	    {
		int const ii(it->i(i));
		int const jj(it->i(j));

		double const r(this->mu()*this->tau()*(s.integrate(generalized_stiffness(e.b(j),e.b(i)))));
		
		Atmp.coeffRef(ii, jj) += r;
	    }
	}
	
	// Force term
	gas_rangeu(i, N)
	{
	    int const ii(it->i(i));

	    double const r1(s.integrate(old_mass(e.b(i))));
	    double r2(0.0);
	    
	    if (this->data().model == 0)
	    {
		r2 = this->lambda()*this->tau()*(s.integrate(source(e.b(i)))); // ChanVeseForce
	    }
	    else
	    {
		r2 = this->tau()*(s.integrate(source(e.b(i)))); // RegionGaussianForce and RegionNonParametricForce
	    }
	    
	    b.coeffRef(ii) += r1;
	    b.coeffRef(ii) += r2;
	}
    }

    // Final matrix
    Eigen::SparseMatrix<double> K(Mtmp+Atmp);
    
    // Compute non-zeros elements
    no_zeros_ = K.nonZeros();
    x_.resize(cdt_.nodes());

    // Solve the algebraic problem using LU factorization for sparse matrix
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::UmfPack> lu(K);
    gas_assert(lu.succeeded());
    lu.solve(b, &x_);

    // Flag problem as solved
    solved=true;
    
    // Export solution to data file
    std::ofstream lset_file("result/level_set.dat", std::ios::app);
    lset_file << x_.size() << std::endl;
    for (iterator_n i(cdt_.nodes_begin()); i != cdt_.nodes_end(); ++i) 
    {
      lset_file << i->x() << "\t" << i->y() << "\t" << x_(i->i()) << std::endl;
    }
    lset_file << std::endl;
    lset_file << std::endl;      
    lset_file.close();
    
    // Clean up memory
    delete force_t;
    force_t = NULL;
    delete myFactory;
    myFactory = NULL;
}


double problem::stopping_res()
{ 
      typedef triangulation::face_iterator_t iterator_f;	  
      typedef triangulation::nodes_iterator_t iterator_n;  
    
      typedef gas::geometry::unit::triangle triangle_t;
      typedef gas::numerical::quadrature::newton_cotes<triangle_t, 6u> method_t;	
      typedef gas::numerical::quadrature::formula<method_t> integrator_t;
    
      // Create force terms associated to the solution at current and previous time step
      ForceFactory* myFactory0 = new ForceFactory(this->data().model, cdt_, x_old_, img_in_, u_0_, this->rho());
      ForceFactory* myFactory1 = new ForceFactory(this->data().model, cdt_, x_, img_in_, u_0_, this->rho());
      Force* f0 = myFactory0->create();
      Force* f1 = myFactory1->create();	  
      
      
      double err(0.0);
      
      if(this->data().model == 0)  // ChanVeseForce
      {
	  Eigen::VectorXd s0;
	  s0.setZero(cdt_.nodes());
	  Eigen::VectorXd s1;
	  s1.setZero(cdt_.nodes());	  	  
	
	  // Compute probability density associated the force terms
	  Eigen::VectorXd P1r1 = f1->p1();
	  Eigen::VectorXd P1r2 = f1->p2();		
	  Eigen::VectorXd P0r1 = f0->p1();		
	  Eigen::VectorXd P0r2 = f0->p2();
	  
	  // Loop over the nodes to set s_h depending on the position of the node with respect to the 
	  // level-set curve (internal or external)
	  for(iterator_n _j(cdt_.nodes_begin()); _j != cdt_.nodes_end(); ++_j)
	  {
	      if(x_(_j->i())>0)
	      {
		  s1(_j->i()) = P1r1(_j->i());
	      }
	      else
	      {
		  s1(_j->i()) = P1r2(_j->i());
	      }
	      
	      if(x_old_(_j->i())>0)
	      {
		  s0(_j->i()) = P0r1(_j->i());
	      }
	      else
	      {
		  s0(_j->i()) = P0r2(_j->i());
	      }    
	  }
	
	  integrator_t s;

	  // Loop over the faces to compute the error
	  for (iterator_f it(cdt_.face_begin()); it != cdt_.face_end(); ++it)
	  {
	    
	    integral_argument S0h(s0, it);
	    integral_argument S1h(s1, it);	
	    
	    s(*it);
	    double const e(s.integrate((S1h.fun()-S0h.fun())*(S1h.fun()-S0h.fun())));
	    err += e;
	  }

	  err = std::sqrt(err);	  
      }
      else  // RegionGaussianForce and RegionNonParametricForce
      {    
	  Eigen::VectorXd ones;
	  ones.setZero(cdt_.faces());
	  for(int k=0; k<cdt_.faces();++k) { ones(k)=1.0; }
	
	  Eigen::VectorXd p1_vec0, p2_vec0;
	  p1_vec0.setZero(cdt_.nodes());
	  p2_vec0.setZero(cdt_.nodes());	  
	  
	  Eigen::VectorXd p1_vec1, p2_vec1;
	  p1_vec1.setZero(cdt_.nodes());
	  p2_vec1.setZero(cdt_.nodes());
	  
	  // Compute probability density associated the force terms
	  p1_vec0 = f0->p1();
	  p2_vec0 = f0->p2();
	  
	  p1_vec1 = f1->p1();
	  p2_vec1 = f1->p2();		      
	  
	  integrator_t s, s0, s1;
	  
	  double A0(0.0);
	  double A1(0.0);
	  
	  Eigen::VectorXd _M0, _m0;
	  _M0.setZero(cdt_.nodes());
	  _m0.setZero(cdt_.nodes());
	  Eigen::VectorXd _M1, _m1;
	  _M1.setZero(cdt_.nodes());
	  _m1.setZero(cdt_.nodes());

	  Eigen::VectorXd diff0, diff1;
	  diff0.setZero(cdt_.nodes());
	  diff1.setZero(cdt_.nodes());		
		    
	  // Loop over the nodes of the triangulation
	  for(iterator_n _j(cdt_.nodes_begin()); _j != cdt_.nodes_end(); ++_j)
	  {   
	      // Compute the functions representing the maximum and the minimum between p1 and p0
	      // at current and previous time step	    
	      _M1(_j->i()) = std::max(p1_vec1(_j->i()),p2_vec1(_j->i()));
	      _m1(_j->i()) = std::min(p1_vec1(_j->i()),p2_vec1(_j->i()));
	      
	      _M0(_j->i()) = std::max(p1_vec0(_j->i()),p2_vec0(_j->i()));
	      _m0(_j->i()) = std::min(p1_vec0(_j->i()),p2_vec0(_j->i()));		    		    

	      // Compute the difference between maximum and minimum
	      diff0(_j->i()) = _M0(_j->i()) - _m0(_j->i());
	      diff1(_j->i()) = _M1(_j->i()) - _m1(_j->i());		    
	  }
	    
	    
	  // Compute the areas described by max and min curves
	  for(iterator_f it(cdt_.face_begin()); it != cdt_.face_end(); ++it)
	  {
	      s0(*it);
	      s1(*it);	    
	    
	      integral_argument M0h(diff0, it);
	      integral_argument M1h(diff1, it);	
	      
	      double const a0(s0.integrate(M0h.fun()));
	      double const a1(s1.integrate(M1h.fun()));
	      
	      A0 += a0;
	      A1 += a1;
	  }

	  // Compute the errors
	  for (iterator_f it(cdt_.face_begin()); it != cdt_.face_end(); ++it)
	  {
	    s(*it);
	    
	    integral_argument Ones(ones, it);
	    
	    double const e((A0-A1)*(A0-A1)*s.integrate(Ones.fun()));
	    
	    err += e;
	  }
	  
	  err = std::sqrt(err);
    }
    
    // Clean up memory
    delete f0;
    f0 = NULL;
    delete f1;
    f1 = NULL;
    delete myFactory0;
    myFactory0 = NULL;
    delete myFactory1;
    myFactory1 = NULL;

    // Return error
    return err;
}


void problem::init_update()
{
    // Get the dimensions of the original image
    double x_scale = curve_.w();
    double y_scale = curve_.h();
    
    // Update initial curve over the new mesh
    img_seg::initial_curve phi_new(cdt_, x_scale, y_scale);
    
    // Update the restriction of initial image to the new triangulation
    u_0_.resize(cdt_.nodes());
    u_0_.setZero(cdt_.nodes());	
    img_in_.Trestrict(cdt_, u_0_); 		
    
    // Update the initial solution to the new triangulation
    x_old_.resize(cdt_.nodes()); 
    x_old_.setZero(cdt_.nodes());	
    phi_new.vectorize(x_old_);
}	


void problem::loop_update()
{
    // Update the restriction to the new mesh
    u_0_.resize(cdt_.nodes()); 
    u_0_.setZero(cdt_.nodes());				
    img_in_.Trestrict(cdt_, u_0_); 
    
    // Update the solution at previous time step
    x_old_.resize(cdt_.nodes()); 
    x_old_.setZero(cdt_.nodes());				
    x_old_=x_;
}	


void problem::updateInterpolator(std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & vec) 
{
    typedef triangulation::nodes_iterator_t iterator_n;
    
    this->InterpolationResize();

    // Loop over the nodes to update the solutions at previous and current time step 
    // For new nodes inserted by the refinement strategy, the value for the solutions 
    // are interpolated 
    for(iterator_n _j(cdt_.nodes_begin()); _j != cdt_.nodes_end(); ++_j)
    {
	for(unsigned int _k(0); _k<vec.size(); ++_k)
	{
	  if( (_j->x() == vec[_k].first(0)) && (_j->y() == vec[_k].first(1)) )
	  {
	      x_old_(_j->i())=vec[_k].second(0);
	      x_(_j->i())=vec[_k].second(1);
	  }
	}
    }
}


void problem::InterpolationResize() 
{
    // Resize the vectors
    x_.resize(cdt_.nodes());	  
    x_old_.resize(cdt_.nodes());
    u_0_.resize(cdt_.nodes());
    
    // Restrict initial image to the mesh
    u_0_.setZero(cdt_.nodes());
    img_in_.Trestrict(cdt_, u_0_); 	      	      
}		


void problem::interpolate(triangulation::face_iterator_t & it, double const & x, double const & y, std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & vec, bool new_point) const 
{
    double r1, r2;
    
    // Store the values for the solution at previous and current time step in a vector 
    // along with the coordinates of the point
    if(new_point) // If the node didn't exist in the previous mesh
    {
	// Interpolate the solutions in the new node
	Eigen::VectorXd vec1(x_old_);
	integral_argument int1(vec1,it);
	r1=int1(x,y);
	Eigen::VectorXd vec2(x_);	      
	integral_argument int2(vec2,it);
	r2=int2(x,y);	      
    }
    else // If node already existing in previous triangulation
    {
	// Retrieve the value from older triangulation
	int j(0);
	
	gas_rangeu(k,3)
	{
	    if((it->x(k)==x) && (it->y(k)==y)) { j=it->i(k); }
	}
	
	r1=x_old_(j);
	r2=x_(j);	      
    }
    
    std::pair<Eigen::Vector2d,Eigen::Vector2d> n;
    n.first(0)=x;
    n.first(1)=y;
    n.second(0)=r1;
    n.second(1)=r2;
    
    vec.push_back(n);
}	
		

} //namespace


#endif 