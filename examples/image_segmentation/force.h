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
 * @file force.h
 * @brief Includes the implementation of the class to manage force term
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_FORCE_H
#define _IMG_SEG_FORCE_H


namespace img_seg {
  

/*! @brief Compute the arctg of a given Eigen::VectorXd
  *  @param a VectorXd of data
  *  @param sol VectorXd to store arctg(a)
  */
void arctg(const Eigen::VectorXd & a, Eigen::VectorXd & sol)
{	
    for(int i_(0); i_<a.size(); ++i_)
    {
	sol(i_) = std::atan(a(i_));
    }	
}


/*! @brief Compute the logarithm of a given Eigen::VectorXd
   *  @param a VectorXd of data
   *  @param sol VectorXd to store log(a)
   */
void ln(const Eigen::VectorXd & a, Eigen::VectorXd & sol)
{
    for(int i_(0); i_<a.size(); ++i_)
    {
	sol(i_)=std::log(a(i_));
    }  
}

/*! @brief Compute the exponential of a given Eigen::VectorXd
   *  @param a VectorXd of data
   *  @param sol VectorXd to store esp(a)
   */
void esp(const Eigen::VectorXd & a, Eigen::VectorXd & sol)
{
    for(int i_(0); i_<a.size(); ++i_)
    {
	sol(i_)=std::exp(a(i_));
    }
}


/*! @brief Compute the square root of a given Eigen::VectorXd
   *  @param a VectorXd of data
   *  @param sol VectorXd to store radq(a)
   */
void radq(const Eigen::VectorXd & a, Eigen::VectorXd & sol)
{
    for(int i_(0); i_<a.size(); ++i_)
    {
	sol(i_)=std::sqrt(a(i_));
    }
}


/*! @brief Compute the power of a given Eigen::VectorXd
   *  @param a VectorXd of data
   *  @param n Exponent
   *  @param sol VectorXd to store pow(a,n)
   */
void pow(const Eigen::VectorXd & a, int n, Eigen::VectorXd & sol)
{	
    for(int i_(0); i_<a.size(); ++i_)
    {
	sol(i_)=std::pow(a(i_),n);
    }	
}  
  
  
/*! @brief Class that builds a linear approximating function that depends on (x,y) starting from a vector of discrete values */
class integral_argument
{

  public:
	  /*! @brief Constructor 
	   *  The constructor builds a linear approximation of a function starting from its nodal values
	   *  @param vec Vector that contains values in discrete points
	   *  @param it Face we are analyzing
	   */
	  integral_argument(Eigen::VectorXd &vec, triangulation::face_iterator_t &it);

	  /*! @brief Operator at, it evaluates the function in the point (x,y) 
	   *  @param x x value of the point
	   *  @param y y value of the point
	   */
	  inline double operator() (double const & x, double const & y) const
	  {
	      return a_coef*x+b_coef*y+d_coef;
	  }
    
  private:
	  // Linear approximation: z = a*x + b*y +d
	  /*! @brief Coefficient for x*/
	  double a_coef;
	  
	  /*! @brief Coefficient for y*/
	  double b_coef;
	  
	  /*! @brief Coefficient of order zero*/
	  double d_coef; 	
  
  public:
	  /*! @brief Function class*/
	  class integral_function : public gas::functional::function<2, integral_function>
	  {
	    public:
		    /*! @brief Constructor */
		    inline integral_function(integral_argument const & F): F_(F) {};
		    
		    /*! @brief Operator at 
		     *  @param x x value of the point
		     *  @param y y value of the point
		     */
		    inline double operator() (double const & x, double const & y) const
		    {
			return F_.a_coef*x+F_.b_coef*y+F_.d_coef;
		    }
		    
	    private:
		    /*! @brief Integral argument object */
		    integral_argument const & F_;	    
	  };    
    
	  /*! @brief Return a function depending on (x,y) from a vector of discrete values */
	  inline integral_function fun() { return integral_function(*this); }

};
    

integral_argument::integral_argument(Eigen::VectorXd &vec, triangulation::face_iterator_t &it)
{  
    // Coordinates of the vertices of the face and value of the solution
    double xA(it->x(0));
    double yA(it->y(0));
    double zA(vec(it->i(0)));
    double xB(it->x(1));
    double yB(it->y(1));
    double zB(vec(it->i(1)));		
    double xC(it->x(2));
    double yC(it->y(2));
    double zC(vec(it->i(2)));
    
    // Impose that (xA,yA,zA), (xB,yB,zB) and (xC,yC,zC) belongs to the approximating function
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
    
    a_coef = coeff(0);
    b_coef = coeff(1);
    d_coef = coeff(2);		
}
    
       
/*! @brief Class that builds a constant approximating function that depends on (x,y) starting from a vector of discrete values */
class int_arg_mod
{
  
  public:
	  /*! @brief Constructor 
	   *  The constructor builds a constant approximation of a function starting from its face values
	   *  @param vec Vector that contains values in discrete points
	   *  @param it Face we are analyzing
	   */
	  inline int_arg_mod(Eigen::VectorXd &vec, triangulation::face_iterator_t &it)
	  {
	      // Since the approximating function is constant we only have its value over the face
	      val=(vec(it->i()));
	  }
    
  private:
	  /*! @brief Value of the function over the face */
	  double val;      

  public:
	  /*! @brief Function class*/
	  class integral_function : public gas::functional::function<2, integral_function>
	  {
	    public:
		    /*! @brief Constructor */
		    inline integral_function(int_arg_mod const & F): F_(F) {};	      
		    
		    /*! @brief Operator at 
		     *  @param x x vlaue of the point
		     *  @param y y value of the point
		     */
		    inline double operator() (double const & x, double const & y) const
		    {
			return F_.val;
		    }
		    
	    private:
		    /*! @brief Integral argument modified object */
		    int_arg_mod const & F_;	  
	  };    
    
	  /*! @brief Return a function depending on (x,y) from a vector of discrete values */
	  inline integral_function fun() { return integral_function(*this); }   
    
};  
  

/*! @brief Abstract class that describes the force term */
class Force
{
  
  public:
	  /*! @brief Constructor 
	   *  @param cdt Triangulation
	   *  @param x_old Solution at previous time step
	   *  @param u0 Initial image restricted to the triangulation
	   *  @param rho Parameter for numerical regularization of dirac delta
	   */
	  Force (triangulation const & cdt, Eigen::VectorXd & x_old, Eigen::VectorXd & u0, double rho): 
	  cdtF_(cdt), x_oldF_(x_old), u0F_(u0), rhoF_(rho) {};
	  
	  /*! @brief Virtual destructor */
	  virtual ~Force() {};
	    
	  /*! @brief Virtual method that returns the force term as a vector */
	  virtual Eigen::VectorXd & vector() = 0;
	  	  
	  /*! @brief Virtual method that returns a criterion to evaluate the segmentation of region 1*/
	  virtual Eigen::VectorXd & p1() = 0;
    
	  /*! @brief Virtual method that returns a criterion to evaluate the segmentation of region 2*/
	  virtual Eigen::VectorXd & p2() = 0;
	    
  protected:
	    /*! @brief Triangulation */
	    triangulation const & cdtF_;    

	    /*! @brief Solution at previous time step */
	    Eigen::VectorXd x_oldF_;
	    
	    /*! @brief Restriction of the initial image to the triangulation */
	    Eigen::VectorXd u0F_;	    
	    
	    /*! @brief Parameter rho */
	    double rhoF_;
	    
  private:
  
};

  
/*! @brief Force term for Chan-Vese model */
class ChanVeseForce : public Force
{
        
  public:
	  /*! @brief Constructor 
	   *  @param cdt Triangulation
	   *  @param x_old Solution at previous time step
	   *  @param u0 Initial image restricted to the triangulation
	   *  @param rho Parameter for numerical regularization of dirac delta
	   */
	  ChanVeseForce (triangulation const & cdt, Eigen::VectorXd & x_old, Eigen::VectorXd & u0, double rho);
    
	  /*! @brief Virtual destructor */
	  virtual ~ChanVeseForce() {};
    
	  /*! @brief Method that returns the force term as a vector */
	  Eigen::VectorXd & vector();
	  	  
	  /*! @brief Returns a criterion to evaluate the segmentation of region 1
	   *  For details, see the treatment of stopping criterion in the related report
	   */
	  inline Eigen::VectorXd & p1() { return p1_; }
    
	  /*! @brief Returns a criterion to evaluate the segmentation of region 2
	   *  For details, see the treatment of stopping criterion in the related report
	   */
	  inline Eigen::VectorXd & p2() { return p2_; }
    
  private:
	  /*! @brief Vector containing force term */
	  Eigen::VectorXd force_;

	  /*! @brief Vector containing the criterion related to region 1 to stop the evolution when the solution is "stationary"*/
	  Eigen::VectorXd p1_;
	  
	  /*! @brief Vector containing the criterion related to region 2 to stop the evolution when the solution is "stationary"*/
	  Eigen::VectorXd p2_;    	    	  
	  
	  /*! @brief Mean value of intensity over region 1*/
	  double mu1_;
	  
	  /*! @brief Mean value of intensity over region 2*/
	  double mu2_;
  
};
  
  
ChanVeseForce::ChanVeseForce (triangulation const & cdt, Eigen::VectorXd & x_old, Eigen::VectorXd & u0, double rho): 
Force(cdt, x_old, u0, rho)
{
    typedef triangulation::face_iterator_t iterator_f;
    typedef triangulation::nodes_iterator_t iterator_n;    
    typedef gas::geometry::unit::triangle triangle_t;
    typedef gas::numerical::quadrature::newton_cotes<triangle_t, 6u> method_t;	
    typedef gas::numerical::quadrature::formula<method_t> integrator_t;

    #define h_rho(v) (0.5*(1+(2*v/M_PI).cwise()))	    
    
    // Criteria related to region 1 and 2 and force term
    p1_.setZero(cdtF_.nodes());
    p2_.setZero(cdtF_.nodes());
    force_.setZero(cdtF_.nodes());
  
    // Compute temporary coefficients vector
    Eigen::VectorXd tmp, arc;
    tmp.setZero(cdtF_.nodes());
    arc.setZero(cdtF_.nodes());
    tmp = x_oldF_/rhoF_;
    arctg(tmp, arc);	
    
    Eigen::VectorXd c0_num_vec, c0_den_vec, c1_num_vec, c1_den_vec;
    c0_num_vec.setZero(cdtF_.nodes());
    c0_den_vec.setZero(cdtF_.nodes());
    c1_num_vec.setZero(cdtF_.nodes());
    c1_den_vec.setZero(cdtF_.nodes());

    c0_num_vec=u0F_*h_rho(arc);	// Heaviside function equals 1 within region 1
    c0_den_vec=h_rho(arc);	
    c1_num_vec=-u0F_*((h_rho(arc)).cwise()-1);	// Heaviside function equals 0 within region 2
    c1_den_vec=-((h_rho(arc)).cwise()-1);	

    // Compute the mean value of intensity in region 1 and region 2
    integrator_t s1, s2, s3, s4, s5;    
    double c0_num(0), c0_den(0), c1_num(0), c1_den(0);
    
    // Loop over the faces to compute the integrals
    for (iterator_f it(cdtF_.face_begin()); it != cdtF_.face_end(); ++it) 
    {    
	// Build the linear approximating functions from the previous vectors
	integral_argument c0_num_fun(c0_num_vec,it);
	integral_argument c0_den_fun(c0_den_vec,it);
	integral_argument c1_num_fun(c1_num_vec,it);	  
	integral_argument c1_den_fun(c1_den_vec,it);	  
	
	s1(*it);
	s2(*it);
	s3(*it);
	s4(*it);

	double const r1(s1.integrate(c0_num_fun.fun()));		
	c0_num += r1;
	
	double const r2(s2.integrate(c0_den_fun.fun()));
	c0_den += r2;
	
	double const r3(s3.integrate(c1_num_fun.fun()));
	c1_num += r3;
	
	double const r4(s4.integrate(c1_den_fun.fun()));
	c1_den += r4;
    }
    
    // Store the mean value of intensity in region 1 and 2
    mu1_=c0_num/c0_den;
    mu2_=c1_num/c1_den;

    // Store the criteria used by stopping_res() function to establish when the solution is "stationary"
    for (iterator_n it(cdtF_.nodes_begin()); it != cdtF_.nodes_end(); ++it) 
    {
	p1_(it->i()) = mu1_;
	p2_(it->i()) = mu2_;
    }    
}  
  
  
Eigen::VectorXd & ChanVeseForce::vector()
{    
    #define right_handCV(a,b) ((((u0F_).cwise()-a)*((u0F_).cwise()-a)) - (((u0F_).cwise()-b)*((u0F_).cwise()-b)))  
    
    force_=right_handCV(mu2_,mu1_);
    
    return force_;
}


  
  
/*! @brief Force term for model based on gaussian regional statistics */
class RegionGaussianForce : public Force
{
        
  public:
	  /*! @brief Constructor 
	   *  @param cdt Triangulation
	   *  @param x_old Solution at previous time step
	   *  @param u0 Initial image restricted to the triangulation
	   *  @param rho Parameter for numerical regularization of dirac delta
	   */
	  RegionGaussianForce (triangulation const & cdt, Eigen::VectorXd & x_old, Eigen::VectorXd & u0, double rho);
	  
	  /*! @brief Virtual destructor*/
	  virtual ~RegionGaussianForce() {};
	  
	  /*! @brief Method that returns the force term as a vector */
	  Eigen::VectorXd & vector();	  
	  
	  /*! @brief Returns a criterion to evaluate the segmentation of region 1
	   *  For details, see the treatment of stopping criterion in the related report
	   */
	  inline Eigen::VectorXd & p1() { return p1_; }
    
	  /*! @brief Returns a criterion to evaluate the segmentation of region 2
	   *  For details, see the treatment of stopping criterion in the related report
	   */
	  inline Eigen::VectorXd & p2() { return p2_; }
       
  private:
	  /*! @brief Vector containing force term */
	  Eigen::VectorXd force_;

	  /*! @brief Vector containing the criterion related to region 1 to stop the evolution when the solution is "stationary"*/
	  Eigen::VectorXd p1_;
	  
	  /*! @brief Vector containing the criterion related to region 2 to stop the evolution when the solution is "stationary"*/
	  Eigen::VectorXd p2_;    	    	  
	  
	  /*! @brief Mean value of intensity over region 1*/
	  double mu1_;
	  
	  /*! @brief Mean value of intensity over region 2*/
	  double mu2_;
	  
	  /*! @brief Standard deviation of intensity over region 1*/
	  double sig1_;
	  
	  /*! @brief Standard deviation of intensity over region 2*/
	  double sig2_;
  
};
  

RegionGaussianForce::RegionGaussianForce (triangulation const & cdt, Eigen::VectorXd & x_old, Eigen::VectorXd & u0, double rho): 
Force(cdt, x_old, u0, rho)
{
    typedef triangulation::face_iterator_t iterator_f;
    typedef gas::geometry::unit::triangle triangle_t;
    typedef gas::numerical::quadrature::newton_cotes<triangle_t, 6u> method_t;	
    typedef gas::numerical::quadrature::formula<method_t> integrator_t;

    #define h_rho(v) (0.5*(1+(2*v/M_PI).cwise()))	    
    #define GP(m,s,z) ((std::exp(-std::pow(z-m,2)/(2*s*s)))/(std::sqrt(2*M_PI)*s))		    
    
    // Criteria related to region 1 and 2 and force term
    p1_.setZero(cdtF_.nodes());
    p2_.setZero(cdtF_.nodes());
    force_.setZero(cdtF_.nodes());
  
    // Compute temporary coefficients vector
    Eigen::VectorXd tmp, arc;
    tmp.setZero(cdtF_.nodes());
    arc.setZero(cdtF_.nodes());
    tmp = x_oldF_/rhoF_;
    arctg(tmp, arc);	    

    Eigen::VectorXd mm1_num_vec, mm1_den_vec, mm2_num_vec, mm2_den_vec;
    mm1_num_vec.setZero(cdtF_.nodes());
    mm1_den_vec.setZero(cdtF_.nodes());
    mm2_num_vec.setZero(cdtF_.nodes());
    mm2_den_vec.setZero(cdtF_.nodes());
    
    Eigen::VectorXd ss1_num_vec, ss1_den_vec, ss2_num_vec, ss2_den_vec;
    ss1_num_vec.setZero(cdtF_.nodes());
    ss1_den_vec.setZero(cdtF_.nodes());
    ss2_num_vec.setZero(cdtF_.nodes());
    ss2_den_vec.setZero(cdtF_.nodes());
    
    mm1_num_vec=u0F_*h_rho(arc); // Heaviside function equals 1 within region 1
    mm1_den_vec=h_rho(arc);	
    mm2_num_vec=-u0F_*((h_rho(arc)).cwise()-1);	// Heaviside function equals 0 within region 2
    mm2_den_vec=-((h_rho(arc)).cwise()-1);	

    // Compute the mean value of intensity in region 1 and region 2   
    integrator_t s1, s2, s3, s4, s5;
    double mm1_num(0), mm1_den(0), mm2_num(0), mm2_den(0);
    double ss1_num(0), ss1_den(0), ss2_num(0), ss2_den(0);    
    
    // Loop over the faces to compute the integrals
    for (iterator_f it(cdtF_.face_begin()); it != cdtF_.face_end(); ++it) 
    {
	// Build the linear approximating functions from the previous vectors
	integral_argument mm1_num_fun(mm1_num_vec,it);
	integral_argument mm1_den_fun(mm1_den_vec,it);
	integral_argument mm2_num_fun(mm2_num_vec,it);	  
	integral_argument mm2_den_fun(mm2_den_vec,it);	  
	
	s1(*it);
	s2(*it);
	s3(*it);
	s4(*it);

	double const r1(s1.integrate(mm1_num_fun.fun()));		
	mm1_num += r1;
	
	double const r2(s2.integrate(mm1_den_fun.fun()));
	mm1_den += r2;
	
	double const r3(s3.integrate(mm2_num_fun.fun()));
	mm2_num += r3;
	
	double const r4(s4.integrate(mm2_den_fun.fun()));
	mm2_den += r4;
    }
    
    // Store the mean value of intensity in region 1 and 2
    mu1_=mm1_num/mm1_den;
    mu2_=mm2_num/mm2_den;
	    
    // Compute temporary coefficients
    Eigen::VectorXd d1, d2;
    d1.setZero(cdtF_.nodes());
    d2.setZero(cdtF_.nodes());	
    
    pow(((u0F_).cwise()-mu1_),2.,d1);
    pow(((u0F_).cwise()-mu2_),2.,d2);	

    ss1_num_vec=d1*h_rho(arc);	// Heaviside function equals 1 within region 1
    ss1_den_vec=h_rho(arc);
    ss2_num_vec=-d2*((h_rho(arc)).cwise()-1); // Heaviside function equals 0 within region 2
    ss2_den_vec=-((h_rho(arc)).cwise()-1);
    
    // Loop over the faces to compute the integrals
    for (iterator_f it(cdtF_.face_begin()); it != cdtF_.face_end(); ++it) 
    {
	// Build the linear approximating functions from the previous vectors
	integral_argument ss1_num_fun(ss1_num_vec,it);
	integral_argument ss1_den_fun(ss1_den_vec,it);
	integral_argument ss2_num_fun(ss2_num_vec,it);	  
	integral_argument ss2_den_fun(ss2_den_vec,it);	  
	
	s1(*it);
	s2(*it);
	s3(*it);
	s4(*it);

	double const r1(s1.integrate(ss1_num_fun.fun()));		
	ss1_num += r1;
	
	double const r2(s2.integrate(ss1_den_fun.fun()));
	ss1_den += r2;
	
	double const r3(s3.integrate(ss2_num_fun.fun()));
	ss2_num += r3;
	
	double const r4(s4.integrate(ss2_den_fun.fun()));
	ss2_den += r4;
    }	
    
    // Store the standard deviation value of intensity in region 1 and 2
    sig1_=ss1_num/ss1_den;
    sig2_=ss2_num/ss2_den;
    
    sig1_ = std::sqrt(sig1_);
    sig2_ = std::sqrt(sig2_);

    // Build probability density of intensities using gaussian distribution
    std::vector<double> _p1(256,0.0);
    std::vector<double> _p2(256,0.0);      
    
    for(int ii=0; ii<256; ii++)
    {
      _p1[ii] = GP(mu1_,sig1_,double(ii)/255.0);
      _p2[ii] = GP(mu2_,sig2_,double(ii)/255.0);	  
    }

    // Export probability densities to file
    std::ofstream density1_file("result/density1.dat", std::ios::app);
    density1_file << _p1.size() << std::endl;
    for(std::vector<double>::iterator _v=_p1.begin(); _v != _p1.end(); _v++)
    {
      density1_file << *_v << std::endl;
    }
    density1_file << std::endl;
    density1_file << std::endl;      
    density1_file.close();
    
    std::ofstream density2_file("result/density2.dat", std::ios::app);
    density2_file << _p2.size() << std::endl;
    for(std::vector<double>::iterator _v=_p2.begin(); _v != _p2.end(); _v++)
    {
      density2_file << *_v << std::endl;
    }
    density2_file << std::endl;
    density2_file << std::endl;      
    density2_file.close();
  
    // Store the criteria used by stopping_res() function to establish when the solution is "stationary"
    Eigen::VectorXd pw1, pw2;
    pw1.setZero(cdtF_.nodes());
    pw2.setZero(cdtF_.nodes());
    
    pow((u0F_.cwise()-mu1_),2.,pw1);      	
    pow((u0F_.cwise()-mu2_),2.,pw2);      	
    
    Eigen::VectorXd e1, e2;
    e1.setZero(cdtF_.nodes());
    e2.setZero(cdtF_.nodes());
    
    esp(-pw1/(2*sig1_*sig1_),e1);
    esp(-pw2/(2*sig2_*sig2_),e2);
    
    p1_ = e1/(std::sqrt(2*M_PI)*sig1_);
    p2_ = e2/(std::sqrt(2*M_PI)*sig2_);	
}
  
  
Eigen::VectorXd & RegionGaussianForce::vector()
{
    // Compute the log-densities from previously determined probability distributions
    Eigen::VectorXd lp1, lp2;
    lp1.setZero(cdtF_.nodes());
    lp2.setZero(cdtF_.nodes());
    ln(p1_,lp1);
    ln(p2_,lp2);
	    
    force_= lp1 - lp2;	
    
    return force_;
}	    
  
  
/*! @brief Force term for model based on regional statistics with probability densities estimated using a non-parametric approach*/
class RegionNonParametricForce : public Force
{
        
  public:    
	  /*! @brief Constructor 
	   *  @param cdt Triangulation
	   *  @param x_old Solution at previous time step
	   *  @param u0 Initial image restricted to the triangulation
	   *  @param rho Parameter for numerical regularization of dirac delta
	   */
	  RegionNonParametricForce (triangulation const & cdt, Eigen::VectorXd & x_old, initial_image & img_in, Eigen::VectorXd & u0, double rho);

	  /*! @brief Virtual destructor*/
	  virtual ~RegionNonParametricForce() {};
	  
	  /*! @brief Method that returns the force term as a vector */
	  Eigen::VectorXd & vector();	  
	  
	  /*! @brief Returns a criterion to evaluate the segmentation of region 1
	   *  For details, see the treatment of stopping criterion in the related report
	   */
	  inline Eigen::VectorXd & p1() { return p1_; }
    
	  /*! @brief Returns a criterion to evaluate the segmentation of region 2
	   *  For details, see the treatment of stopping criterion in the related report
	   */
	  inline Eigen::VectorXd & p2() { return p2_; }   
	  	
  private:
	  /*! @brief Vector containing force term */
	  Eigen::VectorXd force_;
	  
	  /*! @brief Initial image */
	  initial_image img_inF_;
	  
	  /*! @brief Vector containing the criterion related to region 1 to stop the evolution when the solution is "stationary"*/
	  Eigen::VectorXd p1_;
	  
	  /*! @brief Vector containing the criterion related to region 2 to stop the evolution when the solution is "stationary"*/
	  Eigen::VectorXd p2_;    	    	  
	  	  
};
  

RegionNonParametricForce::RegionNonParametricForce (triangulation const & cdt, Eigen::VectorXd & x_old, initial_image & img_in, Eigen::VectorXd & u0, double rho): 
Force(cdt, x_old, u0, rho), img_inF_(img_in)
{
    typedef triangulation::nodes_iterator_t iterator_n;
  
    // Criteria related to region 1 and 2 and force term
    p1_.setZero(cdtF_.nodes());
    p2_.setZero(cdtF_.nodes());
    force_.setZero(cdtF_.nodes());
    
    // Standard deviation of gaussian kernel
    double sigma(8.0);
    
    // Numerical thershold
    double epsilon(std::pow(10,-6));
    
    // Initialize objects needed to estimate the probability densities
    histogram in(0.0,255.0,256);
    histogram out(0.0,255.0,256);      

    std::vector<double> inVal(256,0.0);
    std::vector<double> outVal(256,0.0);
    
    in.get_values(inVal);
    out.get_values(outVal);

    std::fill(inVal.begin(), inVal.end(), 1); 
    std::fill(outVal.begin(), outVal.end(), 1);       
    
    int areaIN(0);
    int areaOUT(0);
      
    // Build the histograms linked to the pixel intensities
    for(iterator_n it(cdtF_.nodes_begin()); it != cdtF_.nodes_end(); ++it)
    {
	// Rescale intensity
	double v = u0F_(it->i())*255.0;
	
	int res = int((v - 0.0) * in.binValue());
	if (res >= in.count_bins()) { --res; }	  
	
	// Identify which region the pixel belongs to
	if (x_oldF_(it->i()) > 0)
	{
	    inVal[res] += 1;
	    areaIN += 1;
	}
	else
	{
	    outVal[res] += 1;
	    areaOUT += 1;
	}	  
    }
    
    double areaIN_inv = areaIN != 0 ? double(1) / double(areaIN) : 1;
    double areaOUT_inv = areaOUT != 0 ? double(1) / double(areaOUT) : 1;
        
    for (int i = 0; i < 256; ++i)
    {
	inVal[i] *= areaIN_inv;
	outVal[i] *= areaOUT_inv;
    }
    
    // Store the vector using an histogram object
    in.set_values(inVal);
    out.set_values(outVal);
    
    // Gaussian smoothing
    in.gaussK_smooth(sigma);
    out.gaussK_smooth(sigma);      
    
    // Build probability densities 
    std::vector<double> _p1(256,0.0);
    std::vector<double> _p2(256,0.0);      

    // Build probabilities based on the information about the bin where the pixel belongs to
    for(iterator_n it(cdtF_.nodes_begin()); it != cdtF_.nodes_end(); ++it)
    {
	double intensity(u0F_(it->i())*255.0);
	
	int int_int = int(intensity);
	
	double res = in.minValue() + (double(int_int) + 0.5)/in.binValue();
		
	double ival = in.get_bin_area(res);
	double oval = out.get_bin_area(res);

	ival = (ival > epsilon) ? ival : epsilon;
	oval = (oval > epsilon) ? oval : epsilon;	  
	
	_p1[int_int] = ival;
	_p2[int_int] = oval;
    }
    
    // Export probability densities to file
    std::ofstream density1_file("result/density1.dat", std::ios::app);
    density1_file << _p1.size() << std::endl;
    for(std::vector<double>::iterator _v=_p1.begin(); _v != _p1.end(); _v++)
    {
      density1_file << *_v << std::endl;
    }
    density1_file << std::endl;
    density1_file << std::endl;      
    density1_file.close();
    
    std::ofstream density2_file("result/density2.dat", std::ios::app);
    density2_file << _p2.size() << std::endl;
    for(std::vector<double>::iterator _v=_p2.begin(); _v != _p2.end(); _v++)
    {
      density2_file << *_v << std::endl;
    }
    density2_file << std::endl;
    density2_file << std::endl;      
    density2_file.close();
    
    // Store the criteria used by stopping_res() function to establish when the solution is "stationary"
    for(iterator_n it(cdtF_.nodes_begin()); it != cdtF_.nodes_end(); ++it)
    {
	double intensity(u0F_(it->i())*255.0);
	
	int res = int((intensity - 0.0) * in.binValue());
	if (res >= in.count_bins()) { --res; }	  	  
	
	p1_(it->i()) = _p1[res];
	p2_(it->i()) = _p2[res];	  
    }      
}


Eigen::VectorXd & RegionNonParametricForce::vector()
{
    // Compute the log-densities from previously determined probability distributions
    Eigen::VectorXd lp1, lp2;
    lp1.setZero(cdtF_.nodes());
    lp2.setZero(cdtF_.nodes());		
    ln(p1_,lp1);
    ln(p2_,lp2);
    
    force_= lp1 - lp2;	

    return force_;
}
  
    
/*! @brief Factory to build force term, based on the model problem in analysis*/    
class ForceFactory
{	
  
    public:
	    /*! @brief Constructor 
	     * 	@param type Type of force term, depending on the model problem
	     *  @param cdt Triangulation
	     *  @param x_old Solution at previous time step
	     *  @param img_in Initial image
	     *  @param u0 Restriction of the initial image to the triangulation
	     *  @param rho Parameter rho for the regularization of dirac delta
	     */
	    ForceFactory (int type, triangulation const & cdt, Eigen::VectorXd & x_old, initial_image & img_in, Eigen::VectorXd & u0, double rho): 
	    type_(type), cdtF_(cdt), x_oldF_(x_old), img_inF_(img_in), u0F_(u0), rhoF_(rho) {};
	
	    /*! @brief Virtual destructor*/
	    virtual ~ForceFactory() {};
	    
	    /*! @brief Create a force term*/
	    Force* create();
	
    private:
	    /*! @brief Parameter to choose which force term has to be considered*/
	    int type_;
      
	    /*! @brief Triangulation */
	    triangulation const & cdtF_;
	    
	    /*! @brief Solution at previous time step*/
	    Eigen::VectorXd x_oldF_;
		  
	    /*! @brief Initial image*/
	    initial_image img_inF_;		
	    
	    /*! @brief Restriction of the initial image to the triangulation */
	    Eigen::VectorXd u0F_;
	    
	    /*! @brief Parameter rho*/
	    double rhoF_;
	    
};
    
    
Force* ForceFactory::create()
{
    // Depending on the parameter type_, calls a different constructor
    if (type_ == 0) { return new ChanVeseForce(cdtF_, x_oldF_, u0F_, rhoF_); }
    else if (type_ == 1) { return new RegionGaussianForce(cdtF_, x_oldF_, u0F_, rhoF_); }
    else if (type_ == 2) { return new RegionNonParametricForce (cdtF_, x_oldF_, img_inF_, u0F_, rhoF_); }
    
    return 0;
}
    

} //namespace


#endif