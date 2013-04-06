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
 * @file initials.h
 * @brief Includes the implementation of the classes to manage initial curves and images
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_INITIALS_H
#define _IMG_SEG_INITIALS_H


namespace img_seg {


/*! @brief Initial image */
class initial_image
{
  
  public:
	  /*! @brief Constructor
	   *  @param img CImg object from CImg library for the treatment of images
	   *  @param w Image width
	   *  @param h Image height
	   */
	  inline initial_image(const cimg_library::CImg<double> & img, double const w, double const h): img_(img), w_(w), h_(h) {} 

	  /*! @brief Destructor */
	  ~initial_image() {};
	  
	  /*! @brief Restrict CImg to the computational mesh
	   *  @param cdt_ Triangulation to restrict original image to
	   *  @param vec Vector to store image restricted to the mesh
	   */	  
	  void Trestrict(triangulation const & cdt_, Eigen::VectorXd & vec) const;	  
	  
	  /*! @brief Return image width */			  
	  inline double w() const { return w_; }

	  /*! @brief Return image height */			  
	  inline double h() const { return h_; }	
	  
  private:	
	  /*! @brief CImge image */
	  const cimg_library::CImg<double> img_;
	  
	  /*! @brief Image width */
	  double const w_;
	  
	  /*! @brief Image height */
	  double const h_;

  public:
	  /*! @brief A wrapper for describing a vector of spatial values as a function of space */
	  class integral_function : public gas::functional::function<2, integral_function>
	  {
	    public:
		    /*! @brief Constructor */
		    inline integral_function(initial_image const & I): I_(I) {};   
	      
		    /*! @brief Destructor */
		    ~integral_function() {};
		    
		    /*! @brief Operator at 
		     *  @param x x-coordinate
		     *  @param y y-coordinate
		     */
		    inline double operator() (double const & x, double const & y) const
		    {
		      return (I_.img_(int((x+I_.w()/2.0)*I_.img_.width()),  int(I_.img_.height() - (y+I_.h()/2.0)*I_.img_.height())))/255.0;		  
		    }
		   
	    private:
		    /*! @brief Initial image class */
		    initial_image const & I_;
	    
	  };  

	  
	  /*! @brief Return a function built from linear approximation of the nodal values of the initial image */
	  inline integral_function fun() { return integral_function(*this); }
	  
};


void initial_image::Trestrict(triangulation const & cdt_, Eigen::VectorXd & vec) const
{	
	typedef triangulation::nodes_iterator_t iterator_n;	

	// Store the values of the image intensity only for points within the triangulation
	for(iterator_n j_(cdt_.nodes_begin()); j_!=cdt_.nodes_end(); ++j_)
	{
	    vec(j_->i())=(img_(int((j_->x()+w_/2.0)*img_.width()),  int(img_.height() - (j_->y()+h_/2.0)*img_.height()))/255.0);
	}
}


/*! @brief Initial level-set curve */
class initial_curve
{

  public:
	  /*! @brief Constructor 
	  *   @param cdt_ Triangulation
	  *   @param w Image width
	  *   @param h Image height
	  */
	  initial_curve(triangulation const &cdt_, double const w, double const h);
	  
	  /*! @brief Destructor */
	  ~initial_curve() {};
		  
	  /*! @brief Export level-set initial vector */	
	  void vectorize(Eigen::VectorXd & vec) const { vec = phi_; }

	  /*! @brief Return image width */		
	  inline double w() const { return w_; }

	  /*! @brief Return image height */		
	  inline double h() const { return h_; }	
	  
  private:
	  /*! @brief Level-set initial vector */
	  Eigen::VectorXd phi_;
	  
	  /*! @brief Image width */
	  double const w_;
	  
	  /*! @brief Image height */
	  double const h_;	
    
};


initial_curve::initial_curve(triangulation const &cdt_, double const w, double const h): w_(w), h_(h)
{
    typedef triangulation::nodes_iterator_t iterator_n;

    phi_.setZero(cdt_.nodes());

    // Possible initial radiuses for level-set curve
    double r(std::min(w_,h_)/2.0);	
//     double r(std::min(w_,h_)/8.0);		
    
    // Initialization of the initial level-set curve
    for(iterator_n i_(cdt_.nodes_begin()); i_!=cdt_.nodes_end(); ++i_)
    {
	// Computing distance from local center
	double b2(std::pow(i_->x(),2.)+std::pow(i_->y(),2.));

	// Setting initial level-set curve as the signed distance 
	// between the point x and the level-set curve
	if(b2-std::pow(r,2.)<0)
	      phi_(i_->i())=r-std::sqrt(b2);
	else
	      phi_(i_->i())=-(std::sqrt(b2)-r);
      
      // Other initialization (Multiple circles inside the domain)
      /*
      double alpha = w_/4.;
      double beta = h_/4.;
      
      double xc;
      double yc;
      
      // Computing the center of the circles - x
      if((i_->x() >= -2.*alpha) && (i_->x() < -alpha))
      {
	xc = -1.5*alpha;	  
      }
      else if((i_->x() >= -alpha) && (i_->x() < 0))
      {
	xc = -0.5*alpha;	  	  
      }
      else if((i_->x() >= 0) && (i_->x() < alpha))
      {
	xc = 0.5*alpha;	  	  
      }
      else
      {
	xc = 1.5*alpha;	  	  
      }
      
      // Computing the center of the circles - y
      if((i_->y() >= -2.*beta) && (i_->y() < -beta))
      {
	yc = -1.5*beta;	  
      }
      else if((i_->y() >= -beta) && (i_->y() < 0))
      {
	yc = -0.5*beta;	  	  
      }
      else if((i_->y() >= 0) && (i_->y() < beta))
      {
	yc = 0.5*beta;	  	  
      }
      else
      {
	yc = 1.5*beta;	  	  
      }
      
      // Computing radius
      double rc(std::min(alpha,beta)/4.0);	    
      
      // Computing distance from local center
      double bc2(std::pow(i_->x()-xc,2.)+std::pow(i_->y()-yc,2.));

      // Setting initial level-set curve as the signed distance 
      // between the point x and the level-set curve
      if(bc2-std::pow(rc,2.)<0)
	    phi_(i_->i())=rc-std::sqrt(bc2);
      else
	    phi_(i_->i())=-(std::sqrt(bc2)-rc);
      */
    }

}


} //namespace


#endif