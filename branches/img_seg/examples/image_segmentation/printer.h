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
 * @file printer.h
 * @brief Includes the implementation of the class to print output to svg
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_PRINTER_H
#define _IMG_SEG_PRINTER_H


namespace img_seg {

  
/*! @brief Class to manage svg output files */	
class svg {

  public:
	  /*! @brief Constructor 
	  *  @param p Problem
	  */
	  svg (problem & p): p_(p), x_(1024u), y_(768u) {};
	  
	  /*! @brief Destructor */
	  ~svg() {};
	  
	  /*! @brief Convert solution color and then to string */	
	  std::string to_string ();

	  /*! @brief Operator << overloading 
	  *   @param out Ostream object
	  *   @param text Object to put in out
	  */	
	  friend std::ostream & operator<< (std::ostream & out, svg & text);	

  private:
	  /*! @brief Problem */
	  problem & p_;

	  /*! @brief Output width */	
	  int x_;

	  /*! @brief Output height */		
	  int y_;
	
};


std::string svg::to_string ()
{  
    typedef triangulation::face_iterator_t iterator_f;
  
    // Style file
    gas::css style;

    // Style for the triangulation 
    style.new_id("grid");
	style.property("stroke", "red");
	style.property("stroke-width", 2);
	style.property("fill", "none");

    // Style for the solution
    style.new_id("solution");
	style.property("fill", "red");

    // Style for the edges
    style.new_id("edges");
	style.property("stroke", "yellow");
	style.property("stroke-width", 4);		

    // Color range
    gas_rangeu(i,256)
    {
	std::ostringstream name;
	std::ostringstream color;

	name << "gray" << i;
	color << "rgb(" << i << ", " << i << ", " << i << ")";

	style.new_class(name.str().c_str());
	style.property("fill", color.str().c_str());
    }

    // Triangulation
    triangulation const & cdt(p_.mesh());
    
    iterator_f i(cdt.face_begin());    

    // Dimensions of the output 
    double _xmin, _xmax, _ymin, _ymax, _rmin, _rmax;

    // Coordinates of the points and values of the solution for face i
    Eigen::Vector3d const x(i->x(0), i->x(1), i->x(2));
    Eigen::Vector3d const y(i->y(0), i->y(1), i->y(2));
    Eigen::Vector3d r;
    
    // Here we set the variable to plot
    // If problem hasn't been solved yet, we plot the initial image or the initial level-set curve
    if(!p_.already_solved())
    {
	// The first time we plot the initial image restricted to the given triangulation
	if(!p_.not_original())
	{
	    r(0)=p_.init_img_vec(i->i(0));
	    r(1)=p_.init_img_vec(i->i(1));
	    r(2)=p_.init_img_vec(i->i(2));	  
	}
	else // We plot the initial level-set curve on the given mesh
	{
	    r(0)=p_.old_solution(i->i(0));
	    r(1)=p_.old_solution(i->i(1));
	    r(2)=p_.old_solution(i->i(2));	      
	}
    }
    else // we plot the solution at current time
    {
	r(0)=p_(i->i(0));
	r(1)=p_(i->i(1));
	r(2)=p_(i->i(2));
    }

    // Set temporary dimension for the output    
    _xmin = x.minCoeff();
    _xmax = x.maxCoeff();

    _ymin = y.minCoeff();
    _ymax = y.maxCoeff();

    _rmin = r.minCoeff();
    _rmax = r.maxCoeff();

    ++i;

    // Loop over the faces to set the variable to plot
    for (; i != cdt.face_end(); ++i)
    {
	// Coordinates of the points and values of the solution for face i      
	Eigen::Vector3d const x(i->x(0), i->x(1), i->x(2));
	Eigen::Vector3d const y(i->y(0), i->y(1), i->y(2));
	Eigen::Vector3d r;
	
	if(!p_.already_solved())
	{
	    // The first time we plot the initial image restricted to the given triangulation
	    if(!p_.not_original())
	    {
		r(0)=p_.init_img_vec(i->i(0));
		r(1)=p_.init_img_vec(i->i(1));
		r(2)=p_.init_img_vec(i->i(2));	  
	    }
	    else // We plot the initial level-set curve on the given mesh
	    {
		r(0)=p_.old_solution(i->i(0));
		r(1)=p_.old_solution(i->i(1));
		r(2)=p_.old_solution(i->i(2));	      
	    }
	}
	else // We plot the solution at current time
	{
	    r(0)=p_(i->i(0));
	    r(1)=p_(i->i(1));
	    r(2)=p_(i->i(2));
	}

	// Update temporary dimension for the output    
	_xmin = std::min(_xmin, x.minCoeff());
	_xmax = std::max(_xmax, x.maxCoeff());

	_ymin = std::min(_ymin, y.minCoeff());
	_ymax = std::max(_ymax, y.maxCoeff());

	_rmin = std::min(_rmin, r.minCoeff());
	_rmax = std::max(_rmax, r.maxCoeff());
    }

    // Resolution and aspect ratio
    int _x, _bx;
    int _y, _by;

    double const _dx(_xmax - _xmin);
    double const _dy(_ymax - _ymin);
    double const _dr(_rmax - _rmin);

    double const _ratio(_dx / _dy);

    if (_ratio < 4./3.)
    {
	_x = 762u * _ratio;
	_y = 762u;
    }
    else
    {
	_x = 1016u;
	_y = 1016u / _ratio;
    }

    _bx = 3u * _ratio;
    _by = 3u;	
    
    // Create SVG file
    gas::svg svg(_x + (2u * _bx), _y + (2u * _by));

    svg.title("Solution");
    svg.description("Level-set based image segmentation");

    svg.style(style);
    
    // Solution //
    svg.open_group("solution");

    for (iterator_f i(cdt.face_begin()); i != cdt.face_end(); ++i)
    {
	// Point to be plotted
	int x[3];
	int y[3];

	x[0] = _bx + _x * (i->x(0) - _xmin) / _dx;
	x[1] = _bx + _x * (i->x(1) - _xmin) / _dx;
	x[2] = _bx + _x * (i->x(2) - _xmin) / _dx;

	y[0] = _by + _y * (_ymax - i->y(0)) / _dy;
	y[1] = _by + _y * (_ymax - i->y(1)) / _dy;
	y[2] = _by + _y * (_ymax - i->y(2)) / _dy;
	
	// Set the variable to plot
	double r;
	
	if(!p_.already_solved())
	{
	      if(!p_.not_original())
	      {
		  r=(p_.init_img_vec(i->i(0)) + p_.init_img_vec(i->i(1)) + p_.init_img_vec(i->i(2))) / 3.;			
	      }
	      else
	      {
		  r=(p_.old_solution(i->i(0)) + p_.old_solution(i->i(1)) + p_.old_solution(i->i(2))) / 3.;
	      }
	}
	else
	{
	      r=(p_(i->i(0)) + p_(i->i(1)) + p_(i->i(2))) / 3.;
	}
	  
	int c(255 * (r - _rmin) / _dr);

	std::ostringstream color;
	color << "gray" << c;

	// Color the triangle
	svg.triangle(x[0], y[0], x[1] ,y[1], x[2], y[2], color.str().c_str());
    }

    svg.close_group();
    
    if(p_.is_refinement()) // We plot the mesh only after refinement
    {
	// Mesh //
	svg.open_group("grid");

	for (iterator_f i(cdt.face_begin()); i != cdt.face_end(); ++i)
	{
	    // Point to be plotted
	    int x[3];
	    int y[3];

	    x[0] = _bx + _x * (i->x(0) - _xmin) / _dx;
	    x[1] = _bx + _x * (i->x(1) - _xmin) / _dx;
	    x[2] = _bx + _x * (i->x(2) - _xmin) / _dx;

	    y[0] = _by + _y * (_ymax - i->y(0)) / _dy;
	    y[1] = _by + _y * (_ymax - i->y(1)) / _dy;
	    y[2] = _by + _y * (_ymax - i->y(2)) / _dy;
	    
	    // Color the triangle
	    svg.triangle(x[0], y[0], x[1] ,y[1], x[2], y[2]);
	}

	svg.close_group();
    }
    else
    {
	// Level-set curve //
	svg.open_group("edges");
	
	for (iterator_f i(cdt.face_begin()); i != cdt.face_end(); ++i)
	{
	    // Points to be plotted
	    int x[3];
	    int y[3];	      

	    x[0] = _bx + _x * (i->x(0) - _xmin) / _dx;
	    x[1] = _bx + _x * (i->x(1) - _xmin) / _dx;
	    x[2] = _bx + _x * (i->x(2) - _xmin) / _dx;
      
	    y[0] = _by + _y * (i->y(0) - _ymin) / _dy;
	    y[1] = _by + _y * (i->y(1) - _ymin) / _dy;
	    y[2] = _by + _y * (i->y(2) - _ymin) / _dy;	  
    
	    // Find which edges are part of the level-set curve
	    if(!p_.already_solved())
	    {
		  if(!p_.not_original()) // using initial image restricted to the triangulation
		  {
		      double m((p_.init_img_vec(i->i(0)) + p_.init_img_vec(i->i(1)) + p_.init_img_vec(i->i(2)))/3.0);
		      
		      // Edge 0-1
		      double v((p_.init_img_vec(i->neighbor_idx(2,0)) + p_.init_img_vec(i->neighbor_idx(2,1)) + p_.init_img_vec(i->neighbor_idx(2,2)))/3.0);
		      if (m*v < 0) { svg.line(x[0], y[0], x[1] ,y[1], "yellow"); }
		      
		      // Edge 1-2
		      v = (p_.init_img_vec(i->neighbor_idx(0,0)) + p_.init_img_vec(i->neighbor_idx(0,1)) + p_.init_img_vec(i->neighbor_idx(0,2)))/3.0;			
		      if (m*v < 0) { svg.line(x[1], y[1], x[2] ,y[2], "yellow"); }	  
		      
		      // Edge 2-0
		      v = (p_.init_img_vec(i->neighbor_idx(1,0)) + p_.init_img_vec(i->neighbor_idx(1,1)) + p_.init_img_vec(i->neighbor_idx(1,2)))/3.0;
		      if (m*v < 0) { svg.line(x[2], y[2], x[0] ,y[0], "yellow"); }		      
		  }
		  else
		  {
		      double m((p_.old_solution(i->i(0)) + p_.old_solution(i->i(1)) + p_.old_solution(i->i(2)))/3.0);
		      
		      // Edge 0-1
		      double v((p_.old_solution(i->neighbor_idx(2,0)) + p_.old_solution(i->neighbor_idx(2,1)) + p_.old_solution(i->neighbor_idx(2,2)))/3.0);
		      if (m*v < 0) { svg.line(x[0], y[0], x[1] ,y[1], "yellow"); }
		      
		      // Edge 1-2
		      v = (p_.old_solution(i->neighbor_idx(0,0)) + p_.old_solution(i->neighbor_idx(0,1)) + p_.old_solution(i->neighbor_idx(0,2)))/3.0;
		      if (m*v < 0) { svg.line(x[1], y[1], x[2] ,y[2], "yellow"); }

		      // Edge 2-0
		      v = (p_.old_solution(i->neighbor_idx(1,0)) + p_.old_solution(i->neighbor_idx(1,1)) + p_.old_solution(i->neighbor_idx(1,2)))/3.0;			
		      if (m*v < 0) { svg.line(x[2], y[2], x[0] ,y[0], "yellow"); }
		  }
	    }
	    else
	    {
		double m((p_(i->i(0)) + p_(i->i(1)) + p_(i->i(2)))/3.0);
		
		// Edge 0-1
		double v((p_(i->neighbor_idx(2,0)) + p_(i->neighbor_idx(2,1)) + p_(i->neighbor_idx(2,2)))/3.0);		  
		if (m*v < 0) { svg.line(x[0], y[0], x[1] ,y[1], "yellow"); }

		// Edge 1-2
		v  =(p_(i->neighbor_idx(0,0)) + p_(i->neighbor_idx(0,1)) + p_(i->neighbor_idx(0,2)))/3.0;
		if (m*v < 0) { svg.line(x[1], y[1], x[2] ,y[2], "yellow"); }		  
		
		// Edge 2-0
		v = (p_(i->neighbor_idx(1,0)) + p_(i->neighbor_idx(1,1)) + p_(i->neighbor_idx(1,2)))/3.0;		  
		if (m*v < 0) { svg.line(x[2], y[2], x[0] ,y[0], "yellow"); }		  
	    }       	
	}
    
	svg.close_group();
    }

    // Return SVG file
    std::ostringstream out;
    out << svg;
    return out.str();
}


std::ostream & operator<< (std::ostream & out, svg & text)
{
    out << text.to_string();
    return out;
}


} //namespace


#endif 