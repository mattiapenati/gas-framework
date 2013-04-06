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
 * @file refinement.h
 * @brief Includes the implementation of the class to perform the refinement of the elements
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_REFINEMENT_H
#define _IMG_SEG_REFINEMENT_H


namespace img_seg {


/*! @brief Abstract class that describes a refinement method */
class RefineMethod
{

  public: 	       
	  typedef triangulation::face_iterator face_iterator_t;    
    
	/*! @brief Constructor 
	  *  @param p Problem
	  *  @param cdt Triangulation
	  *  @param TbR Vector that keeps the information about whether or not elements have to be refined
	  *  @param to_be_inserted Vector that stores the elements to be inserted
	  *  @param to_be_removed Vector that stores the edges to be removed
	  *  @param new_points Vector that stores the values of the solution interpolated on new points in the triangulation
	  *  @param init Boolean variable to flag if refinement is performed at the beginning or duriong evolution loop
	  */
	  inline RefineMethod ( problem const & p, triangulation & cdt, std::vector<bool> & TbR, std::vector<std::list<std::pair<double,double> > > & to_be_inserted, std::vector<std::pair<bool,std::list<int> > > & to_be_removed, std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & new_points, const bool init ): 
	  p_(p), cdt_(cdt), TbR_(TbR), to_be_inserted_(to_be_inserted), to_be_removed_(to_be_removed), new_points_(new_points), init_(init) {};
	  
	  /*! @brief Virtual destructor */
	  virtual ~RefineMethod() {};

	  /*! @brief Virtual method to apply the refinement */
	  virtual void apply() = 0;
    
	  /*! @brief Virtual method to update labelas and reference edges */
	  virtual void update() = 0;	  
	  
    
	  /*! @brief Insert an edge in the vector of edges to be inserted during refinement procedure
	   *  @param left Coordinates of the first vertex of the edge
	   *  @param right Coordinates of the second vertex of the edge
	   */
	  void mark_insertion (const std::pair<double,double> & left, const std::pair<double,double> & right);

	  /*! @brief Mark edge to be removed
	   *  @param it Iterator to the face we are analyzing
	   *  @param i Index of the edge that has to be removed (Local index in face numeration)
	   */
	  void mark_removal (face_iterator_t & it, int const & i);
    
	  /*! @brief Perform one bisection on the triangle
	   *  @param it Iterator to the face to be bisected
	   *  @param ref Index of reference edge to be bisected
	   */
	  void one_bisection(face_iterator_t & it, int const ref);
	  
	  /*! @brief Perform two bisections on the triangle
	   *  @param it Iterator to the face to be bisected
	   *  @param ref Index of reference edge to be bisected
	   *  @param idx Index of the second edge to be refined after reference edge
	   */
	  void two_bisection(face_iterator_t & it, int const ref, int const idx);
	  
	  /*! @brief Perform three bisections on the triangle 
	   *  @param it Iterator to the face to be bisected
	   *  @param ref Index of reference edge to be bisected
	   */
	  void triple_bisection(face_iterator_t & it, int const ref);
       
  protected:
	    /*! @brief Problem */
	    problem const & p_;
	    
	    /*! @brief Triangulation */
	    triangulation & cdt_;	
	    
	    /*! @brief Vector that stores if a triangle has to be refined */
	    std::vector<bool> & TbR_;
	    
	    /*! @brief Vector of elements*/
	    std::vector<std::list<std::pair<double,double> > > & to_be_inserted_;

	    /*! @brief Vector of the edges to be removed*/
	    std::vector<std::pair<bool,std::list<int> > > & to_be_removed_;
	    
	    /*! @brief Vector that interpolates the solution in the new points of the triangulation */
	    std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & new_points_;
	    
	    /*! @brief Boolean variable to identify if the refinement is at the begininng or during evolution */
	    const bool init_;

  private:

};

 
void RefineMethod::mark_insertion (const std::pair<double,double> & left, const std::pair<double,double> & right)
{
    // Extract the information about the vertices of the edge
    std::list<std::pair<double,double> > edge;
    edge.push_back(left);
    edge.push_back(right);
    
    // Add coordinates of the edge in the vector of edges to be inserted
    to_be_inserted_.push_back(edge);
}
  
  
void RefineMethod::mark_removal (face_iterator_t & it, int const & i)
{
    // Extract the list of edges that have to be removed from the triangle
    std::list<int> edge_list;
    edge_list = to_be_removed_[it->i()].second;
    edge_list.push_back(i);

    // Add edges to be removed to the list 
    // Flag triangles that have edges marked for removal
    // Store the list of edges that have to be removed 
    to_be_removed_[it->i()] = std::make_pair(true,edge_list);
}


void RefineMethod::one_bisection(face_iterator_t & it, int const ref)
{
    gas_assert(ref < 3);    
  
    // Get the coordinates of the vertices of the triangle
    Eigen::Vector2d P0(it->x(0), it->y(0));
    Eigen::Vector2d P1(it->x(1), it->y(1));
    Eigen::Vector2d P2(it->x(2), it->y(2));
    
    double Px, Py, F1x, F1y, F2x, F2y;
    double mEx, mEy;
    
    // Index of the edge that will be refined
    int old=it->reference_edge();
    
    // Initialize the points according to the reference edge
    switch (ref)    
    {
      case 0:
      {
	Px = P0(0);
	Py = P0(1);
	F1x = P1(0);
	F1y = P1(1);
	F2x = P2(0);
	F2y = P2(1);	
      }
      break;
      
      case 1:
      {
	Px = P1(0);
	Py = P1(1);
	F1x = P2(0);
	F1y = P2(1);
	F2x = P0(0);
	F2y = P0(1);	
      }
      break;
      
      case 2:
      {
	Px = P2(0);
	Py = P2(1);
	F1x = P0(0);
	F1y = P0(1);
	F2x = P1(0);
	F2y = P1(1);
      }
      break;
    }
    
    // Compute middle point of the reference edge
    mEx = (F1x+F2x)/2;
    mEy = (F1y+F2y)/2;    
    
    // Remove former reference edge F1-F2
    this->mark_removal(it,old);    

    // Add new edges to the vector of edges to be inserted //
    
    // P-m
    if(Px <= mEx)
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(mEx,mEy));
    }
    else
    {
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(Px,Py));      
    }
    
    // F1-m-F2
    if(F1x <= F2x)
    {
      this->mark_insertion(std::make_pair(F1x,F1y),std::make_pair(mEx,mEy));
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(F2x,F2y));
    }
    else
    {
      this->mark_insertion(std::make_pair(F2x,F2y),std::make_pair(mEx,mEy));
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(F1x,F1y));      
    }
    
    // F2-P
    if(F2x <= Px)
    {
      this->mark_insertion(std::make_pair(F2x,F2y),std::make_pair(Px,Py));
    }
    else
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(F2x,F2y));      
    }
    
    // P-F1
    if(Px <= F1x)
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(F1x,F1y));    
    }
    else
    {
      this->mark_insertion(std::make_pair(F1x,F1y),std::make_pair(Px,Py));          
    }

    // Remove former triangle from the vector that stores bisected faces
    cdt_.remove_biFace(it);
    
    // Add newly added triangles to the vector that stores bisected faces
    // First vertex is the last added ones
    std::list<std::pair<double,double> > coord;

    coord.push_back(std::make_pair(mEx,mEy));        
    coord.push_back(std::make_pair(Px,Py));
    coord.push_back(std::make_pair(F1x,F1y));    
    
    cdt_.insert_biFace(coord);
    
    coord.clear();

    coord.push_back(std::make_pair(mEx,mEy));            
    coord.push_back(std::make_pair(F2x,F2y));    
    coord.push_back(std::make_pair(Px,Py));        
    
    cdt_.insert_biFace(coord);
    
    // If refinement during time loop, interpolate current and previous solution on new points
    if(!init_)
    {
      p_.interpolate(it,mEx,mEy,new_points_,1);
    }
    
  }


void RefineMethod::two_bisection(face_iterator_t & it, int const ref, int const idx)
{
    gas_assert(ref < 3);    
    gas_assert(idx < 3);
    
    // Get the coordinates of the vertices of the triangle
    Eigen::Vector2d P0(it->x(0), it->y(0));
    Eigen::Vector2d P1(it->x(1), it->y(1));
    Eigen::Vector2d P2(it->x(2), it->y(2));
    
    double Px, Py, Fx, Fy, Bx, By;
    double mEx, mEy, nEx, nEy;
    
    // Index of the edges that will be refined    
    int* old = new int[2];

    // Boolean value that indicates whether the pending node after first bisection is in 
    // the right triangle or in the left one
    bool right;
    
    // Initialize the points according to the reference edge and the second edge to be refined    
    switch (ref)        
    {
      case 0:
      {
	Px = P0(0);
	Py = P0(1);

	switch (idx)
	{	  
	  case 1:
	  {
	    Fx = P1(0);
	    Fy = P1(1);
	    Bx = P2(0);
	    By = P2(1);	
	    
	    right = true;
	    
	    old[0] = 0;
	    old[1] = 1;
	  }
	  break;
	  
	  case 2:
	  {
	    Fx = P2(0);
	    Fy = P2(1);		    
	    Bx = P1(0);
	    By = P1(1);
	    
	    right = false;
	    
	    old[0] = 0;
	    old[1] = 2;
	  }
	  break;
	}
      }
      break;
      
      case 1:
      {
	Px = P1(0);
	Py = P1(1);

	switch (idx)
	{
	  case 0:
	  {
	    Fx = P0(0);
	    Fy = P0(1);	    
	    Bx = P2(0);
	    By = P2(1);		    

	    right = false;
	    
	    old[0] = 1;
	    old[1] = 0;
	  }
	  break;
	  
	  case 2:
	  {
	    Fx = P2(0);
	    Fy = P2(1);
	    Bx = P0(0);
	    By = P0(1);	
	    
	    right = true;
	    
	    old[0] = 1;
	    old[1] = 2;
	  }
	  break;
	}
      }
      break;
      
      case 2:
      {
	Px = P2(0);
	Py = P2(1);

	switch (idx)
	{
	  case 0:
	  {
	    Fx = P0(0);
	    Fy = P0(1);	    
	    Bx = P1(0);
	    By = P1(1);		    
	    
	    right = true;
	    
	    old[0] = 2;
	    old[1] = 0;
	  }
	  break;
	  
	  case 1:
	  {
	    Fx = P1(0);
	    Fy = P1(1);
	    Bx = P0(0);
	    By = P0(1);	
	    
	    right = false;
	    
	    old[0] = 2;
	    old[1] = 1;
	  }
	  break;
	}
      }
      break;
    }

    // Compute middle point of the reference edge
    mEx = (Fx+Bx)/2;
    mEy = (Fy+By)/2;

    // Compute middle point of the second edge to be refined    
    nEx = (Px+Bx)/2;
    nEy = (Py+By)/2;    
   
    // Remove former reference edge F-B and second edge B-P
    this->mark_removal(it,old[0]);        
    this->mark_removal(it,old[1]);        

    // Add new edges to the vector of edges to be inserted //    
         
    // P-m
    if(Px <= mEx)
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(mEx,mEy));
    }
    else
    {
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(Px,Py));      
    }
    
    // m-n
    if(mEx <= nEx)
    {
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(nEx,nEy));
    }
    else
    {
      this->mark_insertion(std::make_pair(nEx,nEy),std::make_pair(mEx,mEy));      
    }
    
    // P-F
    if(Px <= Fx)
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(Fx,Fy));        
    }
    else
    {
      this->mark_insertion(std::make_pair(Fx,Fy),std::make_pair(Px,Py));              
    }
    
    // F-m-B
    if(Fx <= Bx)
    {
      this->mark_insertion(std::make_pair(Fx,Fy),std::make_pair(mEx,mEy));
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(Bx,By));    
    }
    else
    {
      this->mark_insertion(std::make_pair(Bx,By),std::make_pair(mEx,mEy));
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(Fx,Fy));          
    }
    
    // B-n-P
    if(Bx <= Px)
    {
      this->mark_insertion(std::make_pair(Bx,By),std::make_pair(nEx,nEy));
      this->mark_insertion(std::make_pair(nEx,nEy),std::make_pair(Px,Py));        
    }
    else
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(nEx,nEy));
      this->mark_insertion(std::make_pair(nEx,nEy),std::make_pair(Bx,By));              
    }
    
    // Remove former triangle from the vector that stores bisected faces
    cdt_.remove_biFace(it);
    
    // If right subtriangle
    if(right)
    {
	// Add newly added triangles to the vector that stores bisected faces
	// First vertex is the last added one       
	std::list<std::pair<double,double> > coord;

	coord.push_back(std::make_pair(mEx,mEy));          
	coord.push_back(std::make_pair(Px,Py));
	coord.push_back(std::make_pair(Fx,Fy));    

	cdt_.insert_biFace(coord);
	
	coord.clear();

	coord.push_back(std::make_pair(nEx,nEy));          
	coord.push_back(std::make_pair(Px,Py));        
	coord.push_back(std::make_pair(mEx,mEy));            

	cdt_.insert_biFace(coord);
	
	coord.clear();

	coord.push_back(std::make_pair(nEx,nEy));        
	coord.push_back(std::make_pair(mEx,mEy));            
	coord.push_back(std::make_pair(Bx,By));    

	cdt_.insert_biFace(coord);
    }
    else // left subtriangle
    {
	// Add newly added triangles to the vector that stores bisected faces
	// First vertex is the last added one       
	std::list<std::pair<double,double> > coord;
	
	coord.push_back(std::make_pair(mEx,mEy));          
	coord.push_back(std::make_pair(Fx,Fy));    
	coord.push_back(std::make_pair(Px,Py));
	
	cdt_.insert_biFace(coord);
	
	coord.clear();

	coord.push_back(std::make_pair(nEx,nEy));          
	coord.push_back(std::make_pair(mEx,mEy));            
	coord.push_back(std::make_pair(Px,Py));        
	
	cdt_.insert_biFace(coord);
	
	coord.clear();

	coord.push_back(std::make_pair(nEx,nEy));        
	coord.push_back(std::make_pair(Bx,By));          
	coord.push_back(std::make_pair(mEx,mEy));            

	cdt_.insert_biFace(coord);
    }
    
    delete old;
    
    // If refinement during time loop, interpolate current and previous solution on new points    
    if(!init_)
    {
      p_.interpolate(it,mEx,mEy,new_points_,1);
      p_.interpolate(it,nEx,nEy,new_points_,1);			
    }
  
} 


void RefineMethod::triple_bisection(face_iterator_t & it, int const ref)
{
    gas_assert(ref < 3);    
    
    // Get the coordinates of the vertices of the triangle
    Eigen::Vector2d P0(it->x(0), it->y(0));
    Eigen::Vector2d P1(it->x(1), it->y(1));
    Eigen::Vector2d P2(it->x(2), it->y(2));
    
    double Px, Py;
    double Fx, Fy; // looking at P, F is on the left
    double Bx, By; // looking at P, B is on the right
    
    double mEx, mEy; // F-B
    double nEx, nEy; // P-B
    double oEx, oEy; // P-F 
    
    // Remove all edges of previous triangle
    for(int j=0; j<3; j++) { this->mark_removal(it,j); }

    // Initialize the points according to the reference edge and the second edge to be refined    
    switch (ref)    
    {
      case 0:
      {
	Px = P0(0);
	Py = P0(1);
	Fx = P1(0);
	Fy = P1(1);
	Bx = P2(0);
	By = P2(1);	
      }
      break;
      
      case 1:
      {
	Px = P1(0);
	Py = P1(1);
	Fx = P2(0);
	Fy = P2(1);
	Bx = P0(0);
	By = P0(1);	
      }
      break;
      
      case 2:
      {
	Px = P2(0);
	Py = P2(1);
	Fx = P0(0);
	Fy = P0(1);	    
	Bx = P1(0);
	By = P1(1);		    
      }
      break;
    }

    // Compute middle points of all edges
    mEx = (Fx+Bx)/2;
    mEy = (Fy+By)/2;

    nEx = (Px+Bx)/2;
    nEy = (Py+By)/2;    
   
    oEx = (Px+Fx)/2;
    oEy = (Py+Fy)/2;    
    
    // Add new edges to the vector of edges to be inserted //    
       
    // P-m
    if(Px <= mEx)
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(mEx,mEy));
    }
    else
    {
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(Px,Py));      
    }
    
    // m-n
    if(mEx <= nEx)
    {
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(nEx,nEy));
    }
    else
    {
      this->mark_insertion(std::make_pair(nEx,nEy),std::make_pair(mEx,mEy));      
    }
    
    // m-o
    if(mEx <= oEx)
    {
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(oEx,oEy));
    }
    else
    {
      this->mark_insertion(std::make_pair(oEx,oEy),std::make_pair(mEx,mEy));      
    }
    
    // F-m-B
    if(Fx <= Bx)
    {
      this->mark_insertion(std::make_pair(Fx,Fy),std::make_pair(mEx,mEy));
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(Bx,By));    
    }
    else
    {
      this->mark_insertion(std::make_pair(Bx,By),std::make_pair(mEx,mEy));
      this->mark_insertion(std::make_pair(mEx,mEy),std::make_pair(Fx,Fy));          
    }
    
    // B-n-P
    if(Bx <= Px)
    {
      this->mark_insertion(std::make_pair(Bx,By),std::make_pair(nEx,nEy));
      this->mark_insertion(std::make_pair(nEx,nEy),std::make_pair(Px,Py));        
    }
    else
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(nEx,nEy));
      this->mark_insertion(std::make_pair(nEx,nEy),std::make_pair(Bx,By));              
    }
    
    // P-o-F
    if(Fx <= Px)
    {
      this->mark_insertion(std::make_pair(Fx,Fy),std::make_pair(oEx,oEy));
      this->mark_insertion(std::make_pair(oEx,oEy),std::make_pair(Px,Py));        
    }
    else
    {
      this->mark_insertion(std::make_pair(Px,Py),std::make_pair(oEx,oEy));
      this->mark_insertion(std::make_pair(oEx,oEy),std::make_pair(Fx,Fy));              
    }
   
    // Remove former triangle from the vector that stores bisected faces
    cdt_.remove_biFace(it);
    
    // Add newly added triangles to the vector that stores bisected faces
    // First vertex is the last added one       
    std::list<std::pair<double,double> > coord;

    coord.push_back(std::make_pair(oEx,oEy));          
    coord.push_back(std::make_pair(mEx,mEy));          
    coord.push_back(std::make_pair(Px,Py));

    cdt_.insert_biFace(coord);
    
    coord.clear();

    coord.push_back(std::make_pair(oEx,oEy));          
    coord.push_back(std::make_pair(Fx,Fy));
    coord.push_back(std::make_pair(mEx,mEy));          

    cdt_.insert_biFace(coord);
    
    coord.clear();

    coord.push_back(std::make_pair(nEx,nEy));          
    coord.push_back(std::make_pair(Px,Py));        
    coord.push_back(std::make_pair(mEx,mEy));            

    cdt_.insert_biFace(coord);
    
    coord.clear();

    coord.push_back(std::make_pair(nEx,nEy));        
    coord.push_back(std::make_pair(mEx,mEy));            
    coord.push_back(std::make_pair(Bx,By));    

    cdt_.insert_biFace(coord);

    // If refinement during time loop, interpolate current and previous solution on new points    
    if(!init_)
    {
      p_.interpolate(it,mEx,mEy,new_points_,1);
      p_.interpolate(it,nEx,nEy,new_points_,1);
      p_.interpolate(it,oEx,oEy,new_points_,1);			      
    }
  
}


/*! @brief Regular division refinement method */
class RegularDivisionRefinement : public RefineMethod
{

  public:
	/*! @brief Constructor 
	  *  @param p Problem
	  *  @param cdt Triangulation
	  *  @param TbR Vector that keeps the information about whether or not elements have to be refined
	  *  @param to_be_inserted Vector that stores the elements to be inserted
	  *  @param to_be_removed Vector that stores the edges to be removed
	  *  @param new_points Vector that stores the values of the solution interpolated on new points in the triangulation
	  *  @param init Boolean variable to flag if refinement is performed at the beginning or duriong evolution loop
	  */
	  inline RegularDivisionRefinement ( problem const & p, triangulation & cdt, std::vector<bool> & TbR, std::vector<std::list<std::pair<double,double> > > & to_be_inserted, std::vector<std::pair<bool,std::list<int> > > & to_be_removed, std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & new_points, const bool init ): 
	  RefineMethod(p, cdt, TbR, to_be_inserted, to_be_removed, new_points, init) {};    
    
    	  /*! @brief Virtual destructor */
	  virtual ~RegularDivisionRefinement() {};

	  /*! @brief Method to apply the refinement */
	  void apply();
    
	  /*! @brief Method to update labelas and reference edges */
	  void update();	  
    
	  /*! @brief Performs red refinement over the triangle
	   *  @param it Iterator to the face to be refined
	   */
	  void red_refinement(face_iterator_t & it);
	  
	  /*! @brief Performs green bisection over the face 
	   *  @param it Iterator to the face to be bisected
	   */	  
	  inline void green_refinement(face_iterator_t & it) { this->one_bisection(it, it->reference_edge()); };
	  
	  /*! @brief Performs blue bisection over the face
	   *  @param it Iterator to the face to be bisected
	   *  @param idx Index of the second edge to be refined
	   */	  
	  inline void blue_refinement(face_iterator_t & it, int const idx) { this->two_bisection(it, it->reference_edge(), idx); };

};


void RegularDivisionRefinement::apply()
  {
    std::cout << "Begin refinement" << std::endl;
    
    // For every edge in every face it stores a flag if the corresponding edge has a pending node or 
    // it already has been refined
    int hang_p_nb=0; // number of hanging nodes
    std::vector<std::vector<int> > hanging(cdt_.faces(),std::vector<int>(3,0));

    // Counts triangles that still needs to be checked
    int pending_triangles=0;
    std::vector<int> triangle_status(cdt_.faces(),0);
    
    // Number of refined triangles
    int ref_counter(0);
    
    for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
    {
      // If refinement during time loop, interpolate current and previous solution on existing points      
      if(!init_)
      {
	gas_rangeu(k,3)
	{
	    p_.interpolate(it,it->x(k),it->y(k),new_points_,0);		    
	}
      }      
      
      // If the face has to be refined
      // Mark triangles to be refined
      if (TbR_[it->i()]) 
      {
	  ref_counter++;
	  
	  // Count number of hanging nodes
	  int h_tot=0;
	  gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }

	  // If the triangle has not been refined yet, call red refinement routine
	  if(h_tot>=0)
	  {
	    this->red_refinement(it);
	    hang_p_nb -= h_tot;
	    
	    // Triangle refined
	    triangle_status[it->i()]=-1;

	    // Set neighbors to be checked for the presence of hanging nodes
	    // If so, set the neighbors to be further checked
	    gas_rangeu(k,3)
	    {
	      hanging[it->i()][k]=-1;
	      
	      if(hanging[it->neighbor(k)][it->loc_mirror_vertex(k)] == 0)
	      {
		hanging[it->neighbor(k)][it->loc_mirror_vertex(k)]=1;
		hang_p_nb++;		  		
		triangle_status[it->neighbor(k)]=1;
	      }
	    }
	  }
	  
      } // To be Refined
    } // end face loop
      
    std::cout << "Number of refined triangles: " << ref_counter << std::endl;

    // Checking conformity for the mesh
    // Analyze neighbors to mark edges for refinement to guarantee conformity
    do {  
	  int red_counter=0;
	  
	  // Counts faces with three hanging nodes => Marked for red refinement
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    int h_tot=0;
	    gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }
	    
	    if(h_tot==0) { triangle_status[it->i()]=0; }
	    // Triangles either already refined or with three hanging nodes
	    else if( (h_tot==3) || (h_tot==-3) )
	    {
	      triangle_status[it->i()]=-1;
	      red_counter++;
	    }
	  }

	  // Analyze faces for green and blue refinement by checking if reference edge is marked
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    int h_tot=0;
	    gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }
	    
	    // Triangle with one hanging node
	    if(h_tot == 1)
	    {
	      // If current reference is not marked, mark it
	      if(hanging[it->i()][it->reference_edge()] == 0) 
	      {
		hanging[it->i()][it->reference_edge()]=1;		  
		hang_p_nb++;
		
		triangle_status[it->i()]=0;

		// If the neighbor is not marked for red refinement, mark the edge
		if(triangle_status[it->neighbor(it->reference_edge())] != -1)
		{		
		  if(hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())] == 0)
		  {
		    hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())]=1;		  
		    hang_p_nb++;		
		    triangle_status[it->neighbor(it->reference_edge())]=1;
		  }
		}
	      }
	      else { triangle_status[it->i()]=0; }
	    } 
	    else if(h_tot == 2) // Triangle with two hanging nodes
	    {
	      // If current reference is not marked, mark it
	      if(hanging[it->i()][it->reference_edge()] == 0) 
	      {
		hanging[it->i()][it->reference_edge()]=1;		  
		hang_p_nb++;
		
		triangle_status[it->i()]=-1;
		red_counter++;
		
		// If the neighbor is not marked for red refinement, mark the edge
		if(triangle_status[it->neighbor(it->reference_edge())] != -1)
		{		
		  if(hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())] == 0)
		  {
		    hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())]=1;		  
		    hang_p_nb++;
		    triangle_status[it->neighbor(it->reference_edge())]=1;		    
		  }
		}
	      }
	      else { triangle_status[it->i()]=0; }	      
	    }
	  }
	  
	  // Update number of pending triangles
	  pending_triangles = red_counter;
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    pending_triangles += triangle_status[it->i()];
	  }
	  
    } while (pending_triangles > 0);
	
    std::cout << "Total hanging nodes: " << hang_p_nb << std::endl;
   
    // Number of compatibility triangles
    int comp_counter(0);
    
    for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
    {
  	int h_tot=0;
	int scenario=0;
	
	gas_rangeu(k,3)	{ h_tot+=hanging[it->i()][k]; }
	
	if(h_tot==1) { scenario=1; }
	else if(h_tot==2) { scenario=2; }
	else if(h_tot==3) { scenario=3; }	
	else { scenario=0; }

	// Select refinement method 
	switch(scenario)
	{
	  case 0: // Nothing to do
	  {
	    
	  }
	  break;
	  
	  case 1: // Green refinement
	  {
	    comp_counter++;
	    
	    this->green_refinement(it);
	    hang_p_nb -= h_tot;	    
	    
	    gas_rangeu(k,3) { hanging[it->i()][k]=-1; } // tombstone
	  }
	  break;
	  
	  case 2: // Blue refinement
	  {
	    comp_counter++;

	    int idx;
	    idx = idx;
	    gas_rangeu(k,3)
	    { 
	      if((hanging[it->i()][k]==1) && (k != it->reference_edge())) { idx=k; }
	    }
	    this->blue_refinement(it,idx);
	    hang_p_nb -= h_tot;
	    
	    gas_rangeu(k,3) { hanging[it->i()][k]=-1; } // tombstone
	  }
	  break;
	  
	  case 3: // Red refinement
	  {
	    comp_counter++;	    

	    this->red_refinement(it);
	    hang_p_nb -= h_tot;
	    
	    gas_rangeu(k,3)
	    {
	      hanging[it->i()][k]=-1;
	    }
	    
	  }
	  break;
	}
	
    }
  
    std::cout << "Checking conformity - Triangles to impose compatibility: " << comp_counter << std::endl;
    
    std::cout << "End refinement" << std::endl;

}


void RegularDivisionRefinement::red_refinement(face_iterator_t & it)
{
    // Get the information about the coordinates of the vertices
    Eigen::Vector2d P0(it->x(0), it->y(0));
    Eigen::Vector2d P1(it->x(1), it->y(1));
    Eigen::Vector2d P2(it->x(2), it->y(2));		

    // Middle points of the edges
    double const mE0x((P1(0)+P2(0))/2);
    double const mE0y((P1(1)+P2(1))/2);
    
    double const mE1x((P2(0)+P0(0))/2);
    double const mE1y((P2(1)+P0(1))/2);	
    
    double const mE2x((P0(0)+P1(0))/2);
    double const mE2y((P0(1)+P1(1))/2);
     
    // Remove P0-P1, P1-P2, P2-P0
    for(int j=0; j<3; j++) { this->mark_removal(it,j); }
 
    // Add new edges to the vector of edges to be inserted //    
 
    // m0-m1
    if(mE0x <= mE1x)
    {
      this->mark_insertion(std::make_pair(mE0x,mE0y),std::make_pair(mE1x,mE1y));
    }
    else
    {
      this->mark_insertion(std::make_pair(mE1x,mE1y),std::make_pair(mE0x,mE0y));      
    }
    
    // m1-m2
    if(mE1x <= mE2x)
    {
      this->mark_insertion(std::make_pair(mE1x,mE1y),std::make_pair(mE2x,mE2y));			
    }
    else
    {
      this->mark_insertion(std::make_pair(mE2x,mE2y),std::make_pair(mE1x,mE1y));			      
    }
    
    // m2-m0
    if(mE2x <= mE0x)
    {
      this->mark_insertion(std::make_pair(mE2x,mE2y),std::make_pair(mE0x,mE0y));			    
    }
    else
    {
      this->mark_insertion(std::make_pair(mE0x,mE0y),std::make_pair(mE2x,mE2y));			          
    }
    
    // P0-m2-P1
    if(P0(0) <= P1(0))
    {
      this->mark_insertion(std::make_pair(P0(0),P0(1)),std::make_pair(mE2x,mE2y));			        
      this->mark_insertion(std::make_pair(mE2x,mE2y),std::make_pair(P1(0),P1(1)));
    }
    else
    {
      this->mark_insertion(std::make_pair(P1(0),P1(1)),std::make_pair(mE2x,mE2y));			        
      this->mark_insertion(std::make_pair(mE2x,mE2y),std::make_pair(P0(0),P0(1)));      
    }
    
    // P1-m0-P2    
    if(P1(0) <= P2(0))
    {
      this->mark_insertion(std::make_pair(P1(0),P1(1)),std::make_pair(mE0x,mE0y));
      this->mark_insertion(std::make_pair(mE0x,mE0y),std::make_pair(P2(0),P2(1)));
    }
    else
    {
      this->mark_insertion(std::make_pair(P2(0),P2(1)),std::make_pair(mE0x,mE0y));
      this->mark_insertion(std::make_pair(mE0x,mE0y),std::make_pair(P1(0),P1(1)));      
    }
    
    // P2-m1-P0    
    if(P2(0) <= P0(0))
    {
      this->mark_insertion(std::make_pair(P2(0),P2(1)),std::make_pair(mE1x,mE1y));			
      this->mark_insertion(std::make_pair(mE1x,mE1y),std::make_pair(P0(0),P0(1)));			    
    }
    else
    {
      this->mark_insertion(std::make_pair(P0(0),P0(1)),std::make_pair(mE1x,mE1y));			
      this->mark_insertion(std::make_pair(mE1x,mE1y),std::make_pair(P2(0),P2(1)));			          
    }
    
    // Remove former triangle from the vector that stores red refined faces
    cdt_.remove_redFace(it);    
        
    // Add newly added triangles to the vector that stores red refined faces
    // First vertex is the last added ones   
    std::list<std::pair<double,double> > coord;

    coord.push_back(std::make_pair(mE0x,mE0y));    
    coord.push_back(std::make_pair(mE1x,mE1y));            
    coord.push_back(std::make_pair(mE2x,mE2y));        
    
    cdt_.insert_redFace(coord);
    
    coord.clear();
    
    coord.push_back(std::make_pair(P0(0),P0(1)));    
    coord.push_back(std::make_pair(mE2x,mE2y));        
    coord.push_back(std::make_pair(mE1x,mE1y));        
    
    cdt_.insert_redFace(coord);
    
    coord.clear();

    coord.push_back(std::make_pair(mE2x,mE2y));            
    coord.push_back(std::make_pair(P1(0),P1(1)));    
    coord.push_back(std::make_pair(mE0x,mE0y));        
    
    cdt_.insert_redFace(coord);
    
    coord.clear();

    coord.push_back(std::make_pair(mE1x,mE1y));            
    coord.push_back(std::make_pair(mE0x,mE0y));            
    coord.push_back(std::make_pair(P2(0),P2(1)));    
    
    cdt_.insert_redFace(coord);
    
    coord.clear();
    
    // If refinement during time loop, interpolate current and previous solution on new points    
    if(!init_)
    {
      p_.interpolate(it,mE0x,mE0y,new_points_,1);
      p_.interpolate(it,mE1x,mE1y,new_points_,1);			
      p_.interpolate(it,mE2x,mE2y,new_points_,1);			      
    }
    
}


void RegularDivisionRefinement::update()
{
    // Update labels within triangulation
    cdt_.update_dataRD();    
}


/*! @brief Longest edge bisection refinement method */
class LongestEdgeRefinement : public RefineMethod
{

  public:
	/*! @brief Constructor 
	  *  @param p Problem
	  *  @param cdt Triangulation
	  *  @param TbR Vector that keeps the information about whether or not elements have to be refined
	  *  @param to_be_inserted Vector that stores the elements to be inserted
	  *  @param to_be_removed Vector that stores the edges to be removed
	  *  @param new_points Vector that stores the values of the solution interpolated on new points in the triangulation
	  *  @param init Boolean variable to flag if refinement is performed at the beginning or duriong evolution loop
	  */
	  inline LongestEdgeRefinement ( problem const & p, triangulation & cdt, std::vector<bool> & TbR, std::vector<std::list<std::pair<double,double> > > & to_be_inserted, std::vector<std::pair<bool,std::list<int> > > & to_be_removed, std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & new_points, const bool init ): 
	  RefineMethod(p, cdt, TbR, to_be_inserted, to_be_removed, new_points, init) {};    
    
    	  /*! @brief Virtual destructor */
	  virtual ~LongestEdgeRefinement() {};

	  /*! @brief Method to apply the refinement */
	  void apply();
    
	  /*! @brief Method to update labels and reference edges */
	  void update();	  
    
};  


void LongestEdgeRefinement::apply()
{
    std::cout << "Begin refinement" << std::endl;
    
    // For every edge in every face it stores a flag if the corresponding edge has a pending node or 
    // it already has been refined
    int hang_p_nb=0; // number of hanging nodes
    std::vector<std::vector<int> > hanging(cdt_.faces(),std::vector<int>(3,0));

    // Counts triangles that still needs to be checked
    int pending_triangles=0;
    std::vector<int> triangle_status(cdt_.faces(),0);
    
    // Number of refined triangles
    int ref_counter(0);
    
    for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
    {
      // If refinement during time loop, interpolate current and previous solution on existing points      
      if(!init_)
      {
	gas_rangeu(k,3)
	{
	    p_.interpolate(it,it->x(k),it->y(k),new_points_,0);		    
	}
      }      
      
      // If the face has to be refined
      // Mark triangles to be refined
      if (TbR_[it->i()]) 
      {
	  ref_counter++;
	  // Count number of hanging nodes
	  int h_tot=0;
	  gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }

	  // If the triangle has not been refined yet, call triple bisection routine
	  if(h_tot>=0)
	  {
	    this->triple_bisection(it,it->reference_edge());

	    hang_p_nb -= h_tot;
	    
	    // Triangle refined
	    triangle_status[it->i()]=-1;
	    
	    // Set neighbors to be checked for the presence of hanging nodes
	    // If so, set the neighbors to be further checked
	    gas_rangeu(k,3)
	    {
	      hanging[it->i()][k]=-1;
	      
	      if(hanging[it->neighbor(k)][it->loc_mirror_vertex(k)] == 0)
	      {
		hanging[it->neighbor(k)][it->loc_mirror_vertex(k)]=1;
		hang_p_nb++;		  		
		triangle_status[it->neighbor(k)]=1;
	      }
	    }
	  }

      } // To be Refined
    } // end face loop
        
    std::cout << "Number of refined triangles: " << ref_counter << std::endl;
    
    // Checking conformity for the mesh
    // Analyze neighbors to mark edges for refinement to guarantee conformity
    do {  
	  int red_counter=0;
	  
	  // Counts faces with three hanging nodes => Marked for triple bisection
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    int h_tot=0;
	    gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }
	    
	    if(h_tot==0) { triangle_status[it->i()]=0; }
	    // Triangles either already refined or with three hanging nodes
	    else if( (h_tot==3) || (h_tot==-3) )
	    {
	      triangle_status[it->i()]=-1;
	      red_counter++;
	    }
	  }
      
	  // Analyze faces for green and blue refinement by checking if reference edge is marked
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    int h_tot=0;
	    gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }
	    
	    // Triangle with one hanging node
	    if(h_tot == 1)
	    {
	      // If current reference is not marked, mark it
	      if(hanging[it->i()][it->reference_edge()] == 0) 
	      {
		hanging[it->i()][it->reference_edge()]=1;		  
		hang_p_nb++;
		
		triangle_status[it->i()]=0;
		
		// If the neighbor is not marked for triple bisection, mark the edge
		if(triangle_status[it->neighbor(it->reference_edge())] != -1)
		{		
		  if(hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())] == 0)
		  {
		    hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())]=1;		  
		    hang_p_nb++;		
		    triangle_status[it->neighbor(it->reference_edge())]=1;
		  }
		}
	      }
	      else { triangle_status[it->i()]=0; }
	    } 
	    else if(h_tot == 2) // Triangle with two hanging nodes
	    {
	      // If current reference is not marked, mark it	      
	      if(hanging[it->i()][it->reference_edge()] == 0) 
	      {
		hanging[it->i()][it->reference_edge()]=1;		  
		hang_p_nb++;
		
		triangle_status[it->i()]=-1;
		red_counter++;

		// If the neighbor is not marked for triple bisection, mark the edge		
		if(triangle_status[it->neighbor(it->reference_edge())] != -1)
		{		
		  if(hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())] == 0)
		  {
		    hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())]=1;		  
		    hang_p_nb++;
		    triangle_status[it->neighbor(it->reference_edge())]=1;		    
		  }
		}
	      }
	      else { triangle_status[it->i()]=0; }	      
	    }
	  }
	  	  
	  // Update the number of triangles to be checked for conformity
	  pending_triangles = red_counter;
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    pending_triangles += triangle_status[it->i()];
	  }
	  	  
    } while (pending_triangles > 0);


    std::cout << "Total hanging nodes: " << hang_p_nb << std::endl;
      
    // Number of compatibility triangles
    int comp_counter(0);
  
    for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
    {	
	int h_tot=0;
	int scenario=0;
	
	gas_rangeu(k,3)	{ h_tot+=hanging[it->i()][k]; }
	
	if(h_tot==1) { scenario=1; }
	else if(h_tot==2) { scenario=2; }
	else if(h_tot==3) { scenario=3; }	
	else { scenario=0; }

	// Choose refinement strategy based on the number of hanging nodes
	switch(scenario)
	{
	  case 0: // Nothing to do
	  {

	  }
	  break;
	  
	  case 1: // Green bisection
	  {
	    comp_counter++;

	    this->one_bisection(it,it->reference_edge());
	    hang_p_nb -= h_tot;	    
	    
	    gas_rangeu(k,3) { hanging[it->i()][k]=-1; } // tombstone
	  }
	  break;
	  
	  case 2: // Blue bisection
	  {
	    comp_counter++;

    	    int idx;
	    gas_rangeu(k,3)
	    { 
	      if((hanging[it->i()][k]==1) && (k != it->reference_edge())) { idx=k; }
	    }
	    this->two_bisection(it,it->reference_edge(),idx);
	    hang_p_nb -= h_tot;
	    
	    gas_rangeu(k,3) { hanging[it->i()][k]=-1; } // tombstone
	  }
	  break;
	  
	  case 3: // Triple bisection
	  {
	    comp_counter++;

	    this->triple_bisection(it,it->reference_edge());
	    hang_p_nb -= h_tot;
	    
	    gas_rangeu(k,3)
	    {
	      hanging[it->i()][k]=-1;
	    }	    
	  }
	  break;
	}
      }
  
      std::cout << "Checking conformity - Triangles to impose compatibility: " << comp_counter << std::endl;
      
      std::cout << "End refinement" << std::endl;

}


void LongestEdgeRefinement::update()
{
    // Update labels within triangulation
    cdt_.update_dataLE();
}


/*! @brief Newest Vertex bisection refinement method */
class NewestVertexRefinement : public RefineMethod
{

  public:
	/*! @brief Constructor 
	  *  @param p Problem
	  *  @param cdt Triangulation
	  *  @param TbR Vector that keeps the information about whether or not elements have to be refined
	  *  @param to_be_inserted Vector that stores the elements to be inserted
	  *  @param to_be_removed Vector that stores the edges to be removed
	  *  @param new_points Vector that stores the values of the solution interpolated on new points in the triangulation
	  *  @param init Boolean variable to flag if refinement is performed at the beginning or duriong evolution loop
	  */
	  inline NewestVertexRefinement ( problem const & p, triangulation & cdt, std::vector<bool> & TbR, std::vector<std::list<std::pair<double,double> > > & to_be_inserted, std::vector<std::pair<bool,std::list<int> > > & to_be_removed, std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & new_points, const bool init ): 
	  RefineMethod(p, cdt, TbR, to_be_inserted, to_be_removed, new_points, init) {};    
    
    	  /*! @brief Virtual destructor */
	  virtual ~NewestVertexRefinement() {};

	  /*! @brief Method to apply the refinement */
	  void apply();
    
	  /*! @brief Method to update labelas and reference edges */
	  void update();	  
	  
	  /*! @brief Stores the elements that have been bisected */
	  void store_previousT(face_iterator_t & it);
        
};

  
void NewestVertexRefinement::apply()
{
    std::cout << "Begin refinement" << std::endl;
    
    // For every edge in every face it stores a flag if the corresponding edge has a pending node or 
    // it already has been refined
    int hang_p_nb=0; // number of hanging nodes
    std::vector<std::vector<int> > hanging(cdt_.faces(),std::vector<int>(3,0));

    // Counts triangles that still needs to be checked
    int pending_triangles=0;
    std::vector<int> triangle_status(cdt_.faces(),0);
    
    // Number of refined triangles
    int ref_counter(0);
    
    for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
    {
      // If refinement during time loop, interpolate current and previous solution on existing points      
      if(!init_)
      {
	gas_rangeu(k,3)
	{
	    p_.interpolate(it,it->x(k),it->y(k),new_points_,0);		    
	}
      }      
      
      // If the face has to be refined
      // Mark triangles to be refined
      if (TbR_[it->i()]) 
      {
	  ref_counter++;
	  // Count number of hanging nodes
	  int h_tot=0;
	  gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }

	  // If the triangle has not been refined yet, call triple bisection routine
	  if(h_tot>=0)
	  {
	    this->triple_bisection(it,it->reference_edge());

	    hang_p_nb -= h_tot;
	    
	    // Triangles refined
	    triangle_status[it->i()]=-1;
	    
	    // Set neighbors to be checked for the presence of hanging nodes
	    // If so, set the neighbors to be further checked
	    gas_rangeu(k,3)
	    {
	      hanging[it->i()][k]=-1;
	      
	      if(hanging[it->neighbor(k)][it->loc_mirror_vertex(k)] == 0)
	      {
		hanging[it->neighbor(k)][it->loc_mirror_vertex(k)]=1;
		hang_p_nb++;		  		
		triangle_status[it->neighbor(k)]=1;
	      }
	    }
	  }
	
      }// To be Refined
    } // end face loop
    
    std::cout << "Number of refined triangles: " << ref_counter << std::endl;
    
    // Checking conformity for the mesh
    // Analyze neighbors to mark edges for refinement to guarantee conformity
    do {  
	  int red_counter=0;
	  
	  // Counts faces with three hanging nodes => Marked for triple bisection
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    int h_tot=0;
	    gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }
	    
	    if(h_tot==0) { triangle_status[it->i()]=0; }
	    // Triangles either already refined or with three hanging nodes
	    else if( (h_tot==3) || (h_tot==-3) )
	    {
	      triangle_status[it->i()]=-1;
	      red_counter++;
	    }
	  }

	  // Analyze faces for green and blue refinement by checking if reference edge is marked
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    int h_tot=0;
	    gas_rangeu(k,3) { h_tot+=hanging[it->i()][k]; }
	    
	    // Triangle with one hanging node
	    if(h_tot == 1)
	    {
	      // If current reference is not marked, mark it
	      if(hanging[it->i()][it->reference_edge()] == 0) 
	      {
		hanging[it->i()][it->reference_edge()]=1;		  
		hang_p_nb++;
		
		triangle_status[it->i()]=0;
		
		// If the neighbor is not marked for triple bisection, mark the edge
		if(triangle_status[it->neighbor(it->reference_edge())] != -1)
		{		
		  if(hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())] == 0)
		  {
		    hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())]=1;		  
		    hang_p_nb++;		
		    triangle_status[it->neighbor(it->reference_edge())]=1;
		  }
		}
	      }
	      else { triangle_status[it->i()]=0; }
	    }
	    else if(h_tot == 2) // Triangle with two hanging nodes
	    {
	      // If current reference is not marked, mark it	      
	      if(hanging[it->i()][it->reference_edge()] == 0) 
	      {
		hanging[it->i()][it->reference_edge()]=1;		  
		hang_p_nb++;
		
		triangle_status[it->i()]=-1;
		red_counter++;
		
		// If the neighbor is not marked for triple bisection, mark the edge		
		if(triangle_status[it->neighbor(it->reference_edge())] != -1)
		{		
		  if(hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())] == 0)
		  {
		    hanging[it->neighbor(it->reference_edge())][it->loc_mirror_vertex(it->reference_edge())]=1;		  
		    hang_p_nb++;
		    triangle_status[it->neighbor(it->reference_edge())]=1;		    
		  }
		}
	      }
	      else { triangle_status[it->i()]=0; }	      
	    }
	  }
	  
	  // Update the number of triangles to be checked for conformity	  
	  pending_triangles = red_counter;
	  for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
	  {
	    pending_triangles += triangle_status[it->i()];
	  }

    } while (pending_triangles > 0);

    std::cout << "Total hanging nodes: " << hang_p_nb << std::endl;
   
    // Number of compatibility triangles
    int comp_counter(0);
    
    for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it) 
    {
	int h_tot=0;
	int scenario=0;
	
	gas_rangeu(k,3)	{ h_tot+=hanging[it->i()][k]; }
	
	if(h_tot==1) { scenario=1; }
	else if(h_tot==2) { scenario=2; }
	else if(h_tot==3) { scenario=3; }	
	else { scenario=0; }

	// Choose refinement strategy based on the number of hanging nodes
	switch(scenario)
	{
	  case 0: // Nothing to do
	  {
	    this->store_previousT(it);
	  }
	  break;
	  
	  case 1: // Green bisection
	  {
	    comp_counter++;
	    
	    this->one_bisection(it,it->reference_edge());
	    hang_p_nb -= h_tot;	    
	    
	    gas_rangeu(k,3) { hanging[it->i()][k]=-1; } // tombstone
	  }
	  break;
	  
	  case 2: // Blue bisection
	  {
	    comp_counter++;
	    
    	    int idx;
	    gas_rangeu(k,3)
	    { 
	      if((hanging[it->i()][k]==1) && (k != it->reference_edge())) { idx=k; }
	    }
	    this->two_bisection(it,it->reference_edge(),idx);
	    hang_p_nb -= h_tot;
	    
	    gas_rangeu(k,3) { hanging[it->i()][k]=-1; } // tombstone
	  }
	  break;
	  
	  case 3: // Triple bisection
	  {
	    comp_counter++;

	    this->triple_bisection(it,it->reference_edge());
	    hang_p_nb -= h_tot;
	    
	    gas_rangeu(k,3)
	    {
	      hanging[it->i()][k]=-1;
	    }
	  }
	  break;
	}
	
      }

      std::cout << "Checking conformity - Triangles to impose compatibility: " << comp_counter << std::endl;
      
      std::cout << "End refinement" << std::endl;

}  

  
void NewestVertexRefinement::store_previousT(face_iterator_t & it)
{
    int idx = it->reference_edge();
    
    double Px, Py;
    double Fx, Fy;	  
    double Bx, By;	  
    
    // Based on the reference edge, set values for P, F and B
    switch(idx)
    {
      case 0:
      {
	  Px = it->x(0);
	  Py = it->y(0);	            
	  Fx = it->x(1);
	  Fy = it->y(1);     
	  Bx = it->x(2);
	  By = it->y(2);
      }
      break;
      
      case 1:
      {
	  Px = it->x(1);
	  Py = it->y(1);	            
	  Fx = it->x(2);
	  Fy = it->y(2);     
	  Bx = it->x(0);
	  By = it->y(0);	      
      }
      break;
      
      case 2:
      {
	  Px = it->x(2);
	  Py = it->y(2);	            
	  Fx = it->x(0);
	  Fy = it->y(0);     
	  Bx = it->x(1);
	  By = it->y(1);	      
      }
      break;
    }
    
    // Add current face to the vector of previously bisected triangles
    std::list<std::pair<double,double> > coord;

    coord.push_back(std::make_pair(Px,Py));    
    coord.push_back(std::make_pair(Fx,Fy));            
    coord.push_back(std::make_pair(Bx,By));        
    
    cdt_.insert_biFace(coord);
}
 
 
void NewestVertexRefinement::update()
{
    // Update labels within triangulation
    cdt_.update_dataNV();
} 
   
  
/*! @brief Factory to build refinement method, based on initial data given by user */    
class RefinementFactory
{	
  
    public:
	    /*! @brief Constructor 
	     *  @param type Type of refinement method, depending on given data
	     *  @param p Problem
	     *  @param cdt Triangulation
	     *  @param TbR Vector that keeps the information about whether or not elements have to be refined
	     *  @param to_be_inserted Vector that stores the elements to be inserted
	     *  @param to_be_removed Vector that stores the edges to be removed
	     *  @param new_points Vector that stores the values of the solution interpolated on new points in the triangulation
	     *  @param init Boolean variable to flag if refinement is performed at the beginning or duriong evolution loop
	     */
	    RefinementFactory ( int type, problem const & p, triangulation & cdt, std::vector<bool> & TbR, std::vector<std::list<std::pair<double,double> > > & to_be_inserted, std::vector<std::pair<bool,std::list<int> > > & to_be_removed, std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & new_points, const bool init ):
	    type_(type), p_(p), cdt_(cdt), TbR_(TbR), to_be_inserted_(to_be_inserted), to_be_removed_(to_be_removed), new_points_(new_points), init_(init) {};

	    /*! @brief Virtual destructor */
	    virtual ~RefinementFactory() {};
	    
	    /*! @brief Create a refinement method object*/
	    RefineMethod* create();
	
    private:
	    /*! @brief Parameter that identifies which refinement method has been selected */
	    int type_;
      
	    /*! @brief Problem */
	    problem const & p_;
	    
	    /*! @brief Triangulation */
	    triangulation & cdt_;	
	    
	    /*! @brief Vector that stores if a triangle has to be refined */
	    std::vector<bool> & TbR_;
	    
	    /*! @brief Vector of elements*/
	    std::vector<std::list<std::pair<double,double> > > & to_be_inserted_;

	    /*! @brief Vector of the edges to be removed*/
	    std::vector<std::pair<bool,std::list<int> > > & to_be_removed_;
	    
	    /*! @brief Vector that interpolates the solution in the new points of the triangulation */
	    std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > & new_points_;
	    
	    /*! @brief Boolean variable to identify if the refinement is at the begininng or during evolution */
	    const bool init_;
	    
};
    
    
RefineMethod* RefinementFactory::create()
{
    // Based on the parameter type_ calls a different constructor
    if (type_ == 0) { return new RegularDivisionRefinement( p_, cdt_, TbR_, to_be_inserted_, to_be_removed_, new_points_, init_ ); }
    else if (type_ == 1) { return new LongestEdgeRefinement( p_, cdt_, TbR_, to_be_inserted_, to_be_removed_, new_points_, init_ ); }
    else if (type_ == 2) { return new NewestVertexRefinement( p_, cdt_, TbR_, to_be_inserted_, to_be_removed_, new_points_, init_ ); }
    
    return 0;
}  
  
  
} //namespace


#endif