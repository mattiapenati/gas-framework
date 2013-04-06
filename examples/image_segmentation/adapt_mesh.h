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
 * @brief Includes the implementation of the class to manage mesh adaptivity routines
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_ADAPT_H
#define _IMG_SEG_ADAPT_H


namespace img_seg {

  
/*! @brief Sorting strategy for std::pair < Eigen::Vector2d, Eigen::Vector2d >  */
struct sort_pred
{
    bool operator() (const std::pair<Eigen::Vector2d,Eigen::Vector2d> & n, const std::pair<Eigen::Vector2d,Eigen::Vector2d> & m)
    {
	// Compare x
	if (n.first(0) < m.first(0)) return true;
	else if (n.first(0) > m.first(0)) return false;
	else 
	{
	    // Compare y
	    if (n.first(1) < m.first(1)) return true;
	    else return false;
	}
    }  
};


/*! @brief Comparing strategy for std::pair < Eigen::Vector2d, Eigen::Vector2d > */
struct unique_pred 
{
    bool operator() (const std::pair<Eigen::Vector2d,Eigen::Vector2d> & n, const std::pair<Eigen::Vector2d,Eigen::Vector2d> & m)
    { 
	// Compare x and y
	return ( (n.first(0) == m.first(0)) && (n.first(1) == m.first(1)) );
    }
};


/*! @brief  Sorting strategy for std::list < std::pair < double, double > > */
struct sort_VP
{
    bool operator() (const std::list<std::pair<double,double> > & n, const std::list<std::pair<double,double> > & m)
    {
	// Extract the information from the list 
	// Due to the nature of the algorithm, every list contains exactly two pairs 
	// The first pair contains the coordinates of a point and the second the values of the 
	// solution at previous and current time step
	std::list<std::pair<double,double> >::const_iterator i_n0(n.begin());
	std::list<std::pair<double,double> >::const_iterator i_m0(m.begin());
	std::list<std::pair<double,double> >::const_iterator i_n1(i_n0);
	i_n1++;
	std::list<std::pair<double,double> >::const_iterator i_m1(i_m0);	
	i_m1++;
	
	// Verify if the elements are equal
	if ( ((*i_n0).first == (*i_m0).first) && ((*i_n0).second == (*i_m0).second) && ((*(i_n1)).first == (*(i_m1)).first) && ((*(i_n1)).second == (*(i_m1)).second) )
	{
	    return true;
	}

	// Compare x
	if ((*i_n0).first < (*i_m0).first) return true;
	else if ((*i_n0).first > (*i_m0).first) return false;
	else
	{
	    // Compare y
	    if ((*i_n0).second < (*i_m0).second) return true;
	    else if ((*i_n0).second > (*i_m0).second) return false;
	    else
	    {
		// Compare old solution
		if ((*(i_n1)).first < (*(i_m1)).first) return true;
		else if ((*(i_n1)).first > (*(i_m1)).first) return false;
		else
		{
		    // Compare solution
		    if ((*(i_n1)).second < (*(i_m1)).second) return true;
		    else return false;	      
		}
	    }
	}	
    }
};
  

/*! @brief Comparing strategy for std::list < std::pair < double, double > > */
struct unique_VP
{
    bool operator() (const std::list<std::pair<double,double> > & n, const std::list<std::pair<double,double> > & m)
    { 
	// Extract the information from the list 
	// Due to the nature of the algorithm, every list contains exactly two pairs 
	// The first pair contains the coordinates of a point and the second the values of the 
	// solution at previous and current time step      
	std::list<std::pair<double,double> >::const_iterator i_n0(n.begin());
	std::list<std::pair<double,double> >::const_iterator i_m0(m.begin());
	std::list<std::pair<double,double> >::const_iterator i_n1(i_n0);
	i_n1++;
	std::list<std::pair<double,double> >::const_iterator i_m1(i_m0);	
	i_m1++;
	
	// Verify if the elements are equal
	return ( ((*i_n0).first == (*i_m0).first) && ((*i_n0).second == (*i_m0).second) && ((*(i_n1)).first == (*(i_m1)).first) && ((*(i_n1)).second == (*(i_m1)).second) );
    }
};


/*! @brief Class that manages adaptivity procedures */
class AdaptMesh
{
    public: 
	    typedef triangulation::face_iterator face_iterator_t;
	    
	    /*! @brief Constructor 
	     *  @param cdt Triangulation
	     *  @param stim Posterior class object to compute error estimates and mark elements for refinement
	     */
	    AdaptMesh (triangulation & cdt, posterior & stim);    
	    
	    /*! @brief Destructor */
	    ~AdaptMesh() {};
	    
	    /*! @brief Initial adaptivity routine
	     *  @return A pair whose first elements is the number of inserted nodes and whose second
	     * 		element is the global error
	     */
	    std::pair<int,double> initial_adapt();
	    
	    /*! @brief Adaptivity routine at every time step
	     *  @return A pair whose first elements is the number of inserted nodes and whose second
	     * 		element is the global error
	     */
	    std::pair<int,double> loop_adapt();
    
    private:
	    /*! @brief Triangulation */
	    triangulation & cdt_;
	
	    /*! @brief Posterior object (used for estimate and marking strategy) */
	    posterior & stimator_;
	    
	    /*! @brief Problem */
	    problem & p_;
	
	    /*! @brief Vector that contains the triangles to be inserted */
	    std::vector<std::list<std::pair<double,double> > > to_be_inserted_;

	    /*! @brief Vector that contains the edges to be removed */
	    std::vector<std::pair<bool,std::list<int> > > to_be_removed_;
	    
};
  
  
AdaptMesh::AdaptMesh (triangulation & cdt, posterior & stim): 
cdt_(cdt), stimator_(stim), p_(stim.pb_return()), to_be_inserted_()
{
    for(int i_(0); i_<cdt_.faces(); i_++)
    {
	std::list<int> edge_list;
	to_be_removed_.push_back(std::make_pair(false,edge_list));    
    }       
};


std::pair<int,double> AdaptMesh::initial_adapt()
{    
    // Flag the procedure as "initial adaptivity"
    bool initial_adapt(true);
  
    // Number of nodes present before refinement
    int before(cdt_.nodes());

    // Initial local and global error estimate
    Eigen::VectorXd L2_int_err;
    L2_int_err.setZero(cdt_.faces());
    double cumulative_err(0.0);
    
    for (face_iterator_t i(cdt_.face_begin()); i != cdt_.face_end(); ++i)
    {
	// Compute initial local error estimate
	L2_int_err(i->i()) = stimator_.initial_error(i);
	cumulative_err += std::pow(L2_int_err(i->i()),2);
    }
        
    double L2_glob_err(std::sqrt(cumulative_err));

    // If tolerance is not fulfilled
    if(L2_glob_err > stimator_.init_tol())
    {
	// Mark triangles to be refined using GERS strategy
	std::vector<bool> refine_mark(L2_int_err.size(), false);
	stimator_.mark(L2_int_err,L2_glob_err,refine_mark);

	// Vector to store the results of local interpolation for newly added points
	std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > new_points;		
	
	// Create a strategy for the refinement
	RefinementFactory* myRefinement = new RefinementFactory(stimator_.data().method, p_, cdt_, refine_mark, to_be_inserted_, to_be_removed_, new_points, initial_adapt);
	RefineMethod* strategy = myRefinement->create();
	
	// Apply refinement
	strategy->apply();
	
	// Remove duplicates in the list of faces to be inserted
	std::sort(to_be_inserted_.begin(),to_be_inserted_.end(),sort_VP());
	to_be_inserted_.erase(std::unique(to_be_inserted_.begin(),to_be_inserted_.end(),unique_VP()),to_be_inserted_.end());
	    
	// Remove previously marked edges
	for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it)
	{
	    if (to_be_removed_[it->i()].first)
	    { 
		std::list<int> edge_list = to_be_removed_[it->i()].second;
		
		for(std::list<int>::iterator lt(edge_list.begin()); lt != edge_list.end(); lt++)
		{
		    // Remove edge with index *lt within triangle it
		    cdt_.removeEdge(it,*lt); 
		}
	    }
	}

	// Add previously marked edges
	for(std::vector<std::list<std::pair<double,double> > >::iterator i_(to_be_inserted_.begin()); i_!=to_be_inserted_.end(); i_++)  
	{
	    // Extract information about vertices from the list
	    std::list<std::pair<double,double> >::iterator i_left(i_->begin());
	    std::list<std::pair<double,double> >::iterator i_right = i_left;
	    i_right++;
	    
	    cdt_.insertEdge((*i_left),(*i_right));  	
	}

	// Update labels and reference edges
	strategy->update();

	// Clean up memory
	delete strategy;
	strategy = NULL;
	delete myRefinement;
	myRefinement = NULL;
      
	return std::make_pair((cdt_.nodes()-before),L2_glob_err);
    }
    else // Tolerance fulfilled 
    {
	std::cout<<"Initial tolerance fulfilled -- No need to refine"<<std::endl;	  
	return std::make_pair(0,L2_glob_err);
    }	
}


std::pair<int,double> AdaptMesh::loop_adapt()
{    
    // Flag the procedure as "loop adaptivity"
    bool initial_adapt(false);

    // Number of nodes before refinement
    int before(cdt_.nodes());

    // Local and global error estimates
    Eigen::VectorXd L2_int_err;
    L2_int_err.setZero(cdt_.faces());
    double cumulative_err(0.0);
    
    for (face_iterator_t i(cdt_.face_begin()); i != cdt_.face_end(); ++i)
    {
	// Compute local error estimates
	L2_int_err(i->i()) = stimator_.loop_error(i);
	cumulative_err += std::pow(L2_int_err(i->i()),2);
    }
    
    double L2_glob_err(std::sqrt(cumulative_err));   
    
    // If tolerance is not fulfilled
    if(L2_glob_err > stimator_.loop_tol())
    {
	// Mark triangles to be refined using GERS strategy
	std::vector<bool> refine_mark(L2_int_err.size(), false);
	stimator_.mark(L2_int_err,L2_glob_err,refine_mark);
	
	// Vector to store the results of local interpolation for newly added points
	std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d> > new_points;		
	
	// Create a refinement strategy
	RefinementFactory* myRefinement = new RefinementFactory(stimator_.data().method, p_, cdt_, refine_mark, to_be_inserted_, to_be_removed_, new_points, initial_adapt);
	RefineMethod* strategy = myRefinement->create();
	
	// Apply refinement
	strategy->apply();

	// Remove duplicates in the list of new points of the triangulation
	std::sort(new_points.begin(),new_points.end(),sort_pred());
	new_points.erase(std::unique(new_points.begin(),new_points.end(),unique_pred()),new_points.end());

	// Remove duplicates in the list of faces to be inserted
	std::sort(to_be_inserted_.begin(),to_be_inserted_.end(),sort_VP());
	to_be_inserted_.erase(std::unique(to_be_inserted_.begin(),to_be_inserted_.end(),unique_VP()),to_be_inserted_.end());      
	
	// Remove previously marked edges
	for (face_iterator_t it(cdt_.face_begin()); it != cdt_.face_end(); ++it)
	{
	    if (to_be_removed_[it->i()].first) 
	    { 
		std::list<int> edge_list = to_be_removed_[it->i()].second;
		
		for(std::list<int>::iterator lt(edge_list.begin()); lt != edge_list.end(); lt++)
		{
		    // Remove edge with index *lt within triangle it
		    cdt_.removeEdge(it,*lt); 
		}	  
	    }
	}
	
	// Add previously marked edges
	for(std::vector<std::list<std::pair<double,double> > >::iterator i_(to_be_inserted_.begin()); i_!=to_be_inserted_.end(); i_++)  
	{
	    // Extract information about vertices from the list
	    std::list<std::pair<double,double> >::iterator i_left(i_->begin());
	    std::list<std::pair<double,double> >::iterator i_right = i_left;
	    i_right++;
	    
	    
	    cdt_.insertEdge((*i_left),(*i_right));  	
	}
	
	// Update labels and refinement edges
	strategy->update();
	
	// Update the solution with newly inserted points and their interpolation over the mesh
	p_.updateInterpolator(new_points);
	
	// Clean up memory
	delete strategy;
	strategy = NULL;
	delete myRefinement;
	myRefinement = NULL;
	
	return std::make_pair((cdt_.nodes()-before),L2_glob_err);
    }
    else // Tolerance fulfilled
    {
	std::cout<<"Loop tolerance fulfilled -- No need to refine"<<std::endl;	  
	return std::make_pair(0,L2_glob_err);
    }
}

  
}

#endif