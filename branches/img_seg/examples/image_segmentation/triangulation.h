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
 * @file triangulation.h
 * @brief Includes the implementation of the class that realizes a wrapper to CGAL mesh
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_TRIANGULATION_H
#define _IMG_SEG_TRIANGULATION_H


namespace img_seg {

  
/*! @brief Mesh and geometry information */  
class triangulation {

  private:
	  /*! @brief The information to save in the node of triangulation */
	  class point_info
	  {
	    public:
		    /*! @brief Constructor */
		    inline point_info (): index(0u), removable(true) {};
		    
		    /*! @brief Destructor */
		    ~point_info() {};

		    /*! @brief Index of the node */
		    int index;

		    /*! @brief A flag to check if the node is removable */
		    bool removable;
	  };

	  /*! @brief The information to save in the face of triangulation */
	  class face_info
	  {
	    public:
		    /*! @brief Constructor */
		    inline face_info (): index(0u), ref_idx(0) {};

		    /*! @brief Destructor */
		    ~face_info() {};		    
		    
		    /*! @brief The index of face */
		    int index;

		    /*! @brief Index f the reference edge */
		    int ref_idx;
	  };

	  // CGAL wrapper //
	  // Kernel 
	  struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};

	  // Vertex triangulation
	  typedef CGAL::Triangulation_vertex_base_with_info_2<point_info, K> vb_t;

	  // Face triangulation
	  typedef CGAL::Delaunay_mesh_face_base_2<K,
		      CGAL::Constrained_triangulation_face_base_2<K,
			  CGAL::Triangulation_face_base_with_info_2<face_info , K>
		      >
		  > fb_t;

	  // Data structure
	  typedef CGAL::Triangulation_data_structure_2<vb_t, fb_t> tds_t;
	  typedef tds_t::Face_handle Face_handle;
	  typedef CGAL::Exact_predicates_tag tag_t;

	  // Triangulation
	  typedef CGAL::Constrained_Delaunay_triangulation_2<K, tds_t, tag_t> cdt_t;

private:
	/*! @brief A wrapper to CGAL triangulation face*/
	class face;
	
	/*! @brief An iterator over a face*/	
	class face_iterator;

	/*! @brief Class to manage the nodes of the triangulation */
	class all_nodes;
	
	/*! @brief An iterator over the nodes */	
	class nodes_iterator;

private:
	/*! @brief A wrapper to CGAL triangulation face*/
	class face
	{
	    private:
		    typedef cdt_t::Finite_faces_iterator cgal_face_iterator_t;

		    /*! @brief Original CGAL face iterator */
		    cgal_face_iterator_t it_;
		    
		    /*! @brief Constructor */
		    inline face (): it_() {};
		    inline face (cgal_face_iterator_t const & it): it_(it) {};

	    public:
		    /*! @brief Destructor */
		    ~face() {};
		    
		    /*! @brief Return the index of the face */
		    inline int i () const
		    {
			return it_->info().index;
		    }	      

		    /*! @brief Return the reference edge of the face */
		    inline int reference_edge () const {

			    return it_->info().ref_idx;
		    }		
		    
		    /*! @brief Set the reference edge of the face equal to i*/
		    inline void set_ref_edge (int const & i)
		    {
			gas_assert(i < 3);
			it_->info().ref_idx = i;
		    }		    
		    
		    /*! @brief Return the index of i-th node of the face */
		    inline int i (int const & i) const
		    {
			gas_assert(i < 3);
			return it_->vertex(i)->info().index;
		    }	      
	      
		    /*! @brief Return the first coordinate of i-th node of the face */
		    inline double x (int const & i) const
		    {
			gas_assert(i < 3);
			return it_->vertex(i)->point().x();
		    }

		    /*! @brief Return the second coordinate of i-th node of the face */
		    inline double y (int const & i) const
		    {
			gas_assert(i < 3);
			return it_->vertex(i)->point().y();
		    }

		    /*! @brief Return the index of the neighbor opposite to the i-th edge */
		    inline int neighbor (int const & i) const
		    {
			gas_assert(i < 3);
			return it_->neighbor(i)->info().index;
		    }		

		    /*! @brief Return the index of the j-th node of the neighbor opposite to the i-th edge */
		    inline int neighbor_idx (int const & i, int const & j) const
		    {
			    gas_assert(i < 3);
			    gas_assert(j < 3);
			    return it_->neighbor(i)->vertex(j)->info().index;
		    }		    
		    
		    /*! @brief Return the first coordinate of the j-th node of the neighbor opposite to the i-th edge */
		    inline double neighbor_x (int const & i, int const & j) const
		    {
			gas_assert(i < 3);
			gas_assert(j < 3);
			return it_->neighbor(i)->vertex(j)->point().x();
		    }
		    
		    /*! @brief Return the second coordinate of the j-th node of the neighbor opposite to the i-th edge */
		    inline double neighbor_y (int const & i, int const & j) const
		    {
			gas_assert(i < 3);
			gas_assert(j < 3);
			return it_->neighbor(i)->vertex(j)->point().y();
		    }
		   
		    /*! @brief Return the index of the vertex opposite to the i-th edge within the i-th neighbor*/
		    inline int loc_mirror_vertex (int const & i) const
		    {
			gas_assert(i < 3);
			return it_->mirror_index(i);
		    }

		    /*! @brief Return a face handle for the current face */
		    inline Face_handle face_handle () const 
		    {
			int k=it_->mirror_index(0);
			return it_->neighbor(0)->neighbor(k);
		    }
		    
		    /*! @brief Return the index of the longest edge */
		    int longest_edge() const;
		    
		    // Allow the iterator over the faces to access to the face information
		    friend class face_iterator;
		    
		    // Allow the estimators to access to the face information
		    friend class posterior;

	    };

	    
	    /*! @brief Class to manage the nodes of the triangulation */
	    class all_nodes
	    {
		private:
			typedef cdt_t::Finite_vertices_iterator cgal_nodes_iterator_t;

			/*! @brief Original CGAL node circulator */
			cgal_nodes_iterator_t nodes_;

			/*! @brief Constructor */
			inline all_nodes (): nodes_() {};
			inline all_nodes (cgal_nodes_iterator_t const & nodes): nodes_(nodes) {};
			
		public:
			/*! @brief Destructor */
			~all_nodes() {};		  
		  
			/*! @brief Return the index of the node */
			inline int i () const
			{
			    return nodes_->info().index;
			}
			
			/*! @brief Return the first coordinate of the node */
			inline double x () const
			{
			    return nodes_->point().x();
			}

			/*! @brief Return he second coordinate of the node */
			inline double y () const 
			{
			    return nodes_->point().y();
			}
			
			// Allow the iterator over the nodes to access to the node information
			friend class nodes_iterator;
			
	      };	

public:
	typedef face face_t;
	typedef all_nodes nodes_t;

private:
	/*! @brief An iterator over a face*/	
	class face_iterator
	{

	    private:
		    typedef cdt_t::Finite_faces_iterator cgal_face_iterator_t;
		    typedef face_iterator self_t;
		    
		    /*! @brief Constructor */
		    inline face_iterator (cgal_face_iterator_t const & it): f_(it) {};		    

	    public:
		    /*! @brief Destructor */
		    ~face_iterator() {};		    
		    
		    typedef face_t const * pointer_t;
		    typedef face_t const & reference_t;

		    /*! @brief Constructor */
		    inline face_iterator (): f_() {};
		    inline face_iterator (self_t const & i): f_(i.f_.it_) {};

		    /*! @brief Operator = overloading */
		    inline self_t & operator= (self_t const & i) 
		    {
			f_.it_ = i.f_.it_;
			return *this;
		    }

		    /*! @brief Operator == overloading */
		    inline bool operator== (self_t const & i) const { return (f_.it_ == i.f_.it_); }
		    
		    /*! @brief Operator != overloading */
		    inline bool operator!= (self_t const & i) const { return (f_.it_ != i.f_.it_); }

		    /*! @brief Operator ++ overloading */
		    inline self_t & operator++ () 
		    {
			++(f_.it_);
			return *this;
		    }
		    
		    inline self_t operator++ (int) 
		    {
			self_t _it(*this);
			++(f_.it_);
			return _it;
		    }

		    /*! @brief Operator * overloading */
		    inline reference_t operator* () const { return f_; }
		    
		    /*! @brief Operator -> overloading */
		    inline pointer_t operator-> () const { return &f_; }

	    private:
		    /*! @brief Face class */
		    face f_;

		    // Allow the triangulation and the error estimator to access to the iterator over the face 
		    friend class triangulation;
		    friend class posterior;
		    
	};

		
	/*! @brief An iterator over the nodes */	
	class nodes_iterator
	{

	private:
		typedef cdt_t::Finite_vertices_iterator cgal_nodes_iterator_t;
		typedef nodes_iterator self_t;
		
		/*! @brief Constructor */
		inline nodes_iterator (cgal_nodes_iterator_t const & a): nd_(a) {};

	public:
		/*! @brief Constructor */
		~nodes_iterator() {};
		
		typedef nodes_t const * pointer_t;
		typedef nodes_t const & reference_t;
	
		/*! @brief Constructor */
		inline nodes_iterator (): nd_() {};
		inline nodes_iterator (self_t const & i): nd_(i.nd_.nodes_) {};

		/*! @brief Operator = overloading */
		inline self_t & operator= (self_t const & i) 
		{
		    nd_.nodes_ = i.nd_.nodes_;
		    return *this;
		}

		/*! @brief Operator == overloading */
		inline bool operator== (self_t const & i) const { return (nd_.nodes_ == i.nd_.nodes_); }
		
		/*! @brief Operator != overloading */
		inline bool operator!= (self_t const & i) const { return (nd_.nodes_ != i.nd_.nodes_); }

		/*! @brief Operator ++ overloading */
		inline self_t & operator++ () 
		{
		    ++(nd_.nodes_);
		    return *this;
		}
		
		inline self_t operator++ (int) 
		{
		    self_t _a(*this);
		    ++(nd_.nodes_);
		    return _a;
		}

		/*! @brief Operator * overloading */
		inline reference_t operator* () const { return nd_; }
		
		/*! @brief Operator -> overloading */
		inline pointer_t operator-> () const { return &nd_; }

	private:
		/*! @brief Nodes class */
		all_nodes nd_;		

		// Allow the triangulation to access to the iterator over the nodes
		friend class triangulation;
	};	

public:
	typedef cdt_t::Point point_t;
	typedef face_iterator face_iterator_t;
	typedef nodes_iterator nodes_iterator_t;

	/*! @brief Constructor
	 *  @param begin The iterator to the the first point in the list of nodes
	 *  @param end The iterator to the last point in the list of nodes
	 *  @param h Space step
	 */
	template <typename iterator_>
	triangulation (iterator_ begin, iterator_ end, double h = 0.);

	/*! @brief Destructor */
	~triangulation() {};
	
	/*! @brief Return the number of nodes in the triangulation */
	inline int nodes () const { return n_nodes_; }

	/*! @brief Return the number of faces in the triangulation */
	inline int faces () const { return n_faces_; }

	/*! @brief Return a const iterator to the first face */
	inline face_iterator_t face_begin() const 
	{
	    return face_iterator_t(cdt_.finite_faces_begin());
	}

	/*! @brief Return a const iterator to the last face */
	inline face_iterator_t face_end() const 
	{
	    return face_iterator_t(cdt_.finite_faces_end());
	}
	
	/*! @brief Return a const iterator to the first node */
	inline nodes_iterator_t nodes_begin() const 
	{
	    return nodes_iterator_t(cdt_.finite_vertices_begin());
	}
	
	/*! @brief Return a const iterator to the last node */
	inline nodes_iterator_t nodes_end() const 
	{
	    return nodes_iterator_t(cdt_.finite_vertices_end());
	}
	
	/*! @brief Insert a face in the vector of previously bisected faces
	 *  @param lst A list containing the coordinates of the vertices of the face
	 */
	void insert_biFace (const std::list<std::pair<double,double> > & lst);

	/*! @brief Insert a face in the vector of previously refined faces
	 *  @param lst A list containing the coordinates of the vertices of the face
	 */
	void insert_redFace (const std::list<std::pair<double,double> > & lst);
	
	/*! @brief Remove a face from the vector of previously bisected faces */
	void remove_biFace (face_iterator_t & it);

	/*! @brief Remove a face from the vector of previously refined faces */
	void remove_redFace (face_iterator_t & it);
	
	/*! @brief Verify if a given face is contained in the vector of previously bisected faces
	 *  @param it Face we are looking for
	 *  @param idx If the face is found stores its index
	 */
	std::vector<std::list<std::pair<double,double> > >::iterator is_bisected(face_iterator_t & it, int & idx);

	/*! @brief Verify if a given face is contained in the vector of previously refined faces
	 *  @param it Face we are looking for
	 *  @param idx If the face is found stores its index
	 */
	std::vector<std::list<std::pair<double,double> > >::iterator is_refined(face_iterator_t & it, int & idx);		
		
	/*! @brief Insert an edge as a constraint in the triangulation
	 *  @param left Coordinates of the first vertex of the edge
	 *  @param right Coordinates of the second vertex of the edge
	 */
	void insertEdge (const std::pair<double,double> & left, const std::pair<double,double> & right);

	/*! @brief Remove an edge from the triangulation 
	 *  @param it Face in analysis
	 *  @param i Index of the edge to be removed (using local face numeration)
	 */
	void removeEdge (face_iterator_t & it, const int & i);
	
	/*! @brief Update triangulation after RegularDivision refinement
	 *  New labels for triangles and nodes and new reference edges
	 */
	void update_dataRD();
	
	/*! @brief Update triangulation after LongestEdge refinement
	 *  New labels for triangles and nodes and new reference edges
	 */
	void update_dataLE();
	
	/*! @brief Update triangulation after NewestVertex refinement
	 *  New labels for triangles and nodes and new reference edges
	 */
	void update_dataNV();

private:
	/*! @brief The internal structure for the triangulation */
	cdt_t cdt_;

	/*! @brief Number of nodes */
	int n_nodes_;

	/*! @brief Number of faces */
	int n_faces_;
	
	/*! @brief Container that stores the vertices of bisected triangles */
	std::vector<std::list<std::pair<double,double> > > bisected_;

	/*! @brief Container that stores the vertices of refined triangles */
	std::vector<std::list<std::pair<double,double> > > refined_;

	// Allow the estimators to access the information about the triangulation
	friend class posterior;
	
	// Allow refinement method to access and modify the triangulation
	friend class RefineMethod;
	friend class AdaptMesh;
};


int triangulation::face::longest_edge() const 
{  
    int idx=0;
    
    // Nodes coordinates
    Eigen::Vector2d const P0(this->x(0), this->y(0));
    Eigen::Vector2d const P1(this->x(1), this->y(1));
    Eigen::Vector2d const P2(this->x(2), this->y(2));

    // Edges information
    double const E0d((P1-P2).norm());
    double const E1d((P2-P0).norm());
    double const E2d((P0-P1).norm());
    Eigen::Vector3d edge(E0d,E1d,E2d);

    // Look for the maximum edge
    if(edge(0) >= edge(1))
    {
      if(edge(1) >= edge(2)) { idx=0; }
      else
      {
	if(edge(0) >= edge(2)) { idx=0; }
	else { idx=2; }
      }
    }
    else
    {
      if(edge(0) >= edge(2)) { idx=1; }
      else
      {
	if(edge(1) >= edge(2)) { idx=1; }
	else { idx=2; }
      }			 		  
    }
    
    return idx;
}


template <typename iterator_>
triangulation::triangulation (iterator_ begin, iterator_ end, double criteria): 
bisected_(), refined_() 
{  
    // Insert nodes and edges 
    {
	typedef std::vector<cdt_t::Vertex_handle> vertex_list_t;
	typedef vertex_list_t::const_iterator iterator_t;

	vertex_list_t _vertex_list;

	// Insert nodes
	for (iterator_ _i = begin; _i != end; ++_i) { _vertex_list.push_back(cdt_.insert(*_i)); }

	// Insert edges between nodes as constraints in the triangulation 
	for (iterator_t _i = _vertex_list.begin() + 1; _i != _vertex_list.end(); ++_i)
	{
	    cdt_.insert_constraint(*(_i-1), *_i);
	}
	cdt_.insert_constraint(_vertex_list.front(), _vertex_list.back());
    }

    // Compute minimum length
    {
	double _minimo(0.);

	for (iterator_ i = begin + 1; i != end; ++i) 
	{
	    double const _dx((i->x() - (i-1)->x()));
	    double const _dy((i->y() - (i-1)->y()));
	    double const _dist(_dx * _dx + _dy * _dy);
	    _minimo = std::min(_dist, _minimo);
	}

	double const _dx((begin->x() - (end-1)->x()));
	double const _dy((begin->y() - (end-1)->y()));
	double const _dist(_dx * _dx + _dy * _dy);
	_minimo = std::min(_dist, _minimo);

	if (criteria <= 0. ) criteria = _minimo;
    }

    // Mark nodes which cannot be removed
    {
	typedef cdt_t::Vertex_circulator circulator_t;

	circulator_t _c(cdt_.incident_vertices (cdt_.infinite_vertex()));
	circulator_t _begin(_c);

	do {
	      _c->info().removable = false;
	      ++_c;
	      
	} while (_c != _begin);
    }

    // Refine mesh
    {
	typedef CGAL::Delaunay_mesh_size_criteria_2<cdt_t> Criteria;

	CGAL::refine_Delaunay_mesh_2(cdt_, Criteria(0.125, criteria));
    }

    // Labels for the nodes
    {
	typedef cdt_t::Finite_vertices_iterator iterator_t;

	int n(0);

	for (iterator_t _i(cdt_.finite_vertices_begin()); _i != cdt_.finite_vertices_end(); ++_i)
	{
	    _i->info().index = n++;
	}

	n_nodes_ = n;
    }

    // Labels for the faces and reference edges
    {
	typedef cdt_t::Finite_faces_iterator iterator_t;

	int n(0);

	for (iterator_t i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
	{
	  i->info().index = n++; // label
	  i->info().ref_idx = face_iterator_t(i)->longest_edge();  // reference edge
	}

	n_faces_ = n;
    }
	
}


void triangulation::insert_biFace (const std::list<std::pair<double,double> > & lst)
{
    // Extract the coordinates of the three vertices of the face
    std::list<std::pair<double,double> >::const_iterator l0=lst.begin();
    std::list<std::pair<double,double> >::const_iterator l1=l0;
    l1++;
    std::list<std::pair<double,double> >::const_iterator l2=l1;
    l2++;
	  
    // If the vertices are different, insert the face in the vector of previously bisected faces
    if((*l0 != *l1) && (*l1 != *l2) && (*l2 != *l0))
    {
	bisected_.push_back(lst);  
    }
}


void triangulation::insert_redFace (const std::list<std::pair<double,double> > & lst)
{
    // Extract the coordinates of the three vertices of the face
    std::list<std::pair<double,double> >::const_iterator l0=lst.begin();
    std::list<std::pair<double,double> >::const_iterator l1=l0;
    l1++;
    std::list<std::pair<double,double> >::const_iterator l2=l1;
    l2++;
	
    // If the vertices are different, insert the face in the vector of previously refined faces
    if((*l0 != *l1) && (*l1 != *l2) && (*l2 != *l0))
    {
      refined_.push_back(lst);  
    }
}


void triangulation::remove_biFace (face_iterator_t & it)
{
    // Verify if the face is contained in the vector of previously bisected faces
    int idx;
    std::vector<std::list<std::pair<double,double> > >::iterator pos(this->is_bisected(it,idx));
    
    // If the face is contained in the vector of previously bisected faces remove it
    if(pos!=bisected_.end()) { bisected_.erase(pos); }  
}


void triangulation::remove_redFace (face_iterator_t & it)
{
    // Verify if the face is contained in the vector of previously refined faces
    int idx;
    std::vector<std::list<std::pair<double,double> > >::iterator pos(this->is_refined(it,idx));  
    
    // If the face is contained in the vector of previously refined faces remove it
    if(pos!=refined_.end()) { refined_.erase(pos); } 
}


std::vector<std::list<std::pair<double,double> > >::iterator triangulation::is_bisected(face_iterator_t & it, int & idx)
{
    // If the vector is empty return the end
    if (bisected_.empty()) { return bisected_.end(); }
    else
    {    
	// Loop over the previously bisected elements
	for(std::vector<std::list<std::pair<double,double> > >::iterator i_=bisected_.begin(); i_!=bisected_.end(); i_++)
	{
	    // Extract the vertices of the triangle
	    std::list<std::pair<double,double> >::const_iterator i_0=(*i_).begin();
	    std::list<std::pair<double,double> >::const_iterator i_1=i_0;
	    i_1++;
	    std::list<std::pair<double,double> >::const_iterator i_2=i_1;
	    i_2++;

	    // Test different rotation of the triangles to verify if an element is contained
	    // Vertices are store in clockwise order within a triangle
	    if(it->x(0)==(*i_0).first && it->y(0)==(*i_0).second)
	    {
		if(it->x(1)==(*i_1).first && it->y(1)==(*i_1).second)
		{
		    if(it->x(2)==(*i_2).first && it->y(2)==(*i_2).second)
		    {
			idx = 0;
			return i_;
		    } 
		}
	    }
	    else if(it->x(1)==(*i_0).first && it->y(1)==(*i_0).second)
	    {
		if(it->x(2)==(*i_1).first && it->y(2)==(*i_1).second)
		{
		    if(it->x(0)==(*i_2).first && it->y(0)==(*i_2).second)
		    {
			idx = 1;
			return i_;
		    } 
		}
	    }
	    else if(it->x(2)==(*i_0).first && it->y(2)==(*i_0).second)
	    {
		if(it->x(0)==(*i_1).first && it->y(0)==(*i_1).second)
		{
		    if(it->x(1)==(*i_2).first && it->y(1)==(*i_2).second)
		    {
			idx = 2;
			return i_;
		    } 
		}
	    }	
      }
      
      return bisected_.end(); 
    }
}


std::vector<std::list<std::pair<double,double> > >::iterator triangulation::is_refined(face_iterator_t & it, int & idx)
{
    // If the vector is empty return the end  
    if (refined_.empty()) return refined_.end();
    else
    {
	// Loop over the previously refined elements
	for(std::vector<std::list<std::pair<double,double> > >::iterator i_=refined_.begin(); i_!=refined_.end(); i_++)
	{
	    // Extract the vertices of the triangle
	    std::list<std::pair<double,double> >::const_iterator i_0=(*i_).begin();
	    std::list<std::pair<double,double> >::const_iterator i_1=i_0;
	    i_1++;
	    std::list<std::pair<double,double> >::const_iterator i_2=i_1;
	    i_2++;
	    
	    // Test different rotation of the triangles to verify if an element is contained
	    // Vertices are store in clockwise order within a triangle
	    if(it->x(0)==(*i_0).first && it->y(0)==(*i_0).second)
	    {
		if(it->x(1)==(*i_1).first && it->y(1)==(*i_1).second)
		{
		    if(it->x(2)==(*i_2).first && it->y(2)==(*i_2).second)
		    {
			idx = 0;
			return i_;
		    } 
		}
	    }
	    else if(it->x(1)==(*i_0).first && it->y(1)==(*i_0).second)
	    {
		if(it->x(2)==(*i_1).first && it->y(2)==(*i_1).second)
		{
		    if(it->x(0)==(*i_2).first && it->y(0)==(*i_2).second)
		    {
			idx = 1;
			return i_;
		    } 
		}
	    }
	    else if(it->x(2)==(*i_0).first && it->y(2)==(*i_0).second)
	    {
		if(it->x(0)==(*i_1).first && it->y(0)==(*i_1).second)
		{
		    if(it->x(1)==(*i_2).first && it->y(1)==(*i_2).second)
		    {
			idx = 2;
			return i_;
		    } 
		}
	    }
	}
	    
	return refined_.end();
    }
  
}


void triangulation::insertEdge (const std::pair<double,double> & left, const std::pair<double,double> & right)
{
    // If the vertices of the edge are different, insert the edge as a constraint in the triangulation
   if((left.first != right.first) && (left.second != right.second))
   {
      cdt_.insert_constraint(point_t(left.first,left.second),point_t(right.first,right.second));	
   }
}


void triangulation::removeEdge (face_iterator_t & it, const int & i)
{
    // Remove edge i within face it from the triangulation
    cdt_.remove_constraint(it->face_handle(),i);
}


void triangulation::update_dataRD()
{
    // Labels for the nodes
    typedef cdt_t::Finite_vertices_iterator iterator_n;

    int nn(0);

    for (iterator_n _i(cdt_.finite_vertices_begin()); _i != cdt_.finite_vertices_end(); ++_i)
    {
	_i->info().index = nn++;
    }

    n_nodes_ = nn;

    // Labels for the faces
    typedef cdt_t::Finite_faces_iterator iterator_f;

    int nf(0);

    for (iterator_f i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
    {
	i->info().index = nf++;
    }
    
    n_faces_ = nf;

    // Reference edges
    for (iterator_f i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
    {    
	face_iterator_t it(i);
	
	// Initially we set reference edge equal to longest edge
	int idx = it->longest_edge();
	int tmp;
	
	// If the triangle has previously been bisected or refined we overwrite this information
	if (this->is_bisected(it,tmp)!=bisected_.end()) { idx = tmp; }
	else if(this->is_refined(it,tmp)!=refined_.end()) { idx = tmp; }
	
	i->info().ref_idx = idx;
    }    
}


void triangulation::update_dataLE()
{
    // Labels for the nodes
    typedef cdt_t::Finite_vertices_iterator iterator_n;

    int nn(0);

    for (iterator_n _i(cdt_.finite_vertices_begin()); _i != cdt_.finite_vertices_end(); ++_i)
    {
      _i->info().index = nn++;
    }
    
    n_nodes_ = nn;

    // Labels for the faces
    typedef cdt_t::Finite_faces_iterator iterator_f;

    int nf(0);

    for (iterator_f i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
    {
	i->info().index = nf++;
    }
    
    n_faces_ = nf;
    
    // Reference edges
    for (iterator_f i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
    {     
	face_iterator_t it(i);
      
	int idx = it->longest_edge();
	
	i->info().ref_idx = idx;
    }
}


void triangulation::update_dataNV()
{
    // Labels for the nodes
    typedef cdt_t::Finite_vertices_iterator iterator_n;

    int nn(0);

    for (iterator_n _i(cdt_.finite_vertices_begin()); _i != cdt_.finite_vertices_end(); ++_i)
    {
	_i->info().index = nn++;
    }
    
    n_nodes_ = nn;

    // Labels for the faces
    typedef cdt_t::Finite_faces_iterator iterator_f;

    int nf(0);

    for (iterator_f i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
    {
	i->info().index = nf++;
    }
    
    n_faces_ = nf;
    
    // Reference edges
    for (iterator_f i(cdt_.finite_faces_begin()); i != cdt_.finite_faces_end(); ++i)
    {     
	face_iterator_t it(i);
      
	// Initially we set reference edge equal to longest edge
	int idx = it->longest_edge();
	int tmp;
	
	// If the triangle has previously been bisected we overwrite this information
	if (this->is_bisected(it,tmp)!=bisected_.end()) { idx = tmp; }
	
	i->info().ref_idx = idx;
    }
    
}


} //namespace


#endif 