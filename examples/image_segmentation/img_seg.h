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
 * @file img_seg.h
 * @brief Includes the headers for the image segmentation problem
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_H
#define _IMG_SEG_H


#include <gas>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>
#include <string>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <Eigen/Array>

#include "CImg.h"

#include <libconfig.h++>


/*!
 * @namespace img_seg
 * @brief Namespace for image segmentation
 */
namespace img_seg {


/*! @brief Enumeration to set the type of the force term */  
enum force_type
{
  ChanVese,
  RegionGaussian,
  RegionNonParametric
};
  

/*! @brief Enumeration to set the refinement method */  
enum refinement_type
{
  RegularDivision,
  LongestEdge,
  NewestVertex
};
  

/*! @brief Struct to store problem data
 *  @param model Model to build force term: ChanVese, RegionGaussian, RegionNonParametric
 *  @param h_step Initial space length to build a new mesh
 *  @param lambda Parameter of the model
 *  @param epsilon Parameter for numerical regularization
 *  @param rho Parameter for numerical regularization
 *  @param mu Parameter of the model
 *  @param tau Time step
 *  @param stop_criterion Tolerance for stopping evolutionary problem
 */
struct pb_data
{ 
  enum force_type model;
  
  double h_step;
  
  double lambda;
  
  double epsilon;
  
  double rho;
  
  double mu;
  
  double tau;
  
  double stop_criterion;
  
};


/*! @brief Struct to store refinement data
 *  @param initial_refinement Boolean value - Initial refinement? Yes/No
 *  @param initial_tolerance Tolerance for initial mesh refinement
 *  @param max_iteration Maximum number of iteration for initial refinement
 *  @param loop_refinement Boolean value - Refinement at every time step? Yes/No
 *  @param loop_tolerance Tolerance for mesh refinement at every time step
 *  @param nu Parameter for GERS marking strategy
 *  @param theta_star Parameter for GERS marking strategy
 *  @param method Refinement method: RegularDivision, LongestEdge, NewestVertex
 */
struct ref_data
{
  bool initial_refinement;
  
  double initial_tolerance;
  
  int max_iteration;
  
  bool loop_refinement;
  
  double loop_tolerance;
  
  double nu;
  
  double theta_star;
  
  enum refinement_type method;
  
};
  

/*! @brief Function that reads data from file and stores them 
 *  @param inputpath Path for data file
 *  @param path Path for image file
 *  @param res Path to save results
 *  @param pb Struct to store problem data
 *  @param rf Struct to store refinement data 
 */
void readData(const std::string & inputpath, std::string & path, std::string & res, struct pb_data & pb, struct ref_data & rf);


} //namespace


#include "readData.h"
#include "triangulation.h"
#include "initials.h"
#include "histogram.h"
#include "force.h"
#include "problem.h"
#include "posterior.h"
#include "refinement.h"
#include "adapt_mesh.h"
#include "printer.h"


#endif 