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
 * @file main.cpp
 * @brief Main program for image segmentation
 * @author Matteo Giacomini
 * @date 2013
 */


#include "img_seg.h"


int main (int argc, char * argv[])
{
    typedef img_seg::triangulation::point_t point_t;

    // Read the parameter passed to the executable file
    const char *input=argv[1];

    std::stringstream ss;
    std::string path;
    std::string res;
    struct img_seg::pb_data data;
    struct img_seg::ref_data data_r;
    
    ss << input;
    
    // Read data from input file and stores in the data structs
    readData(ss.str(), path, res, data, data_r);
  
    // Create a CImg image to be processed
    cimg_library::CImg<double> img(path.c_str());

    int width(img.width());
    int height(img.height());
    double ratio(double(width)/double(height));

    std::cout<<"Dimension of image to be processed: "<<width<<"x"<<height<<" pixels"<<std::endl;

    // Start timer
    gas::chrono timer;
    timer.start();

    // Boundary definition
    std::vector<point_t> boundary;
    double nw, nh;
    double x_scale, y_scale;

    if(width<100 || height<100)
    {
	std::cerr<<"Image too small. Minimum dimensions: 100x100 pixels"<<std::endl;
	return 0;
    }
    else
    {
	if (ratio >=1)
	{
	    nw = 90;
	    nh = nw/ratio;
	    x_scale = 1.0;
	    y_scale = x_scale/ratio;
	}
	else
	{
	    nh = 90;
	    nw = nh*ratio;
	    y_scale = 1.0;
	    x_scale = y_scale*ratio;
	}
    }

    // Geometry information
    int nelem_w(width/nw), nelem_h(height/nh);
    int const N(2*nelem_w+2*nelem_h);
    double H(data.h_step);
    
    // Build boundary
     boundary.reserve(N);
    
    for (int i(0); i < nelem_w; ++i)
	    boundary.push_back(point_t(-x_scale/2.0+(i*nw/width)*x_scale,-y_scale/2.0));
    for(int j(0); j< nelem_h-1; ++j)
	    boundary.push_back(point_t(-x_scale/2.0+(width/width)*x_scale,-y_scale/2.0+(j*nh/height)*y_scale));
    for (int i(0); i < nelem_w; ++i)
	    boundary.push_back(point_t(-x_scale/2.0+((width-i*nw)/width)*x_scale,-y_scale/2.0+(height/height)*y_scale));
    for(int j(0); j< nelem_h; ++j)
	    boundary.push_back(point_t(-x_scale/2.0,-y_scale/2.0+((height-j*nh)/height)*y_scale));

    // Build the triangulation
    img_seg::triangulation mesh(boundary.begin(), boundary.end(), H);

    // Build the initial level-set curve
    img_seg::initial_curve phi0(mesh, x_scale, y_scale);

    // Build the restriction of the initial image to the triangulation 
    img_seg::initial_image u0(img, x_scale, y_scale);

    // Build the problem
    img_seg::problem problem(mesh, phi0, u0, data);

    // Plot initial image
    {	
	std::string rp(res);
	rp += "/initial_image.svg";

	img_seg::svg out(problem);
	std::ofstream out_file;
	out_file.open(rp.c_str());
	out_file << out;
	out_file.close();
    }

    // Refinement information
    std::pair<int, double> n = std::make_pair(0, 0.0);

    // Initial refinement script
    if(data_r.initial_refinement)
    {
	int i(0);
	
	do {
		// Running a refinement step
		problem.update_refinement(true);
	  
		// Build error estimate
		img_seg::posterior stimator(problem, data_r);		
		
		// Run initial adaptivity procedure
		img_seg::AdaptMesh adapt_mesh1(mesh, stimator);
		n = adapt_mesh1.initial_adapt();		

		// Update data
		problem.init_update();

		// Plot refined solution
		{
		    std::stringstream ssI;
		    ssI << res << "/initial_image-refined_" << i << ".svg";

		    img_seg::svg out(problem);
		    std::ofstream out_file(ssI.str().c_str());
		    out_file << out;
		}

		// Refinement log
		std::cout << "Adaptivity iteration # " << i << std::endl;
		std::cout << "Global error: " << n.second << std::endl;		
		std::cout << "Number of added nodes: " << n.first << std::endl;		
		std::cout << "Global number of nodes: " << mesh.nodes() << std::endl;
		std::cout << "Global number of faces: " << mesh.faces() << std::endl;	
		
		++i;
		
	} while( (i < data_r.max_iteration) && (n.second > data_r.initial_tolerance) );
      }

      // We are not running the first step
      if(!problem.not_original()) { problem.setNotOriginal(); }

      // Running a solver step
      problem.update_refinement(false);	
	
      // Plot initial solution 
      {
	  std::string rp(res);
	  rp += "/initial_solution.svg";

	  img_seg::svg out(problem);
	  std::ofstream out_file;
	  out_file.open(rp.c_str());
	  out_file << out;
	  out_file.close();		
      }

      // Initialize global error, evolutionary error and iteration counter for loop execution
      double stop_e(100.0);
      double glob_e(100.0);
      int _m(0);

      // Loop to solve the problem and iteratively refine the triangulation
      while(stop_e > data.stop_criterion)
      {
	  std::cout<<"Iteration #"<<_m<<std::endl;

	  // Solver
	  problem.solve();

	  // Running a solver step
	  problem.update_refinement(false);

	  // Plot the solution
	  {
	      std::stringstream ss1;
	      ss1 << res << "/solution_" << _m << ".svg";

	      img_seg::svg out1(problem);
	      std::ofstream out_file1(ss1.str().c_str());
	      out_file1 << out1;
	  }
			
	  // Update stop criterion	    
	  stop_e = problem.stopping_res();		

	  // Loop refinement script
	  if(data_r.loop_refinement)
	  {
	      // If tolerance isn't fulfilled
	      if(glob_e > data_r.loop_tolerance)
	      {
		// Running adaptivity step
		problem.update_refinement(true);
		
		// Build loop error estimate
		img_seg::posterior stimator(problem, data_r);		
		
		// Run loop adaptivity procedure
		img_seg::AdaptMesh adapt_mesh2(mesh, stimator);
		n = adapt_mesh2.loop_adapt();
		
		// Compute global error
		glob_e = n.second;

		// Update data
		problem.loop_update();
		
		// Plot refined solution
		{	    
		    std::stringstream ss2;
		    ss2 << res << "/solution_" << _m << "-refined.svg";

		    img_seg::svg out2(problem);
		    std::ofstream out_file2(ss2.str().c_str());
		    out_file2 << out2;
		}
		
		// Refinement log
		std::cout << "Global error: " << n.second << std::endl;		
		std::cout << "Number of added nodes: " << n.first << std::endl;		
	      }
	      else
	      {
		// Update data
		problem.loop_update();
		
		// Log
		std::cout << "Global error: " << n.second << std::endl;
		std::cout << "Tolerance fulfilled" << std::endl;
	      }
	  }
	  else
	  {
	      // Update data
	      problem.loop_update();	      
	  }
			
	  // Log
	  std::cout << "Global number of nodes: " << mesh.nodes() << std::endl;
	  std::cout << "Global number of faces: " << mesh.faces() << std::endl;	
	  std::cout << "Non-zeros elements: " << problem.no_zeros() << std::endl;
	  std::cout << "Evolutionary error: " << stop_e << std::endl;				
		  
	  // Export error to data file
	  std::ofstream errT_file("result/error_evolution.dat", std::ios::app);
	  errT_file << _m << "\t" << mesh.faces() << "\t" << stop_e << std::endl;
	  errT_file.close();

	  // Increase counter
	  _m++;

	  // Timer stop
	  timer.stop();

	  // Final time
	  std::cout << "-- Elapsed time: " << timer << std::endl;
      }

      return 0;

}