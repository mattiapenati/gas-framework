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
 * @file readData.h
 * @brief Routine to read data from file
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_READ_H
#define _IMG_SEG_READ_H


namespace img_seg
{


/*! @brief Function that reads a file and stores the data for the execution of the program
 *  @param inputpath Path for the file containing the data
 *  @param path Path for the image to be processed
 *  @param res Path to the folder to store the results
 *  @param pb Struct containing the problem data
 *  @param rf Struct containing the refinement data
 */  
void readData(const std::string & inputpath, std::string & path, std::string & res, struct pb_data & pb, struct ref_data & rf)
{
    libconfig::Config cfg;

    // Read file
    try
    {
        cfg.readFile(inputpath.c_str());
    }
    catch (const libconfig::FileIOException &fioex)
    {
        throw "I/O error while reading file.";
    }
    catch (const libconfig::ParseException &pex)
    {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
        << " - " << pex.getError() << std::endl;
        throw "Abort";
    }

    const libconfig::Setting& root = cfg.getRoot();
    
    // Look for "file_path"
    if (root.lookupValue("file_path",path))
    {
        std::cout<< "File path: " << path << std::endl << std::endl;
    }
    else
    {
        throw "No 'file_path' setting in configuration file.";
    }
     
    // Look for "result_path"
    if (root.lookupValue("result_path",res))
    {
        std::cout<< "Result path: " << res << std::endl << std::endl;
    }
    else
    {
        throw "No 'result_path' setting in configuration file.";
    }
   
    // Problem parameters //
    if (!(root.exists("problem")))
    {
        throw "No 'problem' setting in configuration file.";
    }
    const libconfig::Setting& problem = root["problem"];

    std::cout << "-- Problem general settings --" << std::endl;
     
    // Look for  "model" setting
    std::string model_string;
    std::map<std::string,enum force_type> strtomod;
    
    strtomod["chanvese"] = ChanVese;
    strtomod["ChanVese"] = ChanVese;
    strtomod["CHANVESE"] = ChanVese;    
    strtomod["cv"] = ChanVese;
    strtomod["CV"] = ChanVese;
    
    strtomod["regiongaussian"] = RegionGaussian;
    strtomod["RegionGaussian"] = RegionGaussian;
    strtomod["REGIONGAUSSIAN"] = RegionGaussian;    
    strtomod["gaussian"] = RegionGaussian;
    strtomod["Gaussian"] = RegionGaussian;    
    strtomod["GAUSSIAN"] = RegionGaussian;        
    strtomod["g"] = RegionGaussian;
    strtomod["G"] = RegionGaussian;        
    
    strtomod["regionnonparametric"] = RegionNonParametric;
    strtomod["RegionNonParametric"] = RegionNonParametric;
    strtomod["REGIONNONPARAMETRIC"] = RegionNonParametric;    
    strtomod["nonparametric"] = RegionNonParametric;
    strtomod["NonParametric"] = RegionNonParametric;    
    strtomod["NONPARAMETRIC"] = RegionNonParametric;
    strtomod["parzen"] = RegionNonParametric;
    strtomod["Parzen"] = RegionNonParametric;            
    strtomod["PARZEN"] = RegionNonParametric;            
    strtomod["p"] = RegionNonParametric;
    strtomod["P"] = RegionNonParametric;            
 
    if (problem.lookupValue("model", model_string))
    {
        pb.model = strtomod[model_string];
        std::cout<< "Mathematical model: " << model_string << std::endl;
    }
    else
    {
        throw "No 'model' setting in configuration file.";
    }
    
    // Look for "h_step" setting
    if (problem.lookupValue("h_step", pb.h_step))
    {
        std::cout<< "h_step: " << pb.h_step << std::endl;
    }
    else
    {
        throw "No 'h_step' setting in configuration file.";
    }
    
    // Look for "lambda" setting
    if (problem.lookupValue("lambda", pb.lambda))
    {
        std::cout<< "lambda: " << pb.lambda << std::endl;
    }
    else
    {
        throw "No 'lambda' setting in configuration file.";
    }

    // Look for "epsilon" setting
    if (problem.lookupValue("epsilon", pb.epsilon))
    {
        std::cout<< "epsilon: " << pb.epsilon << std::endl;
    }
    else
    {
        throw "No 'epsilon' setting in configuration file.";
    }
   
    // Look for "rho" setting
    if (problem.lookupValue("rho", pb.rho))
    {
        std::cout<< "rho: " << pb.rho << std::endl;
    }
    else
    {
        throw "No 'rho' setting in configuration file.";
    }   
   
    // Look for "mu" setting
    if (problem.lookupValue("mu", pb.mu))
    {
        std::cout<< "mu: " << pb.mu << std::endl;
    }
    else
    {
        throw "No 'mu' setting in configuration file.";
    }   
  
    // Look for "tau" setting
    if (problem.lookupValue("tau", pb.tau))
    {
        std::cout<< "tau: " << pb.tau << std::endl;
    }
    else
    {
        throw "No 'tau' setting in configuration file.";
    }     
  
    // Look for "stop_criterion" setting
    if (problem.lookupValue("stop_criterion", pb.stop_criterion))
    {
        std::cout<< "Stopping criterion: " << pb.stop_criterion << std::endl;
    }
    else
    {
        throw "No 'stop_criterion' setting in configuration file.";
    }      
    
    std::cout << std::endl;
    
    // Refinement parameters //
    if (!(root.exists("refinement")))
    {
        throw "No 'refinement' setting in configuration file.";
    }
    const libconfig::Setting& refinement = root["refinement"];

    std::cout << "-- Mesh adaptivity settings --" << std::endl;
        
    // Look for "initial_refinement" setting
    if (refinement.lookupValue("initial_refinement", rf.initial_refinement))
    {
	std::string s;
	if(rf.initial_refinement) { s = "yes"; }
	else { s = "no"; }
	
        std::cout<< "Initial refinement: " << s << std::endl;
    }
    else
    {
        throw "No 'initial_refinement' setting in configuration file.";
    }
    
    // Look for "initial_tolerance" setting
    if (refinement.lookupValue("initial_tolerance", rf.initial_tolerance))
    {
	if(rf.initial_refinement)
	{
	  std::cout<< "Initial tolerance: " << rf.initial_tolerance << std::endl;
	}
    }
    else
    {
        throw "No 'initial_tolerance' setting in configuration file.";
    }
    
    // Look for "max_iteration" setting
    if (refinement.lookupValue("max_iteration", rf.max_iteration))
    {
	if(rf.initial_refinement)
	{      
	  std::cout<< "Max number of initial refinement iterations: " << rf.max_iteration << std::endl;
	}
    }
    else
    {
        throw "No 'max_iteration' setting in configuration file.";
    }   
    
    // Look for "loop_refinement" setting
    if (refinement.lookupValue("loop_refinement", rf.loop_refinement))
    {
      	std::string s;
	if(rf.loop_refinement) { s = "yes"; }
	else { s = "no"; }
      
        std::cout<< "Loop refinement: " << s << std::endl;
    }
    else
    {
        throw "No 'loop_refinement' setting in configuration file.";
    }     
  
    // Look for "loop_tolerance" setting
    if (refinement.lookupValue("loop_tolerance", rf.loop_tolerance))
    {
	if(rf.loop_refinement)
	{      
	  std::cout<< "Loop tolerance: " << rf.loop_tolerance << std::endl;
	}
    }
    else
    {
        throw "No 'loop_tolerance' setting in configuration file.";
    }   
  
     // Look for "nu" setting
    if (refinement.lookupValue("nu", rf.nu))
    {
	if( (rf.initial_refinement) || (rf.loop_refinement) )
	{      
	  std::cout<< "nu: " << rf.nu << std::endl;
	}
    }
    else
    {
        throw "No 'nu' setting in configuration file.";
    }     
      
    // Look for "theta_star" setting
    if (refinement.lookupValue("theta_star", rf.theta_star))
    {
	if( (rf.initial_refinement) || (rf.loop_refinement) )
	{            
	  std::cout<< "theta_star: " << rf.theta_star << std::endl;
	}
    }
    else
    {
        throw "No 'theta_star' setting in configuration file.";
    }          
    
    // Look for "method" setting 
    std::string method_string;
    std::map<std::string,enum refinement_type> strtoref;
    
    strtoref["regulardivision"] = RegularDivision;
    strtoref["RegularDivision"] = RegularDivision;
    strtoref["REGULARDIVISION"] = RegularDivision;    
    strtoref["redgreenblue"] = RegularDivision;
    strtoref["RedGreenBlue"] = RegularDivision;
    strtoref["REDGEENBLUE"] = RegularDivision;
    strtoref["rd"] = RegularDivision;
    strtoref["RD"] = RegularDivision;
    strtoref["rgb"] = RegularDivision;
    strtoref["RGB"] = RegularDivision;        
    
    strtoref["longestedge"] = LongestEdge;
    strtoref["LongestEdge"] = LongestEdge;
    strtoref["LONGESTEDGE"] = LongestEdge;    
    strtoref["le"] = LongestEdge;
    strtoref["LE"] = LongestEdge;    
   
    strtoref["newestvertex"] = NewestVertex;
    strtoref["NewestVertex"] = NewestVertex;
    strtoref["NEWESTVERTEX"] = NewestVertex;    
    strtoref["nv"] = NewestVertex;
    strtoref["NV"] = NewestVertex;    
 
    if (refinement.lookupValue("method", method_string))
    {
        rf.method = strtoref[method_string];
	if( (rf.initial_refinement) || (rf.loop_refinement) )
	{      
	  std::cout<< "Refinement strategy: " << method_string << std::endl;
	}
    }
    else
    {
        throw "No 'method' setting in configuration file.";
    }
    
    std::cout << std::endl;
}


} //namespace


#endif