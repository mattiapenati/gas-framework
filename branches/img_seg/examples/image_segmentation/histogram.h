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
 * @file histogram.h
 * @brief Includes the implementation of the histogram for Parzen estimate
 * @author Matteo Giacomini
 * @date 2013
 */


#ifndef _IMG_SEG_HISTOGRAM_H
#define _IMG_SEG_HISTOGRAM_H


namespace img_seg {

  
/*! @brief Class to build the histogram of the intensities starting from an image */  
class histogram
{
  public:
	  /*! @brief Constructor
	   *  @param m Minimum value within the histogram
	   *  @param M Maximum value within the histogram
	   *  @param bin_nbr Number of bins of the histogram
	   */
	  histogram(double const m, double const M, int const bin_nbr);
	  
	  /*! @brief Destructor */
	  ~histogram() {};
	  
	  /*! @brief Return the number of bins */
	  inline int count_bins() const { return values_.size(); }

	  /*! @brief Return the minimum value within the histogram */
	  inline double minValue() const { return min_; }

	  /*! @brief Return the maximum value within the histogram */
	  inline double maxValue() const { return max_; }
	  
	  /*! @brief Return bin size*/
	  inline double binValue() const { return bin_dim_; }
	  
	  /*! @brief Return bin area associated to a given value
	   *  @param val Value to look for and for which computing the area 	   
	   */
	  double get_bin_area(double const val) const;

	  /*! @brief Return true if a value is contained in a bin 
	   *  @param val Value to look for
	   *  @param idx Index of the bin
	   */
	  bool is_contained(double const val, int const idx) const;
	  
	  /*! @brief Locate the bin a value belongs to 
	   *  @param Value to look for
	   */
	  int which_bin(double const val) const;	  
	  
	  /*! @brief Copy the values of the histogram into a vector 
	   *  @param vec Vector to store the histogram in
	   */
	  void get_values(std::vector<double>& vec) const { vec = values_; }
	  
	  /*! @brief Store the values of a vector into the histogram 
	   *  @param vec Vector from which building the histogram
	   */
	  void set_values(const std::vector<double>& vec) { values_ = vec; }
	  
	  /*! @brief Gaussian smoothing 
	   *  @param sigma Standard deviation of the Gaussian kernel
	   */
	  void gaussK_smooth(double const sigma);
	        
  private:
	  /*! @brief Minimum value */
	  double const min_;

	  /*! @brief Maximum value */	  
	  double const max_;

	  /*! @brief Bin dimension */	  
	  double bin_dim_;

	  /*! @brief Histogram values stored in a vector */	  
	  std::vector<double> values_;
};


histogram::histogram(double const m, double const M, int const bin_nbr): min_(m), max_(M)
{
    values_.resize(bin_nbr, 0.0);
    bin_dim_ = double(bin_nbr) / (max_ - min_);
}


double histogram::get_bin_area(double const val) const
{
    int idx = this->which_bin(val);
    
    if (idx < 0 || idx >= this->count_bins()) { return 0; }
    
    return values_[idx] * bin_dim_;	
}


bool histogram::is_contained(double const val, int const idx) const 
{
    if ((val > idx*bin_dim_) && (val < (idx+1)*bin_dim_)) 
      return true;
    else
      return false;
}


int histogram::which_bin(double const val) const 
{
    int k=0;
    
    // Loop over the bins
    for(; k<this->count_bins(); k++)
    {
      if(this->is_contained(val,k)) { break; }
    }
    
    if(k < this->count_bins())
      return k;
    else
    {
      std::cerr << "Value not present in any bin." << std::endl;
      return 0;
    }
}


void histogram::gaussK_smooth(double const sigma)
{
    int dim = this->count_bins();
    
    if (dim > 1)
    {
      // Build a temporary CImg image and copy the histogram values in
      cimg_library::CImg<double> temp(dim, 1, 1, 1);
      std::copy(&values_[0], &values_[0] + dim, temp.data());
      // Gaussian smoothing
      temp.blur(sigma);
      // Copy smoothed values back to the histogram vector
      std::copy(temp.data(), temp.data() + dim, &values_[0]);
    }
    else
    {
      std::cerr << "Just one intensity value present, no gaussian smoothing possible." << std::endl;
      return;
    }
}


} // namespace


#endif 