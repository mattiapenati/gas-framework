/*
 * Copyright (c) 2009, Politecnico di Milano
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
 * @file chrono.h
 * @brief A chronometer to profile your code
 */

#ifndef _gas_chrono_
#define _gas_chrono_

#include <ctime>
#include <iostream>

namespace gas {

/*! @brief A chronometer to profile your applications */
class chrono {

public:
	/*! @brief Create a new chronometer */
	inline chrono (): time_start_(0), time_elapsed_(0) {
	}

	/*! @brief Start the chronometer */
	inline chrono & start () {
		time_elapsed_ = 0;
		time_start_ = std::clock();
		return * this;
	}

	/*! @brief Stop the chronometer */
	inline chrono & stop () {
		time_elapsed_ += (std::clock()-time_start_);
		return *this;
	}

	/*! @brief Pause the chronometer */
	inline chrono & pause () {
		return stop();
	}

	/*! @brief Restart the chronometer */
	inline chrono & restart () {
		time_start_ = std::clock();
		return *this;
	}

	/*! @brief Give the elapsed time */
	inline std::clock_t elapsed () const {
		return time_elapsed_;
	}

	/*! @brief Give the elapsed time in seconds  */
	inline double elapsed_s () const {
		return time_elapsed_ / CLOCKS_PER_SEC;
	}

private:
	/*! @brief The start time */
	std::clock_t time_start_;

	/*! @brief The elapsed time */
	std::clock_t time_elapsed_;

};

/*! @brief The overloading of operator<< to print a chrono objec */
inline std::ostream & operator<< (std::ostream & out, chrono const & c) {
	return (out << c.elapsed_s() << "s");
}

}

#endif // _gas_chrono_
