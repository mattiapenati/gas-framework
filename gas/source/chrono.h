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

#ifndef GAS_CHRONO_H
#define GAS_CHRONO_H

#include <ctime>
#include <iostream>

namespace gas {

/*! @brief A chronometer to profile your applications
 *
 * @code
 * gas::chrono timer;
 * timer.start();
 * // ...
 * timer.pause();
 * // ...
 * timer.start();
 * // ...
 * timer.stop();
 * std::cout << "Elapsed time: " << timer << std::endl;
 * @endcode
 */
class chrono {

public:
	/*! @brief Create a new chronometer */
	inline chrono (): m_time_start(0), m_time_elapsed(0) {
	}

	/*! @brief Start the chronometer */
	inline chrono & start () {
		m_time_elapsed = 0;
		m_time_start = std::clock();
		return * this;
	}

	/*! @brief Stop the chronometer */
	inline chrono & stop () {
		m_time_elapsed += (std::clock()-m_time_start);
		return *this;
	}

	/*! @brief Pause the chronometer */
	inline chrono & pause () {
		return stop();
	}

	/*! @brief Restart the chronometer */
	inline chrono & restart () {
		m_time_start = std::clock();
		return *this;
	}

	/*! @brief Give the elapsed time */
	inline std::clock_t elapsed () const {
		return m_time_elapsed;
	}

	/*! @brief Give the elapsed time in seconds  */
	inline double elapsed_s () const {
		return m_time_elapsed / CLOCKS_PER_SEC;
	}

private:
	/*! @brief The start time */
	std::clock_t m_time_start;

	/*! @brief The elapsed time */
	std::clock_t m_time_elapsed;

};

/*! @brief The overloading of operator<< to print a chrono objec */
inline std::ostream & operator<< (std::ostream & out, chrono const & c) {
	return (out << c.elapsed_s() << "s");
}

}

#endif // GAS_CHRONO_H
