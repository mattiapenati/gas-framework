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
 * @file method.h
 * @brief Some generical implementation of quadrature formulae
 */

#ifndef GAS_NUMERICAL_QUADRATURE_METHOD_H
#define GAS_NUMERICAL_QUADRATURE_METHOD_H

#include <Eigen/Core>

#include <cmath>
#include "../../gas"
#include "../../geometry/unit/unit"

namespace gas { namespace numerical { namespace quadrature {

/*!
 * @brief A general implementation for 1d formulae, with a list of all nodes
 *        and weights
 * @param unit_ The unit geometry
 * @param nodes_ The number of nodes
 * @param data_ The struct that contains the nodes and weights
 */
template <typename unit_, int nodes_, typename data_>
class method_1 {

private:
	/*! @brief The dimension of unit_t */
	static int const d_ = gas::geometry::unit::info<unit_>::d;

public:
	/*!
	 * @brief Initialize the method, copying the value of nodes and weights
	 */
	method_1 () {
		gas_static_assert((d_ == 1), This_is_a_1d_quadrature_formula);
		gas_rangeu(i, nodes_)
			m_x(i, 0) = data_::x_[i];
		gas_rangeu(i, nodes_)
			m_w(i) = data_::w_[i];
	}

private:
	/*!
	 * @brief Map the nodes and weight
	 * @param m The map used to recalculate nodes and weights
	 */
	template <typename map_>
	void map (map_ const & m) {
		gas_rangeu(i, nodes_) {
			m_x(i, 0) = m.x(data_::x_[i]);
			m_w(i) = std::abs(m.det(data_::x_[i])) * data_::w_[i];
		}
	}

	/*!
	 * @brief Apply the quadrature formula with current nodes and weights
	 * @param f Any kind of object with the possibility to call f(...)
	 * @return The computed numerical integral
	 */
	template <typename function_>
	double apply (function_ const & f) {
		Eigen::Matrix<double, nodes_, 1> m_f;
		gas_rangeu(i, nodes_)
			m_f(i) = f(m_x(i, 0));
		return m_f.dot(m_w);
	}
	
	
private:
	/*! @brief The list of nodes of quadrature formula */
	Eigen::Matrix<double, nodes_, d_> m_x;

	/*! @brief The list of weights of quadrature formula */
	Eigen::Matrix<double, nodes_, 1> m_w;

	template <typename method__, typename map__> friend class formula;

};


/*!
 * @brief A general implementation for 2d formulae, with a list of all nodes
 *        and weights
 * @param unit_ The unit geometry
 * @param nodes_ The number of nodes
 * @param data_ The struct that contains the nodes and weights
 */
template <typename unit_, int nodes_, typename data_>
class method_2 {

private:
	/*! @brief The dimension of unit_t */
	static int const d_ = gas::geometry::unit::info<unit_>::d;

public:
	/*!
	 * @brief Initialize the method, copying the value of nodes and weights
	 */
	method_2 () {
		gas_static_assert((d_ == 2), This_is_a_2d_quadrature_formula);
		gas_rangeu(i, nodes_) {
			m_x(i, 0) = data_::x_[i];
			m_x(i, 1) = data_::y_[i];
		}
		gas_rangeu(i, nodes_)
			m_w(i) = data_::w_[i];
	}

private:
	/*!
	 * @brief Map the nodes and weight
	 * @param m The map used to recalculate nodes and weights
	 */
	template <typename map_>
	void map (map_ const & m) {
		gas_rangeu(i, nodes_) {
			m_x(i, 0) = m.x(data_::x_[i], data_::y_[i]);
			m_x(i, 1) = m.y(data_::x_[i], data_::y_[i]);
			m_w(i) = std::abs(m.det(data_::x_[i], data_::y_[i])) * data_::w_[i];
		}
	}

	/*!
	 * @brief Apply the quadrature formula with current nodes and weights
	 * @param f Any kind of object with the possibility to call f(...)
	 * @return The computed numerical integral
	 */
	template <typename function_>
	double apply (function_ const & f) {
		Eigen::Matrix<double, nodes_, 1> m_f;
		gas_rangeu(i, nodes_)
			m_f(i) = f(m_x(i, 0), m_x(i, 1));
		return m_f.dot(m_w);
	}

private:
	/*! @brief The list of nodes of quadrature formula */
	Eigen::Matrix<double, nodes_, d_> m_x;

	/*! @brief The list of weights of quadrature formula */
	Eigen::Matrix<double, nodes_, 1> m_w;

	template <typename method__, typename map__> friend class formula;

};


/*!
 * @brief A general implementation for 3d formulae, with a list of all nodes
 *        and weights
 * @param unit_ The unit geometry
 * @param nodes_ The number of nodes
 * @param data_ The struct that contains the nodes and weights
 */
template <typename unit_, int nodes_, typename data_>
class method_3 {

private:
	/*! @brief The dimension of unit_t */
	static int const d_ = gas::geometry::unit::info<unit_>::d;

public:
	/*!
	 * @brief Initialize the method, copying the value of nodes and weights
	 */
	method_3 () {
		gas_static_assert((d_ == 3), This_is_a_3d_quadrature_formula);
		gas_rangeu(i, nodes_) {
			m_x(i, 0) = data_::x_[i];
			m_x(i, 1) = data_::y_[i];
			m_x(i, 2) = data_::z_[i];
		}
		gas_rangeu(i, nodes_)
			m_w(i) = data_::w_[i];
	}

private:
	/*!
	 * @brief Map the nodes and weight
	 * @param m The map used to recalculate nodes and weights
	 */
	template <typename map_>
	void map (map_ const & m) {
		gas_rangeu(i, nodes_) {
			m_x(i, 0) = m.x(data_::x_[i], data_::y_[i], data_::z_[i]);
			m_x(i, 1) = m.y(data_::x_[i], data_::y_[i], data_::z_[i]);
			m_x(i, 2) = m.z(data_::x_[i], data_::y_[i], data_::z_[i]);
			m_w(i) = std::abs(m.det(data_::x_[i], data_::y_[i], data_::z_[i])) * data_::w_[i];
		}
	}

	/*!
	 * @brief Apply the quadrature formula with current nodes and weights
	 * @param f Any kind of object with the possibility to call f(...)
	 * @return The computed numerical integral
	 */
	template <typename function_>
	double apply (function_ const & f) {
		Eigen::Matrix<double, nodes_ ,1> m_f;
		gas_rangeu(i, nodes_)
			m_f(i) = f(m_x(i, 0), m_x(i, 1), m_x(i, 2));
		return m_f.dot(m_w);
	}

private:
	/*! @brief The list of nodes of quadrature formula */
	Eigen::Matrix<double, nodes_, d_> m_x;

	/*! @brief The list of weights of quadrature formula */
	Eigen::Matrix<double, nodes_, 1> m_w;

	template <typename method__, typename map__> friend class formula;

};

} } }

#endif // GAS_NUMERICAL_QUADRATURE_METHOD_H
