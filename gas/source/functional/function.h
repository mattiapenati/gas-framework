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
 * @file function.h
 * @brief Classes and functions for expression templates on functions
 */

#ifndef _gas_functional_function_
#define _gas_functional_function_

#include "../gas"

namespace gas { namespace functional {

template <unsigned int d_, typename derived_>
class function;

template <typename derived_>
class function<1u, derived_> {

public:
	inline double operator() (double const & x) const {
		return static_cast<derived_ const *>(this)->operator()(x);
	}

};

template <typename derived_>
class function<2u, derived_> {

public:
	inline double operator() (double const & x, double const & y) const {
		return static_cast<derived_ const *>(this)->operator()(x, y);
	}

};

template <typename derived_>
class function<3u, derived_> {

public:
	inline double operator() (double const & x, double const & y, double const & z) const {
		return static_cast<derived_ const *>(this)->operator()(x, y, z);
	}

};

template <typename operand_> class base_function_expression;
template <unsigned int d_, typename operand_, typename operator_> class function_expression;
template <unsigned int d_, typename left_, typename right_, typename operator_> class function_binary;


template <typename operand_, typename operator_>
class function_expression<1u, operand_, operator_>: public function<1u, function_expression<1u, operand_, operator_> >  {

public:
	inline function_expression (operand_ const & o): o_(o) {
	}

public:
	inline double operator() (double const & x) const {
		return operator_::eval(o_(x));
	}

private:
	operand_ const & o_;

};

template <typename operand_, typename operator_>
class function_expression<2u, operand_, operator_>: public function<2u, function_expression<2u, operand_, operator_> >  {

public:
	inline function_expression (operand_ const & o): o_(o) {
	}

public:
	inline double operator() (double const & x, double const & y) const {
		return operator_::eval(o_(x, y));
	}

private:
	operand_ const & o_;

};

template <typename operand_, typename operator_>
class function_expression<3u, operand_, operator_>: public function<3u, function_expression<3u, operand_, operator_> > {

public:
	inline function_expression (operand_ const & o): o_(o) {
	}

public:
	inline double operator() (double const & x, double const & y, double const & z) const {
		return operator_::eval(o_(x, y, z));
	}

private:
	operand_ const & o_;

};

template <typename left_, typename right_, typename operator_>
class function_binary<1u, left_, right_, operator_> {

public:
	inline function_binary (left_ const & l, right_ const & r): l_(l), r_(r) {
	}

public:
	inline double operator() (double const & x) const {
		return operator_::eval(l_(x), r_(x));
	}

private:
	left_ const & l_;
	right_ const & r_;

};

template <typename left_, typename right_, typename operator_>
class function_binary<2u, left_, right_, operator_> {

public:
	inline function_binary (left_ const & l, right_ const & r): l_(l), r_(r) {
	}

public:
	inline double operator() (double const & x, double const & y) const {
		return operator_::eval(l_(x, y), r_(x, y));
	}

private:
	left_ const & l_;
	right_ const & r_;

};

template <typename left_, typename right_, typename operator_>
class function_binary<3u, left_, right_, operator_> {

public:
	inline function_binary (left_ const & l, right_ const & r): l_(l), r_(r) {
	}

public:
	inline double operator() (double const & x, double const & y, double const & z) const {
		return operator_::eval(l_(x, y, z), r_(x, y, z));
	}

private:
	left_ const & l_;
	right_ const & r_;

};

template <unsigned int d_, typename derived_left_, typename derived_right_>
inline function_expression<d_, function_binary<d_, function<d_, derived_left_>, function<d_, derived_right_>, gas::add>, gas::id>
operator+ (function<d_, derived_left_> const & l, function<d_, derived_right_> const & r) {
	typedef function_binary<d_, function<d_, derived_left_>, function<d_, derived_right_>, gas::add> binary_t;
	typedef function_expression<d_, binary_t, gas::id> expression_t;

	return expression_t(binary_t(l, r));
}

template <unsigned int d_, typename derived_left_, typename derived_right_>
inline function_expression<d_, function_binary<d_, function<d_, derived_left_>, function<d_, derived_right_>, gas::sub>, gas::id>
operator- (function<d_, derived_left_> const & l, function<d_, derived_right_> const & r) {
	typedef function_binary<d_, function<d_, derived_left_>, function<d_, derived_right_>, gas::sub> binary_t;
	typedef function_expression<d_, binary_t, gas::id> expression_t;

	return expression_t(binary_t(l, r));
}

template <unsigned int d_, typename derived_left_, typename derived_right_>
inline function_expression<d_, function_binary<d_, function<d_, derived_left_>, function<d_, derived_right_>, gas::mul>, gas::id>
operator* (function<d_, derived_left_> const & l, function<d_, derived_right_> const & r) {
	typedef function_binary<d_, function<d_, derived_left_>, function<d_, derived_right_>, gas::mul> binary_t;
	typedef function_expression<d_, binary_t, gas::id> expression_t;

	return expression_t(binary_t(l, r));
}

template <unsigned int d_, typename derived_left_, typename derived_right_>
inline function_expression<d_, function_binary<d_, function<d_, derived_left_>, function<d_, derived_right_>, gas::div>, gas::id>
operator/ (function<d_, derived_left_> const & l, function<d_, derived_right_> const & r) {
	typedef function_binary<d_, function<d_, derived_left_>, function<d_, derived_right_>, gas::div> binary_t;
	typedef function_expression<d_, binary_t, gas::id> expression_t;

	return expression_t(binary_t(l, r));
}

// TODO moltiplicazione per una costante

} }

#endif // _gas_functional_function_
