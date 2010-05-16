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

#ifndef GAS_FUNCTIONAL_FUNCTION_H
#define GAS_FUNCTIONAL_FUNCTION_H

#include "../gas"

namespace gas { namespace functional {

// Tipo della funzione
template <int d_, typename derived_> class function;

/*! @brief Funzione di una sola variabile */
template <typename derived_>
class function<1, derived_>
{

public:
	inline double operator() (double const x) const {
		return (*static_cast<derived_ const *>(this))(x);
	}

};

/*! @brief Funzione di due variabili */
template <typename derived_>
class function<2, derived_>
{

public:
	inline double operator() (double const x, double const y) const {
		return (*static_cast<derived_ const *>(this))(x, y);
	}

};

/*! @brief Funzione di tre variabili */
template <typename derived_>
class function<3, derived_>
{

public:
	inline double operator() (double const x, double const y, double const z) const {
		return (*static_cast<derived_ const *>(this))(x, y, z);
	}

};

// Tipi delle espressioni tra funzioni
template <int d_, typename operand_, typename operator_> class function_expression;
template <int d_, typename left_, typename right_, typename operator_> class function_binary;


/*! @brief Espressione di una variabile */
template <typename operand_, typename operator_>
class function_expression<1, operand_, operator_>
	: public function<1, function_expression<1u, operand_, operator_> >
{

public:
	inline function_expression (operand_ const & o): o_(o) {
	}

public:
	inline double operator() (double const x) const {
		return operator_::eval(o_(x));
	}

private:
	operand_ const & o_;

};

/*! @brief Espressione di due variabili */
template <typename operand_, typename operator_>
class function_expression<2, operand_, operator_>
	: public function<2, function_expression<2u, operand_, operator_> >
{

public:
	inline function_expression (operand_ const & o): o_(o) {
	}

public:
	inline double operator() (double const x, double const y) const {
		return operator_::eval(o_(x, y));
	}

private:
	operand_ const & o_;

};

/*! @brief Espressione di tre variabili */
template <typename operand_, typename operator_>
class function_expression<3, operand_, operator_>
	: public function<3, function_expression<3u, operand_, operator_> >
{

public:
	inline function_expression (operand_ const & o): o_(o) {
	}

public:
	inline double operator() (double const x, double const y, double const z) const {
		return operator_::eval(o_(x, y, z));
	}

private:
	operand_ const & o_;

};

/*! @brief Espressione binaria di una variabile */
template <typename left_, typename right_, typename operator_>
class function_binary<1, left_, right_, operator_>
	: public function<1, function_binary<1, left_, right_, operator_> >
{

public:
	inline function_binary (left_ const & l, right_ const & r): l_(l), r_(r) {
	}

public:
	inline double operator() (double const x) const {
		return operator_::eval(l_(x), r_(x));
	}

private:
	left_ const & l_;
	right_ const & r_;

};

/*! @brief Espressione binaria di due variabili */
template <typename left_, typename right_, typename operator_>
class function_binary<2, left_, right_, operator_>
	: public function<2, function_binary<2, left_, right_, operator_> >
{

public:
	inline function_binary (left_ const & l, right_ const & r): l_(l), r_(r) {
	}

public:
	inline double operator() (double const x, double const y) const {
		return operator_::eval(l_(x, y), r_(x, y));
	}

private:
	left_ const & l_;
	right_ const & r_;

};

/*! @brief Espressione binaria di tre variabili */
template <typename left_, typename right_, typename operator_>
class function_binary<3, left_, right_, operator_>
	: public function<3, function_binary<3, left_, right_, operator_> >
{

public:
	inline function_binary (left_ const & l, right_ const & r): l_(l), r_(r) {
	}

public:
	inline double operator() (double const x, double const y, double const z) const {
		return operator_::eval(l_(x, y, z), r_(x, y, z));
	}

private:
	left_ const & l_;
	right_ const & r_;

};

// definizione degli operatori

template <int d_, typename left_, typename right_>
inline function_binary<d_, function<d_, left_>, function<d_, right_>, gas::add>
operator+ (function<d_, left_> const & l, function<d_, right_> const & r)
{
	return function_binary<d_, function<d_, left_>, function<d_, right_>, gas::add>(l, r);
}


template <int d_, typename left_, typename right_>
inline function_binary<d_, function<d_, left_>, function<d_, right_>, gas::sub>
operator- (function<d_, left_> const & l, function<d_, right_> const & r)
{
	return function_binary<d_, function<d_, left_>, function<d_, right_>, gas::sub>(l, r);
}


template <int d_, typename left_, typename right_>
inline function_binary<d_, function<d_, left_>, function<d_, right_>, gas::mul>
operator* (function<d_, left_> const & l, function<d_, right_> const & r)
{
	return function_binary<d_, function<d_, left_>, function<d_, right_>, gas::mul>(l, r);
}


template <int d_, typename left_, typename right_>
inline function_binary<d_, function<d_, left_>, function<d_, right_>, gas::div>
operator/ (function<d_, left_> const & l, function<d_, right_> const & r)
{
	return function_binary<d_, function<d_, left_>, function<d_, right_>, gas::div>(l, r);
}

// TODO moltiplicazione per una costante

} }

#endif // GAS_FUNCTIONAL_FUNCTION_H
