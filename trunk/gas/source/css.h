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

#ifndef GAS_CSS_H
#define GAS_CSS_H

#include <iostream>
#include <sstream>

namespace gas {

/*! @brief A class that help the creation of CSS file */
class css {

public:
	/*! @brief The constructor */
	css ();
	/*! @brief The destructor */
	~css ();

	/*!
	 * @brief Define a new class
	 * @param c The class name
	 * @return The reference to the current object
	 */
	css & new_class (char const * c);
	/*!
	 * @brief Define a new id
	 * @param i The id name
	 * @return The reference to the current object
	 */
	css & new_id (char const * i);
	/*!
	 * @brief Define a new tag
	 * @param i The tag name
	 * @return The reference to the current object
	 */
	css & newTag (char const * t);

	/*!
	 * @brief Define a new property in the current selector
	 * @param name The property name
	 * @param value The property value
	 * @return The reference to the current object
	 */
	css & property (char const * name, char const * value);
	/*!
	 * @brief Define a new property in the current selector
	 * @param name The property name
	 * @param value The property value
	 * @return The reference to the current object
	 */
	css & property (char const * name, int const value);

	/*!
	 * @brief Create a one line comment in the CSS
	 * @param c The comment
	 * @return The reference to the current object
	 */
	css & comment (char const * c);

private:
	/* check if this is the first selector */
	bool m_first;
	/* the stile content */
	std::ostringstream m_buff;

	/* close the current selector */
	void end();

	/* print out the style to an output stream */
	friend std::ostream & operator<< (std::ostream & out, css & style);

};

css::css (): m_first(true), m_buff() {
}

css::~css () {
	end();
}

css & css::new_class (char const * c) {
	if (m_first)
		m_first = false;
	else
		m_buff << "}" << std::endl;
	m_buff << "." << c << " {" << std::endl;
	return *this;
}

css & css::new_id (char const * i) {
	if (m_first)
		m_first = false;
	else
		m_buff << "}" << std::endl;
	m_buff << "#" << i << " {" << std::endl;
	return *this;
}

css & css::newTag (char const * t) {
	if (m_first)
		m_first = false;
	else
		m_buff << "}" << std::endl;
	m_buff << t << " {" << std::endl;
	return *this;
}

css & css::property (char const * name, char const * value) {
	m_buff << "\t" << name << ": " << value << ";" << std::endl;
	return *this;
}

css & css::property (char const * name, int const value) {
	m_buff << "\t" << name << ": " << value << ";" << std::endl;
	return *this;
}

css & css::comment (char const * c) {
	m_buff << "// " << c << std::endl;
	return *this;
}

void css::end () {
	if (!m_first) {
		m_buff << "}" << std::endl;
		m_first = true;
	}
}

std::ostream & operator<< (std::ostream & out, css & style) {
	style.end();
	out << style.m_buff.str();
	return out;
}

}

#endif // GAS_CSS_H
