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

#ifndef _gas_css_
#define _gas_css_

#include <iostream>
#include <sstream>

namespace gas {

class css {

public:
	css ();
	~css ();

	css & new_class (char const * c);
	css & new_id (char const * i);
	css & newTag (char const * t);

	css & property (char const * name, char const * value);
	css & property (char const * name, int const & value);

	css & comment (char const * c);

private:
	bool first_;
	std::ostringstream buff_;

	void end();

	friend std::ostream & operator<< (std::ostream & out, css & style);

};

css::css (): first_(true), buff_() {
}

css::~css () {
	end();
}

css & css::new_class (char const * c) {
	if (first_)
		first_ = false;
	else
		buff_ << "}" << std::endl;
	buff_ << "." << c << " {" << std::endl;
	return *this;
}

css & css::new_id (char const * i) {
	if (first_)
		first_ = false;
	else
		buff_ << "}" << std::endl;
	buff_ << "#" << i << " {" << std::endl;
	return *this;
}

css & css::newTag (char const * t) {
	if (first_)
		first_ = false;
	else
		buff_ << "}" << std::endl;
	buff_ << t << " {" << std::endl;
	return *this;
}

css & css::property (char const * name, char const * value) {
	buff_ << "\t" << name << ": " << value << ";" << std::endl;
	return *this;
}

css & css::property (char const * name, int const & value) {
	buff_ << "\t" << name << ": " << value << ";" << std::endl;
	return *this;
}

css & css::comment (char const * c) {
	buff_ << "// " << c << std::endl;
	return *this;
}

void css::end () {
	if (!first_) {
		buff_ << "}" << std::endl;
		first_ = true;
	}
}

std::ostream & operator<< (std::ostream & out, css & style) {
	style.end();
	out << style.buff_.str();
	return out;
}

}

#endif // _gas_css_
