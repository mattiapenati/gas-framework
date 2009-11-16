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

#ifndef _gas_svg_
#define _gas_svg_

#include <iostream>
#include <sstream>
#include <stack>
#include <string>

#include "css.h"

namespace gas {

class svg {

public:
	svg (unsigned int const & width, unsigned int const & height);
	~svg ();

	svg & title (char const * tit);
	svg & description (char const * desc);
	svg & style (css & style);

	svg & open_group (char const * id, char const * clas = "");
	svg & close_group ();

	svg & line (int const & x1, int const & y1, int const & x2, int const & y2, char const * clas = "");
	svg & rectangle (int const & x1, int const & y1, int const & x2, int const & y2, char const * clas = "");
	svg & circle (int const & xc, int const & yc, int const & r, char const * clas = "");
	svg & triangle (int const & x1, int const & y1, int const & x2, int const & y2, int const & x3, int const & y3, char const * clas = "");

	svg & text (int const & x, int const & y, char const * text, char const * clas = "");

private:
	std::stack<std::string> tags_;
	std::ostringstream buff_;
	unsigned int width_;
	unsigned int height_;

	void end ();

	friend std::ostream & operator<< (std::ostream & out, svg & image);

};

svg::svg (unsigned int const & width, unsigned int const & height): buff_(), tags_(), width_(width), height_(height) {
	buff_ << "<?xml version=\"1.0\" standalone=\"yes\"?>" << std::endl;
	buff_ << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << std::endl;
	buff_ << "\t\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << std::endl;
	buff_ << "<svg width=\"" << width << "\" height=\"" << height << "\" " <<
		"version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;
	tags_.push(std::string("svg"));
}

svg::~svg () {
	end();
}

svg & svg::title (char const * tit) {
	buff_ << "<title>" << tit << "</title>" << std::endl;
	return *this;
}

svg & svg::description (char const * desc) {
	buff_ << "<desc>" << desc << "</desc>" << std::endl;
	return *this;
}

svg & svg::style (css & style) {
	buff_ << "<defs><style type=\"text/css\"><![CDATA[" << std::endl;
	buff_ << style << std::endl;
	buff_ << "]]></style></defs>" << std::endl;
	return *this;
}

svg & svg::open_group (char const * id, char const * clas){
	buff_ << "<g id=\"" << id << "\"";
	if (clas != "")
		buff_ << " class=\"" << clas << "\"";
	buff_ << ">" << std::endl;
	tags_.push(std::string("g"));
	return *this;
}

svg & svg::close_group () {
	if (!tags_.empty()) {
		buff_ << "</" << tags_.top() << ">" << std::endl;
		tags_.pop();
	}
	return *this;
}

svg & svg::line (int const & x1, int const & y1, int const & x2, int const & y2, char const * clas) {
	int X1(x1); int Y1(height_-y1);
	int X2(x2); int Y2(height_-y2);
	buff_ << "<path d=\"M" << X1 << "," << Y1 << " L" << X2 << "," << Y2 << "\"";
	if (clas != "")
		buff_ << " class=\"" << clas << "\"";
	buff_ << "/>" << std::endl;
	return *this;
}

svg & svg::rectangle (int const & x1, int const & y1, int const & x2, int const & y2, char const * clas) {
	int X(x1); int Y(height_-y2);
	int W(x2-x1); int H(y2-y1);
	buff_ << "<rect x=\"" << X << "\" y=\"" << Y << "\" " <<
		"width=\"" << W << "\" height=\"" << H << "\"";
	if (clas != "")
		buff_ << " class=\"" << clas << "\"";
	buff_ << "/>" << std::endl;
	return *this;
}
svg & svg::circle (int const & xc, int const & yc, int const & r, char const * clas) {
	int X(xc); int Y(height_-yc);
	buff_ << "<circle cx=\"" << X << "\" cy=\"" << Y << "\" " <<
		"r=\"" << r << "\"";
	if (clas != "")
		buff_ << " class=\"" << clas << "\"";
	buff_ << "/>" << std::endl;
	return *this;
}

svg & svg::triangle (int const & x1, int const & y1, int const & x2, int const & y2, int const & x3, int const & y3, char const * clas) {
	buff_ << "<polygon points=\"" <<
			x1 << "," << y1 << " " <<
			x2 << "," << y2 << " " <<
			x3 << "," << y3 << "\"";
	if (clas != "")
		buff_ << " class=\"" << clas << "\"";
	buff_ << " />" << std::endl;
	return *this;
}

svg & svg::text (int const & x, int const & y, char const * text, char const * clas) {
	int X(x); int Y(height_-y);
	buff_ << "<text x=\"" << X << "\" y=\"" << Y << "\"";
	if (clas != "")
		buff_ << " class=\"" << clas << "\"";
	buff_ << ">" << text << "</text>" << std::endl;
	return *this;
}

void svg::end () {
	while (!tags_.empty()) {
		buff_ << "</" << tags_.top() << ">" << std::endl;
		tags_.pop();
	}
}

std::ostream & operator<< (std::ostream & out, svg & image) {
	image.end();
	out << image.buff_.str();
	return out;
}

}

#endif // _gas_svg_
