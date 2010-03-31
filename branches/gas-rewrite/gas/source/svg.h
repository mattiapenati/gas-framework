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

#ifndef GAS_SVG_H
#define GAS_SVG_H

#include <iostream>
#include <sstream>
#include <stack>
#include <string>
#include <cstring>

#include "css.h"

namespace gas {

/*! @brief A class that help the creation of SVG image file */
class svg {

public:
	/*!
	 * @brief The constructor
	 * @param width The width of the image
	 * @param height The height of the image
	 */
	svg (int const width, int const height);
	/*! @brief The destructor */
	~svg ();

	/*!
	 * Set the title of the document (not the name of the file)
	 * @param tit The title of the document
	 * @return The reference to the current object
	 */
	svg & title (char const * tit);
	/*!
	 * Set the description of the document (not the name of the file)
	 * @param desc The description of the document
	 * @return The reference to the current object
	 */
	svg & description (char const * desc);
	/*!
	 * Set the style of the document (not the name of the file)
	 * @param style The style of the document
	 * @return The reference to the current object
	 */
	svg & style (css & style);

	/*!
	 * Open a new group (the tag <g>)
	 * @param id The id of the tag
	 * @param clas The class of the tag
	 * @return The reference to the current object
	 */
	svg & open_group (char const * id, char const * clas = "\0");
	/*!
	 * Close the current opened group
	 * @return The reference to the current object
	 */
	svg & close_group ();

	svg & line (int const x1, int const y1, int const x2, int const y2, char const * clas = "\0");
	svg & rectangle (int const x1, int const y1, int const x2, int const y2, char const * clas = "\0");
	svg & circle (int const xc, int const yc, int const r, char const * clas = "\0");
	svg & triangle (int const x1, int const y1, int const x2, int const y2, int const x3, int const y3, char const * clas = "\0");

	/*!
	 * Insert a text in the image
	 * @param x The first coordinate of text position
	 * @param y The second coordinate of text position
	 * @param text The text
	 * @param clas The class of the tag
	 * @return The reference to the current object
	 */
	svg & text (int const x, int const y, char const * text, char const * clas = "\0");

private:
	/* the stack containing the opened tag */
	std::stack<std::string> m_tags;
	/* the buffer that store the file */
	std::ostringstream m_buff;

	/* the dimension of image */
	int const m_width;
	int const m_height;

	/* close all opened tag */
	void end ();

	friend std::ostream & operator<< (std::ostream & out, svg & image);

};

svg::svg (int const width, int const height)
	: m_tags(),
	  m_buff(),
	  m_width(width),
	  m_height(height)
{
	GAS_PRE((width > 0) and (height > 0));

	m_buff << "<?xml version=\"1.0\" standalone=\"yes\"?>" << std::endl;
	m_buff << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << std::endl;
	m_buff << "\t\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << std::endl;
	m_buff << "<svg width=\"" << m_width << "\" height=\"" << m_height << "\" " <<
		"version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;
	m_tags.push(std::string("svg"));
}

svg::~svg () {
	end();
}

svg & svg::title (char const * tit) {
	m_buff << "<title>" << tit << "</title>" << std::endl;
	return *this;
}

svg & svg::description (char const * desc) {
	m_buff << "<desc>" << desc << "</desc>" << std::endl;
	return *this;
}

svg & svg::style (css & style) {
	m_buff << "<defs><style type=\"text/css\"><![CDATA[" << std::endl;
	m_buff << style << std::endl;
	m_buff << "]]></style></defs>" << std::endl;
	return *this;
}

svg & svg::open_group (char const * id, char const * clas){
	m_buff << "<g id=\"" << id << "\"";
	if (std::strlen(clas))
		m_buff << " class=\"" << clas << "\"";
	m_buff << ">" << std::endl;
	m_tags.push(std::string("g"));
	return *this;
}

svg & svg::close_group () {
	if (!m_tags.empty()) {
		m_buff << "</" << m_tags.top() << ">" << std::endl;
		m_tags.pop();
	}
	return *this;
}

svg & svg::line (int const x1, int const y1, int const x2, int const y2, char const * clas) {
	int const X1(x1); int const Y1(m_height-y1);
	int const X2(x2); int const Y2(m_height-y2);
	m_buff << "<path d=\"M" << X1 << "," << Y1 << " L" << X2 << "," << Y2 << "\"";
	if (std::strlen(clas))
		m_buff << " class=\"" << clas << "\"";
	m_buff << "/>" << std::endl;
	return *this;
}

svg & svg::rectangle (int const x1, int const y1, int const x2, int const y2, char const * clas) {
	int const X(x1); int const Y(m_height-y2);
	int const W(x2-x1); int const H(y2-y1);
	m_buff << "<rect x=\"" << X << "\" y=\"" << Y << "\" " <<
		"width=\"" << W << "\" height=\"" << H << "\"";
	if (std::strlen(clas))
		m_buff << " class=\"" << clas << "\"";
	m_buff << "/>" << std::endl;
	return *this;
}
svg & svg::circle (int const xc, int const yc, int const r, char const * clas) {
	int const X(xc); int const Y(m_height-yc);
	m_buff << "<circle cx=\"" << X << "\" cy=\"" << Y << "\" " <<
		"r=\"" << r << "\"";
	if (std::strlen(clas))
		m_buff << " class=\"" << clas << "\"";
	m_buff << "/>" << std::endl;
	return *this;
}

svg & svg::triangle (int const x1, int const y1, int const x2, int const y2, int const x3, int const y3, char const * clas) {
	m_buff << "<polygon points=\"" <<
			x1 << "," << y1 << " " <<
			x2 << "," << y2 << " " <<
			x3 << "," << y3 << "\"";
	if (std::strlen(clas))
		m_buff << " class=\"" << clas << "\"";
	m_buff << " />" << std::endl;
	return *this;
}

svg & svg::text (int const x, int const y, char const * text, char const * clas) {
	int const X(x); int const Y(m_height-y);
	m_buff << "<text x=\"" << X << "\" y=\"" << Y << "\"";
	if (std::strlen(clas))
		m_buff << " class=\"" << clas << "\"";
	m_buff << ">" << text << "</text>" << std::endl;
	return *this;
}

void svg::end () {
	while (!m_tags.empty()) {
		m_buff << "</" << m_tags.top() << ">" << std::endl;
		m_tags.pop();
	}
}

std::ostream & operator<< (std::ostream & out, svg & image) {
	image.end();
	out << image.m_buff.str();
	return out;
}

}

#endif // GAS_SVG_H
