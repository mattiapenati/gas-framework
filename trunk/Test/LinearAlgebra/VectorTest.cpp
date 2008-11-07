/*
 * Copyright (c) 2008, Davide Ferrarese & Mattia Penati
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
 * 3. Neither the name of the <ORGANIZATION> nor the names of its
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
 

#include <iostream>
#include <string>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestFixture.h>
#include <cppunit/ui/text/TestRunner.h>

//#define GAS_USE_META
#include <Gas/Gas.h>

using namespace LinearAlgebra;

namespace Test {
	class VectorTest: public CppUnit::TestFixture {
		private:
			Vector<100, double> v1, v2, v3, v4, v5;
		public:
			void setUp() {
				range(i, 0, 100) {
					v1[i] = double(i)/12.;
					v2[i] = double(i)/12.;
					v3[i] = double(i)/3.;
					v4[i] = double(i)/12. + double(i)/3.;
				}
				v5 = double(0);
			}
			void tearDown() {
			}

			/**
			 * All the test
			 */
			void testEquality() {
				CPPUNIT_ASSERT(v1 == v2);
				CPPUNIT_ASSERT(!(v2 == v3));
			}
			void testCopy() {
				v5 = v4;
				CPPUNIT_ASSERT(v4 == v5);
			}
			void testSelfAddition() {
				v2 += v3;
				CPPUNIT_ASSERT(v2 == v4);
			}
			void testSelfSubtraction() {
				v4 -= v3;
				CPPUNIT_ASSERT(v4 == v2);
			}
			void testAddition() {
				CPPUNIT_ASSERT(v4 == (v2 + v3));
			}

			static CppUnit::Test *suite() {
				CppUnit::TestSuite *s = new CppUnit::TestSuite("VectorTest");
				s->addTest(
					new CppUnit::TestCaller<Test::VectorTest>(
						"testEquality", &Test::VectorTest::testEquality
					)
				);
				s->addTest(
					new CppUnit::TestCaller<Test::VectorTest>(
						"testCopy", &Test::VectorTest::testCopy
					)
				);
				s->addTest(
					new CppUnit::TestCaller<Test::VectorTest>(
						"testSelfAddition", &Test::VectorTest::testSelfAddition
					)
				);
				s->addTest(
					new CppUnit::TestCaller<Test::VectorTest>(
						"testSelfSubtraction", &Test::VectorTest::testSelfSubtraction
					)
				);
				s->addTest(
					new CppUnit::TestCaller<Test::VectorTest>(
						"testAddition", &Test::VectorTest::testAddition
					)
				);
				return s;
			}
	};
}

int main( int argc, char **argv) {
	CppUnit::TextUi::TestRunner runner;
	runner.addTest(Test::VectorTest::suite());
	runner.run();
	return 0;
}
