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
 *    this software without specific prior written permision.
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

#include <Gas/Gas.h>
#include <iostream>
#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>

using namespace Common;
using namespace std;

#define N_SMALL 300000000
#define N_BIG   3000
#define SMALL     10
#define BIG       100000

int main (int argc, char **argv) {

	// Timing
	struct tms start, stop;
	long r1, r2;

	// Constant
	double a = 3.0;
	double t;

	// Standard array
	double vStd[SMALL], wStd[SMALL], xStd[SMALL];
	double *pStd = new double[BIG];
	double *qStd = new double[BIG];
	double *rStd = new double[BIG];
	// Array<N, T>
	Array<double, SMALL> v, w, x;
	Array<double, BIG> p, q, r;

	range(i, 0, SMALL) {
		vStd[i] = double(i);
		wStd[i] = double(i);
		xStd[i] = 0.;
		v[i] = double(i);
		x[i] = 0.;
	}
	w = v;

	range(i, 0, BIG) {
		pStd[i] = double(i);
		qStd[i] = double(i);
		rStd[i] = 0.;
		p[i] = double(i);
		r[i] = 0.;
	}
	q = p;

	cout << "The output must be the same!" << endl << endl;

	cout << "Small array test" << endl;
	cout << "================" << endl;

	// axpy
	r1 = times(&start);
	range(i, 0, N_SMALL) {
		range(j, 0, SMALL) {
			xStd[j] = a * vStd[j] + wStd[j];
		}
		t = vStd[0];
		vStd[0] = vStd[1];
		vStd[1] = vStd[2];
		vStd[2] = t;
	}
	r2 = times(&stop);
	cout << "Output: " << xStd[0] << " ~ " << xStd[1] << " ~ " << xStd[2] << " =time=> ";
	cout << "C-array axpy: " << float(stop.tms_utime - start.tms_utime)/HZ << "s" << endl;

	// axpy
	r1 = times(&start);
	range(i, 0, N_SMALL) {
		range(j, 0, SMALL) {
			x[j] = a * v[j] + w[j];
		}
		t = v[0];
		v[0] = v[1];
		v[1] = v[2];
		v[2] = t;
	}
	r2 = times(&stop);
	cout << "Output: "<< x[0] << " ~ " << x[1] << " ~ " << x[2] << " =time=> ";
	cout << "Array<N, T> axpy: " << float(stop.tms_utime - start.tms_utime)/HZ << "s" << endl;

	cout << endl;
	cout << "Big array test" << endl;
	cout << "==============" << endl;

	// axpy
	r1 = times(&start);
	range(i, 0, N_BIG) {
		range(j, 0, BIG) {
			pStd[j] = a * qStd[j] + rStd[j];
		}
		t = pStd[0];
		pStd[0] = pStd[1];
		pStd[1] = pStd[2];
		pStd[2] = t;
	}
	r2 = times(&stop);
	cout << "Output: " << pStd[0] << " ~ " << pStd[1] << " ~ " << pStd[2] << " =time=> ";
	cout << "C-array axpy: " << float(stop.tms_utime - start.tms_utime)/HZ << "s" << endl;

	// axpy
	r1 = times(&start);
	range(i, 0, N_BIG) {
		range(j, 0, BIG) {
			p[j] = a * q[j] + r[j];
		}
		t = p[0];
		p[0] = p[1];
		p[1] = p[2];
		p[2] = t;
	}
	r2 = times(&stop);
	cout << "Output: "<< p[0] << " ~ " << p[1] << " ~ " << p[2] << " =time=> ";
	cout << "Array<N, T> axpy: " << float(stop.tms_utime - start.tms_utime)/HZ << "s" << endl;

	delete[] pStd;
	delete[] qStd;
	delete[] rStd;

	return 0;
}
