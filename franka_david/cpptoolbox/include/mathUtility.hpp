/***************************************************************************** 
 *  \author 
 *  	Dr. Fares Abu-Dakka
 *
 *  \version 
 *	V1.0
 *
 *  \date $Date: 2010-2015$
 *
 *  \contact: fabudakk@ing.uc3m.es
 *
 *  \file   mathUtility.hpp
 *     1) provides some general functions and macro definitions
 *     2) randome generation functions
 *     3) Timing class: provides a simple tool for computation time calculations
 *        time can be calculated in nano-, micro-, mili- seconds and in seconds
 *        the posibility to use tic toc functions as in matlab.
 *        tic; ....your code.....; toc;
 */
/**************************************************************************
*   This library is free software; you can redistribute it and/or         *
*   modify it under the terms of the GNU Lesser General Public            *
*   License as published by the Free Software Foundation; either          *
*   version 2.1 of the License, or (at your option) any later version.    *
*                                                                         *
*   This library is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
*   Lesser General Public License for more details.                       *
*                                                                         *
*   You should have received a copy of the GNU Lesser General Public      *
*   License along with this library; if not, write to the Free Software   *
*   Foundation, Inc., 59 Temple Place,                                    *
*   Suite 330, Boston, MA  02111-1307  USA                                *
*                                                                         *
***************************************************************************/
#ifndef __MATH_UTILITY__
#define __MATH_UTILITY__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include "include.h"
#include <chrono>
#include <thread>
#include <algorithm>
#include <vector>

namespace math{

#ifndef Pi
#define Pi		3.1415926535897932384626433832795
#endif
#ifndef TwoPi
#define TwoPi	6.28318530717959
#endif
#ifndef Pi_2	//Pi/2
#define Pi_2	1.57079632679489661923
#endif

#ifndef PI
#define PI			3.1415926535897932384626433832795
#endif
#ifndef TwoPI	//2*Pi
#define TwoPI		6.283185307179586476925286766559
#endif
#ifndef PI_2
#define PI_2		1.5707963267948966192313216916398
#endif

#ifndef M_PI
#define M_PI		3.1415926535897932384626433832795
#endif
#ifndef M_2_PI
#define M_2_PI		6.283185307179586476925286766559
#endif
#ifndef M_PI_2
#define M_PI_2		1.5707963267948966192313216916398
#endif
#ifndef M_3_PI_4
#define M_3_PI_4	2.3561944901923449288469825374596
#endif
#ifndef EPSilon
#define EPSilon 2.2204460492503131e-016	/* smallest such that 1.0+DBL_EPSILON != 1.0 */
#endif

#ifndef epsiLON
#define epsiLON 0.0000001
#endif

#ifndef STRING_MAX
#define STRING_MAX		256
#endif

#define GRAVITY 9.81

#ifndef SWAP_
#define SWAP_
inline void Swap(Real &A, Real &B){	Real C;	C = A; A = B; B = C;}
#endif

#ifndef SQR_
#define SQR_
inline Real SQR(Real const &A){	return A * A;}
#endif

#ifndef EQUAL_
#define EQUAL_
inline bool EQUAL(Real a,Real b,Real eps=EPSilon){
    Real tmp=(a-b);
    return ((eps>tmp)&& (tmp>-eps) );
}
#endif

#ifndef MAX_
#define MAX_
template<class T>
inline T MAX(T const &A, T const &B){	return (A>B) ? A : B; }
#endif

#ifndef MIN_
#define MIN_
inline Real MIN(Real const &A, Real const &B){	return (A<B) ? A : B; }
#endif

// angle conversion
#ifndef DEG2RAD_
#define DEG2RAD_
inline Real deg2rad(const Real angle_deg){ return angle_deg*M_PI/180; }
#endif

#ifndef RAD2DEG_
#define RAD2DEG_

inline Real rad2deg(const Real angle_rad){ return angle_rad*180/M_PI; }
#endif

#ifndef ISZERO_
#define ISZERO_
inline bool isZero(const Real x,Real eps=EPSilon){ return fabs(x) < eps ? true : false; }
#endif

#ifndef SIGN
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

#ifndef ERROR
#define ERROR        0
#endif

#ifndef SUCCESS
#define SUCCESS      1
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef ZERO
#define ZERO static_cast<Real>(0.0)
#endif

#ifndef HALF
#define HALF static_cast<Real>(0.5)
#endif

#ifndef ONE
#define ONE static_cast<Real>(1.0)
#endif

#ifndef TWO
#define TWO static_cast<Real>(2.0)
#endif

#ifdef _MSC_VER
#ifndef __PRETTY_FUNCTION__
#define __PRETTY_FUNCTION__ __FUNCDNAME__
#endif
#endif
#define IK2PI  ((Real)6.28318530717959)
#define IKPI  ((Real)3.14159265358979)
#define IKPI_2  ((Real)1.57079632679490)
#define IKFAST_ASSERT(b) { if( !(b) ) { std::stringstream ss; ss << "ikfast exception: " << __FILE__ << ":" << __LINE__ << ": " <<__PRETTY_FUNCTION__ << ": Assertion '" << #b << "' failed"; throw std::runtime_error(ss.str()); } }
#ifdef _MSC_VER
#ifndef isnan
#define isnan _isnan
#endif
#endif // _MSC_VER
#ifndef MATH_SINCOS_THRESH
#define MATH_SINCOS_THRESH ((Real)0.000001)
#endif

#ifdef __cplusplus_
extern "C" {
#endif





Real drand (Real min, Real max, int seed=-1, int n=4);
Real Simpledrand(Real a, Real b, int seed=-1);
inline void random(Real& a) {	a = 1.98*rand()/(Real)RAND_MAX -0.99;}
Real signof(Real a);
inline bool is_even(int x) { return ((x&1) == 0); }// check this!!!!
template<class T>
inline T round(T d){return (d > 0.0) ? ::floor(d + 0.5) : ::ceil(d - 0.5);}
// Calculates log2 of number.
inline Real Log2( Real n ){  return ::log( n ) / ::log( 2.0 );  }
void log_quat(double q1[], double q2[], double log_q[]);
void mult_quat(double q1[], double q2[], double q[]);
void conjugate_quat(double q1[], double q2[]);

#ifdef __cplusplus_
}
#endif
template<class T>
vector<T> linspace(T X1, T X2, int n)
/*
 * Template linspace function (similar to MatLab)
 * generates N points between X1 and X2.
*/
{
	vector<T> result;
	// vector iterator
	int iterator = 0;

	for (int i = 0; i <= n - 2; i++)
	{
		T temp = X1 + i*(X2 - X1) / (floor((T)n) - 1);
		result.insert(result.begin() + iterator, temp);
		iterator += 1;
	}
	//iterator += 1;
	result.insert(result.begin() + iterator, X2);
	return result;
}


};

/// Simple timing
class Timing{
	typedef std::chrono::high_resolution_clock my_clock;
	typedef std::chrono::microseconds MICROseconds;
	typedef std::chrono::milliseconds MILLIseconds;
	typedef std::chrono::nanoseconds  NANOseconds;
	typedef std::chrono::duration<Real> Seconds;

	my_clock::time_point _st;
	my_clock::time_point get_st() {return _st;}
	my_clock::time_point gettime(){return my_clock::now();}

	static NANOseconds timediffNs(const my_clock::time_point& t1, const my_clock::time_point& t0){
		return std::chrono::duration_cast<NANOseconds>(t1 - t0);
	}
	static MICROseconds timediffUs(const my_clock::time_point& t1, const my_clock::time_point& t0){
		return std::chrono::duration_cast<MICROseconds>(t1 - t0);
	}
	static MILLIseconds timediffMs(const my_clock::time_point& t1,const my_clock::time_point& t0){
		return std::chrono::duration_cast<MILLIseconds>(t1 - t0);
	}
	static Seconds timediffs(const my_clock::time_point& t1,const my_clock::time_point& t0){
		return (t1 - t0);
	}
public:
	Timing() : _st(my_clock::now()){}
	my_clock::time_point reset() { _st = my_clock::now(); return _st;}
	int time_Ns()	{ return static_cast<int>(timediffNs(gettime(), _st).count());}
	int time_Us()	{ return static_cast<int>(timediffUs(gettime(), _st).count());}
	int time_Ms()	{ return static_cast<int>(timediffMs(gettime(), _st).count()); }
	int time_s()	{ return static_cast<int>(timediffs (gettime(), _st).count()); }
	static void SLEEP_Ms(unsigned int ms){std::this_thread::sleep_for(MILLIseconds(ms));}
	static void SLEEP_Us(unsigned int ms){std::this_thread::sleep_for(MICROseconds(ms));}
	static int timediff_Ns(Timing t0, Timing t1){ return static_cast<int>(timediffNs(t0.get_st(), t1.get_st()).count()); }
	static int timediff_Us(Timing t0, Timing t1){ return static_cast<int>(timediffUs(t0.get_st(), t1.get_st()).count()); }
	static int timediff_Ms(Timing t0, Timing t1){ return static_cast<int>(timediffMs(t0.get_st(), t1.get_st()).count()); }
	static Real timediff_s(Timing t0, Timing t1){ return static_cast<Real>(timediffs (t0.get_st(), t1.get_st()).count()); }
	void print_Ns()	{ std::cout<<"\nElapssed time  is: "<< timediffNs(gettime(), _st).count() <<" micro s.\n";}
	void print_Us()	{ std::cout<<"\nElapssed time  is: "<< timediffUs(gettime(), _st).count() <<" micro s.\n";}
	void print_Ms()	{ std::cout<<"\nElapssed time  is: "<< timediffMs(gettime(), _st).count() <<" ms.\n";}
	void print_s()	{ std::cout<<"\nElapssed time  is: "<< timediffs (gettime(), _st).count() <<" s.\n";}
	void print(){
		std::cout<<"\nElapssed time  is: "<< timediffNs(gettime(), _st).count() <<" nano s.\n";
		std::cout<<"\nElapssed time  is: "<< timediffUs(gettime(), _st).count() <<" micro s.\n";
		std::cout<<"Elapssed time  is: "<< timediffMs(gettime(), _st).count() <<" ms.\n";
		std::cout<<"Elapssed time  is: "<< timediffs (gettime(), _st).count() <<" s.\n";
	}
	static my_clock::time_point _tic;
};

#ifndef tic
#define tic Timing::_tic = std::chrono::high_resolution_clock::now()
#endif
#ifndef toc
#define toc std::cout<<"\nElapssed time  is: "<< std::chrono::duration<Real>(std::chrono::high_resolution_clock::now() - Timing::_tic).count() <<" s.\n"
#endif

/*int time_Us()	{ return std::chrono::duration_cast<MICROseconds>(my_clock::now()-st).count();}
int time_Ms()	{ return std::chrono::duration_cast<MILLIseconds>(my_clock::now()-st).count();}
int time_s()	{ return std::chrono::duration_cast<Seconds>(my_clock::now()-st).count();}*/

#endif
