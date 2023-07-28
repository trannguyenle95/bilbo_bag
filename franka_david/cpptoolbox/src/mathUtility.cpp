
#include "mathUtility.hpp"

std::chrono::high_resolution_clock::time_point Timing::_tic = std::chrono::high_resolution_clock::now();

namespace math{

//-------------------------------------------------------------------------------------------
Real drand (Real min, Real max, int seed, int n) {
	if ( seed >= 0 )
		srand(seed);
	long k=1;
	long r = ((long)rand() << 15) | rand();
	for (int i=0;i<n; i++)
		k*=10;
	return ((r % (int)(max*k - min*k + 1)) + (int)(min*k)) / (Real) k;
}
//-------------------------------------------------------------------------------------------
Real Simpledrand(Real a, Real b, int seed) {
	if ( seed >= 0 )
		srand(seed);
	Real random = ((Real) rand()) / (Real) RAND_MAX;
	Real diff = b - a;
	Real r = random * diff;
	return a + r;
}
//-------------------------------------------------------------------------------------------
Real signof(Real a){
	//calculates the SIGNOF function
	//return	: SIGNOF(a)
	if (a > 0)		return  1.0;
	else
		if (a < 0)	return -1.0;
		else		return  0.0;
}
void conjugate_quat(double q1[], double q2[])
{
	// Quaternion conjugation
	q2[0] = q1[0];
	for (int i = 1; i <= 3; i++)
		q2[i] = -q1[i];
}
void mult_quat(double q1[], double q2[], double q[])
{
	// Quaternion multiplication
	q[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];

	for (int i = 1; i <= 3; i++)
		q[i] = q1[0] * q2[i] + q2[0] * q1[i];
	q[1] += q1[2] * q2[3] - q1[3] * q2[2];
	q[2] += q1[3] * q2[1] - q1[1] * q2[3];
	q[3] += q1[1] * q2[2] - q1[2] * q2[1];
}
void log_quat(double q1[], double q2[], double log_q[])
/*!
@brief Calculates logarithm of orientation difference between quaternions
*/
{
	double q[4], q2c[4], tmp, theta;
	int i;

	conjugate_quat(q2, q2c);
	mult_quat(q1, q2c, q);

	// Normalization to prevent any crashes
	tmp = 0;
	for (i = 0; i < 4; i++)
		tmp += SQR(q[i]);
	tmp = sqrt(tmp);
	for (i = 0; i < 4; i++)
		q[i] = q[i] / tmp;

	tmp = 0;
	for (i = 1; i <= 3; i++)
		tmp += SQR(q[i]);
	tmp = sqrt(tmp);
	if (tmp > EPSilon) {
		theta = acos(q[0]);
		for (i = 1; i <= 3; i++)
			log_q[i] = theta * q[i] / tmp;
	}
	else if (q[0] > 0)
		log_q[1] = log_q[2] = log_q[3] = 0;
	else {
		tmp = 0;
		for (i = 1; i <= 3; i++)
			tmp += SQR(log_q[i]);
		tmp = sqrt(tmp);

		if (tmp > EPSilon)
			for (i = 1; i <= 3; i++)
				log_q[i] = M_PI * log_q[i] / tmp;
		else
			log_q[1] = log_q[2] = log_q[3] = 0;
	}
	log_q[0] = 0;
}

};
