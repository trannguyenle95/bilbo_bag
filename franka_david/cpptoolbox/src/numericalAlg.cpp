/*
 * numericalAlg.cpp
 *
 *  Created on: Nov 10, 2015
 *      Author: fares
 */


#include "numericalAlg.h"


using namespace std;

//****************************************************************************80
int dchdc ( double a[], int lda, int p, double work[], int ipvt[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DCHDC computes the Cholesky decomposition of a positive definite matrix.
//
//  Discussion:
//
//    A pivoting option allows the user to estimate the condition of a
//    positive definite matrix or determine the rank of a positive
//    semidefinite matrix.
//
//    For positive definite matrices, INFO = P is the normal return.
//
//    For pivoting with positive semidefinite matrices, INFO will
//    in general be less than P.  However, INFO may be greater than
//    the rank of A, since rounding error can cause an otherwise zero
//    element to be positive.  Indefinite systems will always cause
//    INFO to be less than P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*P].
//    On input, A contains the matrix whose decomposition is to
//    be computed.  Only the upper half of A need be stored.
//    The lower part of the array a is not referenced.
//    On output, A contains in its upper half the Cholesky factor
//    of the input matrix, as it has been permuted by pivoting.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int P, the order of the matrix.
//
//    Input, double WORK[P] is a work array.
//
//    Input/output, int IPVT[P].
//    On input, IPVT contains integers that control the selection
//    of the pivot elements, if pivoting has been requested.
//    Each diagonal element A(K,K) is placed in one of three classes
//    according to the value of IPVT(K).
//
//      > 0, then X(K) is an initial element.
//      = 0, then X(K) is a free element.
//      < 0, then X(K) is a final element.
//
//    Before the decomposition is computed, initial elements are moved by
//    symmetric row and column interchanges to the beginning of the array A
//    and final elements to the end.  Both initial and final elements are
//    frozen in place during the computation and only free elements are moved.
//    At the K-th stage of the reduction, if A(K,K) is occupied by a free
//    element, it is interchanged with the largest free element A(L,L) with
//    K <= L.  IPVT is not referenced if JOB is 0.
//
//    On output, IPVT(J) contains the index of the diagonal element
//    of A that was moved into the J-th position, if pivoting was requested.
//
//    Input, int JOB, initiates column pivoting.
//    0, no pivoting is done.
//    nonzero, pivoting is done.
//
//    Output, int DCHDC, contains the index of the last positive diagonal
//    element of the Cholesky factor.
//
{
	int info = p, j, jp, jt, k, l, maxl, pl = 1, pu = 0;
	double maxdia;
	bool negk;
	bool swapk;
	double temp;

	//!<  Pivoting has been requested.
	//!<  Rearrange the the elements according to IPVT.
	if ( job != 0 ){
		for ( k = 1; k <= p; k++ ){
			swapk = ( 0 < ipvt[k-1] );
			negk = ( ipvt[k-1] < 0 );

			if ( negk )	{	ipvt[k-1] = -k;		}
			else		{	ipvt[k-1] = k;		}

			if ( swapk ){
				if ( k != pl ){
					dswap ( pl-1, a+0+(k-1)*lda, 1, a+0+(pl-1)*lda, 1 );

					temp = a[k-1+(k-1)*lda];
					a[k-1+(k-1)*lda] = a[pl-1+(pl-1)*lda];
					a[pl-1+(pl-1)*lda] = temp;

					for ( j = pl+1; j <= p; j++ ){
						if ( j < k ){
							temp = a[pl-1+(j-1)*lda];
							a[pl-1+(j-1)*lda] = a[j-1+(k-1)*lda];
							a[j-1+(k-1)*lda] = temp;
						}else if ( k < j ){
							temp = a[k-1+(j-1)*lda];
							a[k-1+(j-1)*lda] = a[pl-1+(j-1)*lda];
							a[pl-1+(j-1)*lda] = temp;
						}
					}
					ipvt[k-1] = ipvt[pl-1];
					ipvt[pl-1] = k;
				}
				pl = pl + 1;
			}
		}
		pu = p;

		for ( k = p; pl <= k; k-- ){
			if ( ipvt[k-1] < 0 ){
				ipvt[k-1] = -ipvt[k-1];

				if ( pu != k ){
					dswap ( k-1, a+0+(k-1)*lda, 1, a+0+(pu-1)*lda, 1 );

					temp = a[k-1+(k-1)*lda];
					a[k-1+(k-1)*lda] = a[pu-1+(pu-1)*lda];
					a[pu-1+(pu-1)*lda] = temp;

					for ( j = k+1; j <= p; j++ ){
						if ( j < pu ){
							temp = a[k-1+(j-1)*lda];
							a[k-1+(j-1)*lda] = a[j-1+(pu-1)*lda];
							a[j-1+(pu-1)*lda] = temp;
						}else if ( pu < j ){
							temp = a[k-1+(j-1)*lda];
							a[k-1+(j-1)*lda] = a[pu-1+(j-1)*lda];
							a[pu-1+(j-1)*lda] = temp;
						}
					}
					jt = ipvt[k-1];
					ipvt[k-1] = ipvt[pu-1];
					ipvt[pu-1] = jt;
				}
				pu = pu - 1;
			}
		}
	}

	for ( k = 1; k <= p; k++ ){
		//  Reduction loop.
		maxdia = a[k-1+(k-1)*lda];
		maxl = k;
		//  Determine the pivot element.
		if ( pl <= k && k < pu ){
			for ( l = k+1; l <= pu; l++ ){
				if ( maxdia < a[l-1+(l-1)*lda] ){
					maxdia = a[l-1+(l-1)*lda];
					maxl = l;
				}
			}
		}
		//  Quit if the pivot element is not positive.
		if ( maxdia <= ZERO ){	info = k - 1;		return info;		}
		//  Start the pivoting and update IPVT.
		if ( k != maxl ){
			dswap ( k-1, a+0+(k-1)*lda, 1, a+0+(maxl-1)*lda, 1 );
			a[maxl-1+(maxl-1)*lda] = a[k-1+(k-1)*lda];
			a[k-1+(k-1)*lda] = maxdia;
			jp = ipvt[maxl-1];
			ipvt[maxl-1] = ipvt[k-1];
			ipvt[k-1] = jp;
		}
		//  Reduction step.
		//  Pivoting is contained across the rows.
		work[k-1] = sqrt ( a[k-1+(k-1)*lda] );
		a[k-1+(k-1)*lda] = work[k-1];

		for ( j = k+1; j <= p; j++ ){
			if ( k != maxl ){
				if ( j < maxl ){
					temp = a[k-1+(j-1)*lda];
					a[k-1+(j-1)*lda] = a[j-1+(maxl-1)*lda];
					a[j-1+(maxl-1)*lda] = temp;
				}else if ( maxl < j ){
					temp = a[k-1+(j-1)*lda];
					a[k-1+(j-1)*lda] = a[maxl-1+(j-1)*lda];
					a[maxl-1+(j-1)*lda] = temp;
				}
			}
			a[k-1+(j-1)*lda] = a[k-1+(j-1)*lda] / work[k-1];
			work[j-1] = a[k-1+(j-1)*lda];
			temp = -a[k-1+(j-1)*lda];
			daxpy ( j-k, temp, work+k, 1, a+k+(j-1)*lda, 1 );
		}
	}
	return info;
}
//****************************************************************************80
int dchdd ( double r[], int ldr, int p, double x[], double z[], int ldz,
		int nz, double y[], double rho[], double c[], double s[] )
//****************************************************************************80
//
//  Purpose:
//
//    DCHDD downdates an augmented Cholesky decomposition.
//
//  Discussion:
//
//    DCHDD can also downdate the triangular factor of an augmented QR
//    decomposition.
//
//    Specifically, given an upper triangular matrix R of order P, a
//    row vector X, a column vector Z, and a scalar Y, DCHDD
//    determines an orthogonal matrix U and a scalar ZETA such that
//
//          (R   Z )     (RR  ZZ)
//      U * (      )  =  (      ),
//          (0 ZETA)     ( X   Y)
//
//    where RR is upper triangular.
//
//    If R and Z have been obtained from the factorization of a least squares
//    problem, then RR and ZZ are the factors corresponding to the problem
//    with the observation (X,Y) removed.  In this case, if RHO
//    is the norm of the residual vector, then the norm of
//    the residual vector of the downdated problem is
//    sqrt ( RHO * RHO - ZETA * ZETA ). DCHDD will simultaneously downdate
//    several triplets (Z, Y, RHO) along with R.
//
//    For a less terse description of what DCHDD does and how
//    it may be applied, see the LINPACK guide.
//
//    The matrix U is determined as the product U(1)*...*U(P)
//    where U(I) is a rotation in the (P+1,I)-plane of the form
//
//      ( C(I)      -S(I)    )
//      (                    ).
//      ( S(I)       C(I)    )
//
//    The rotations are chosen so that C(I) is real.
//
//    The user is warned that a given downdating problem may be impossible
//    to accomplish or may produce inaccurate results.  For example, this
//    can happen if X is near a vector whose removal will reduce the
//    rank of R.  Beware.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double R[LDR*P], the upper triangular matrix that
//    is to be  downdated.  The part of R below the diagonal is not referenced.
//
//    Input, int LDR, the leading dimension of the array R.
//    LDR must be at least P.
//
//    Input, int P, the order of the matrix R.
//
//    Input, double X[P], the row vector that is to be removed from R.
//
//    Input/output, double Z[LDZ*NZ], an array of NZ P-vectors
//    which are to be downdated along with R.
//
//    Input, int LDZ, the leading dimension of the array Z.
//    LDZ must be at least P.
//
//    Input, int NZ, the number of vectors to be downdated.
//    NZ may be zero, in which case Z, Y, and RHO are not referenced.
//
//    Input, double Y[NZ], the scalars for the downdating of
//    the vectors Z.
//
//    Input/output, double RHO[NZ], the norms of the residual vectors.
//    On output these have been changed along with R and Z.
//
//    Output, double C[P], S[P], the cosines and sines of the
//    transforming rotations.
//
//    Output, int DCHDD, return flag.
//     0, the entire downdating was successful.
//    -1, if R could not be downdated.  In this case, all quantities
//        are left unaltered.
//     1, if some RHO could not be downdated.  The offending RHO's are
//        set to -1.
//
{
	double a, alpha, azeta, b, norm, scale, t, xx, zeta;
	int i, ii, info = 0, j;
	//  Solve R' * A = X, placing the result in the array S.
	s[0] = x[0] / r[0+0*ldr];

	for ( j = 2; j <= p; j++ ){
		s[j-1] = x[j-1] - ddot ( j-1, r+0+(j-1)*ldr, 1, s, 1 );
		s[j-1] = s[j-1] / r[j-1+(j-1)*ldr];
	}

	norm = dnrm2 ( p, s, 1 );

	if ( ONE <= norm ){		info = -1;		return info;	}

	alpha = sqrt ( ONE - norm * norm );
	//  Determine the transformations.
	for ( ii = 1; ii <= p; ii++ ){
		i = p - ii + 1;
		scale = alpha + r8_abs ( s[i-1] );
		a = alpha / scale;
		b = s[i-1] / scale;
		norm = sqrt ( a * a + b * b );
		c[i-1] = a / norm;
		s[i-1] = b / norm;
		alpha = scale * norm;
	}
	//  Apply the transformations to R.
	for ( j = 1; j <= p; j++ ){
		xx = ZERO;
		for ( ii = 1; ii <= j; ii++ ){
			i = j - ii + 1;
			t = c[i-1] * xx + s[i-1] * r[i-1+(j-1)*ldr];
			r[i-1+(j-1)*ldr] = c[i-1] * r[i-1+(j-1)*ldr] - s[i-1] * xx;
			xx = t;
		}
	}
	//  If required, downdate Z and RHO.
	for ( j = 1; j <= nz; j++ ){
		zeta = y[j-1];
		for ( i = 1; i <= p; i++ ){
			z[i-1+(j-1)*ldz] = ( z[i-1+(j-1)*ldz] - s[i-1] * zeta ) / c[i-1];
			zeta = c[i-1] * zeta - s[i-1] * z[i-1+(j-1)*ldz];
		}

		azeta = r8_abs ( zeta );

		if ( rho[j-1] < azeta ){	info = 1;		rho[j-1] = -ONE;		}
		else{			rho[j-1] = rho[j-1] * sqrt ( ONE - pow ( azeta / rho[j-1], 2 ) );		}
	}
	return info;
}
//****************************************************************************80
void dchex ( double r[], int ldr, int p, int k, int l, double z[], int ldz,
		int nz, double c[], double s[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DCHEX updates the Cholesky factorization of a positive definite matrix.
//
//  Discussion:
//
//    The factorization has the form
//
//      A = R' * R
//
//    where A is a positive definite matrix of order P.
//
//    The updating involves diagonal permutations of the form
//
//      E' * A * E
//
//    where E is a permutation matrix.  Specifically, given
//    an upper triangular matrix R and a permutation matrix
//    E (which is specified by K, L, and JOB), DCHEX determines
//    an orthogonal matrix U such that
//
//      U * R * E = RR,
//
//    where RR is upper triangular.  At the user's option, the
//    transformation U will be multiplied into the array Z.
//    If A = X'*X, so that R is the triangular part of the
//    QR factorization of X, then RR is the triangular part of the
//    QR factorization of X*E, that is, X with its columns permuted.
//
//    For a less terse description of what DCHEX does and how
//    it may be applied, see the LINPACK guide.
//
//    The matrix Q is determined as the product U(L-K)*...*U(1)
//    of plane rotations of the form
//
//      (    C(I)       S(I) )
//      (                    ),
//      (   -S(I)       C(I) )
//
//    where C(I) is real, the rows these rotations operate on
//    are described below.
//
//    There are two types of permutations, which are determined
//    by the value of JOB.
//
//    1, right circular shift.  The columns are rearranged in the order:
//
//         1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.
//
//       U is the product of L-K rotations U(I), where U(I)
//       acts in the (L-I,L-I+1)-plane.
//
//    2, left circular shift: the columns are rearranged in the order
//
//         1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.
//
//       U is the product of L-K rotations U(I), where U(I)
//       acts in the (K+I-1,K+I)-plane.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double R[LDR*P].  On input, the upper
//    triangular factor that is to be updated.  Elements of R below the
//    diagonal are not referenced.  On output, R has been updated.
//
//    Input, int LDR, the leading dimension of the array R.
//    LDR must be at least P.
//
//    Input, int P, the order of the matrix R.
//
//    Input, int K, the first column to be permuted.
//
//    Input, int L, the last column to be permuted.
//    L must be strictly greater than K.
//
//    Input/output double Z[LDZ*NZ], an array of NZ P-vectors into
//    which the transformation U is multiplied.  Z is not referenced if NZ = 0.
//    On output, Z has been updated.
//
//    Input, int LDZ, the leading dimension of the array Z.
//    LDZ must be at least P.
//
//    Input, int NZ, the number of columns of the matrix Z.
//
//    Output, double C[P], S[P], the cosines and sines of the
//    transforming rotations.
//
//    Input, int JOB, determines the type of permutation.
//    1, right circular shift.
//    2, left circular shift.
//
{
	int i, ii, il, iu, j, jj, lm1 = l - 1, lmk = l - k;
	double t;
	//  Right circular shift.
	if ( job == 1 ){
		//  Reorder the columns.
		for ( i = 1; i <= l; i++ ){
			ii = l - i + 1;
			s[i-1] = r[ii-1+(l-1)*ldr];
		}

		for ( jj = k; jj <= lm1; jj++ ){
			j = lm1 - jj + k;
			for ( i = 1; i <= j; i++ ){		r[i-1+(j)*ldr] = r[i-1+(j-1)*ldr];		}
			r[j+(j)*ldr] = ZERO;
		}

		for ( i = 1; i <= k-1; i++ ){
			ii = l - i + 1;
			r[i-1+(k-1)*ldr] = s[ii-1];
		}
		//  Calculate the rotations.
		t = s[0];
		for ( i = 1; i <= lmk; i++ ){
			drotg ( s+i, &t, c+i-1, s+i-1 );
			t = s[i];
		}

		r[k-1+(k-1)*ldr] = t;

		for ( j = k+1; j <= p; j++ ){
			il = i4_max ( 1, l-j+1 );
			for ( ii = il; ii <= lmk; ii++ ){
				i = l - ii;
				t = c[ii-1] * r[i-1+(j-1)*ldr] + s[ii-1] * r[i+(j-1)*ldr];
				r[i+(j-1)*ldr] = c[ii-1] * r[i+(j-1)*ldr] - s[ii-1] * r[i-1+(j-1)*ldr];
				r[i-1+(j-1)*ldr] = t;
			}
		}
		//  If required, apply the transformations to Z.
		for ( j = 1; j <= nz; j++ ){
			for ( ii = 1; ii <= lmk; ii++ ){
				i = l - ii;
				t = c[ii-1] * z[i-1+(j-1)*ldr] + s[ii-1] * z[i+(j-1)*ldr];
				z[i+(j-1)*ldr] = c[ii-1] * z[i+(j-1)*ldr] - s[ii-1] * z[i-1+(j-1)*ldr];
				z[i-1+(j-1)*ldr] = t;
			}
		}
	}else{//  Left circular shift.
		//  Reorder the columns.
		for ( i = 1; i <= k; i++ ){
			ii = lmk + i;
			s[ii-1] = r[i-1+(k-1)*ldr];
		}

		for ( j = k; j <= lm1; j++ ){
			for ( i = 1; i <= j; i++ ){		r[i-1+(j-1)*ldr] = r[i-1+(j)*ldr];		}
			jj = j - k + 1;
			s[jj-1] = r[j+(j)*ldr];
		}

		for ( i = 1; i <= k; i++ ){
			ii = lmk + i;
			r[i-1+(l-1)*ldr] = s[ii-1];
		}

		for ( i = k+1; i <= l; i++ ){		r[i-1+(l-1)*ldr] = ZERO;		}
		//  Reduction loop.
		for ( j = k; j <= p; j++ ){
			//  Apply the rotations.
			if ( j != k ){
				iu = i4_min ( j-1, l-1 );

				for ( i = k; i <= iu; i++ ){
					ii = i - k + 1;
					t = c[ii-1] * r[i-1+(j-1)*ldr] + s[ii-1] * r[i+(j-1)*ldr];
					r[i+(j-1)*ldr] = c[ii-1] * r[i+(j-1)*ldr]
												 - s[ii-1] * r[i-1+(j-1)*ldr];
					r[i-1+(j-1)*ldr] = t;
				}
			}

			if ( j < l ){
				jj = j - k + 1;
				t = s[jj-1];
				drotg ( r+j-1+(j-1)*ldr, &t, c+jj-1, s+jj-1 );
			}
		}
		//  Apply the rotations to Z.
		for ( j = 1; j <= nz; j++ ){
			for ( i = k; i <= lm1; i++ ){
				ii = i - k + 1;
				t = c[ii-1] * z[i-1+(j-1)*ldr] + s[ii-1] * z[i+(j-1)*ldr];
				z[i+(j-1)*ldr] = c[ii-1] * z[i+(j-1)*ldr] - s[ii-1] * z[i-1+(j-1)*ldr];
				z[i-1+(j-1)*ldr] = t;
			}
		}
	}

	return;
}
//****************************************************************************80
void dchud ( double r[], int ldr, int p, double x[], double z[], int ldz,
		int nz, double y[], double rho[], double c[], double s[] )
//****************************************************************************80
//
//  Purpose:
//
//    DCHUD updates an augmented Cholesky decomposition.
//
//  Discussion:
//
//    DCHUD can also update the triangular part of an augmented QR
//    decomposition.
//
//    Specifically, given an upper triangular matrix R of order P, a row vector
//    X, a column vector Z, and a scalar Y, DCHUD determines a unitary matrix
//    U and a scalar ZETA such that
//
//           (R  Z)     (RR   ZZ )
//      U  * (    )  =  (        ),
//           (X  Y)     ( 0  ZETA)
//
//    where RR is upper triangular.
//
//    If R and Z have been obtained from the factorization of a least squares
//    problem, then RR and ZZ are the factors corresponding to the problem
//    with the observation (X,Y) appended.  In this case, if RHO is the
//    norm of the residual vector, then the norm of the residual vector of
//    the updated problem is sqrt ( RHO * RHO + ZETA * ZETA ).  DCHUD will
//    simultaneously update several triplets (Z, Y, RHO).
//
//    For a less terse description of what DCHUD does and how
//    it may be applied, see the LINPACK guide.
//
//    The matrix U is determined as the product U(P)*...*U(1),
//    where U(I) is a rotation in the (I,P+1) plane of the form
//
//      (     C(I)      S(I) )
//      (                    ).
//      (    -S(I)      C(I) )
//
//    The rotations are chosen so that C(I) is real.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double R[LDR*P], the upper triangular matrix to be
//    updated.  The part of R below the diagonal is not referenced.
//    On output, the matrix has been updated.
//
//    Input, int LDR, the leading dimension of the array R.
//    LDR must be at least equal to P.
//
//    Input, int P, the order of the matrix R.
//
//    Input, double X[P], the row to be added to R.
//
//    Input/output, double Z[LDZ*NZ], contains NZ P-vectors
//    to be updated with R.
//
//    Input, int LDZ, the leading dimension of the array Z.
//    LDZ must be at least P.
//
//    Input, int NZ, the number of vectors to be updated.  NZ may be
//    zero, in which case Z, Y, and RHO are not referenced.
//
//    Input, double Y[NZ], the scalars for updating the vectors Z.
//
//    Input/output, double RHO[NZ].  On input, the norms of the
//    residual vectors to be updated.  If RHO(J) is negative, it is left
//    unaltered.
//
//    Output, double C[P], S[P], the cosines and sines of the
//    transforming rotations.
//
{
	double azeta, scale, t, xj, zeta;
	int i, j;
	//  Update R.
	for ( j = 1; j <= p; j++ ){
		xj = x[j-1];
		//  Apply the previous rotations.
		for ( i = 1; i <= j-1; i++ ){
			t = c[i-1] * r[i-1+(j-1)*ldz] + s[i-1] * xj;
			xj = c[i-1] * xj - s[i-1] * r[i-1+(j-1)*ldz];
			r[i-1+(j-1)*ldz] = t;
		}
		//  Compute the next rotation.
		drotg ( r+j-1+(j-1)*ldr, &xj, c+j-1, s+j-1 );
	}
	//  If required, update Z and RHO.
	for ( j = 1; j <= nz; j++ ){
		zeta = y[j-1];
		for ( i = 1; i <= p; i++ ){
			t =    c[i-1] * z[i-1+(j-1)*ldz] + s[i-1] * zeta;
			zeta = c[i-1] * zeta   - s[i-1] * z[i-1+(j-1)*ldz];
			z[i-1+(j-1)*ldz] = t;
		}

		azeta = r8_abs ( zeta );

		if ( azeta != ZERO && ZERO <= rho[j-1] ){
			scale = azeta + rho[j-1];
			rho[j-1] = scale * sqrt ( pow ( azeta / scale, 2 ) + pow ( rho[j-1] / scale, 2 ) );
		}
	}
	return;
}
//****************************************************************************80
double dgbco ( double abd[], int lda, int n, int ml, int mu, int ipvt[], double z[] )
//****************************************************************************80
//
//  Purpose:
//
//    DGBCO factors a real band matrix and estimates its condition.
//
//  Discussion:
//
//    If RCOND is not needed, DGBFA is slightly faster.
//
//    To solve A*X = B, follow DGBCO by DGBSL.
//
//    To compute inverse(A)*C, follow DGBCO by DGBSL.
//
//    To compute determinant(A), follow DGBCO by DGBDI.
//
//  Example:
//
//    If the original matrix is
//
//      11 12 13  0  0  0
//      21 22 23 24  0  0
//       0 32 33 34 35  0
//       0  0 43 44 45 46
//       0  0  0 54 55 56
//       0  0  0  0 65 66
//
//    then for proper band storage,
//
//      N = 6, ML = 1, MU = 2, 5 <= LDA and ABD should contain
//
//       *  *  *  +  +  +      * = not used
//       *  * 13 24 35 46      + = used for pivoting
//       * 12 23 34 45 56
//      11 22 33 44 55 66
//      21 32 43 54 65  *
//
//  Band storage:
//
//    If A is a band matrix, the following program segment
//    will set up the input.
//
//      ml = (band width below the diagonal)
//      mu = (band width above the diagonal)
//      m = ml + mu + 1
//
//      do j = 1, n
//        i1 = max ( 1, j-mu )
//        i2 = min ( n, j+ml )
//        do i = i1, i2
//          k = i - j + m
//          abd(k,j) = a(i,j)
//        }
//      }
//
//    This uses rows ML+1 through 2*ML+MU+1 of ABD.  In addition, the first
//    ML rows in ABD are used for elements generated during the
//    triangularization.  The total number of rows needed in ABD is
//    2*ML+MU+1.  The ML+MU by ML+MU upper left triangle and the ML by ML
//    lower right triangle are not referenced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double ABD[LDA*N].  On input, the matrix in band
//    storage.  The columns of the matrix are stored in the columns of ABD and
//    the diagonals of the matrix are stored in rows ML+1 through 2*ML+MU+1
//    of ABD.  On output, an upper triangular matrix in band storage and
//    the multipliers which were used to obtain it.  The factorization can
//    be written A = L*U where L is a product of permutation and unit lower
//    triangular matrices and U is upper triangular.
//
//    Input, int LDA, the leading dimension of the array ABD.
//    2*ML + MU + 1 <= LDA is required.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, MU, the number of diagonals below and above the
//    main diagonal.  0 <= ML < N, 0 <= MU < N.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Workspace, double Z[N], a work vector whose contents are
//    usually unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
//    Output, double DGBCO, an estimate of the reciprocal condition number RCOND
//    of A.  For the system A*X = B, relative perturbations in A and B of size
//    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
//    If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0D+00
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
{
	int i, info, is, j, ju, k, l = ml + 1, la, lm, lz, m, mm;
	double anorm = ZERO, ek, rcond, s, sm, t, wk, wkm, ynorm;
	//  Compute the 1-norm of A.
	is = l + mu;

	for ( j = 1; j <= n; j++ ){
		anorm = r8_max ( anorm, dasum ( l, abd+is-1+(j-1)*lda, 1 ) );
		if ( ml + 1 < is )	{		is = is - 1;	}
		if ( j <= mu )		{		l = l + 1;		}
		if ( n - ml <= j )	{		l = l - 1;		}
	}
	//  Factor.
	info = dgbfa ( abd, lda, n, ml, mu, ipvt );
	//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
	//
	//  Estimate = norm(Z)/norm(Y) where  a*z = y  and A'*Y = E.
	//
	//  A' is the transpose of A.  The components of E are
	//  chosen to cause maximum local growth in the elements of W where
	//  U'*W = E.  The vectors are frequently rescaled to avoid
	//  overflow.
	//
	//  Solve U' * W = E.
	ek = ONE;
	for ( i = 1; i <= n; i++ ){		z[i-1] = ZERO;	}
	m = ml + mu + 1;
	ju = 0;

	for ( k = 1; k <= n; k++ ){
		if ( z[k-1] != ZERO ){		ek = ek * r8_sign ( -z[k-1] );		}
		if ( r8_abs ( abd[m-1+(k-1)*lda] ) < r8_abs ( ek - z[k-1] ) ){
			s = r8_abs ( abd[m-1+(k-1)*lda] ) / r8_abs ( ek - z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];		}
			ek = s * ek;
		}

		wk = ek - z[k-1];
		wkm = -ek - z[k-1];
		s = r8_abs ( wk );
		sm = r8_abs ( wkm );

		if ( abd[m-1+(k-1)*lda] != ZERO ){
			wk = wk / abd[m-1+(k-1)*lda];
			wkm = wkm / abd[m-1+(k-1)*lda];
		}else{
			wk = ONE;
			wkm = ONE;
		}

		ju = i4_min ( i4_max ( ju, mu+ipvt[k-1] ), n );
		mm = m;

		if ( k+1 <= ju ){
			for ( j = k+1; j <= ju; j++ ){
				mm = mm - 1;
				sm = sm + r8_abs ( z[j-1] + wkm * abd[mm-1+(j-1)*lda] );
				z[j-1] = z[j-1] + wk * abd[mm-1+(j-1)*lda];
				s = s + r8_abs ( z[j-1] );
			}

			if ( s < sm ){
				t = wkm - wk;
				wk = wkm;
				mm = m;
				for ( j = k+1; j <= ju; ju++ ){
					mm = mm - 1;
					z[j-1] = z[j-1] + t * abd[mm-1+(j-1)*lda];
				}
			}
		}
		z[k-1] = wk;
	}

	s = dasum ( n, z, 1 );

	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	//  Solve L' * Y = W.
	for ( k = n; 1 <= k; k-- ){
		lm = i4_min ( ml, n-k );

		if ( k < m ){		z[k-1] = z[k-1] + ddot ( lm, abd+m+(k-1)*lda, 1, z+k, 1 );		}

		if ( ONE < r8_abs ( z[k-1] ) ){
			s = ONE / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ )	{		z[i-1] = s * z[i-1];			}
		}
		l = ipvt[k-1];
		t = z[l-1];
		z[l-1] = z[k-1];
		z[k-1] = t;
	}

	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	ynorm = ONE;
	//  Solve L * V = Y.
	for ( k = 1; k <= n; k++ ){
		l = ipvt[k-1];
		t = z[l-1];
		z[l-1] = z[k-1];
		z[k-1] = t;
		lm = i4_min ( ml, n-k );

		if ( k < n ){	daxpy ( lm, t, abd+m+(k-1)*lda, 1, z+k, 1 );	}

		if ( ONE < r8_abs ( z[k-1] ) ){
			s = ONE / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];		}
			ynorm = s * ynorm;
		}
	}
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;
	//  Solve U * Z = W.
	for ( k = n; 1 <= k; k-- ){
		if ( r8_abs ( abd[m-1+(k-1)*lda] ) < r8_abs ( z[k-1] ) ){
			s = r8_abs ( abd[m-1+(k-1)*lda] ) / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ynorm = s * ynorm;
		}

		if ( abd[m-1+(k-1)*lda] != ZERO ){	z[k-1] = z[k-1] / abd[m-1+(k-1)*lda];	}
		else{			z[k-1] = ONE;		}

		lm = i4_min ( k, m ) - 1;
		la = m - lm;
		lz = k - lm;
		t = -z[k-1];
		daxpy ( lm, t, abd+la-1+(k-1)*lda, 1, z+lz-1, 1 );
	}
	//  Make ZNORM = 1.0.
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){
		z[i-1] = s * z[i-1];
	}
	ynorm = s * ynorm;

	if ( anorm != ZERO ){	rcond = ynorm / anorm;	}
	else{		rcond = ZERO;	}

	return rcond;
}
//****************************************************************************80
void dgbdi ( double abd[], int lda, int n, int ml, int mu, int ipvt[], double det[2] )
//****************************************************************************80
//
//  Purpose:
//
//    DGBDI computes the determinant of a band matrix factored by DGBCO or DGBFA.
//
//  Discussion:
//
//    If the inverse is needed, use DGBSL N times.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double ABD[LDA*N], the LU factor information from DGBCO
//    or DGBFA.
//
//    Input, int LDA, the leading dimension of the array ABD.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, MU, the number of diagonals below and above the
//    main diagonal.  0 <= ML < N, 0 <= MU < N.
//
//    Input, int IPVT[N], the pivot vector from DGBCO or DGBFA.
//
//    Output, double DET[2], the determinant of the original matrix,
//    if requested.
//      determinant = DET[0] * 10.0**DET[1]
//    with  1.0D+00 <= abs ( DET[0] ) < 10.0D+00 or DET[0] = 0.0D+00.
//
{
	int i, m = ml + mu + 1;
	det[0] = ONE;
	det[1] = ZERO;

	for ( i = 1; i <= n; i++ ){
		if ( ipvt[i-1] != i )	{		det[0] = -det[0];		}

		det[0] = det[0] * abd[m-1+(i-1)*lda];

		if ( det[0] == ZERO )	{		return;		}

		while ( r8_abs ( det[0] ) < ONE ){
			det[0] = det[0] * static_cast<Real>(10.0);
			det[1] = det[1] - ONE;
		}

		while ( static_cast<Real>(10.0) <= r8_abs ( det[0] ) ){
			det[0] = det[0] / static_cast<Real>(10.0);
			det[1] = det[1] + ONE;
		}
	}
	return;
}
//****************************************************************************80
int dgbfa ( double abd[], int lda, int n, int ml, int mu, int ipvt[] )
//****************************************************************************80
//
//  Purpose:
//
//    DGBFA factors a real band matrix by elimination.
//
//  Discussion:
//
//    DGBFA is usually called by DGBCO, but it can be called
//    directly with a saving in time if RCOND is not needed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double ABD[LDA*N].  On input, the matrix in band
//    storage.  The columns of the matrix are stored in the columns of ABD
//    and the diagonals of the matrix are stored in rows ML+1 through
//    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
//    and the multipliers which were used to obtain it.  The factorization
//    can be written A = L*U where L is a product of permutation and unit lower
//    triangular matrices and U is upper triangular.
//
//    Input, int LDA, the leading dimension of the array ABD.
//    2*ML + MU + 1 <= LDA is required.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, MU, the number of diagonals below and above the
//    main diagonal.  0 <= ML < N, 0 <= MU < N.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, integer DGBFA, error flag.
//    0, normal value.
//    K, if U(K,K) == 0.0D+00.  This is not an error condition for this
//      subroutine, but it does indicate that DGBSL will divide by zero if
//      called.  Use RCOND in DGBCO for a reliable indication of singularity.
//
{
	int i, i0, info = 0, j, j0, j1, ju, jz, k, l, lm, m = ml + mu + 1, mm;
	double t;
	//  Zero initial fill-in columns.
	j0 = mu + 2;
	j1 = i4_min ( n, m ) - 1;

	for ( jz = j0; jz <= j1; jz++ ){
		i0 = m + 1 - jz;
		for ( i = i0; i <= ml; i++ ){		abd[i-1+(jz-1)*lda] = ZERO;		}
	}

	jz = j1;
	ju = 0;
	//  Gaussian elimination with partial pivoting.
	for ( k = 1; k <= n-1; k++ ){
		//  Zero out the next fill-in column.
		jz = jz + 1;
		if ( jz <= n ){		for ( i = 1; i <= ml; i++ )	{	abd[i-1+(jz-1)*lda] = ZERO;	}		}
		//  Find L = pivot index.
		lm = i4_min ( ml, n-k );
		l = idamax ( lm+1, abd+m-1+(k-1)*lda, 1 ) + m - 1;
		ipvt[k-1] = l + k - m;
		//  Zero pivot implies this column already triangularized.
		if ( abd[l-1+(k-1)*lda] == ZERO )	{		info = k;		}
		//  Interchange if necessary.
		else{
			if ( l != m ){
				t = abd[l-1+(k-1)*lda];
				abd[l-1+(k-1)*lda] = abd[m-1+(k-1)*lda];
				abd[m-1+(k-1)*lda] = t;
			}
			//  Compute multipliers.
			t = -ONE / abd[m-1+(k-1)*lda];
			dscal ( lm, t, abd+m+(k-1)*lda, 1 );
			//  Row elimination with column indexing.
			ju = i4_min ( i4_max ( ju, mu+ipvt[k-1] ), n );
			mm = m;

			for ( j = k+1; j <= ju; j++ ){
				l = l - 1;
				mm = mm - 1;
				t = abd[l-1+(j-1)*lda];
				if ( l != mm ){
					abd[l-1+(j-1)*lda] = abd[mm-1+(j-1)*lda];
					abd[mm-1+(j-1)*lda] = t;
				}
				daxpy ( lm, t, abd+m+(k-1)*lda, 1, abd+mm+(j-1)*lda, 1 );
			}
		}
	}

	ipvt[n-1] = n;

	if ( abd[m-1+(n-1)*lda] == ZERO ){	info = n;	}

	return info;
}
//****************************************************************************80
void dgbsl ( double abd[], int lda, int n, int ml, int mu, int ipvt[],
		double b[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DGBSL solves a real banded system factored by DGBCO or DGBFA.
//
//  Discussion:
//
//    DGBSL can solve either A * X = B  or  A' * X = B.
//
//    A division by zero will occur if the input factor contains a
//    zero on the diagonal.  Technically this indicates singularity
//    but it is often caused by improper arguments or improper
//    setting of LDA.  It will not occur if the subroutines are
//    called correctly and if DGBCO has set 0.0 < RCOND
//    or DGBFA has set INFO == 0.
//
//    To compute inverse(A) * C  where C is a matrix with P columns:
//
//      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
//
//      if ( rcond is too small ) then
//        exit
//      end if
//
//      do j = 1, p
//        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double ABD[LDA*N], the output from DGBCO or DGBFA.
//
//    Input, integer LDA, the leading dimension of the array ABD.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, MU, the number of diagonals below and above the
//    main diagonal.  0 <= ML < N, 0 <= MU < N.
//
//    Input, int IPVT[N], the pivot vector from DGBCO or DGBFA.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
//    Input, int JOB, job choice.
//    0, solve A*X=B.
//    nonzero, solve A'*X=B.
//
{
	int k, l, la, lb, lm, m = mu + ml + 1;
	double t;
	//  JOB = 0, Solve A * x = b.
	//  First solve L * y = b.
	if ( job == 0 ){
		if ( 0 < ml ){
			for ( k = 1; k <= n-1; k++ ){
				lm = i4_min ( ml, n-k );
				l = ipvt[k-1];
				t = b[l-1];
				if ( l != k ){
					b[l-1] = b[k-1];
					b[k-1] = t;
				}
				daxpy ( lm, t, abd+m+(k-1)*lda, 1, b+k, 1 );
			}
		}
		//  Now solve U * x = y.
		for ( k = n; 1 <= k; k-- ){
			b[k-1] = b[k-1] / abd[m-1+(k-1)*lda];
			lm = i4_min ( k, m ) - 1;
			la = m - lm;
			lb = k - lm;
			t = -b[k-1];
			daxpy ( lm, t, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
		}
	}
	//  JOB nonzero, solve A' * x = b.
	//  First solve U' * y = b.
	else{
		for ( k = 1; k <= n; k++ ){
			lm = i4_min ( k, m ) - 1;
			la = m - lm;
			lb = k - lm;
			t = ddot ( lm, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
			b[k-1] = ( b[k-1] - t ) / abd[m-1+(k-1)*lda];
		}
		//  Now solve L' * x = y.
		if ( 0 < ml ){
			for ( k = n-1; 1 <= k; k-- ){
				lm = i4_min ( ml, n-k );
				b[k-1] = b[k-1] + ddot ( lm, abd+m+(k-1)*lda, 1, b+k, 1 );
				l = ipvt[k-1];
				if ( l != k ){
					t = b[l-1];
					b[l-1] = b[k-1];
					b[k-1] = t;
				}
			}
		}
	}
	return;
}
//****************************************************************************80
double dgeco ( double a[], int lda, int n, int ipvt[], double z[] )
//****************************************************************************80
//
//  Purpose:
//
//    DGECO factors a real matrix and estimates its condition number.
//
//  Discussion:
//
//    If RCOND is not needed, DGEFA is slightly faster.
//
//    To solve A * X = B, follow DGECO by DGESL.
//
//    To compute inverse ( A ) * C, follow DGECO by DGESL.
//
//    To compute determinant ( A ), follow DGECO by DGEDI.
//
//    To compute inverse ( A ), follow DGECO by DGEDI.
//
//    For the system A * X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    EPSILON/RCOND.
//
//    If RCOND is so small that the logical expression
//      1.0D+00 + RCOND == 1.0D+00
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate
//    underflows.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, a matrix to be
//    factored.  On output, the LU factorization of the matrix.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix A.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, double Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
//
//    Output, double DGECO, the value of RCOND, an estimate
//    of the reciprocal condition number of A.
//
{
	double anorm, ek, rcond, s, sm, t, wk, wkm, ynorm;
	int i, info, j, k, l;
	//  Compute the L1 norm of A.
	anorm = ZERO;
	for ( j = 1; j <= n; j++ ){		anorm = r8_max ( anorm, dasum ( n, a+0+(j-1)*lda, 1 ) );	}
	//  Compute the LU factorization.
	info = dgefa ( a, lda, n, ipvt );
	/*!
	 * RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )
	 *
	 * estimate of norm(inverse(A)) = norm(Z) / norm(Y)
	 *
	 * where
	 * 		A * Z = Y
	 * and
	 * 		A' * Y = E
	 *
	 * The components of E are chosen to cause maximum local growth in the
	 * elements of W, where U'*W = E.  The vectors are frequently rescaled
	 * to avoid overflow.
	 *
	 * Solve U' * W = E.
	 */
	ek = ONE;
	for ( i = 1; i <= n; i++ ){		z[i-1] = ZERO;	}

	for ( k = 1; k <= n; k++ ){
		if ( z[k-1] != ZERO ){		ek = ek * r8_sign ( -z[k-1] );		}

		if ( r8_abs ( a[k-1+(k-1)*lda] ) < r8_abs ( ek - z[k-1] ) ){
			s = r8_abs ( a[k-1+(k-1)*lda] ) / r8_abs ( ek - z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];		}
			ek = s * ek;
		}

		wk = ek - z[k-1];
		wkm = -ek - z[k-1];
		s = r8_abs ( wk );
		sm = r8_abs ( wkm );

		if ( a[k-1+(k-1)*lda] != ZERO ){
			wk = wk / a[k-1+(k-1)*lda];
			wkm = wkm / a[k-1+(k-1)*lda];
		}else{
			wk = ONE;
			wkm = ONE;
		}

		if ( k+1 <= n ){
			for ( j = k+1; j <= n; j++ ){
				sm = sm + r8_abs ( z[j-1] + wkm * a[k-1+(j-1)*lda] );
				z[j-1] = z[j-1] + wk * a[k-1+(j-1)*lda];
				s = s + r8_abs ( z[j-1] );
			}

			if ( s < sm ){
				t = wkm - wk;
				wk = wkm;
				for ( i = k+1; i <= n; i++ ){	z[i-1] = z[i-1] + t * a[k-1+(i-1)*lda];		}
			}
		}
		z[k-1] = wk;
	}

	t = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / t;	}
	//  Solve L' * Y = W
	for ( k = n; 1 <= k; k-- ){
		z[k-1] = z[k-1] + ddot ( n - k, a+k+(k-1)*lda, 1, z+k, 1 );

		if ( ONE < r8_abs ( z[k-1] ) ){
			t = r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / t;	}
		}

		l = ipvt[k-1];

		t = z[l-1];
		z[l-1] = z[k-1];
		z[k-1] = t;
	}
	t = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / t;	}

	ynorm = ONE;
	//  Solve L * V = Y.
	for ( k = 1; k <= n; k++ ){
		l = ipvt[k-1];

		t = z[l-1];
		z[l-1] = z[k-1];
		z[k-1] = t;

		for ( i = k+1; i <= n; i++ ){	z[i-1] = z[i-1] + t * a[i-1+(k-1)*lda];		}

		if ( ONE < r8_abs ( z[k-1] ) ){
			ynorm = ynorm / r8_abs ( z[k-1] );
			t = r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / t;	}
		}
	}
	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	ynorm = ynorm / s;
	//  Solve U * Z = V.
	for ( k = n; 1 <= k; k-- ){
		if ( r8_abs ( a[k-1+(k-1)*lda] ) < r8_abs ( z[k-1] ) ){
			s = r8_abs ( a[k-1+(k-1)*lda] ) / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ynorm = s * ynorm;
		}

		if ( a[k-1+(k-1)*lda] != ZERO ){		z[k-1] = z[k-1] / a[k-1+(k-1)*lda];		}
		else{			z[k-1] = ONE;		}
		for ( i = 1; i <= k-1; i++ ){	z[i-1] = z[i-1] - z[k-1] * a[i-1+(k-1)*lda];	}
	}
	//  Normalize Z in the L1 norm.
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;

	if ( anorm != ZERO ){	rcond = ynorm / anorm;	}
	else{		rcond = ZERO;	}

	return rcond;
}
//****************************************************************************80
void dgedi ( double a[], int lda, int n, int ipvt[], double det[],
		double work[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DGEDI computes the determinant and inverse of a matrix factored by DGECO or DGEFA.
//
//  Discussion:
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal and the inverse is requested.
//    It will not occur if the subroutines are called correctly
//    and if DGECO has set 0.0 < RCOND or DGEFA has set INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N], on input, the LU factor information,
//    as output by DGECO or DGEFA.  On output, the inverse
//    matrix if requested.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix A.
//
//    Input, int IPVT[N], the pivot vector from DGECO or DGEFA.
//
//    Workspace, double WORK[N].
//
//    Output, double DET[2], the determinant of original matrix if
//    requested.  The determinant = DET[0] * pow ( 10.0, DET[1] )
//    with  1.0 <= abs ( DET[0] ) < 10.0 or DET[0] == 0.0.
//
//    Input, int JOB, specifies what is to be computed.
//    11, both determinant and inverse.
//    01, inverse only.
//    10, determinant only.
//
{
	int i, j, k, l;
	double t;
	//  Compute the determinant.
	if ( job / 10 != 0 ){
		det[0] = ONE;
		det[1] = ZERO;

		for ( i = 1; i <= n; i++ ){
			if ( ipvt[i-1] != i ){		det[0] = -det[0];		}
			det[0] = det[0] * a[i-1+(i-1)*lda];

			if ( det[0] == ZERO ){		break;		}

			while ( r8_abs ( det[0] ) < ONE ){
				det[0] = det[0] * static_cast<Real>(10.0);
				det[1] = det[1] - ONE;
			}
			while ( static_cast<Real>(10.0) <= r8_abs ( det[0] ) ){
				det[0] = det[0] / static_cast<Real>(10.0);
				det[1] = det[1] + ONE;
			}
		}
	}
	//  Compute inverse(U).
	if ( ( job % 10 ) != 0 ){
		for ( k = 1; k <= n; k++ ){
			a[k-1+(k-1)*lda] = ONE / a[k-1+(k-1)*lda];
			t = -a[k-1+(k-1)*lda];
			dscal ( k-1, t, a+0+(k-1)*lda, 1 );

			for ( j = k+1; j <= n; j++ ){
				t = a[k-1+(j-1)*lda];
				a[k-1+(j-1)*lda] = ZERO;
				daxpy ( k, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
			}
		}
		//  Form inverse(U) * inverse(L).
		for ( k = n-1; 1 <= k; k-- ){
			for ( i = k+1; i <= n; i++ ){
				work[i-1] = a[i-1+(k-1)*lda];
				a[i-1+(k-1)*lda] = ZERO;
			}

			for ( j = k+1; j <= n; j++ ){
				t = work[j-1];
				daxpy ( n, t, a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
			}

			l = ipvt[k-1];
			if ( l != k ){		dswap ( n, a+0+(k-1)*lda, 1, a+0+(l-1)*lda, 1 );		}
		}
	}
	return;
}
//****************************************************************************80
int dgefa ( double a[], int lda, int n, int ipvt[] )
//****************************************************************************80
//
//  Purpose:
//
//    DGEFA factors a real general matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].
//    On intput, the matrix to be factored.
//    On output, an upper triangular matrix and the multipliers used to obtain
//    it.  The factorization can be written A=L*U, where L is a product of
//    permutation and unit lower triangular matrices, and U is upper triangular.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix A.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int DGEFA, singularity indicator.
//    0, normal value.
//    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
//    but it does indicate that DGESL or DGEDI will divide by zero if called.
//    Use RCOND in DGECO for a reliable indication of singularity.
//
{
	int info = 0, j, k, l;
	double t;
	//  Gaussian elimination with partial pivoting.
	for ( k = 1; k <= n-1; k++ ){
		//  Find L = pivot index.
		l = idamax ( n-k+1, a+(k-1)+(k-1)*lda, 1 ) + k - 1;
		ipvt[k-1] = l;
		//  Zero pivot implies this column already triangularized.
		if ( a[l-1+(k-1)*lda] == ZERO ){
			info = k;
			continue;
		}
		//  Interchange if necessary.
		if ( l != k ){
			t = a[l-1+(k-1)*lda];
			a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
			a[k-1+(k-1)*lda] = t;
		}
		//  Compute multipliers.
		t = -ONE / a[k-1+(k-1)*lda];

		dscal ( n-k, t, a+k+(k-1)*lda, 1 );
		//  Row elimination with column indexing.
		for ( j = k+1; j <= n; j++ ){
			t = a[l-1+(j-1)*lda];
			if ( l != k ){
				a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
				a[k-1+(j-1)*lda] = t;
			}
			daxpy ( n-k, t, a+k+(k-1)*lda, 1, a+k+(j-1)*lda, 1 );
		}
	}

	ipvt[n-1] = n;

	if ( a[n-1+(n-1)*lda] == ZERO ){		info = n;	}

	return info;
}
//****************************************************************************80
void dgesl ( double a[], int lda, int n, int ipvt[], double b[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DGESL solves a real general linear system A * X = B.
//
//  Discussion:
//
//    DGESL can solve either of the systems A * X = B or A' * X = B.
//
//    The system matrix must have been factored by DGECO or DGEFA.
//
//    A division by zero will occur if the input factor contains a
//    zero on the diagonal.  Technically this indicates singularity
//    but it is often caused by improper arguments or improper
//    setting of LDA.  It will not occur if the subroutines are
//    called correctly and if DGECO has set 0.0 < RCOND
//    or DGEFA has set INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double A[LDA*N], the output from DGECO or DGEFA.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix A.
//
//    Input, int IPVT[N], the pivot vector from DGECO or DGEFA.
//
//    Input/output, double B[N].
//    On input, the right hand side vector.
//    On output, the solution vector.
//
//    Input, int JOB.
//    0, solve A * X = B;
//    nonzero, solve A' * X = B.
//
{
	int k, l;
	double t;
	//  Solve A * X = B.
	if ( job == 0 ){
		for ( k = 1; k <= n-1; k++ ){
			l = ipvt[k-1];
			t = b[l-1];

			if ( l != k ){
				b[l-1] = b[k-1];
				b[k-1] = t;
			}
			daxpy ( n-k, t, a+k+(k-1)*lda, 1, b+k, 1 );
		}

		for ( k = n; 1 <= k; k-- ){
			b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
			t = -b[k-1];
			daxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
		}
	}else{	//  Solve A' * X = B.
		for ( k = 1; k <= n; k++ ){
			t = ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
			b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
		}

		for ( k = n-1; 1 <= k; k-- ){
			b[k-1] = b[k-1] + ddot ( n-k, a+k+(k-1)*lda, 1, b+k, 1 );
			l = ipvt[k-1];

			if ( l != k ){
				t = b[l-1];
				b[l-1] = b[k-1];
				b[k-1] = t;
			}
		}
	}
	return;
}
//****************************************************************************80
int dgtsl ( int n, double c[], double d[], double e[], double b[] )
//****************************************************************************80
//
//  Purpose:
//
//    DGTSL solves a general tridiagonal linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, int N, the order of the tridiagonal matrix.
//
//    Input/output, double C[N], contains the subdiagonal of the
//    tridiagonal matrix in entries C(2:N).  On output, C is destroyed.
//
//    Input/output, double D[N].  On input, the diagonal of the
//    matrix.  On output, D is destroyed.
//
//    Input/output, double E[N], contains the superdiagonal of the
//    tridiagonal matrix in entries E(1:N-1).  On output E is destroyed.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
//    Output, int DGTSL, error flag.
//    0, normal value.
//    K, the K-th element of the diagonal becomes exactly zero.  The
//       subroutine returns if this error condition is detected.
//
{
	int info = 0, k;
	double t;

	c[0] = d[0];

	if ( 2 <= n ){
		d[0] = e[0];
		e[0] = e[n-1] = ZERO;

		for ( k = 1; k <= n - 1; k++ ){
			//  Find the larger of the two rows.
			if ( r8_abs ( c[k-1] ) <= r8_abs ( c[k] ) ){
				//  Interchange rows.
				t = c[k];		c[k] = c[k-1];			c[k-1] = t;
				t = d[k];		d[k] = d[k-1];			d[k-1] = t;
				t = e[k];		e[k] = e[k-1];			e[k-1] = t;
				t = b[k];		b[k] = b[k-1];			b[k-1] = t;
			}
			//  Zero elements.
			if ( c[k-1] == ZERO ){
				info = k;
				return info;
			}

			t = -c[k] / c[k-1];
			c[k] = d[k] + t * d[k-1];
			d[k] = e[k] + t * e[k-1];
			e[k] = ZERO;
			b[k] = b[k] + t * b[k-1];
		}
	}

	if ( c[n-1] == ZERO ){
		info = n;
		return info;
	}
	//  Back solve.
	b[n-1] = b[n-1] / c[n-1];

	if ( 1 < n ){
		b[n-2] = ( b[n-2] - d[n-2] * b[n-1] ) / c[n-2];

		for ( k = n-2; 1 <= k; k-- ){	b[k-1] = ( b[k-1] - d[k-1] * b[k] - e[k-1] * b[k+1] ) / c[k-1];		}

	}
	return info;
}
//****************************************************************************80
double dpbco ( double abd[], int lda, int n, int m, double z[] )
//****************************************************************************80
//
//  Purpose:
//
//    DPBCO factors a real symmetric positive definite banded matrix.
//
//  Discussion:
//
//    DPBCO also estimates the condition of the matrix.
//
//    If RCOND is not needed, DPBFA is slightly faster.
//
//    To solve A*X = B, follow DPBCO by DPBSL.
//
//    To compute inverse(A)*C, follow DPBCO by DPBSL.
//
//    To compute determinant(A), follow DPBCO by DPBDI.
//
//  Band storage:
//
//    If A is a symmetric positive definite band matrix, the following
//    program segment will set up the input.
//
//      m = (band width above diagonal)
//      do j = 1, n
//        i1 = max (1, j-m)
//        do i = i1, j
//          k = i-j+m+1
//          abd(k,j) = a(i,j)
//        }
//      }
//
//    This uses M + 1 rows of A, except for the M by M upper left triangle,
//    which is ignored.
//
//    For example, if the original matrix is
//
//      11 12 13  0  0  0
//      12 22 23 24  0  0
//      13 23 33 34 35  0
//       0 24 34 44 45 46
//       0  0 35 45 55 56
//       0  0  0 46 56 66
//
//    then N = 6, M = 2  and ABD should contain
//
//       *  * 13 24 35 46
//       * 12 23 34 45 56
//      11 22 33 44 55 66
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double ABD[LDA*N].  On input, the matrix to be
//    factored.  The columns of the upper triangle are stored in the columns
//    of ABD and the diagonals of the upper triangle are stored in the rows
//    of ABD.  On output, an upper triangular matrix R, stored in band form,
//    so that A = R'*R.  If INFO /= 0, the factorization is not complete.
//
//    Input, int LDA, the leading dimension of the array ABD.
//    M+1 <= LDA is required.
//
//    Input, int N, the order of the matrix.
//
//    Input, int M, the number of diagonals above the main diagonal.
//
//    Output, double Z[N], a work vector whose contents are usually
//    unimportant.  If A is singular to working precision, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//    If INFO /= 0, Z is unchanged.
//
//    Output, double DPBCO, an estimate of the reciprocal condition number
//    RCOND.  For the system A*X = B, relative perturbations in A and B of size
//    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
//    If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0D+00
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
{
	int i, info, j, j2, k, l, la, lb, lm, mu;
	double anorm, ek, rcond, s, sm, t, wk, wkm, ynorm;
	//  Find the norm of A.
	for ( j = 1; j <= n; j++ ){
		l = i4_min ( j, m+1 );
		mu = i4_max ( m+2-j, 1 );
		z[j-1] = dasum ( l, abd+mu-1+(j-1)*lda, 1 );
		k = j - l;
		for ( i = mu; i <= m; i++ ){
			k = k + 1;
			z[k-1] = z[k-1] + r8_abs ( abd[i-1+(j-1)*lda] );
		}
	}
	anorm = ZERO;
	for ( i = 1; i <= n; i++ ){	anorm = r8_max ( anorm, z[i-1] );	}
	//  Factor.
	info = dpbfa ( abd, lda, n, m );

	if ( info != 0 ){
		rcond = ZERO;
		return rcond;
	}
	/*!
	 * RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
	 *
	 * Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
	 *
	 * The components of E are chosen to cause maximum local
	 * growth in the elements of W where R'*W = E.
	 *
	 * The vectors are frequently rescaled to avoid overflow.
	 *
	 * Solve R' * W = E.
	 */
	ek = ONE;
	for ( i = 1; i <= n; i++ ){		z[i-1] = ZERO;	}

	for ( k = 1; k <= n; k++ ){
		if ( z[k-1] != ZERO ){		ek = ek * r8_sign ( -z[k-1] );	}

		if ( abd[m+(k-1)*lda] < r8_abs ( ek - z[k-1] ) ){
			s = abd[m+(k-1)*lda] / r8_abs ( ek - z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ek = s * ek;
		}
		wk = ek - z[k-1];
		wkm = -ek - z[k-1];
		s = r8_abs ( wk );
		sm = r8_abs ( wkm );
		wk = wk / abd[m+(k-1)*lda];
		wkm = wkm / abd[m+(k-1)*lda];
		j2 = i4_min ( k+m, n );
		i = m + 1;

		if ( k+1 <= j2 ){
			for ( j = k+1; j <= j2; j++ ){
				i = i - 1;
				sm = sm + r8_abs ( z[j-1] + wkm * abd[i-1+(j-1)*lda] );
				z[j-1] = z[j-1] + wk * abd[i-1+(j-1)*lda];
				s = s + r8_abs ( z[j-1] );
			}

			if ( s < sm ){
				t = wkm - wk;
				wk = wkm;
				i = m + 1;

				for ( j = k+1; j <= j2; j++ ){
					i = i - 1;
					z[j-1] = z[j-1] + t * abd[i-1+(j-1)*lda];
				}
			}
		}
		z[k-1] = wk;
	}
	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){	z[i-1] = z[i-1] / s;	}
	//  Solve R * Y = W.
	for ( k = n; 1 <= k; k-- ){
		if ( abd[m+(k-1)*lda] < r8_abs ( z[k-1] ) )	{
			s = abd[m+(k-1)*lda] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
		}
		z[k-1] = z[k-1] / abd[m+(k-1)*lda];
		lm = i4_min ( k-1, m );
		la = m + 1 - lm;
		lb = k - lm;
		t = -z[k-1];
		daxpy ( lm, t, abd+la-1+(k-1)*lda, 1, z+lb-1, 1 );
	}
	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	ynorm = ONE;
	//  Solve R' * V = Y.
	for ( k = 1; k <= n; k++ ){
		lm = i4_min ( k-1, m );
		la = m + 1 - lm;
		lb = k - lm;

		z[k-1] = z[k-1] - ddot ( lm, abd+la-1+(k-1)*lda, 1, z+lb-1, 1 );

		if ( abd[m+(k-1)*lda] < r8_abs ( z[k-1] ) ){
			s = abd[m+(k-1)*lda] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ynorm = s * ynorm;
		}
		z[k-1] = z[k-1] / abd[m+(k-1)*lda];
	}
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;
	//  Solve R * Z = W.
	for ( k = n; 1 <= k; k-- ){
		if ( abd[m+(k-1)*lda] < r8_abs ( z[k-1] ) ){
			s = abd[m+(k-1)*lda] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){
				z[i-1] = s * z[i-1];
			}
			ynorm = s * ynorm;
		}
		z[k-1] = z[k-1] / abd[m+(k-1)*lda];
		lm = i4_min ( k-1, m );
		la = m + 1 - lm;
		lb = k - lm;
		t = -z[k-1];
		daxpy ( lm, t, abd+la-1+(k-1)*lda, 1, z+lb-1, 1 );
	}
	//  Make ZNORM = 1.0.
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;

	if ( anorm != ZERO ){	rcond = ynorm / anorm;	}
	else{	rcond = ZERO;	}

	return rcond;
}
//****************************************************************************80
void dpbdi ( double abd[], int lda, int n, int m, double det[] )
//****************************************************************************80
//
//  Purpose:
//
//    DPBDI computes the determinant of a matrix factored by DPBCO or DPBFA.
//
//  Discussion:
//
//    If the inverse is needed, use DPBSL N times.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double ABD[LDA*N], the output from DPBCO or DPBFA.
//
//    Input, int LDA, the leading dimension of the array ABD.
//
//    Input, int N, the order of the matrix.
//
//    Input, int M, the number of diagonals above the main diagonal.
//
//    Output, double DET[2], the determinant of the original
//    matrix in the form
//      determinant = DET[0] * 10.0**DET[1]
//    with 1.0D+00 <= DET[0] < 10.0D+00 or DET[0] == 0.0D+00.
//
{
	int i;
	double s;
	//  Compute the determinant.
	det[0] = ONE;
	det[1] = ZERO;
	s = static_cast<Real>(10.0);

	for ( i = 1; i <= n; i++ ){
		det[0] = det[0] * abd[m+(i-1)*lda] * abd[m+(i-1)*lda];

		if ( det[0] == ZERO ){		return;		}

		while ( det[0] < ONE ){
			det[0] = det[0] * s;
			det[1] = det[1] - ONE;
		}

		while ( s <= det[0] ){
			det[0] = det[0] / s;
			det[1] = det[1] + ONE;
		}
	}
	return;
}
//****************************************************************************80
int dpbfa ( double abd[], int lda, int n, int m )
//****************************************************************************80
//
//  Purpose:
//
//    DPBFA factors a symmetric positive definite matrix stored in band form.
//
//  Discussion:
//
//    DPBFA is usually called by DPBCO, but it can be called
//    directly with a saving in time if RCOND is not needed.
//
//    If A is a symmetric positive definite band matrix,
//    the following program segment will set up the input.
//
//      m = (band width above diagonal)
//      do j = 1, n
//        i1 = max ( 1, j-m )
//        do i = i1, j
//          k = i-j+m+1
//          abd(k,j) = a(i,j)
//        end do
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double ABD[LDA*N].  On input, the matrix to be
//    factored.  The columns of the upper triangle are stored in the columns
//    of ABD and the diagonals of the upper triangle are stored in the
//    rows of ABD.  On output, an upper triangular matrix R, stored in band
//    form, so that A = R' * R.
//
//    Input, int LDA, the leading dimension of the array ABD.
//    M+1 <= LDA is required.
//
//    Input, int N, the order of the matrix.
//
//    Input, int M, the number of diagonals above the main diagonal.
//
//    Output, int DPBFA, error indicator.
//    0, for normal return.
//    K, if the leading minor of order K is not positive definite.
//
{
	int ik, info, j, jk, k, mu;
	double s, t;

	for ( j = 1; j <= n; j++ ){
		s = ZERO;
		ik = m + 1;
		jk = i4_max ( j - m, 1 );
		mu = i4_max ( m + 2 - j, 1 );

		for ( k = mu; k <= m; k++ ){
			t = abd[k-1+(j-1)*lda]
					- ddot ( k-mu, abd+ik-1+(jk-1)*lda, 1, abd+mu-1+(j-1)*lda, 1 );
			t = t / abd[m+(jk-1)*lda];
			abd[k-1+(j-1)*lda] = t;
			s = s + t * t;
			ik = ik - 1;
			jk = jk + 1;
		}
		s = abd[m+(j-1)*lda] - s;

		if ( s <= ZERO ){
			info = j;
			return info;
		}
		abd[m+(j-1)*lda] = sqrt ( s );
	}
	info = 0;

	return info;
}
//****************************************************************************80
void dpbsl ( double abd[], int lda, int n, int m, double b[] )
//****************************************************************************80
//
//  Purpose:
//
//    DPBSL solves a real SPD band system factored by DPBCO or DPBFA.
//
//  Discussion:
//
//    The matrix is assumed to be a symmetric positive definite (SPD)
//    band matrix.
//
//    To compute inverse(A) * C  where C is a matrix with P columns:
//
//      call dpbco ( abd, lda, n, rcond, z, info )
//
//      if ( rcond is too small .or. info /= 0) go to ...
//
//      do j = 1, p
//        call dpbsl ( abd, lda, n, c(1,j) )
//      end do
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal.  Technically this indicates
//    singularity but it is usually caused by improper subroutine
//    arguments.  It will not occur if the subroutines are called
//    correctly and INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double ABD[LDA*N], the output from DPBCO or DPBFA.
//
//    Input, int LDA, the leading dimension of the array ABD.
//
//    Input, int N, the order of the matrix.
//
//    Input, int M, the number of diagonals above the main diagonal.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
{
	int k, la, lb, lm;
	double t;
	//  Solve R'*Y = B.
	for ( k = 1; k <= n; k++ ){
		lm = i4_min ( k-1, m );
		la = m + 1 - lm;
		lb = k - lm;
		t = ddot ( lm, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
		b[k-1] = ( b[k-1] - t ) / abd[m+(k-1)*lda];
	}
	//  Solve R*X = Y.
	for ( k = n; 1 <= k; k-- ){
		lm = i4_min ( k-1, m );
		la = m + 1 - lm;
		lb = k - lm;
		b[k-1] = b[k-1] / abd[m+(k-1)*lda];
		t = -b[k-1];
		daxpy ( lm, t, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
	}
	return;
}
//****************************************************************************80
double dpoco ( double a[], int lda, int n, double z[] )
//****************************************************************************80
//
//  Purpose:
//
//    DPOCO factors a real symmetric positive definite matrix and estimates its condition.
//
//  Discussion:
//
//    If RCOND is not needed, DPOFA is slightly faster.
//
//    To solve A*X = B, follow DPOCO by DPOSL.
//
//    To compute inverse(A)*C, follow DPOCO by DPOSL.
//
//    To compute determinant(A), follow DPOCO by DPODI.
//
//    To compute inverse(A), follow DPOCO by DPODI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the symmetric
//    matrix to be factored.  Only the diagonal and upper triangle are used.
//    On output, an upper triangular matrix R so that A = R'*R where R'
//    is the transpose.  The strict lower triangle is unaltered.
//    If INFO /= 0, the factorization is not complete.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Output, double Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//    If INFO /= 0, Z is unchanged.
//
//    Output, double DPOCO, an estimate of the reciprocal
//    condition of A.  For the system A*X = B, relative perturbations in
//    A and B of size EPSILON may cause relative perturbations in X of
//    size EPSILON/RCOND.  If RCOND is so small that the logical expression
//      1.0D+00 + RCOND == 1.0D+00
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
{
	double anorm, ek, rcond, s, sm, t, wk, wkm, ynorm;
	int i, info, j, k;
	//  Find norm of A using only upper half.
	for ( j = 1; j <= n; j++ ){
		z[j-1] = dasum ( j, a+0+(j-1)*lda, 1 );
		for ( i = 1; i <= j-1; i++ ){		z[i-1] = z[i-1] + r8_abs ( a[i-1+(j-1)*lda] );		}
	}

	anorm = ZERO;
	for ( i = 1; i <= n; i++ ){		anorm = r8_max ( anorm, z[i-1] );	}
	//  Factor.
	info = dpofa ( a, lda, n );

	if ( info != 0 ){
		rcond = ZERO;
		return rcond;
	}
	/*!
	 * RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
	 *
	 * Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
	 *
	 * The components of E are chosen to cause maximum local
	 * growth in the elements of W where R'*W = E.
	 *
	 * The vectors are frequently rescaled to avoid overflow.
	 *
	 * Solve R' * W = E.
	 */
	ek = ONE;
	for ( i = 1; i <= n; i++ ){		z[i-1] = ZERO;	}

	for ( k = 1; k <= n; k++ ){
		if ( z[k-1] != ZERO ){		ek = ek * r8_sign ( -z[k-1] );		}

		if ( a[k-1+(k-1)*lda] < r8_abs ( ek - z[k-1] ) ){
			s = a[k-1+(k-1)*lda] / r8_abs ( ek - z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ek = s * ek;
		}

		wk = ek - z[k-1];
		wkm = -ek - z[k-1];
		s = r8_abs ( wk );
		sm = r8_abs ( wkm );
		wk = wk / a[k-1+(k-1)*lda];
		wkm = wkm / a[k-1+(k-1)*lda];

		if ( k + 1 <= n ){
			for ( j = k+1; j <= n; j++ ){
				sm = sm + r8_abs ( z[j-1] + wkm * a[k-1+(j-1)*lda] );
				z[j-1] = z[j-1] + wk * a[k-1+(j-1)*lda];
				s = s + r8_abs ( z[j-1] );
			}

			if ( s < sm ){
				t = wkm - wk;
				wk = wkm;
				for ( j = k+1; j <= n; j++ ){	z[j-1] = z[j-1] + t * a[k-1+(j-1)*lda];		}
			}
		}
		z[k-1] = wk;
	}
	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	//  Solve R * Y = W.
	for ( k = n; 1 <= k; k-- ){
		if ( a[k-1+(k-1)*lda] < r8_abs ( z[k-1] ) )	{
			s = a[k-1+(k-1)*lda] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
		}
		z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
		t = -z[k-1];
		daxpy ( k-1, t, a+0+(k-1)*lda, 1, z, 1 );
	}
	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	ynorm = ONE;
	//  Solve R' * V = Y.
	for ( k = 1; k <= n; k++ ){
		z[k-1] = z[k-1] - ddot ( k-1, a+0+(k-1)*lda, 1, z, 1 );

		if ( a[k-1+(k-1)*lda] < r8_abs ( z[k-1] ) ){
			s = a[k-1+(k-1)*lda] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ynorm = s * ynorm;
		}
		z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
	}
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;
	//  Solve R * Z = V.
	for ( k = n; 1 <= k; k-- ){
		if ( a[k-1+(k-1)*lda] < r8_abs ( z[k-1] ) ){
			s = a[k-1+(k-1)*lda] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ynorm = s * ynorm;
		}
		z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
		t = -z[k-1];
		daxpy ( k-1, t, a+0+(k-1)*lda, 1, z, 1 );
	}
	//  Make ZNORM = 1.0.
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;

	if ( anorm != ZERO ){	rcond = ynorm / anorm;	}
	else{		rcond = ZERO;	}

	return rcond;
}
//****************************************************************************80
void dpodi ( double a[], int lda, int n, double det[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DPODI computes the determinant and inverse of a certain matrix.
//
//  Discussion:
//
//    The matrix is real symmetric positive definite.
//    DPODI uses the factors computed by DPOCO, DPOFA or DQRDC.
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal and the inverse is requested.
//    It will not occur if the subroutines are called correctly
//    and if DPOCO or DPOFA has set INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the output A from
//    DPOCO or DPOFA, or the output X from DQRDC.  On output, if DPOCO or
//    DPOFA was used to factor A then DPODI produces the upper half of
//    inverse(A).  If DQRDC was used to decompose X then DPODI produces
//    the upper half of inverse(X'*X) where X' is the transpose.
//    Elements of A below the diagonal are unchanged.  If the units digit
//    of JOB is zero, A is unchanged.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix A.
//
//    Input, int JOB, specifies the task.
//    11, both determinant and inverse.
//    01, inverse only.
//    10, determinant only.
//
//    Output, double DET[2], the determinant of A or of X'*X
//    if requested.
//      determinant = DET[0] * 10.0**DET[1]
//    with 1.0D+00 <= DET[0] < 10.0D+00 or DET[0] == 0.0D+00.
//
{
	int i, j, k;
	double s, t;
	//  Compute the determinant.
	if ( job / 10 != 0 ){
		det[0] = ONE;
		det[1] = ZERO;
		s = static_cast<Real>(10.0);

		for ( i = 1; i <= n; i++ ){
			det[0] = det[0] * a[i-1+(i-1)*lda] * a[i-1+(i-1)*lda];

			if ( det[0] == ZERO ){		break;		}

			while ( det[0] < ONE ){
				det[0] = det[0] * s;
				det[1] = det[1] - ONE;
			}

			while ( s <= det[0] ){
				det[0] = det[0] / s;
				det[1] = det[1] + ONE;
			}
		}
	}
	//  Compute inverse(R).
	if ( ( job % 10 ) != 0 ){
		for ( k = 1; k <= n; k++ ){
			a[k-1+(k-1)*lda] = ONE / a[k-1+(k-1)*lda];
			t = -a[k-1+(k-1)*lda];
			dscal ( k-1, t, a+0+(k-1)*lda, 1 );

			for ( j = k+1; j <= n; j++ )
			{
				t = a[k-1+(j-1)*lda];
				a[k-1+(j-1)*lda] = ZERO;
				daxpy ( k, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
			}
		}
		//  Form inverse(R) * (inverse(R))'.
		for ( j = 1; j <= n; j++ ){
			for ( k = 1; k <= j-1; k++ ){
				t = a[k-1+(j-1)*lda];
				daxpy ( k, t, a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
			}
			t = a[j-1+(j-1)*lda];
			dscal ( j, t, a+0+(j-1)*lda, 1 );
		}
	}
	return;
}
//****************************************************************************80
int dpofa ( double a[], int lda, int n )
//****************************************************************************80
//
//  Purpose:
//
//    DPOFA factors a real symmetric positive definite matrix.
//
//  Discussion:
//
//    DPOFA is usually called by DPOCO, but it can be called
//    directly with a saving in time if RCOND is not needed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the symmetric matrix
//    to be  factored.  Only the diagonal and upper triangle are used.
//    On output, an upper triangular matrix R so that A = R'*R
//    where R' is the transpose.  The strict lower triangle is unaltered.
//    If INFO /= 0, the factorization is not complete.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int DPOFA, error flag.
//    0, for normal return.
//    K, signals an error condition.  The leading minor of order K is not
//    positive definite.
//
{
	int info, j, k;
	double s, t;

	for ( j = 1; j <= n; j++ ){
		s = ZERO;

		for ( k = 1; k <= j-1; k++ ){
			t = a[k-1+(j-1)*lda] - ddot ( k-1, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
			t = t / a[k-1+(k-1)*lda];
			a[k-1+(j-1)*lda] = t;
			s = s + t * t;
		}

		s = a[j-1+(j-1)*lda] - s;

		if ( s <= ZERO ){	info = j;		return info;		}

		a[j-1+(j-1)*lda] = sqrt ( s );
	}

	info = 0;

	return info;
}
//****************************************************************************80
void dposl ( double a[], int lda, int n, double b[] )
//****************************************************************************80
//
//  Purpose:
//
//    DPOSL solves a linear system factored by DPOCO or DPOFA.
//
//  Discussion:
//
//    To compute inverse(A) * C where C is a matrix with P columns:
//
//      call dpoco ( a, lda, n, rcond, z, info )
//
//      if ( rcond is not too small .and. info == 0 ) then
//        do j = 1, p
//          call dposl ( a, lda, n, c(1,j) )
//        end do
//      end if
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal.  Technically this indicates
//    singularity but it is usually caused by improper subroutine
//    arguments.  It will not occur if the subroutines are called
//    correctly and INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double A[LDA*N], the output from DPOCO or DPOFA.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
{
	int k;
	double t;
	//  Solve R' * Y = B.
	for ( k = 1; k <= n; k++ ){
		t = ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
		b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
	}
	//  Solve R * X = Y.
	for ( k = n; 1 <= k; k-- ){
		b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
		t = -b[k-1];
		daxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
	}

	return;
}
//****************************************************************************80
double dppco ( double ap[], int n, double z[] )
//****************************************************************************80
//
//  Purpose:
//
//    DPPCO factors a real symmetric positive definite matrix in packed form.
//
//  Discussion:
//
//    DPPCO also estimates the condition of the matrix.
//
//    If RCOND is not needed, DPPFA is slightly faster.
//
//    To solve A*X = B, follow DPPCO by DPPSL.
//
//    To compute inverse(A)*C, follow DPPCO by DPPSL.
//
//    To compute determinant(A), follow DPPCO by DPPDI.
//
//    To compute inverse(A), follow DPPCO by DPPDI.
//
//  Packed storage:
//
//    The following program segment will pack the upper triangle of
//    a symmetric matrix.
//
//      k = 0
//      do j = 1, n
//        do i = 1, j
//          k = k + 1
//          ap[k-1] = a(i,j)
//        }
//      }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double AP[N*(N+1)/2].  On input, the packed
//    form of a symmetric matrix A.  The columns of the upper triangle are
//    stored sequentially in a one-dimensional array.  On output, an upper
//    triangular matrix R, stored in packed form, so that A = R'*R.
//    If INFO /= 0, the factorization is not complete.
//
//    Input, int N, the order of the matrix.
//
//    Output, double Z[N], a work vector whose contents are usually
//    unimportant.  If A is singular to working precision, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//    If INFO /= 0, Z is unchanged.
//
//    Output, double DPPCO, an estimate of the reciprocal condition number RCOND
//    of A.  For the system A*X = B, relative perturbations in A and B of size
//    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
//    If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0D+00
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
{
	double anorm, ek, rcond, s, sm, t, wk, wkm, ynorm;
	int i, ij, info, j, j1, k, kj, kk;
	//  Find the norm of A.
	j1 = 1;
	for ( j = 1; j <= n; j++ ){
		z[j-1] = dasum ( j, ap+j1-1, 1 );
		ij = j1;
		j1 = j1 + j;
		for ( i = 1; i <= j-1; i++ ){
			z[i-1] = z[i-1] + r8_abs ( ap[ij-1] );
			ij = ij + 1;
		}
	}
	anorm = ZERO;
	for ( i = 1; i <= n; i++ ){		anorm = r8_max ( anorm, z[i-1] );	}
	//  Factor.
	info = dppfa ( ap, n );

	if ( info != 0 ){		rcond = ZERO;		return rcond;	}
	/*!
	 * RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
	 *
	 * Estimate = norm(Z)/norm(Y) where A * Z = Y and A * Y = E.
	 *
	 * The components of E are chosen to cause maximum local
	 * growth in the elements of W where R'*W = E.
	 *
	 * The vectors are frequently rescaled to avoid overflow.
	 *
	 * Solve R' * W = E.
	 */
	ek = ONE;
	for ( i = 1; i <= n; i++ ){		z[i-1] = ZERO;	}

	kk = 0;

	for ( k = 1; k <= n; k++ ){
		kk = kk + k;

		if ( z[k-1] != ZERO ){	ek = ek * r8_sign ( -z[k-1] );	}

		if ( ap[kk-1] < r8_abs ( ek - z[k-1] ) ){
			s = ap[kk-1] / r8_abs ( ek - z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ek = s * ek;
		}
		wk = ek - z[k-1];
		wkm = -ek - z[k-1];
		s = r8_abs ( wk );
		sm = r8_abs ( wkm );
		wk = wk / ap[kk-1];
		wkm = wkm / ap[kk-1];
		kj = kk + k;

		if ( k + 1 <= n ){
			for ( j = k + 1; j <= n; j++ ){
				sm = sm + r8_abs ( z[j-1] + wkm * ap[kj-1] );
				z[j-1] = z[j-1] + wk * ap[kj-1];
				s = s + r8_abs ( z[j-1] );
				kj = kj + j;
			}

			if ( s < sm ){
				t = wkm - wk;
				wk = wkm;
				kj = kk + k;

				for ( j = k+1; j <= n; j++ ){
					z[j-1] = z[j-1] + t * ap[kj-1];
					kj = kj + j;
				}
			}
		}
		z[k-1] = wk;
	}

	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	//  Solve R * Y = W.
	for ( k = n; 1 <= k; k-- ){
		if ( ap[kk-1] < r8_abs ( z[k-1] ) ){
			s = ap[kk-1] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
		}
		z[k-1] = z[k-1] / ap[kk-1];
		kk = kk - k;
		t = -z[k-1];
		daxpy ( k-1, t, ap+kk, 1, z, 1 );
	}
	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	ynorm = ONE;
	//  Solve R' * V = Y.
	for ( k = 1; k <= n; k++ ){
		z[k-1] = z[k-1] - ddot ( k-1, ap+kk, 1, z, 1 );
		kk = kk + k;

		if ( ap[kk-1] < r8_abs ( z[k-1] ) ){
			s = ap[kk-1] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ynorm = s * ynorm;
		}
		z[k-1] = z[k-1] / ap[kk-1];
	}
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;
	//  Solve R * Z = V.
	for ( k = n; 1 <= k; k-- ){
		if ( ap[kk-1] < r8_abs ( z[k-1] ) ){
			s = ap[kk-1] / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
			ynorm = s * ynorm;
		}
		z[k-1] = z[k-1] / ap[kk-1];
		kk = kk - k;
		t = -z[k-1];
		daxpy ( k-1, t, ap+kk, 1, z, 1 );
	}
	//  Make ZNORM = 1.0.
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;

	if ( anorm != ZERO ){	rcond = ynorm / anorm;	}
	else{		rcond = ZERO;	}

	return rcond;
}
//****************************************************************************80
void dppdi ( double ap[], int n, double det[2], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DPPDI computes the determinant and inverse of a matrix factored by DPPCO or DPPFA.
//
//  Discussion:
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal and the inverse is requested.
//    It will not occur if the subroutines are called correctly
//    and if DPOCO or DPOFA has set INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double AP[N*(N+1)/2].  On input, the output from
//    DPPCO or DPPFA.  On output, the upper triangular half of the
//    inverse, if requested.
//
//    Input, int N, the order of the matrix.
//
//    Output, double DET[2], the determinant of the original matrix
//    if requested.
//      determinant = DET[0] * 10.0**DET[1]
//    with  1.0D+00 <= DET[0] < 10.0D+00 or DET[0] == 0.0D+00.
//
//    Input, int JOB, job request.
//    11, both determinant and inverse.
//    01, inverse only.
//    10, determinant only.
//
{
	int i, ii, j, j1, jj, k, k1, kj, kk;
	double s, t;
	//  Compute the determinant.
	if ( job / 10 != 0 ){
		det[0] = ONE;
		det[1] = ZERO;
		s = static_cast<Real>(10.0);
		ii = 0;

		for ( i = 1; i <= n; i++ ){
			ii = ii + i;

			det[0] = det[0] * ap[ii-1] * ap[ii-1];

			if ( det[0] == ZERO ){	break;	}

			while ( det[0] < ONE ){
				det[0] = det[0] * s;
				det[1] = det[1] - ONE;
			}

			while ( s <= det[0] ){
				det[0] = det[0] / s;
				det[1] = det[1] + ONE;
			}
		}
	}
	//  Compute inverse(R).
	if ( ( job % 10 ) != 0 ){
		kk = 0;

		for ( k = 1; k <= n; k++ ){
			k1 = kk + 1;
			kk = kk + k;
			ap[kk-1] = ONE / ap[kk-1];
			t = -ap[kk-1];
			dscal ( k-1, t, ap+k1-1, 1 );
			j1 = kk + 1;
			kj = kk + k;

			for ( j = k + 1; j <= n; j++ ){
				t = ap[kj-1];
				ap[kj-1] = ZERO;
				daxpy ( k, t, ap+k1-1, 1, ap+j1-1, 1 );
				j1 = j1 + j;
				kj = kj + j;
			}
		}
		//  Form inverse(R) * (inverse(R))'.
		jj = 0;
		for ( j = 1; j <= n; j++ ){
			j1 = jj + 1;
			jj = jj + j;
			k1 = 1;
			kj = j1;

			for ( k = 1; k <= j-1; k++ ){
				t = ap[kj-1];
				daxpy ( k, t, ap+j1-1, 1, ap+k1-1, 1 );
				k1 = k1 + k;
				kj = kj + 1;
			}
			t = ap[jj-1];
			dscal ( j, t, ap+j1-1, 1 );
		}
	}
	return;
}
//****************************************************************************80
int dppfa ( double ap[], int n )
//****************************************************************************80
//
//  Purpose:
//
//    DPPFA factors a real symmetric positive definite matrix in packed form.
//
//  Discussion:
//
//    DPPFA is usually called by DPPCO, but it can be called
//    directly with a saving in time if RCOND is not needed.
//
//  Packed storage:
//
//    The following program segment will pack the upper
//    triangle of a symmetric matrix.
//
//      k = 0
//      do j = 1, n
//        do i = 1, j
//          k = k + 1
//          ap(k) = a(i,j)
//        end do
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double AP[N*(N+1)/2].  On input, the packed
//    form of a symmetric matrix A.  The columns of the upper triangle are
//    stored sequentially in a one-dimensional array.  On output, an upper
//    triangular matrix R, stored in packed form, so that A = R'*R.
//
//    Input, int N, the order of the matrix.
//
//    Output, int DPPFA, error flag.
//    0, for normal return.
//    K, if the leading minor of order K is not positive definite.
//
{
	int info = 0, j, jj = 0, k, kj, kk;
	double s, t;

	for ( j = 1; j <= n; j++ ){
		s = ZERO;
		kj = jj;
		kk = 0;

		for ( k = 1; k <= j-1; k++ ){
			kj = kj + 1;
			t = ap[kj-1] - ddot ( k-1, ap+kk, 1, ap+jj, 1 );
			kk = kk + k;
			t = t / ap[kk-1];
			ap[kj-1] = t;
			s = s + t * t;
		}

		jj = jj + j;
		s = ap[jj-1] - s;

		if ( s <= ZERO ){	info = j;		return info;	}

		ap[jj-1] = sqrt ( s );
	}
	return info;
}
//****************************************************************************80
void dppsl ( double ap[], int n, double b[] )
//****************************************************************************80
//
//  Purpose:
//
//    DPPSL solves a real symmetric positive definite system factored by DPPCO or DPPFA.
//
//  Discussion:
//
//    To compute inverse(A) * C where C is a matrix with P columns
//
//      call dppco ( ap, n, rcond, z, info )
//
//      if ( rcond is too small .or. info /= 0 ) then
//        exit
//      end if
//
//      do j = 1, p
//        call dppsl ( ap, n, c(1,j) )
//      end do
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal.  Technically this indicates
//    singularity but it is usually caused by improper subroutine
//    arguments.  It will not occur if the subroutines are called
//    correctly and INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double AP[N*(N+1)/2], the output from DPPCO or DPPFA.
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
{
	int k, kk = 0;
	double t;

	for ( k = 1; k <= n; k++ ){
		t = ddot ( k-1, ap+kk, 1, b, 1 );
		kk = kk + k;
		b[k-1] = ( b[k-1] - t ) / ap[kk-1];
	}

	for ( k = n; 1 <= k; k-- ){
		b[k-1] = b[k-1] / ap[kk-1];
		kk = kk - k;
		t = -b[k-1];
		daxpy ( k-1, t, ap+kk, 1, b, 1 );
	}
	return;
}
//****************************************************************************80
void dptsl ( int n, double d[], double e[], double b[] )
//****************************************************************************80
//
//  Purpose:
//
//    DPTSL solves a positive definite tridiagonal linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D[N], on input the diagonal of the
//    tridiagonal matrix.  On output, D is destroyed.
//
//    Input, double E[N], the offdiagonal of the tridiagonal matrix in
//    entries E(1:N-1).
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
{
	int k, kbm1, ke, kf, kp1, nm1d2;
	double t1, t2;
	//  Check for 1 x 1 case.
	if ( n == 1 ){
		b[0] = b[0] / d[0];
		return;
	}

	nm1d2 = ( n - 1 ) / 2;

	if ( 2 < n ){
		kbm1 = n - 1;
		//  Zero top half of subdiagonal and bottom half of superdiagonal.
		for ( k = 1; k <= nm1d2; k++ ){
			t1 = e[k-1] / d[k-1];
			d[k] = d[k] - t1 * e[k-1];
			b[k] = b[k] - t1 * b[k-1];
			t2 = e[kbm1-1] / d[kbm1];
			d[kbm1-1] = d[kbm1-1] - t2 * e[kbm1-1];
			b[kbm1-1] = b[kbm1-1] - t2 * b[kbm1];
			kbm1 = kbm1 - 1;
		}
	}

	kp1 = nm1d2 + 1;
	//  Clean up for possible 2 x 2 block at center.
	if ( ( n % 2 ) == 0 ){
		t1 = e[kp1-1] / d[kp1-1];
		d[kp1] = d[kp1] - t1 * e[kp1-1];
		b[kp1] = b[kp1] - t1 * b[kp1-1];
		kp1 = kp1 + 1;
	}
	//  Back solve starting at the center, going towards the top and bottom.
	b[kp1-1] = b[kp1-1] / d[kp1-1];

	if ( 2 < n ){
		k = kp1 - 1;
		ke = kp1 + nm1d2 - 1;

		for ( kf = kp1; kf <= ke; kf++ ){
			b[k-1] = ( b[k-1] - e[k-1] * b[k] ) / d[k-1];
			b[kf] = ( b[kf] - e[kf-1] * b[kf-1] ) / d[kf];
			k = k - 1;
		}
	}

	if ( ( n % 2 ) == 0 ){	b[0] = ( b[0] - e[0] * b[1] ) / d[0];	}

	return;
}
//****************************************************************************80
void dqrdc ( double a[], int lda, int n, int p, double qraux[], int jpvt[],
		double work[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DQRDC computes the QR factorization of a real rectangular matrix.
//
//  Discussion:
//
//    DQRDC uses Householder transformations.
//
//    Column pivoting based on the 2-norms of the reduced columns may be
//    performed at the user's option.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A(LDA,P).  On input, the N by P matrix
//    whose decomposition is to be computed.  On output, A contains in
//    its upper triangle the upper triangular matrix R of the QR
//    factorization.  Below its diagonal A contains information from
//    which the orthogonal part of the decomposition can be recovered.
//    Note that if pivoting has been requested, the decomposition is not that
//    of the original matrix A but that of A with its columns permuted
//    as described by JPVT.
//
//    Input, int LDA, the leading dimension of the array A.  LDA must
//    be at least N.
//
//    Input, int N, the number of rows of the matrix A.
//
//    Input, int P, the number of columns of the matrix A.
//
//    Output, double QRAUX[P], contains further information required
//    to recover the orthogonal part of the decomposition.
//
//    Input/output, integer JPVT[P].  On input, JPVT contains integers that
//    control the selection of the pivot columns.  The K-th column A(*,K) of A
//    is placed in one of three classes according to the value of JPVT(K).
//      > 0, then A(K) is an initial column.
//      = 0, then A(K) is a free column.
//      < 0, then A(K) is a final column.
//    Before the decomposition is computed, initial columns are moved to
//    the beginning of the array A and final columns to the end.  Both
//    initial and final columns are frozen in place during the computation
//    and only free columns are moved.  At the K-th stage of the
//    reduction, if A(*,K) is occupied by a free column it is interchanged
//    with the free column of largest reduced norm.  JPVT is not referenced
//    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
//    original matrix that has been interchanged into the K-th column, if
//    pivoting was requested.
//
//    Workspace, double WORK[P].  WORK is not referenced if JOB == 0.
//
//    Input, int JOB, initiates column pivoting.
//    0, no pivoting is done.
//    nonzero, pivoting is done.
//
{
	int j, jp, l, lup, maxj, pl = 1, pu = 0;
	bool swapj;
	double maxnrm, nrmxl, t, tt;
	//  If pivoting is requested, rearrange the columns.
	if ( job != 0 ){
		for ( j = 1; j <= p; j++ ){
			swapj = ( 0 < jpvt[j-1] );

			if ( jpvt[j-1] < 0 ){	jpvt[j-1] = -j;		}
			else{			jpvt[j-1] = j;		}

			if ( swapj ){
				if ( j != pl )	{	dswap ( n, a+0+(pl-1)*lda, 1, a+0+(j-1), 1 );	}
				jpvt[j-1] = jpvt[pl-1];
				jpvt[pl-1] = j;
				pl = pl + 1;
			}
		}
		pu = p;

		for ( j = p; 1 <= j; j-- ){
			if ( jpvt[j-1] < 0 ){
				jpvt[j-1] = -jpvt[j-1];

				if ( j != pu ){
					dswap ( n, a+0+(pu-1)*lda, 1, a+0+(j-1)*lda, 1 );
					jp = jpvt[pu-1];
					jpvt[pu-1] = jpvt[j-1];
					jpvt[j-1] = jp;
				}
				pu = pu - 1;
			}
		}
	}
	//  Compute the norms of the free columns.
	for ( j = pl; j <= pu; j++ ){	qraux[j-1] = dnrm2 ( n, a+0+(j-1)*lda, 1 );	}

	for ( j = pl; j <= pu; j++ ){	work[j-1] = qraux[j-1];	}
	//  Perform the Householder reduction of A.
	lup = i4_min ( n, p );

	for ( l = 1; l <= lup; l++ ){
		//  Bring the column of largest norm into the pivot position.
		if ( pl <= l && l < pu ){
			maxnrm = ZERO;
			maxj = l;
			for ( j = l; j <= pu; j++ ){
				if ( maxnrm < qraux[j-1] ){
					maxnrm = qraux[j-1];
					maxj = j;
				}
			}

			if ( maxj != l ){
				dswap ( n, a+0+(l-1)*lda, 1, a+0+(maxj-1)*lda, 1 );
				qraux[maxj-1] = qraux[l-1];
				work[maxj-1] = work[l-1];
				jp = jpvt[maxj-1];
				jpvt[maxj-1] = jpvt[l-1];
				jpvt[l-1] = jp;
			}
		}
		//  Compute the Householder transformation for column L.
		qraux[l-1] = ZERO;

		if ( l != n ){
			nrmxl = dnrm2 ( n-l+1, a+l-1+(l-1)*lda, 1 );

			if ( nrmxl != ZERO ){
				if ( a[l-1+(l-1)*lda] != ZERO ){nrmxl = nrmxl * r8_sign ( a[l-1+(l-1)*lda] );	}

				dscal ( n-l+1, ONE / nrmxl, a+l-1+(l-1)*lda, 1 );
				a[l-1+(l-1)*lda] = ONE + a[l-1+(l-1)*lda];
				//  Apply the transformation to the remaining columns, updating the norms.
				for ( j = l + 1; j <= p; j++ ){
					t = -ddot ( n-l+1, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 )
            												/ a[l-1+(l-1)*lda];
					daxpy ( n-l+1, t, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 );

					if ( pl <= j && j <= pu ){
						if ( qraux[j-1] != ZERO ){
							tt = ONE - pow ( r8_abs ( a[l-1+(j-1)*lda] ) / qraux[j-1], 2 );
							tt = r8_max ( tt, ZERO );
							t = tt;
							tt = ONE + static_cast<Real>(0.05) * tt * pow ( qraux[j-1] / work[j-1], 2 );

							if ( tt != ONE ){	qraux[j-1] = qraux[j-1] * sqrt ( t );	}
							else{
								qraux[j-1] = dnrm2 ( n-l, a+l+(j-1)*lda, 1 );
								work[j-1] = qraux[j-1];
							}
						}
					}
				}
				//  Save the transformation.
				qraux[l-1] = a[l-1+(l-1)*lda];
				a[l-1+(l-1)*lda] = -nrmxl;
			}
		}
	}
	return;
}
//****************************************************************************80
int dqrsl ( double a[], int lda, int n, int k, double qraux[], double y[],
		double qy[], double qty[], double b[], double rsd[], double ab[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DQRSL computes transformations, projections, and least squares solutions.
//
//  Discussion:
//
//    DQRSL requires the output of DQRDC.
//
//    For K <= min(N,P), let AK be the matrix
//
//      AK = ( A(JPVT[0]), A(JPVT(2)), ..., A(JPVT(K)) )
//
//    formed from columns JPVT[0], ..., JPVT(K) of the original
//    N by P matrix A that was input to DQRDC.  If no pivoting was
//    done, AK consists of the first K columns of A in their
//    original order.  DQRDC produces a factored orthogonal matrix Q
//    and an upper triangular matrix R such that
//
//      AK = Q * (R)
//               (0)
//
//    This information is contained in coded form in the arrays
//    A and QRAUX.
//
//    The parameters QY, QTY, B, RSD, and AB are not referenced
//    if their computation is not requested and in this case
//    can be replaced by dummy variables in the calling program.
//    To save storage, the user may in some cases use the same
//    array for different parameters in the calling sequence.  A
//    frequently occuring example is when one wishes to compute
//    any of B, RSD, or AB and does not need Y or QTY.  In this
//    case one may identify Y, QTY, and one of B, RSD, or AB, while
//    providing separate arrays for anything else that is to be
//    computed.
//
//    Thus the calling sequence
//
//      dqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
//
//    will result in the computation of B and RSD, with RSD
//    overwriting Y.  More generally, each item in the following
//    list contains groups of permissible identifications for
//    a single calling sequence.
//
//      1. (Y,QTY,B) (RSD) (AB) (QY)
//
//      2. (Y,QTY,RSD) (B) (AB) (QY)
//
//      3. (Y,QTY,AB) (B) (RSD) (QY)
//
//      4. (Y,QY) (QTY,B) (RSD) (AB)
//
//      5. (Y,QY) (QTY,RSD) (B) (AB)
//
//      6. (Y,QY) (QTY,AB) (B) (RSD)
//
//    In any group the value returned in the array allocated to
//    the group corresponds to the last member of the group.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double A[LDA*P], contains the output of DQRDC.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the number of rows of the matrix AK.  It must
//    have the same value as N in DQRDC.
//
//    Input, int K, the number of columns of the matrix AK.  K
//    must not be greater than min(N,P), where P is the same as in the
//    calling sequence to DQRDC.
//
//    Input, double QRAUX[P], the auxiliary output from DQRDC.
//
//    Input, double Y[N], a vector to be manipulated by DQRSL.
//
//    Output, double QY[N], contains Q * Y, if requested.
//
//    Output, double QTY[N], contains Q' * Y, if requested.
//
//    Output, double B[K], the solution of the least squares problem
//      minimize norm2 ( Y - AK * B),
//    if its computation has been requested.  Note that if pivoting was
//    requested in DQRDC, the J-th component of B will be associated with
//    column JPVT(J) of the original matrix A that was input into DQRDC.
//
//    Output, double RSD[N], the least squares residual Y - AK * B,
//    if its computation has been requested.  RSD is also the orthogonal
//    projection of Y onto the orthogonal complement of the column space
//    of AK.
//
//    Output, double AB[N], the least squares approximation Ak * B,
//    if its computation has been requested.  AB is also the orthogonal
//    projection of Y onto the column space of A.
//
//    Input, integer JOB, specifies what is to be computed.  JOB has
//    the decimal expansion ABCDE, with the following meaning:
//
//      if A != 0, compute QY.
//      if B != 0, compute QTY.
//      if C != 0, compute QTY and B.
//      if D != 0, compute QTY and RSD.
//      if E != 0, compute QTY and AB.
//
//    Note that a request to compute B, RSD, or AB automatically triggers
//    the computation of QTY, for which an array must be provided in the
//    calling sequence.
//
//    Output, int DQRSL, is zero unless the computation of B has
//    been requested and R is exactly singular.  In this case, INFO is the
//    index of the first zero diagonal element of R, and B is left unaltered.
//
{
	bool cab, cb, cqty, cqy, cr;
	int i, info = 0, j, jj, ju;
	double t, temp;
	//  Determine what is to be computed.
	cqy =  (   job / 10000          != 0 );
	cqty = ( ( job %  10000 )       != 0 );
	cb =   ( ( job %   1000 ) / 100 != 0 );
	cr =   ( ( job %    100 ) /  10 != 0 );
	cab =  ( ( job %     10 )       != 0 );

	ju = i4_min ( k, n-1 );
	//  Special action when N = 1.
	if ( ju == 0 ){
		if ( cqy )	{		qy[0] = y[0];		}
		if ( cqty )	{		qty[0] = y[0];		}
		if ( cab )	{		ab[0] = y[0];		}
		if ( cb ){
			if ( a[0+0*lda] == ZERO ){	info = 1;	}
			else{		b[0] = y[0] / a[0+0*lda];	}
		}

		if ( cr )	{		rsd[0] = ZERO;		}
		return info;
	}
	//  Set up to compute QY or QTY.
	if ( cqy )	{	for ( i = 1; i <= n; i++ ){		qy[i-1] = y[i-1];	}	}
	if ( cqty )	{	for ( i = 1; i <= n; i++ ){		qty[i-1] = y[i-1];	}	}
	//  Compute QY.
	if ( cqy ){
		for ( jj = 1; jj <= ju; jj++ ){
			j = ju - jj + 1;

			if ( qraux[j-1] != ZERO ){
				temp = a[j-1+(j-1)*lda];
				a[j-1+(j-1)*lda] = qraux[j-1];
				t = -ddot ( n-j+1, a+j-1+(j-1)*lda, 1, qy+j-1, 1 ) / a[j-1+(j-1)*lda];
				daxpy ( n-j+1, t, a+j-1+(j-1)*lda, 1, qy+j-1, 1 );
				a[j-1+(j-1)*lda] = temp;
			}
		}
	}
	//  Compute Q'*Y.
	if ( cqty ){
		for ( j = 1; j <= ju; j++ ){
			if ( qraux[j-1] != ZERO ){
				temp = a[j-1+(j-1)*lda];
				a[j-1+(j-1)*lda] = qraux[j-1];
				t = -ddot ( n-j+1, a+j-1+(j-1)*lda, 1, qty+j-1, 1 ) / a[j-1+(j-1)*lda];
				daxpy ( n-j+1, t, a+j-1+(j-1)*lda, 1, qty+j-1, 1 );
				a[j-1+(j-1)*lda] = temp;
			}
		}
	}
	//  Set up to compute B, RSD, or AB.
	if ( cb )				{	for ( i = 1; i <= k; i++ )	{	b[i-1] = qty[i-1];	}	}
	if ( cab )				{	for ( i = 1; i <= k; i++ )	{	ab[i-1] = qty[i-1];	}	}
	if ( cr && k < n )		{	for ( i = k+1; i <= n; i++ ){	rsd[i-1] = qty[i-1];}	}
	if ( cab && k+1 <= n )	{	for ( i = k+1; i <= n; i++ ){	ab[i-1] = ZERO;		}	}
	if ( cr )				{	for ( i = 1; i <= k; i++ )	{	rsd[i-1] = ZERO;		}	}
	//  Compute B.
	if ( cb ){
		for ( jj = 1; jj <= k; jj++ ){
			j = k - jj + 1;

			if ( a[j-1+(j-1)*lda] == ZERO ){	info = j;	break;		}

			b[j-1] = b[j-1] / a[j-1+(j-1)*lda];

			if ( j != 1 ){
				t = -b[j-1];
				daxpy ( j-1, t, a+0+(j-1)*lda, 1, b, 1 );
			}
		}
	}
	//  Compute RSD or AB as required.
	if ( cr || cab ){
		for ( jj = 1; jj <= ju; jj++ ){
			j = ju - jj + 1;

			if ( qraux[j-1] != ZERO ){
				temp = a[j-1+(j-1)*lda];
				a[j-1+(j-1)*lda] = qraux[j-1];

				if ( cr ){
					t = -ddot ( n-j+1, a+j-1+(j-1)*lda, 1, rsd+j-1, 1 )
            												/ a[j-1+(j-1)*lda];
					daxpy ( n-j+1, t, a+j-1+(j-1)*lda, 1, rsd+j-1, 1 );
				}

				if ( cab ){
					t = -ddot ( n-j+1, a+j-1+(j-1)*lda, 1, ab+j-1, 1 )
            												/ a[j-1+(j-1)*lda];
					daxpy ( n-j+1, t, a+j-1+(j-1)*lda, 1, ab+j-1, 1 );
				}
				a[j-1+(j-1)*lda] = temp;
			}
		}
	}
	return info;
}
//****************************************************************************80
double dsico ( double a[], int lda, int n, int kpvt[], double z[] )
//****************************************************************************80
//
//  Purpose:
//
//    DSICO factors a real symmetric matrix and estimates its condition.
//
//  Discussion:
//
//    If RCOND is not needed, DSIFA is slightly faster.
//
//    To solve A * X = B, follow DSICO by DSISL.
//
//    To compute inverse(A)*C, follow DSICO by DSISL.
//
//    To compute inverse(A), follow DSICO by DSIDI.
//
//    To compute determinant(A), follow DSICO by DSIDI.
//
//    To compute inertia(A), follow DSICO by DSIDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the symmetric
//    matrix to be factored.  Only the diagonal and upper triangle are used.
//    On output, a block diagonal matrix and the multipliers which
//    were used to obtain it.  The factorization can be written A = U*D*U'
//    where U is a product of permutation and unit upper triangular
//    matrices, U' is the transpose of U, and D is block diagonal
//    with 1 by 1 and 2 by 2 blocks.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int KPVT[N], pivot indices.
//
//    Output, double Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
//    Output, double DSICO, an estimate of the reciprocal condition number RCOND
//    of A.  For the system A*X = B, relative perturbations in A and B of size
//    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
//    If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0D+00
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
{
	double rcond, s, t, ynorm, ak, akm1, anorm, bk, bkm1, denom, ek;
	int i, info, j, k, kp, kps, ks;
	//  Find the norm of A, using only entries in the upper half of the matrix.
	for ( j = 1; j <= n; j++ ){
		z[j-1] = dasum ( j, a+0+(j-1)*lda, 1 );
		for ( i = 1; i <= j-1; i++ ){	z[i-1] = z[i-1] + r8_abs ( a[i-1+(j-1)*lda] );	}
	}

	anorm = ZERO;
	for ( i = 1; i <= n; i++ ){	anorm = r8_max ( anorm, z[i-1] );	}
	//  Factor.
	info = dsifa ( a, lda, n, kpvt );
	/*!
	 * RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
	 *
	 * Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
	 *
	 * The components of E are chosen to cause maximum local
	 * growth in the elements of W where U*D*W = E.
	 *
	 * The vectors are frequently rescaled to avoid overflow.
	 *
	 * Solve U * D * W = E.
	 */
	ek = ONE;
	for ( i = 1; i <= n; i++ ){		z[i-1] = ZERO;	}
	k = n;

	while ( k != 0 ){
		if ( kpvt[k-1] < 0 ){	ks = 2;		}
		else{			ks = 1;		}
		kp = abs ( kpvt[k-1] );
		kps = k + 1 - ks;

		if ( kp != kps ){
			t = z[kps-1];
			z[kps-1] = z[kp-1];
			z[kp-1] = t;
		}

		if ( z[k-1] != ZERO ){	ek = ek * r8_sign ( z[k-1] );	}

		z[k-1] = z[k-1] + ek;
		daxpy ( k-ks, z[k-2], a+0+(k-1)*lda, 1, z, 1 );

		if ( ks != 1 ){
			if ( z[k-2] != ZERO ){
				ek = ek * r8_sign ( z[k-2] );
			}
			z[k-2] = z[k-2] + ek;
			daxpy ( k-ks, z[k-2], a+0+(k-2)*lda, 1, z, 1 );
		}

		if ( ks != 2 ){
			if ( r8_abs ( a[k-1+(k-1)*lda] ) < r8_abs ( z[k-1] ) ){
				s = r8_abs ( a[k-1+(k-1)*lda] ) / r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
				ek = s * ek;
			}

			if ( a[k-1+(k-1)*lda] != ZERO ){		z[k-1] = z[k-1] / a[k-1+(k-1)*lda];		}
			else{		z[k-1] = ONE;		}
		}else{
			ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
			akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
			bk = z[k-1] / a[k-2+(k-1)*lda];
			bkm1 = z[k-2] / a[k-2+(k-1)*lda];
			denom = ak * akm1 - ONE;
			z[k-1] = ( akm1 * bk - bkm1 ) / denom;
			z[k-2] = ( ak * bkm1 - bk ) / denom;
		}
		k = k - ks;
	}
	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	//  Solve U' * Y = W.
	k = 1;

	while ( k <= n ){
		if ( kpvt[k-1] < 0 ){	ks = 2;		}
		else{		ks = 1;		}

		if ( k != 1 ){
			z[k-1] = z[k-1] + ddot ( k-1, a+0+(k-1)*lda, 1, z, 1 );

			if ( ks == 2 ){		z[k] = z[k] + ddot ( k-1, a+0+k*lda, 1, z, 1 );		}

			kp = abs ( kpvt[k-1] );

			if ( kp != k ){
				t = z[k-1];
				z[k-1] = z[kp-1];
				z[kp-1] = t;
			}
		}
		k = k + ks;
	}
	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = z[i-1] / s;	}
	ynorm = ONE;
	//  Solve U * D * V = Y.
	k = n;

	while ( k != 0 ){
		if ( kpvt[k-1] < 0 ){	ks = 2;		}
		else{		ks = 1;		}

		if ( k != ks ){
			kp = abs ( kpvt[k-1] );
			kps = k + 1 - ks;

			if ( kp != kps ){
				t = z[kps-1];
				z[kps-1] = z[kp-1];
				z[kp-1] = t;
			}

			daxpy ( k-ks, z[k-1], a+0+(k-1)*lda, 1, z, 1 );

			if ( ks == 2 ){	daxpy ( k-ks, z[k-2], a+0+(k-2)*lda, 1, z, 1 );	}
		}

		if ( ks != 2 ){
			if ( r8_abs ( a[k-1+(k-1)*lda] ) < r8_abs ( z[k-1] ) ){
				s = r8_abs ( a[k-1+(k-1)*lda] ) / r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
				ynorm = s * ynorm;
			}

			if ( a[k-1+(k-1)*lda] != ZERO ){		z[k-1] = z[k-1] / a[k-1+(k-1)*lda];		}
			else{	z[k-1] = ONE;	}
		}else{
			ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
			akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
			bk = z[k-1] / a[k-2+(k-1)*lda];
			bkm1 = z[k-2] / a[k-2+(k-1)*lda];
			denom = ak * akm1 - ONE;
			z[k-1] = ( akm1 * bk - bkm1 ) / denom;
			z[k-2] = ( ak * bkm1 - bk ) / denom;
		}
		k = k - ks;
	}

	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;
	//  Solve U' * Z = V.
	k = 1;

	while ( k <= n ){
		if ( kpvt[k-1] < 0 ){	ks = 2;		}
		else{	ks = 1;		}

		if ( k != 1 ){
			z[k-1] = z[k-1] + ddot ( k-1, a+0+(k-1)*lda, 1, z, 1 );
			if ( ks == 2 ){	z[k] = z[k] + ddot ( k-1, a+0+k*lda, 1, z, 1 );		}
			kp = abs ( kpvt[k-1] );

			if ( kp != k ){
				t = z[k-1];
				z[k-1] = z[kp-1];
				z[kp-1] = t;
			}
		}
		k = k + ks;
	}
	//  Make ZNORM = 1.0.
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;

	if ( anorm != ZERO ){	rcond = ynorm / anorm;	}
	else{	rcond = ZERO;	}

	return rcond;
}
//****************************************************************************80
void dsidi ( double a[], int lda, int n, int kpvt[], double det[2],
		int inert[3], double work[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DSIDI computes the determinant, inertia and inverse of a real symmetric matrix.
//
//  Discussion:
//
//    DSIDI uses the factors from DSIFA.
//
//    A division by zero may occur if the inverse is requested
//    and DSICO has set RCOND == 0.0D+00 or DSIFA has set INFO /= 0.
//
//    Variables not requested by JOB are not used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A(LDA,N).  On input, the output from DSIFA.
//    On output, the upper triangle of the inverse of the original matrix,
//    if requested.  The strict lower triangle is never referenced.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int KPVT[N], the pivot vector from DSIFA.
//
//    Output, double DET[2], the determinant of the original matrix,
//    if requested.
//      determinant = DET[0] * 10.0**DET[1]
//    with 1.0D+00 <= abs ( DET[0] ) < 10.0D+00 or DET[0] = 0.0.
//
//    Output, int INERT(3), the inertia of the original matrix,
//    if requested.
//    INERT(1) = number of positive eigenvalues.
//    INERT(2) = number of negative eigenvalues.
//    INERT(3) = number of zero eigenvalues.
//
//    Workspace, double WORK[N].
//
//    Input, int JOB, specifies the tasks.
//    JOB has the decimal expansion ABC where
//    If C /= 0, the inverse is computed,
//    If B /= 0, the determinant is computed,
//    If A /= 0, the inertia is computed.
//    For example, JOB = 111 gives all three.
//
{
	double ak, akkp1, akp1, d, t, temp;
	bool dodet, doert, doinv;
	int j, jb, k, ks, kstep;

	doinv = ( job %   10 )       != 0;
	dodet = ( job %  100 ) /  10 != 0;
	doert = ( job % 1000 ) / 100 != 0;

	if ( dodet || doert ){
		if ( doert ){
			inert[0] = 0;		inert[1] = 0;		inert[2] = 0;
		}

		if ( dodet ){
			det[0] = ONE;		det[1] = ZERO;
		}

		t = ZERO;

		for ( k = 1; k <= n; k++ ){
			d = a[k-1+(k-1)*lda];
			/*!
			 * 2 by 2 block.
			 * use det (d  s)  =  (d/t * c - t) * t,  t = abs ( s )
			 *         (s  c)
			 * to avoid underflow/overflow troubles.
			 *
			 * Take two passes through scaling.  Use T for flag.
			 */
			if ( kpvt[k-1] <= 0 ){
				if ( t == ZERO ){
					t = r8_abs ( a[k-1+k*lda] );
					d = ( d / t ) * a[k+k*lda] - t;
				}else{
					d = t;		t = ZERO;
				}
			}

			if ( doert ){
				if ( ZERO < d )			{	inert[0] = inert[0] + 1;	}
				else if ( d < ZERO )	{	inert[1] = inert[1] + 1;	}
				else if ( d == ZERO )	{	inert[2] = inert[2] + 1;	}
			}

			if ( dodet ){
				det[0] = det[0] * d;

				if ( det[0] != ZERO ){
					while ( r8_abs ( det[0] ) < ONE ){
						det[0] = det[0] * static_cast<Real>(10.0);
						det[1] = det[1] - ONE;
					}

					while ( static_cast<Real>(10.0) <= r8_abs ( det[0] ) ){
						det[0] = det[0] / static_cast<Real>(10.0);
						det[1] = det[1] + ONE;
					}
				}
			}
		}
	}
	//  Compute inverse(A).
	if ( doinv ){
		k = 1;

		while ( k <= n ){
			if ( 0 <= kpvt[k-1] ){
				//  1 by 1.
				a[k-1+(k-1)*lda] = ONE / a[k-1+(k-1)*lda];

				if ( 2 <= k ){
					dcopy ( k-1, a+0+(k-1)*lda, 1, work, 1 );

					for ( j = 1; j <= k-1; j++ ){
						a[j-1+(k-1)*lda] = ddot ( j, a+0+(j-1)*lda, 1, work, 1 );
						daxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
					}
					a[k-1+(k-1)*lda] = a[k-1+(k-1)*lda]
										 + ddot ( k-1, work, 1, a+0+(k-1)*lda, 1 );
				}
				kstep = 1;
			}else{	//  2 by 2.
				t = r8_abs ( a[k-1+k*lda] );
				ak = a[k-1+(k-1)*lda] / t;
				akp1 = a[k+k*lda] / t;
				akkp1 = a[k-1+k*lda] / t;
				d = t * ( ak * akp1 - ONE );
				a[k-1+(k-1)*lda] = akp1 / d;
				a[k+k*lda] = ak / d;
				a[k-1+k*lda] = -akkp1 / d;

				if ( 2 <= k ){
					dcopy ( k-1, a+0+k*lda, 1, work, 1 );

					for ( j = 1; j <= k-1; j++ ){
						a[j-1+k*lda] = ddot ( j, a+0+(j-1)*lda, 1, work, 1 );
						daxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+k*lda, 1 );
					}
					a[k+k*lda] = a[k+k*lda] + ddot ( k-1, work, 1, a+0+k*lda, 1 );
					a[k-1+k*lda] = a[k-1+k*lda]
									 + ddot ( k-1, a+0+(k-1)*lda, 1, a+0+k*lda, 1 );
					dcopy ( k-1, a+0+(k-1)*lda, 1, work, 1 );

					for ( j = 1; j <= k-1; j++ ){
						a[j-1+(k-1)*lda] = ddot ( j, a+0+(j-1)*lda, 1, work, 1 );
						daxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
					}
					a[k-1+(k-1)*lda] = a[k-1+(k-1)*lda]
										 + ddot ( k-1, work, 1, a+0+(k-1)*lda, 1 );
				}
				kstep = 2;
			}
			//  Swap.
			ks = abs ( kpvt[k-1] );

			if ( ks != k ){
				dswap ( ks, a+0+(ks-1)*lda, 1, a+0+(k-1)*lda, 1 );

				for ( jb = ks; jb <= k; jb++ ){
					j = k + ks - jb;
					temp = a[j-1+(k-1)*lda];
					a[j-1+(k-1)*lda] = a[ks-1+(j-1)*lda];
					a[ks-1+(j-1)*lda] = temp;
				}

				if ( kstep != 1 ){
					temp = a[ks-1+k*lda];
					a[ks-1+k*lda] = a[k-1+k*lda];
					a[k-1+k*lda] = temp;
				}
			}
			k = k + kstep;
		}
	}
	return;
}
//****************************************************************************80
int dsifa ( double a[], int lda, int n, int kpvt[] )
//****************************************************************************80
//
//  Purpose:
//
//    DSIFA factors a real symmetric matrix.
//
//  Discussion:
//
//    To solve A*X = B, follow DSIFA by DSISL.
//
//    To compute inverse(A)*C, follow DSIFA by DSISL.
//
//    To compute determinant(A), follow DSIFA by DSIDI.
//
//    To compute inertia(A), follow DSIFA by DSIDI.
//
//    To compute inverse(A), follow DSIFA by DSIDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the symmetric matrix
//    to be factored.  Only the diagonal and upper triangle are used.
//    On output, a block diagonal matrix and the multipliers which
//    were used to obtain it.  The factorization can be written A = U*D*U'
//    where U is a product of permutation and unit upper triangular
//    matrices, U' is the transpose of U, and D is block diagonal
//    with 1 by 1 and 2 by 2 blocks.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int KPVT[N], the pivot indices.
//
//    Output, integer DSIFA, error flag.
//    0, normal value.
//    K, if the K-th pivot block is singular.  This is not an error
//    condition for this subroutine, but it does indicate that DSISL
//    or DSIDI may divide by zero if called.
//
{
	double absakk, ak, akm1, alpha, bk, bkm1, colmax, denom, mulkm1, mulk, rowmax, t;
	int imax, imaxp1, info = 0, j, jj, jmax, k, km1, kstep;
	bool swap;
	//  ALPHA is used in choosing pivot block size.
	alpha = ( ONE + sqrt ( 17.0 ) ) / 8.0;
	//  Main loop on K, which goes from N to 1.
	k = n;

	while ( 0 < k ){
		if ( k == 1 ){
			kpvt[0] = 1;
			if ( a[0+0*lda] == ZERO ){	info = 1;	}
			return info;
		}
		/*
		 * This section of code determines the kind of
		 * elimination to be performed.  When it is completed,
		 * KSTEP will be set to the size of the pivot block, and
		 * SWAP will be set to .true. if an interchange is required.
		 */
		km1 = k - 1;
		absakk = r8_abs ( a[k-1+(k-1)*lda] );
		//  Determine the largest off-diagonal element in column K.
		imax = idamax ( k-1, a+0+(k-1)*lda, 1 );
		colmax = r8_abs ( a[imax-1+(k-1)*lda] );

		if ( alpha * colmax <= absakk ){
			kstep = 1;
			swap = false;
		}else{	//  Determine the largest off-diagonal element in row IMAX.
			rowmax = ZERO;
			imaxp1 = imax + 1;
			for ( j = imaxp1; j <= k; j++ ){rowmax = r8_max ( rowmax, r8_abs ( a[imax-1+(j-1)*lda] ) );	}

			if ( imax != 1 ){
				jmax = idamax ( imax-1, a+0+(imax-1)*lda, 1 );
				rowmax = r8_max ( rowmax, r8_abs ( a[jmax-1+(imax-1)*lda] ) );
			}

			if ( alpha * rowmax <= r8_abs ( a[imax-1+(imax-1)*lda] ) ){
				kstep = 1;
				swap = true;
			}else if ( alpha * colmax * ( colmax / rowmax ) <= absakk ){
				kstep = 1;
				swap = false;
			}else{
				kstep = 2;
				swap = ( imax != k-1 );
			}
		}
		/*!
		 * Column K is zero.
		 * Set INFO and iterate the loop.
		 */
		if ( r8_max ( absakk, colmax ) == ZERO ){
			kpvt[k-1] = k;
			info = k;
		}
		/*!
		 * 1 x 1 pivot block.
		 *
		 * Perform an interchange.
		 */
		else if ( kstep != 2 ){
			if ( swap ){
				dswap ( imax, a+0+(imax-1)*lda, 1, a+0+(k-1)*lda, 1 );

				for ( jj = imax; jj <= k; jj++ ){
					j = k + imax - jj;
					t = a[j-1+(k-1)*lda];
					a[j-1+(k-1)*lda] = a[imax-1+(j-1)*lda];
					a[imax-1+(j-1)*lda] = t;
				}
			}
			//  Perform the elimination.
			for ( jj = 1; jj <= k-1; jj++ ){
				j = k - jj;
				mulk = -a[j-1+(k-1)*lda] / a[k-1+(k-1)*lda];
				t = mulk;
				daxpy ( j, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
				a[j-1+(k-1)*lda] = mulk;
			}
			//  Set the pivot array.
			if ( swap ){	kpvt[k-1] = imax;	}
			else{			kpvt[k-1] = k;		}
		}
		/*!
		 * 2 x 2 pivot block.
		 * Perform an interchange.
		 */
		else{
			if ( swap ){
				dswap ( imax, a+0+(imax-1)*lda, 1, a+0+(k-2)*lda, 1 );

				for ( jj = imax; jj <= k-1; jj++ ){
					j = k-1 + imax - jj;
					t = a[j-1+(k-1)*lda];
					a[j-1+(k-1)*lda] = a[imax-1+(j-1)*lda];
					a[imax-1+(j-1)*lda] = t;
				}

				t = a[k-2+(k-1)*lda];
				a[k-2+(k-1)*lda] = a[imax-1+(k-1)*lda];
				a[imax-1+(k-1)*lda] = t;
			}
			//  Perform the elimination.
			if ( k-2 != 0 ){
				ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
				akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
				denom = ONE - ak * akm1;

				for ( jj = 1; jj <= k-2; jj++ ){
					j = k-1 - jj;
					bk = a[j-1+(k-1)*lda] / a[k-2+(k-1)*lda];
					bkm1 = a[j-1+(k-2)*lda] / a[k-2+(k-1)*lda];
					mulk = ( akm1 * bk - bkm1 ) / denom;
					mulkm1 = ( ak * bkm1 - bk ) / denom;
					t = mulk;
					daxpy ( j, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
					t = mulkm1;
					daxpy ( j, t, a+0+(k-2)*lda, 1, a+0+(j-1)*lda, 1 );
					a[j-1+(k-1)*lda] = mulk;
					a[j-1+(k-2)*lda] = mulkm1;
				}
			}
			//  Set the pivot array.
			if ( swap ){	kpvt[k-1] = -imax;	}
			else	{		kpvt[k-1] = 1 - k;	}
			kpvt[k-2] = kpvt[k-1];
		}
		k = k - kstep;
	}
	return info;
}
//****************************************************************************80
void dsisl ( double a[], int lda, int n, int kpvt[], double b[] )
//****************************************************************************80
//
//  Purpose:
//
//    DSISL solves a real symmetric system factored by DSIFA.
//
//  Discussion:
//
//    To compute inverse(A) * C where C is a matrix with P columns
//
//      call dsifa ( a, lda, n, kpvt, info )
//
//      if ( info == 0 ) then
//        do j = 1, p
//          call dsisl ( a, lda, n, kpvt, c(1,j) )
//        end do
//      end if
//
//    A division by zero may occur if the inverse is requested
//    and DSICO has set RCOND == 0.0D+00 or DSIFA has set INFO /= 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double A[LDA*N], the output from DSIFA.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int KPVT[N], the pivot vector from DSIFA.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
{
	double ak, akm1, bk, bkm1, denom, temp;
	int k = n, kp;
	//  Loop backward applying the transformations and D inverse to B.
	while ( 0 < k ){
		if ( 0 <= kpvt[k-1] ){
			//  1 x 1 pivot block.
			if ( k != 1 ){
				kp = kpvt[k-1];
				//  Interchange.
				if ( kp != k ){
					temp = b[k-1];
					b[k-1] = b[kp-1];
					b[kp-1] = temp;
				}
				//  Apply the transformation.
				daxpy ( k-1, b[k-1], a+0+(k-1)*lda, 1, b, 1 );
			}
			//  Apply D inverse.
			b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
			k = k - 1;
		}else{
			//  2 x 2 pivot block.
			if ( k != 2 ){
				kp = abs ( kpvt[k-1] );
				//  Interchange.
				if ( kp != k-1 ){
					temp = b[k-2];
					b[k-2] = b[kp-1];
					b[kp-1] = temp;
				}
				//  Apply the transformation.
				daxpy ( k-2, b[k-1], a+0+(k-1)*lda, 1, b, 1 );
				daxpy ( k-2, b[k-2], a+0+(k-2)*lda, 1, b, 1 );
			}
			//  Apply D inverse.
			ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
			akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
			bk = b[k-1] / a[k-2+(k-1)*lda];
			bkm1 = b[k-2] / a[k-2+(k-1)*lda];
			denom = ak * akm1 - ONE;
			b[k-1] = ( akm1 * bk - bkm1 ) / denom;
			b[k-2] = ( ak * bkm1 - bk ) / denom;
			k = k - 2;
		}
	}
	//  Loop forward applying the transformations.
	k = 1;

	while ( k <= n ){
		if ( 0 <= kpvt[k-1] ){
			//  1 x 1 pivot block.
			if ( k != 1 ){
				//  Apply the transformation.
				b[k-1] = b[k-1] + ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
				kp = kpvt[k-1];
				//  Interchange.
				if ( kp != k ){
					temp = b[k-1];
					b[k-1] = b[kp-1];
					b[kp-1] = temp;
				}
			}
			k = k + 1;
		}else{
			//  2 x 2 pivot block.
			if ( k != 1 ){
				//  Apply the transformation.
				b[k-1] = b[k-1] + ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
				b[k] = b[k] + ddot ( k-1, a+0+k*lda, 1, b, 1 );
				kp = abs ( kpvt[k-1] );
				//  Interchange.
				if ( kp != k ){
					temp = b[k-1];
					b[k-1] = b[kp-1];
					b[kp-1] = temp;
				}
			}
			k = k + 2;
		}
	}
	return;
}
//****************************************************************************80
double dspco ( double ap[], int n, int kpvt[], double z[] )
//****************************************************************************80
//
//  Purpose:
//
//    DSPCO factors a real symmetric matrix stored in packed form.
//
//  Discussion:
//
//    DSPCO uses elimination with symmetric pivoting and estimates
//    the condition of the matrix.
//
//    If RCOND is not needed, DSPFA is slightly faster.
//
//    To solve A*X = B, follow DSPCO by DSPSL.
//
//    To compute inverse(A)*C, follow DSPCO by DSPSL.
//
//    To compute inverse(A), follow DSPCO by DSPDI.
//
//    To compute determinant(A), follow DSPCO by DSPDI.
//
//    To compute inertia(A), follow DSPCO by DSPDI.
//
//  Packed storage:
//
//    The following program segment will pack the upper triangle of a
//    symmetric matrix.
//
//      k = 0
//      do j = 1, n
//        do i = 1, j
//          k = k + 1
//          ap[k-1] = a(i,j)
//        }
//      }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double AP[N*(N+1)/2].  On input, the packed form
//    of a symmetric matrix A.  The columns of the upper triangle are stored
//    sequentially in a one-dimensional array.  On output, a block diagonal
//    matrix and the multipliers which were used to obtain it, stored in
//    packed form.  The factorization can be written A = U*D*U'
//    where U is a product of permutation and unit upper triangular
//    matrices, U' is the transpose of U, and D is block diagonal
//    with 1 by 1 and 2 by 2 blocks.
//
//    Input, int N, the order of the matrix.
//
//    Output, int KPVT[N], the pivot indices.
//
//    Output, double Z[N] a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
//    Output, double DSPCO, an estimate of the reciprocal condition number RCOND
//    of A.  For the system A*X = B, relative perturbations in A and B of size
//    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
//    If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0D+00
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
{
	double ak, akm1, anorm, bk, bkm1, denom, ek, rcond, s, t, ynorm;
	int i, ij, ik, ikm1, ikp1, info, j, j1, k, kk, km1k, km1km1, kp, kps, ks;
	//  Find norm of A using only upper half.
	j1 = 1;
	for ( j = 1; j <= n; j++ ){
		z[j-1] = dasum ( j, ap+j1-1, 1 );
		ij = j1;
		j1 = j1 + j;
		for ( i = 1; i <= j-1; i++ ){
			z[i-1] = z[i-1] + r8_abs ( ap[ij-1] );
			ij = ij + 1;
		}
	}
	anorm = ZERO;
	for ( i = 1; i <= n; i++ ){	anorm = r8_max ( anorm, z[i-1] );	}
	//  Factor.
	info = dspfa ( ap, n, kpvt );
	/*!
	 * RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
	 *
	 * Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
	 *
	 * The components of E are chosen to cause maximum local
	 * growth in the elements of W where U*D*W = E.
	 *
	 * The vectors are frequently rescaled to avoid overflow.
	 *
	 * Solve U * D * W = E.
	 */
	ek = ONE;
	for ( i = 1; i <= n; i++ ){		z[i-1] = ZERO;	}

	k = n;
	ik = ( n * ( n - 1 ) ) / 2;

	while ( k != 0 ){
		kk = ik + k;
		ikm1 = ik - ( k - 1 );

		if ( kpvt[k-1] < 0 ){	ks = 2;	}
		else{	ks = 1;	}

		kp = abs ( kpvt[k-1] );
		kps = k + 1 - ks;

		if ( kp != kps ){
			t = z[kps-1];
			z[kps-1] = z[kp-1];
			z[kp-1] = t;
		}

		if ( z[k-1] != ZERO ){	ek = ek * r8_sign ( z[k-1] );	}

		z[k-1] = z[k-1] + ek;
		daxpy ( k-ks, z[k-1], ap+ik, 1, z, 1 );

		if ( ks != 1 ){
			if ( z[k-2] != ZERO ){	ek = ek * r8_sign ( z[k-2] );	}
			z[k-2] = z[k-2] + ek;
			daxpy ( k-ks, z[k-2], ap+ikm1, 1, z, 1 );
		}

		if ( ks != 2 ){
			if ( r8_abs ( ap[kk-1] ) < r8_abs ( z[k-1] ) ){
				s = r8_abs ( ap[kk-1] ) / r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ ){	z[i-1] = s * z[i-1];	}
				ek = s * ek;
			}

			if ( ap[kk-1] != ZERO ){	z[k-1] = z[k-1] / ap[kk-1];	}
			else{	z[k-1] = ONE;	}
		}else{
			km1k = ik + k - 1;
			km1km1 = ikm1 + k - 1;
			ak = ap[kk-1] / ap[km1k-1];
			akm1 = ap[km1km1-1] / ap[km1k-1];
			bk = z[k-1] / ap[km1k-1];
			bkm1 = z[k-2] / ap[km1k-1];
			denom = ak * akm1 - ONE;
			z[k-1] = ( akm1 * bk - bkm1 ) / denom;
			z[k-2] = ( ak * bkm1 - bk ) / denom;
		}
		k = k - ks;
		ik = ik - k;
		if ( ks == 2 ){	ik = ik - ( k + 1 );	}
	}

	s = dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){	z[i-1] = z[i-1] / s;	}
	//  Solve U' * Y = W.
	k = 1;
	ik = 0;

	while ( k <= n ){
		if ( kpvt[k-1] < 0 ){	ks = 2;	}	else{		ks = 1;		}

		if ( k != 1 ){
			z[k-1] = z[k-1] + ddot ( k-1, ap+ik, 1, z, 1 );
			ikp1 = ik + k;

			if ( ks == 2 ){	z[k] = z[k] + ddot ( k-1, ap+ikp1, 1, z, 1 );		}

			kp = abs ( kpvt[k-1] );

			if ( kp != k ){
				t = z[k-1];
				z[k-1] = z[kp-1];
				z[kp-1] = t;
			}
		}

		ik = ik + k;
		if ( ks == 2 ){		ik = ik + ( k + 1 );	}
		k = k + ks;
	}
	for ( i = 1; i <= n; i++ ){	z[i-1] = s * z[i-1];}
	ynorm = ONE;
	//  Solve U * D * V = Y.
	k = n;

	ik = ( n * ( n - 1 ) ) / 2;

	while ( 0 < k ){
		kk = ik + k;
		ikm1 = ik - ( k - 1 );

		if ( kpvt[k-1] < 0 ){	ks = 2;	}	else{	ks = 1;	}

		if ( k != ks ){
			kp = abs ( kpvt[k-1] );
			kps = k + 1 - ks;

			if ( kp != kps ){
				t = z[kps-1];
				z[kps-1] = z[kp-1];
				z[kp-1] = t;
			}

			daxpy ( k-ks, z[k-1], ap+ik, 1, z, 1 );

			if ( ks == 2 ){		daxpy ( k-ks, z[k-2], ap+ikm1, 1, z, 1 );	}
		}

		if ( ks != 2 ){
			if ( r8_abs ( ap[kk-1] ) < r8_abs ( z[k-1] ) ){
				s = r8_abs ( ap[kk-1] ) / r8_abs ( z[k-1] );
				for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
				ynorm = s * ynorm;
			}

			if ( ap[kk-1] != ZERO ){	z[k-1] = z[k-1] / ap[kk-1];	}
			else{	z[k-1] = ONE;	}
		}else{
			km1k = ik + k - 1;
			km1km1 = ikm1 + k - 1;
			ak = ap[kk-1] / ap[km1k-1];
			akm1 = ap[km1km1-1] / ap[km1k-1];
			bk = z[k-1] / ap[km1k-1];
			bkm1 = z[k-2] / ap[km1k-1];
			denom = ak * akm1 - ONE;
			z[k-1] = ( akm1 * bk - bkm1 ) / denom;
			z[k-2] = ( ak * bkm1 - bk ) / denom;
		}
		k = k - ks;
		ik = ik - k;
		if ( ks == 2 ){		ik = ik - ( k + 1 );	}
	}
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){		z[i-1] = s * z[i-1];	}
	ynorm = s * ynorm;
	//  Solve U' * Z = V.
	k = 1;
	ik = 0;

	while ( k <= n ){
		if ( kpvt[k-1] < 0 ){	ks = 2;	}	else{	ks = 1;	}

		if ( k != 1 ){
			z[k-1] = z[k-1] + ddot ( k-1, ap+ik, 1, z, 1 );
			ikp1 = ik + k;

			if ( ks == 2 ){	z[k] = z[k] + ddot ( k-1, ap+ikp1, 1, z, 1 );	}

			kp = abs ( kpvt[k-1] );

			if ( kp != k ){
				t = z[k-1];
				z[k-1] = z[kp-1];
				z[kp-1] = t;
			}
		}

		ik = ik + k;
		if ( ks == 2 ){	ik = ik + ( k + 1 );	}
		k = k + ks;
	}
	//  Make ZNORM = 1.0.
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){	z[i-1] = s * z[i-1];}
	ynorm = s * ynorm;

	if ( anorm != ZERO ){	rcond = ynorm / anorm;}	else{	rcond = ZERO;}

	return rcond;
}
//****************************************************************************80
void dspdi ( double ap[], int n, int kpvt[], double det[2], int inert[3],
		double work[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DSPDI computes the determinant, inertia and inverse of a real symmetric matrix.
//
//  Discussion:
//
//    DSPDI uses the factors from DSPFA, where the matrix is stored in
//    packed form.
//
//    A division by zero will occur if the inverse is requested
//    and DSPCO has set RCOND == 0.0D+00 or DSPFA has set INFO /= 0.
//
//    Variables not requested by JOB are not used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double AP[(N*(N+1))/2].  On input, the output from
//    DSPFA.  On output, the upper triangle of the inverse of the original
//    matrix, stored in packed form, if requested.  The columns of the upper
//    triangle are stored sequentially in a one-dimensional array.
//
//    Input, int N, the order of the matrix.
//
//    Input, int KPVT[N], the pivot vector from DSPFA.
//
//    Output, double DET[2], the determinant of the original matrix,
//    if requested.
//      determinant = DET[0] * 10.0**DET[1]
//    with 1.0D+00 <= abs ( DET[0] ) < 10.0D+00 or DET[0] = 0.0.
//
//    Output, int INERT[3], the inertia of the original matrix, if requested.
//    INERT(1) = number of positive eigenvalues.
//    INERT(2) = number of negative eigenvalues.
//    INERT(3) = number of zero eigenvalues.
//
//    Workspace, double WORK[N].
//
//    Input, int JOB, has the decimal expansion ABC where:
//      if A /= 0, the inertia is computed,
//      if B /= 0, the determinant is computed,
//      if C /= 0, the inverse is computed.
//    For example, JOB = 111  gives all three.
//
{
	double ak, akkp1, akp1, d, t, temp;
	bool dodet, doert, doinv;
	int ij, ik, ikp1, iks, j,jb, jk, jkp1, k, kk, kkp1, km1, ks, ksj, kskp1, kstep;

	doinv = ( job %   10 )       != 0;
	dodet = ( job %  100 ) /  10 != 0;
	doert = ( job % 1000 ) / 100 != 0;

	if ( dodet || doert ){
		if ( doert ){
			inert[0] = 0;		inert[1] = 0;		inert[2] = 0;
		}

		if ( dodet ){
			det[0] = ONE;		det[1] = ZERO;
		}

		t = ZERO;
		ik = 0;

		for ( k = 1; k <= n; k++ ){
			kk = ik + k;
			d = ap[kk-1];
			/*!
			 * 2 by 2 block
			 * use det (d  s)  =  (d/t * c - t) * t,  t = abs ( s )
			 *         (s  c)
			 * to avoid underflow/overflow troubles.
			 *
			 * Take two passes through scaling.  Use T for flag.
			 */
			if ( kpvt[k-1] <= 0 ){
				if ( t == ZERO ){
					ikp1 = ik + k;
					kkp1 = ikp1 + k;
					t = r8_abs ( ap[kkp1-1] );
					d = ( d / t ) * ap[kkp1] - t;
				}else{
					d = t;
					t = ZERO;
				}
			}

			if ( doert ){
				if ( ZERO < d )			{	inert[0] = inert[0] + 1;	}
				else if ( d < ZERO )	{	inert[1] = inert[1] + 1;	}
				else if ( d == ZERO )	{	inert[2] = inert[2] + 1;	}
			}

			if ( dodet ){
				det[0] = det[0] * d;

				if ( det[0] != ZERO ){
					while ( r8_abs ( det[0] ) < ONE ){
						det[0] = det[0] * static_cast<Real>(10.0);
						det[1] = det[1] - ONE;
					}

					while ( static_cast<Real>(10.0) <= r8_abs ( det[0] ) ){
						det[0] = det[0] / static_cast<Real>(10.0);
						det[1] = det[1] + ONE;
					}
				}
			}
			ik = ik + k;
		}
	}
	//  Compute inverse(A).
	if ( doinv ){
		k = 1;		ik = 0;

		while ( k <= n ){
			km1 = k - 1;
			kk = ik + k;
			ikp1 = ik + k;
			kkp1 = ikp1 + k;

			if ( 0 <= kpvt[k-1] ){
				//  1 by 1.
				ap[kk-1] = ONE / ap[kk-1];

				if ( 2 <= k ){
					dcopy ( k-1, ap+ik, 1, work, 1 );
					ij = 0;

					for ( j = 1; j <= k-1; j++ ){
						jk = ik + j;
						ap[jk-1] = ddot ( j, ap+ij, 1, work, 1 );
						daxpy ( j-1, work[j-1], ap+ij, 1, ap+ik, 1 );
						ij = ij + j;
					}
					ap[kk-1] = ap[kk-1] + ddot ( k-1, work, 1, ap+ik, 1 );
				}
				kstep = 1;
			}else{
				//  2 by 2.
				t = r8_abs ( ap[kkp1-1] );
				ak = ap[kk-1] / t;
				akp1 = ap[kkp1] / t;
				akkp1 = ap[kkp1-1] / t;
				d = t * ( ak * akp1 - ONE );
				ap[kk-1] = akp1 / d;
				ap[kkp1] = ak / d;
				ap[kkp1-1] = -akkp1 / d;

				if ( 1 <= km1 ){
					dcopy ( km1, ap+ikp1, 1, work, 1 );
					ij = 0;

					for ( j = 1; j <= km1; j++ ){
						jkp1 = ikp1 + j;
						ap[jkp1-1] = ddot ( j, ap+ij, 1, work, 1 );
						daxpy ( j-1, work[j-1], ap+ij, 1, ap+ikp1, 1 );
						ij = ij + j;
					}

					ap[kkp1] = ap[kkp1] + ddot ( km1, work, 1, ap+ikp1, 1 );
					ap[kkp1-1] = ap[kkp1-1] + ddot ( km1, ap+ik, 1, ap+ikp1, 1 );
					dcopy ( km1, ap+ik, 1, work, 1 );
					ij = 0;

					for ( j = 1; j <= km1; j++ ){
						jk = ik + j;
						ap[jk-1] = ddot ( j, ap+ij, 1, work, 1 );
						daxpy ( j-1, work[j-1], ap+ij, 1, ap+ik, 1 );
						ij = ij + j;
					}
					ap[kk-1] = ap[kk-1] + ddot ( km1, work, 1, ap+ik, 1 );
				}
				kstep = 2;
			}
			//  Swap.
			ks = abs ( kpvt[k-1] );

			if ( ks != k ){
				iks = ( ks * ( ks - 1 ) ) / 2;
				dswap ( ks, ap+iks, 1, ap+ik, 1 );
				ksj = ik + ks;

				for ( jb = ks; jb <= k; jb++ ){
					j = k + ks - jb;
					jk = ik + j;
					temp = ap[jk-1];
					ap[jk-1] = ap[ksj-1];
					ap[ksj-1] = temp;
					ksj = ksj - ( j - 1 );
				}

				if ( kstep != 1 ){
					kskp1 = ikp1 + ks;
					temp = ap[kskp1-1];
					ap[kskp1-1] = ap[kkp1-1];
					ap[kkp1-1] = temp;
				}
			}
			ik = ik + k;
			if ( kstep == 2 ){	ik = ik + k + 1;	}
			k = k + kstep;
		}
	}
	return;
}
//****************************************************************************80
int dspfa ( double ap[], int n, int kpvt[] )
//****************************************************************************80
//
//  Purpose:
//
//    DSPFA factors a real symmetric matrix stored in packed form.
//
//  Discussion:
//
//    To solve A*X = B, follow DSPFA by DSPSL.
//
//    To compute inverse(A)*C, follow DSPFA by DSPSL.
//
//    To compute determinant(A), follow DSPFA by DSPDI.
//
//    To compute inertia(A), follow DSPFA by DSPDI.
//
//    To compute inverse(A), follow DSPFA by DSPDI.
//
//  Packed storage:
//
//    The following program segment will pack the upper triangle of a
//    symmetric matrix.
//
//      k = 0
//      do j = 1, n
//        do i = 1, j
//          k = k + 1
//          ap(k) = a(i,j)
//        end do
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double AP[(N*(N+1))/2].  On input, the packed form of a
//    symmetric matrix A.  The columns of the upper triangle are stored
//    sequentially in a one-dimensional array.  On output, a block diagonal
//    matrix and the multipliers which were used to obtain it stored in
//    packed form.  The factorization can be written A = U*D*U' where U
//    is a product of permutation and unit upper triangular matrices, U'
//    is the transpose of U, and D is block diagonal with 1 by 1 and 2
//    by 2 blocks.
//
//    Input, int N, the order of the matrix.
//
//    Output, int KPVT[N], the pivot indices.
//
//    Output, int DSPFA, error flag.
//    0, normal value.
//    K, if the K-th pivot block is singular.  This is not an error
//    condition for this subroutine, but it does indicate that DSPSL or
//    DSPDI may divide by zero if called.
//
{
	double absakk, ak, akm1, alpha, bk, bkm1, colmax, denom, mulk, mulkm1, rowmax, t;
	int ij, ijj, ik, ikm1, im, imax, imaxp1, imim, imj, imk, info = 0, j, jj;
	int jk, jkm1, jmax, jmim, k, kk, km1, km1k, km1km1, kstep;
	bool swap;
	//  ALPHA is used in choosing pivot block size.
	alpha = ( ONE + sqrt ( 17.0 ) ) / 8.0;

	//  Main loop on K, which goes from N to 1.
	k = n;
	ik = ( n * ( n - 1 ) ) / 2;

	for ( ; ; ){
		//  Leave the loop if K = 0 or K = 1.
		if ( k == 0 ){	break;	}

		if ( k == 1 ){
			kpvt[0] = 1;
			if ( ap[0] == ZERO ){	info = 1;	}
			break;
		}
		/*!
		 * This section of code determines the kind of elimination to be performed.
		 * When it is completed, KSTEP will be set to the size of the pivot block,
		 * and SWAP will be set to .true. if an interchange is required.
		 */
		km1 = k - 1;
		kk = ik + k;
		absakk = r8_abs ( ap[kk-1] );
		//  Determine the largest off-diagonal element in column K.
		imax = idamax ( k-1, ap+ik, 1 );
		imk = ik + imax;
		colmax = r8_abs ( ap[imk-1] );

		if ( alpha * colmax <= absakk ){
			kstep = 1;
			swap = false;
		}else{	//  Determine the largest off-diagonal element in row IMAX.
			rowmax = ZERO;
			imaxp1 = imax + 1;
			im = ( imax * ( imax - 1 ) ) / 2;
			imj = im + 2 * imax;

			for ( j = imaxp1; j <= k; j++ ){
				rowmax = r8_max ( rowmax, r8_abs ( ap[imj-1] ) );
				imj = imj + j;
			}

			if ( imax != 1 ){
				jmax = idamax ( imax-1, ap+im, 1 );
				jmim = jmax + im;
				rowmax = r8_max ( rowmax, r8_abs ( ap[jmim-1] ) );
			}

			imim = imax + im;

			if ( alpha * rowmax <= r8_abs ( ap[imim-1] ) ){
				kstep = 1;
				swap = true;
			}
			else if ( alpha * colmax * ( colmax / rowmax ) <= absakk ){
				kstep = 1;
				swap = false;
			}else{
				kstep = 2;
				swap = imax != km1;
			}
		}
		//  Column K is zero.  Set INFO and iterate the loop.
		if ( r8_max ( absakk, colmax ) == ZERO ){
			kpvt[k-1] = k;
			info = k;
		}else{
			if ( kstep != 2 ){
				//  1 x 1 pivot block.
				if ( swap ){
					//  Perform an interchange.
					dswap ( imax, ap+im, 1, ap+ik, 1 );
					imj = ik + imax;

					for ( jj = imax; jj <= k; jj++ ){
						j = k + imax - jj;
						jk = ik + j;
						t = ap[jk-1];
						ap[jk-1] = ap[imj-1];
						ap[imj-1] = t;
						imj = imj - ( j - 1 );
					}
				}
				//  Perform the elimination.
				ij = ik - ( k - 1 );

				for ( jj = 1; jj <= km1; jj++ ){
					j = k - jj;
					jk = ik + j;
					mulk = -ap[jk-1] / ap[kk-1];
					t = mulk;
					daxpy ( j, t, ap+ik, 1, ap+ij, 1 );
					ijj = ij + j;
					ap[jk-1] = mulk;
					ij = ij - ( j - 1 );
				}
				//  Set the pivot array.
				if ( swap ){	kpvt[k-1] = imax;	}	else{	kpvt[k-1] = k;	}
			}else{
				//  2 x 2 pivot block.
				km1k = ik + k - 1;
				ikm1 = ik - ( k - 1 );
				//  Perform an interchange.
				if ( swap ){
					dswap ( imax, ap+im, 1, ap+ikm1, 1 );
					imj = ikm1 + imax;

					for ( jj = imax; jj <= km1; jj++ ){
						j = km1 + imax - jj;
						jkm1 = ikm1 + j;
						t = ap[jkm1-1];
						ap[jkm1-1] = ap[imj-1];
						ap[imj-1] = t;
						imj = imj - ( j - 1 );
					}
					t = ap[km1k-1];
					ap[km1k-1] = ap[imk-1];
					ap[imk-1] = t;
				}
				//  Perform the elimination.
				if ( k-2 != 0 ){
					ak = ap[kk-1] / ap[km1k-1];
					km1km1 = ikm1 + k - 1;
					akm1 = ap[km1km1-1] / ap[km1k-1];
					denom = ONE - ak * akm1;
					ij = ik - ( k - 1 ) - ( k - 2 );

					for ( jj = 1; jj <= k-2; jj++ ){
						j = km1 - jj;
						jk = ik + j;
						bk = ap[jk-1] / ap[km1k-1];
						jkm1 = ikm1 + j;
						bkm1 = ap[jkm1-1] / ap[km1k-1];
						mulk = ( akm1 * bk - bkm1 ) / denom;
						mulkm1 = ( ak * bkm1 - bk ) / denom;
						t = mulk;
						daxpy ( j, t, ap+ik, 1, ap+ij, 1 );
						t = mulkm1;
						daxpy ( j, t, ap+ikm1, 1, ap+ij, 1 );
						ap[jk-1] = mulk;
						ap[jkm1-1] = mulkm1;
						ijj = ij + j;
						ij = ij - ( j - 1 );
					}
				}
				//  Set the pivot array.
				if ( swap ){kpvt[k-1] = -imax;}	else{	kpvt[k-1] = 1 - k;}
				kpvt[k-2] = kpvt[k-1];
			}
		}

		ik = ik - ( k - 1 );
		if ( kstep == 2 ){	ik = ik - ( k - 2 );}

		k = k - kstep;
	}
	return info;
}
//****************************************************************************80
void dspsl ( double ap[], int n, int kpvt[], double b[] )
//****************************************************************************80
//
//  Purpose:
//
//    DSPSL solves the real symmetric system factored by DSPFA.
//
//  Discussion:
//
//    To compute inverse(A) * C where C is a matrix with P columns:
//
//      call dspfa ( ap, n, kpvt, info )
//
//      if ( info /= 0 ) go to ...
//
//      do j = 1, p
//        call dspsl ( ap, n, kpvt, c(1,j) )
//      end do
//
//    A division by zero may occur if DSPCO has set RCOND == 0.0D+00
//    or DSPFA has set INFO /= 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double AP[(N*(N+1))/2], the output from DSPFA.
//
//    Input, int N, the order of the matrix.
//
//    Input, int KPVT[N], the pivot vector from DSPFA.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
{
	double ak, akm1, bk, bkm1, denom, temp;
	int ik, ikm1, ikp1, k, kk, km1k, km1km1, kp;
	//  Loop backward applying the transformations and D inverse to B.
	k = n;
	ik = ( n * ( n - 1 ) ) / 2;

	while ( 0 < k ){
		kk = ik + k;

		if ( 0 <= kpvt[k-1] ){
			//  1 x 1 pivot block.
			if ( k != 1 ){
				kp = kpvt[k-1];
				//  Interchange.
				if ( kp != k ){
					temp = b[k-1];
					b[k-1] = b[kp-1];
					b[kp-1] = temp;
				}
				//  Apply the transformation.
				daxpy ( k-1, b[k-1], ap+ik, 1, b, 1 );
			}
			//  Apply D inverse.
			b[k-1] = b[k-1] / ap[kk-1];
			k = k - 1;
			ik = ik - k;
		}else{
			//  2 x 2 pivot block.
			ikm1 = ik - ( k - 1 );

			if ( k != 2 ){
				kp = abs ( kpvt[k-1] );
				//  Interchange.
				if ( kp != k-1 ){
					temp = b[k-2];
					b[k-2] = b[kp-1];
					b[kp-1] = temp;
				}
				//  Apply the transformation.
				daxpy ( k-2, b[k-1], ap+ik, 1, b, 1 );
				daxpy ( k-2, b[k-2], ap+ikm1, 1, b, 1 );
			}
			//  Apply D inverse.
			km1k = ik + k - 1;
			kk = ik + k;
			ak = ap[kk-1] / ap[km1k-1];
			km1km1 = ikm1 + k - 1;
			akm1 = ap[km1km1-1] / ap[km1k-1];
			bk = b[k-1] / ap[km1k-1];
			bkm1 = b[k-2] / ap[km1k-1];
			denom = ak * akm1 - ONE;
			b[k-1] = ( akm1 * bk - bkm1 ) / denom;
			b[k-2] = ( ak * bkm1 - bk ) / denom;
			k = k - 2;
			ik = ik - ( k + 1 ) - k;
		}
	}
	//  Loop forward applying the transformations.
	k = 1;
	ik = 0;

	while ( k <= n ){
		if ( 0 <= kpvt[k-1] ){
			//  1 x 1 pivot block.
			if ( k != 1 ){
				//  Apply the transformation.
				b[k-1] = b[k-1] + ddot ( k-1, ap+ik, 1, b, 1 );
				kp = kpvt[k-1];
				//  Interchange.
				if ( kp != k ){
					temp = b[k-1];
					b[k-1] = b[kp-1];
					b[kp-1] = temp;
				}
			}
			ik = ik + k;
			k = k + 1;
		}else{
			//  2 x 2 pivot block.
			if ( k != 1 ){
				//  Apply the transformation.
				b[k-1] = b[k-1] + ddot ( k-1, ap+ik, 1, b, 1 );
				ikp1 = ik + k;
				b[k] = b[k] + ddot ( k-1, ap+ikp1, 1, b, 1 );
				kp = abs ( kpvt[k-1] );
				//  Interchange.
				if ( kp != k ){
					temp = b[k-1];
					b[k-1] = b[kp-1];
					b[kp-1] = temp;
				}
			}
			ik = ik + k + k + 1;
			k = k + 2;
		}
	}
	return;
}
//****************************************************************************80
int dsvdc ( double a[], int lda, int m, int n, double s[], double e[],
		double u[], int ldu, double v[], int ldv, double work[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DSVDC computes the singular value decomposition of a real rectangular matrix.
//
//  Discussion:
//
//    This routine reduces an M by N matrix A to diagonal form by orthogonal
//    transformations U and V.  The diagonal elements S(I) are the singular
//    values of A.  The columns of U are the corresponding left singular
//    vectors, and the columns of V the right singular vectors.
//
//    The form of the singular value decomposition is then
//
//      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the M by N matrix whose
//    singular value decomposition is to be computed.  On output, the matrix
//    has been destroyed.  Depending on the user's requests, the matrix may
//    contain other useful information.
//
//    Input, int LDA, the leading dimension of the array A.
//    LDA must be at least M.
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix A.
//
//    Output, double S[MM], where MM = min(M+1,N).  The first
//    min(M,N) entries of S contain the singular values of A arranged in
//    descending order of magnitude.
//
//    Output, double E[MM], where MM = min(M+1,N), ordinarily contains zeros.
//    However see the discussion of INFO for exceptions.
//
//    Output, double U[LDU*K].  If JOBA = 1 then K = M;
//    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of left singular
//    vectors.  U is not referenced if JOBA = 0.  If M <= N or if JOBA = 2, then
//    U may be identified with A in the subroutine call.
//
//    Input, int LDU, the leading dimension of the array U.
//    LDU must be at least M.
//
//    Output, double V[LDV*N], the N by N matrix of right singular vectors.
//    V is not referenced if JOB is 0.  If N <= M, then V may be identified
//    with A in the subroutine call.
//
//    Input, int LDV, the leading dimension of the array V.
//    LDV must be at least N.
//
//    Workspace, double WORK[M].
//
//    Input, int JOB, controls the computation of the singular
//    vectors.  It has the decimal expansion AB with the following meaning:
//      A =  0, do not compute the left singular vectors.
//      A =  1, return the M left singular vectors in U.
//      A >= 2, return the first min(M,N) singular vectors in U.
//      B =  0, do not compute the right singular vectors.
//      B =  1, return the right singular vectors in V.
//
//    Output, int *DSVDC, status indicator INFO.
//    The singular values (and their corresponding singular vectors)
//    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
//    Thus if *INFO is 0, all the singular values and their vectors are
//    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
//    matrix with the elements of S on its diagonal and the elements of E on
//    its superdiagonal.  Thus the singular values of A and B are the same.
//
{
	double b, c, cs, el, emm1, f, g;
	double scale, shift, sl, sm, smm1, sn, t, t1, test, ztest;
	int i, info = 0, iter, j, jobu, k, kase, kk, l, ll, lls, ls, lu;
	int maxit = 30, mm, mm1, mn, mp1, nct, nctp1, ncu, nrt, nrtp1;
	bool wantu = false, wantv = false;
	//  Determine what is to be computed.
	jobu = ( job % 100 ) / 10;

	if ( 1 < jobu )	{	ncu = i4_min ( m, n );	}else{	ncu = m;	}
	if ( jobu != 0 ){	wantu = true;			}
	if ( ( job % 10 ) != 0 ){	wantv = true;	}
	/*!
	 * Reduce A to bidiagonal form, storing the diagonal elements
	 * in S and the super-diagonal elements in E.
	 */
	nct = i4_min ( m-1, n );
	nrt = i4_max ( 0, i4_min ( m, n-2 ) );
	lu = i4_max ( nct, nrt );

	for ( l = 1; l <= lu; l++ ){
		/*!
		 * Compute the transformation for the L-th column and
		 * place the L-th diagonal in S(L).
		 */
		if ( l <= nct ){
			s[l-1] = dnrm2 ( m-l+1, a+l-1+(l-1)*lda, 1 );

			if ( s[l-1] != ZERO ){
				if ( a[l-1+(l-1)*lda] != ZERO ){
					s[l-1] = r8_sign ( a[l-1+(l-1)*lda] ) * r8_abs ( s[l-1] );
				}
				dscal ( m-l+1, ONE / s[l-1], a+l-1+(l-1)*lda, 1 );
				a[l-1+(l-1)*lda] = ONE + a[l-1+(l-1)*lda];
			}
			s[l-1] = -s[l-1];
		}

		for ( j = l+1; j <= n; j++ ){
			//  Apply the transformation.
			if ( l <= nct && s[l-1] != ZERO ){
				t = - ddot ( m-l+1, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 )
        												  / a[l-1+(l-1)*lda];
				daxpy ( m-l+1, t, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 );
			}
			/*!
			 * Place the L-th row of A into E for the
			 * subsequent calculation of the row transformation.
			 */
			e[j-1] = a[l-1+(j-1)*lda];
		}
		//  Place the transformation in U for subsequent back multiplication.
		if ( wantu && l <= nct ){
			for ( i = l; i <= m; i++ ){	u[i-1+(l-1)*ldu] = a[i-1+(l-1)*lda]; }
		}

		if ( l <= nrt ){
			/*!
			 * Compute the L-th row transformation and place the
			 * L-th superdiagonal in E(L).
			 */
			e[l-1] = dnrm2 ( n-l, e+l, 1 );

			if ( e[l-1] != ZERO ){
				if ( e[l] != ZERO ){	e[l-1] = r8_sign ( e[l] ) * r8_abs ( e[l-1] );}
				dscal ( n-l, ONE / e[l-1], e+l, 1 );
				e[l] = ONE + e[l];
			}

			e[l-1] = -e[l-1];
			//  Apply the transformation.
			if ( l+1 <= m && e[l-1] != ZERO ){
				for ( j = l+1; j <= m; j++ ){	work[j-1] = ZERO;	}
				for ( j = l+1; j <= n; j++ ){	daxpy ( m-l, e[j-1], a+l+(j-1)*lda, 1, work+l, 1 );	}
				for ( j = l+1; j <= n; j++ ){	daxpy ( m-l, -e[j-1]/e[l], work+l, 1, a+l+(j-1)*lda, 1 );}
			}
			//  Place the transformation in V for subsequent back multiplication.
			if ( wantv ){	for ( j = l+1; j <= n; j++ ){	v[j-1+(l-1)*ldv] = e[j-1];	}	}
		}
	}
	//  Set up the final bidiagonal matrix of order MN.
	mn = i4_min ( m + 1, n );
	nctp1 = nct + 1;
	nrtp1 = nrt + 1;

	if ( nct < n )		{	s[nctp1-1] = a[nctp1-1+(nctp1-1)*lda];	}
	if ( m < mn )		{	s[mn-1] = ZERO;							}
	if ( nrtp1 < mn )	{	e[nrtp1-1] = a[nrtp1-1+(mn-1)*lda];		}

	e[mn-1] = ZERO;
	//  If required, generate U.
	if ( wantu ){
		for ( i = 1; i <= m; i++ ){
			for ( j = nctp1; j <= ncu; j++ ){	u[(i-1)+(j-1)*ldu] = ZERO;	}
		}

		for ( j = nctp1; j <= ncu; j++ )	{	u[j-1+(j-1)*ldu] = ONE;		}

		for ( ll = 1; ll <= nct; ll++ ){
			l = nct - ll + 1;

			if ( s[l-1] != ZERO ){
				for ( j = l+1; j <= ncu; j++ ){
					t = - ddot ( m-l+1, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 )
            												/ u[l-1+(l-1)*ldu];
					daxpy ( m-l+1, t, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 );
				}

				dscal ( m-l+1, -ONE, u+(l-1)+(l-1)*ldu, 1 );
				u[l-1+(l-1)*ldu] = ONE + u[l-1+(l-1)*ldu];
				for ( i = 1; i <= l-1; i++ ){	u[i-1+(l-1)*ldu] = ZERO;		}
			}else{
				for ( i = 1; i <= m; i++ )	{	u[i-1+(l-1)*ldu] = ZERO;		}
				u[l-1+(l-1)*ldu] = ONE;
			}
		}
	}
	//  If it is required, generate V.
	if ( wantv ){
		for ( ll = 1; ll <= n; ll++ ){
			l = n - ll + 1;

			if ( l <= nrt && e[l-1] != ZERO ){
				for ( j = l+1; j <= n; j++ ){
					t = - ddot ( n-l, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 )
            												/ v[l+(l-1)*ldv];
					daxpy ( n-l, t, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 );
				}

			}
			for ( i = 1; i <= n; i++ )	{	v[i-1+(l-1)*ldv] = ZERO;		}
			v[l-1+(l-1)*ldv] = ONE;
		}
	}
	//  Main iteration loop for the singular values.
	mm = mn;
	iter = 0;

	while ( 0 < mn ){
		//  If too many iterations have been performed, set flag and return.
		if ( maxit <= iter ){	info = mn;		return info;		}
		/*!
		 * This section of the program inspects for
		 * negligible elements in the S and E arrays.
		 *
		 * On completion the variables KASE and L are set as follows:
		 *
		 * KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
		 * KASE = 2     if S(L) is negligible and L < MN
		 * KASE = 3     if E(L-1) is negligible, L < MN, and
		 * 				   S(L), ..., S(MN) are not negligible (QR step).
		 * 				   KASE = 4     if E(MN-1) is negligible (convergence).
		 */
		for ( ll = 1; ll <= mn; ll++ ){
			l = mn - ll;

			if ( l == 0 ){		break;		}

			test = r8_abs ( s[l-1] ) + r8_abs ( s[l] );
			ztest = test + r8_abs ( e[l-1] );

			if ( ztest == test ){	e[l-1] = ZERO;		break;		}
		}

		if ( l == mn - 1 ){	kase = 4;	}
		else{
			mp1 = mn + 1;

			for ( lls = l+1; lls <= mn+1; lls++ ){
				ls = mn - lls + l + 1;

				if ( ls == l )		{	break;		}

				test = ZERO;
				if ( ls != mn )		{	test = test + r8_abs ( e[ls-1] );	}

				if ( ls != l + 1 )	{	test = test + r8_abs ( e[ls-2] );	}

				ztest = test + r8_abs ( s[ls-1] );

				if ( ztest == test ){	s[ls-1] = ZERO;		break;			}

			}

			if ( ls == l )			{	kase = 3;	}
			else if ( ls == mn )	{	kase = 1;	}
			else					{	kase = 2;		l = ls;				}
		}

		l = l + 1;
		//  Deflate negligible S(MN).
		if ( kase == 1 ){
			mm1 = mn - 1;
			f = e[mn-2];
			e[mn-2] = ZERO;

			for ( kk = 1; kk <= mm1; kk++ ){
				k = mm1 - kk + l;
				t1 = s[k-1];
				drotg ( &t1, &f, &cs, &sn );
				s[k-1] = t1;

				if ( k != l ){
					f = -sn * e[k-2];
					e[k-2] = cs * e[k-2];
				}

				if ( wantv ){drot ( n, v+0+(k-1)*ldv, 1, v+0+(mn-1)*ldv, 1, cs, sn );}
			}
		}else if ( kase == 2 ){		//  Split at negligible S(L).
			f = e[l-2];
			e[l-2] = ZERO;

			for ( k = l; k <= mn; k++ ){
				t1 = s[k-1];
				drotg ( &t1, &f, &cs, &sn );
				s[k-1] = t1;
				f = - sn * e[k-1];
				e[k-1] = cs * e[k-1];
				if ( wantu ){
					drot ( m, u+0+(k-1)*ldu, 1, u+0+(l-2)*ldu, 1, cs, sn );
				}
			}
		}else if ( kase == 3 ){	//  Perform one QR step.
			//  Calculate the shift.
			scale = r8_max ( r8_abs ( s[mn-1] ),
					r8_max ( r8_abs ( s[mn-2] ),
							r8_max ( r8_abs ( e[mn-2] ),
									r8_max ( r8_abs ( s[l-1] ), r8_abs ( e[l-1] ) ) ) ) );

			sm = s[mn-1] / scale;
			smm1 = s[mn-2] / scale;
			emm1 = e[mn-2] / scale;
			sl = s[l-1] / scale;
			el = e[l-1] / scale;
			b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0;
			c = ( sm * emm1 ) * ( sm * emm1 );
			shift = ZERO;

			if ( b != ZERO || c != ZERO ){
				shift = sqrt ( b * b + c );
				if ( b < ZERO ){
					shift = -shift;
				}
				shift = c / ( b + shift );
			}

			f = ( sl + sm ) * ( sl - sm ) - shift;
			g = sl * el;
			//  Chase zeros.
			mm1 = mn - 1;

			for ( k = l; k <= mm1; k++ ){
				drotg ( &f, &g, &cs, &sn );

				if ( k != l ){	e[k-2] = f;		}

				f = cs * s[k-1] + sn * e[k-1];
				e[k-1] = cs * e[k-1] - sn * s[k-1];
				g = sn * s[k];
				s[k] = cs * s[k];

				if ( wantv ){	drot ( n, v+0+(k-1)*ldv, 1, v+0+k*ldv, 1, cs, sn );	}

				drotg ( &f, &g, &cs, &sn );
				s[k-1] = f;
				f = cs * e[k-1] + sn * s[k];
				s[k] = -sn * e[k-1] + cs * s[k];
				g = sn * e[k];
				e[k] = cs * e[k];

				if ( wantu && k < m ){	drot ( m, u+0+(k-1)*ldu, 1, u+0+k*ldu, 1, cs, sn );	}
			}
			e[mn-2] = f;
			iter = iter + 1;
		}else if ( kase == 4 ){	//  Convergence.
			//  Make the singular value nonnegative.
			if ( s[l-1] < ZERO ){
				s[l-1] = -s[l-1];
				if ( wantv ){	dscal ( n, -ONE, v+0+(l-1)*ldv, 1 );	}
			}
			//  Order the singular value.
			for ( ; ; ){
				if ( l == mm )			{	break;	}
				if ( s[l] <= s[l-1] )	{	break;	}

				t = s[l-1];
				s[l-1] = s[l];
				s[l] = t;

				if ( wantv && l < n ){	dswap ( n, v+0+(l-1)*ldv, 1, v+0+l*ldv, 1 );}
				if ( wantu && l < m ){	dswap ( m, u+0+(l-1)*ldu, 1, u+0+l*ldu, 1 );}

				l = l + 1;
			}
			iter = 0;
			mn = mn - 1;
		}
	}
	return info;
}
//****************************************************************************80
double dtrco ( double t[], int ldt, int n, double z[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DTRCO estimates the condition of a real triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double T[LDT*N], the triangular matrix.  The zero
//    elements of the matrix are not referenced, and the corresponding
//    elements of the array can be used to store other information.
//
//    Input, int LDT, the leading dimension of the array T.
//
//    Input, int N, the order of the matrix.
//
//    Output, double Z[N] a work vector whose contents are usually
//    unimportant.  If T is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
//    Input, int JOB, indicates the shape of T:
//    0, T is lower triangular.
//    nonzero, T is upper triangular.
//
//    Output, double DTRCO, an estimate of the reciprocal condition RCOND
//    of T.  For the system T*X = B, relative perturbations in T and B of size
//    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
//    If RCOND is so small that the logical expression
//      1.0D+00 + RCOND == 1.0D+00
//    is true, then T may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
{
	int i, i1, j, j1, j2, k, kk, l;
	bool lower = ( job == 0 );
	double ek, rcond, s, sm, temp, tnorm = 0, w, wk, wkm, ynorm;
	//  Compute the 1-norm of T.
	for ( j = 1; j <= n; j++ ){
		if ( lower ){
			l = n + 1 - j;
			i1 = j;
		}else{
			l = j;		i1 = 1;
		}
		tnorm = r8_max ( tnorm, dasum ( l, t+i1-1+(j-1)*ldt, 1 ) );
	}
	/*!
	 * RCOND = 1/(norm(T)*(estimate of norm(inverse(T)))).
	 *
	 * Estimate = norm(Z)/norm(Y) where T * Z = Y and T' * Y = E.
	 *
	 * T' is the transpose of T.
	 *
	 * The components of E are chosen to cause maximum local
	 * growth in the elements of Y.
	 *
	 * The vectors are frequently rescaled to avoid overflow.
	 *
	 * Solve T' * Y = E.
	 */
	ek = ONE;
	for ( i = 1; i <= n; i++ ){	z[i-1] = ZERO;	}

	for ( kk = 1; kk <= n; kk++ ){
		if ( lower ){	k = n + 1 - kk;	}
		else		{	k = kk;			}

		if ( z[k-1] != ZERO ){ ek = r8_sign ( -z[k-1] ) * ek;}

		if ( r8_abs ( t[k-1+(k-1)*ldt] ) < r8_abs ( ek - z[k-1] ) ){
			s = r8_abs ( t[k-1+(k-1)*ldt] ) / r8_abs ( ek - z[k-1] );
			for ( i = 1; i <= n; i++ ){	z[i-1] = s * z[i-1];}
			ek = s * ek;
		}

		wk = ek - z[k-1];
		wkm = -ek - z[k-1];
		s = r8_abs ( wk );
		sm = r8_abs ( wkm );

		if ( t[k-1+(k-1)*ldt] != ZERO ){
			wk = wk / t[k-1+(k-1)*ldt];
			wkm = wkm / t[k-1+(k-1)*ldt];
		}else{
			wk = ONE;		wkm = ONE;
		}

		if ( kk != n ){
			if ( lower ){
				j1 = 1;			j2 = k - 1;
			}else{
				j1 = k + 1;		j2 = n;
			}

			for ( j = j1; j <= j2; j++ ){
				sm = sm + r8_abs ( z[j-1] + wkm * t[k-1+(j-1)*ldt] );
				z[j-1] = z[j-1] + wk * t[k-1+(j-1)*ldt];
				s = s + r8_abs ( z[j-1] );
			}

			if ( s < sm ){
				w = wkm - wk;
				wk = wkm;
				for ( j = j1; j <= j2; j++ ){z[j-1] = z[j-1] + w * t[k-1+(j-1)*ldt];}
			}
		}
		z[k-1] = wk;
	}

	temp = dasum ( n, z, 1 );

	for ( i = 1; i <= n; i++ ){	z[i-1] = z[i-1] / temp; }

	ynorm = ONE;
	//  Solve T * Z = Y.
	for ( kk = 1; kk <= n; kk++ ){
		if ( lower ){	k = kk;			}
		else		{	k = n + 1 - kk;	}

		if ( r8_abs ( t[k-1+(k-1)*ldt] ) < r8_abs ( z[k-1] ) ){
			s = r8_abs ( t[k-1+(k-1)*ldt] ) / r8_abs ( z[k-1] );
			for ( i = 1; i <= n; i++ ){	z[i-1] = s * z[i-1];}
			ynorm = s * ynorm;
		}

		if ( t[k-1+(k-1)*ldt] != ZERO ){	z[k-1] = z[k-1] / t[k-1+(k-1)*ldt];	}
		else{	z[k-1] = ONE;	}

		if ( lower ){	i1 = k + 1;	}
		else		{	i1 = 1;		}

		if ( kk < n ){
			w = -z[k-1];
			daxpy ( n-kk, w, t+i1-1+(k-1)*ldt, 1, z+i1-1, 1 );
		}
	}
	//  Make ZNORM = 1.0.
	s = ONE / dasum ( n, z, 1 );
	for ( i = 1; i <= n; i++ ){	z[i-1] = s * z[i-1];}
	ynorm = s * ynorm;

	if ( tnorm != ZERO )	{	rcond = ynorm / tnorm;	}
	else	{	rcond = ZERO;	}

	return rcond;
}
//****************************************************************************80
int dtrdi ( double t[], int ldt, int n, double det[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DTRDI computes the determinant and inverse of a real triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double T[LDT*N].
//    On input, T contains the triangular matrix.  The zero elements of the
//    matrix are not referenced, and the corresponding elements of the array
//    can be used to store other information.
//    On output, T contains the inverse matrix, if it was requested.
//
//    Input, int LDT, the leading dimension of T.
//
//    Input, int N, the order of the matrix.
//
//    Output, double DET[2], the determinant of the matrix, if
//    requested.  The determinant = DET[0] * 10.0**DET[1], with
//    1.0 <= abs ( DET[0] ) < 10.0, or DET[0] == 0.
//
//    Input, int JOB, specifies the shape of T, and the task.
//    010, inverse of lower triangular matrix.
//    011, inverse of upper triangular matrix.
//    100, determinant only.
//    110, determinant and inverse of lower triangular.
//    111, determinant and inverse of upper triangular.
//
//    Output, int DTRDI.
//    If the inverse was requested, then
//    0, if the system was nonsingular;
//    nonzero, if the system was singular.
//
{
	int i, info, j, k;
	double temp;
	//  Determinant.
	info = 0;

	if ( job / 100 != 0 ){
		det[0] = ONE;
		det[1] = ZERO;

		for ( i = 1; i <= n; i++ ){
			det[0] = det[0] * t[i-1+(i-1)*ldt];

			if ( det[0] == ZERO ){	break;	}

			while ( r8_abs ( det[0] ) < ONE ){
				det[0] = det[0] * static_cast<Real>(10.0);
				det[1] = det[1] - ONE;
			}

			while ( static_cast<Real>(10.0) <= r8_abs ( det[0] ) ){
				det[0] = det[0] / static_cast<Real>(10.0);
				det[1] = det[1] + ONE;
			}
		}
	}

	if ( ( ( job / 10 ) % 10 ) == 0 ){	return info;	}
	//  Inverse of an upper triangular matrix.
	if ( ( job % 10 ) != 0 ){
		info = 0;

		for ( k = 1; k <= n; k++ ){
			if ( t[k-1+(k-1)*ldt] == ZERO ){
				info = k;
				break;
			}

			t[k-1+(k-1)*ldt] = ONE / t[k-1+(k-1)*ldt];
			temp = -t[k-1+(k-1)*ldt];
			dscal ( k-1, temp, t+0+(k-1)*ldt, 1 );

			for ( j = k + 1; j <= n; j++ ){
				temp = t[k-1+(j-1)*ldt];
				t[k-1+(j-1)*ldt] = ZERO;
				daxpy ( k, temp, t+0+(k-1)*ldt, 1, t+0+(j-1)*ldt, 1 );
			}
		}
	}else{	//  Inverse of a lower triangular matrix.
		info = 0;

		for ( k = n; 1 <= k; k-- ){
			if ( t[k-1+(k-1)*ldt] == ZERO ){
				info = k;
				break;
			}

			t[k-1+(k-1)*ldt] = ONE / t[k-1+(k-1)*ldt];
			temp = -t[k-1+(k-1)*ldt];

			if ( k != n ){	dscal ( n-k, temp, t+k+(k-1)*ldt, 1 );	}

			for ( j = 1; j <= k-1; j++ ){
				temp = t[k-1+(j-1)*ldt];
				t[k-1+(j-1)*ldt] = ZERO;
				daxpy ( n-k+1, temp, t+k-1+(k-1)*ldt, 1, t+k-1+(j-1)*ldt, 1 );
			}
		}
	}
	return info;
}
//****************************************************************************80
int dtrsl ( double t[], int ldt, int n, double b[], int job )
//****************************************************************************80
//
//  Purpose:
//
//    DTRSL solves triangular linear systems.
//
//  Discussion:
//
//    DTRSL can solve T * X = B or T' * X = B where T is a triangular
//    matrix of order N.
//
//    Here T' denotes the transpose of the matrix T.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double T[LDT*N], the matrix of the system.  The zero
//    elements of the matrix are not referenced, and the corresponding
//    elements of the array can be used to store other information.
//
//    Input, int LDT, the leading dimension of the array T.
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
//    Input, int JOB, specifies what kind of system is to be solved:
//    00, solve T * X = B, T lower triangular,
//    01, solve T * X = B, T upper triangular,
//    10, solve T'* X = B, T lower triangular,
//    11, solve T'* X = B, T upper triangular.
//
//    Output, int DTRSL, singularity indicator.
//    0, the system is nonsingular.
//    nonzero, the index of the first zero diagonal element of T.
//
{
	int info, j, jj, kase;
	double temp;
	//  Check for zero diagonal elements.
	for ( j = 1; j <= n; j++ ){
		if ( t[j-1+(j-1)*ldt] == ZERO ){
			info = j;		return info;
		}
	}

	info = 0;
	//  Determine the task and go to it.
	if ( ( job % 10 ) == 0 ){	kase = 1;	}
	else{	kase = 2;	}

	if ( ( job % 100 ) / 10 != 0 ){	kase = kase + 2;}
	//  Solve T * X = B for T lower triangular.
	if ( kase == 1 ){
		b[0] = b[0] / t[0+0*ldt];
		for ( j = 2; j <= n; j++ ){
			temp = -b[j-2];
			daxpy ( n-j+1, temp, t+(j-1)+(j-2)*ldt, 1, b+j-1, 1 );
			b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
		}
	}else if ( kase == 2 ){	//  Solve T * X = B for T upper triangular.
		b[n-1] = b[n-1] / t[n-1+(n-1)*ldt];
		for ( jj = 2; jj <= n; jj++ ){
			j = n - jj + 1;
			temp = -b[j];
			daxpy ( j, temp, t+0+j*ldt, 1, b, 1 );
			b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
		}
	}else if ( kase == 3 ){	//  Solve T' * X = B for T lower triangular.
		b[n-1] = b[n-1] / t[n-1+(n-1)*ldt];
		for ( jj = 2; jj <= n; jj++ ){
			j = n - jj + 1;
			b[j-1] = b[j-1] - ddot ( jj-1, t+j+(j-1)*ldt, 1, b+j, 1 );
			b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
		}
	}else if ( kase == 4 ){	//  Solve T' * X = B for T upper triangular.
		b[0] = b[0] / t[0+0*ldt];
		for ( j = 2; j <= n; j++ ){
			b[j-1] = b[j-1] - ddot ( j-1, t+0+(j-1)*ldt, 1, b, 1 );
			b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
		}
	}

	return info;
}
//****************************************************************************80
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )
//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
	int i, ix, iy, m;

	if ( n <= 0 )	{	return;	}
	if ( da == ZERO ){	return;	}
	//  Code for unequal increments or equal increments
	//  not equal to 1.
	if ( incx != 1 || incy != 1 ){
		if ( 0 <= incx ){	ix = 0;	}
		else{	ix = ( - n + 1 ) * incx;	}

		if ( 0 <= incy ){	iy = 0;		}
		else{	iy = ( - n + 1 ) * incy;}

		for ( i = 0; i < n; i++ ){
			dy[iy] = dy[iy] + da * dx[ix];
			ix = ix + incx;
			iy = iy + incy;
		}
	}else{	//  Code for both increments equal to 1.
		m = n % 4;

		for ( i = 0; i < m; i++ ){	dy[i] = dy[i] + da * dx[i];	}

		for ( i = m; i < n; i = i + 4 ){
			dy[i  ] = dy[i  ] + da * dx[i  ];
			dy[i+1] = dy[i+1] + da * dx[i+1];
			dy[i+2] = dy[i+2] + da * dx[i+2];
			dy[i+3] = dy[i+3] + da * dx[i+3];
		}
	}
	return;
}
//****************************************************************************80
void dscal ( int n, double sa, double x[], int incx )
//****************************************************************************80
//
//  Purpose:
//
//    DSCAL scales a vector by a constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double SA, the multiplier.
//
//    Input/output, double X[*], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of X.
//
{
	int i;
	int ix;
	int m;

	if ( n <= 0 ){	}
	else if ( incx == 1 ){
		m = n % 5;

		for ( i = 0; i < m; i++ ){	x[i] = sa * x[i];	}
		for ( i = m; i < n; i = i + 5 ){
			x[i]   = sa * x[i];
			x[i+1] = sa * x[i+1];
			x[i+2] = sa * x[i+2];
			x[i+3] = sa * x[i+3];
			x[i+4] = sa * x[i+4];
		}
	}else{
		if ( 0 <= incx ){ix = 0;}
		else{	ix = ( - n + 1 ) * incx;	}

		for ( i = 0; i < n; i++ ){
			x[ix] = sa * x[ix];
			ix = ix + incx;
		}
	}
	return;
}
//****************************************************************************80
void drot ( int n, double x[], int incx, double y[], int incy, double c, double s )
//****************************************************************************80
//
//  Purpose:
//
//    DROT applies a plane rotation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to be rotated.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to be rotated.
//
//    Input, int INCY, the increment between successive elements of Y.
//
//    Input, double C, S, parameters (presumably the cosine and
//    sine of some angle) that define a plane rotation.
//
{
	int i, ix, iy;
	double stemp;

	if ( n <= 0 ){	}
	else if ( incx == 1 && incy == 1 ){
		for ( i = 0; i < n; i++ ){
			stemp = c * x[i] + s * y[i];
			y[i]  = c * y[i] - s * x[i];
			x[i]  = stemp;
		}
	}else{
		if ( 0 <= incx ){	ix = 0;						}
		else			{	ix = ( - n + 1 ) * incx;	}

		if ( 0 <= incy ){	iy = 0;						}
		else			{	iy = ( - n + 1 ) * incy;	}

		for ( i = 0; i < n; i++ ){
			stemp = c * x[ix] + s * y[iy];
			y[iy] = c * y[iy] - s * x[ix];
			x[ix] = stemp;
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return;
}
//****************************************************************************80
void drotg ( double *sa, double *sb, double *c, double *s )
//****************************************************************************80
//
//  Purpose:
//
//    DROTG constructs a Givens plane rotation.
//
//  Discussion:
//
//    Given values A and B, this routine computes
//
//    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
//          = sign ( B ) if abs ( A ) <= abs ( B );
//
//    R     = SIGMA * ( A * A + B * B );
//
//    C = A / R if R is not 0
//      = 1     if R is 0;
//
//    S = B / R if R is not 0,
//        0     if R is 0.
//
//    The computed numbers then satisfy the equation
//
//    (  C  S ) ( A ) = ( R )
//    ( -S  C ) ( B ) = ( 0 )
//
//    The routine also computes
//
//    Z = S     if abs ( A ) > abs ( B ),
//      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
//      = 1     if C is 0.
//
//    The single value Z encodes C and S, and hence the rotation:
//
//    If Z = 1, set C = 0 and S = 1;
//    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
//    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input/output, double *SA, *SB,  On input, SA and SB are the values
//    A and B.  On output, SA is overwritten with R, and SB is
//    overwritten with Z.
//
//    Output, double *C, *S, the cosine and sine of the Givens rotation.
//
{
	double r, roe, scale, z;

	if ( fabs ( *sb ) < fabs ( *sa ) ){	roe = *sa;	}
	else{	roe = *sb;	}

	scale = fabs ( *sa ) + fabs ( *sb );

	if ( scale == ZERO ){
		*c = ONE;		*s = ZERO;		r = ZERO;
	}else{
		r = scale * sqrt ( ( *sa / scale ) * ( *sa / scale )
				+ ( *sb / scale ) * ( *sb / scale ) );
		r = r8_sign ( roe ) * r;
		*c = *sa / r;
		*s = *sb / r;
	}

	if ( ZERO < fabs ( *c ) && fabs ( *c ) <= *s ){	z = ONE / *c;}
	else{	z = *s;	}

	*sa = r;	*sb = z;

	return;
}
//****************************************************************************80
void dswap ( int n, double x[], int incx, double y[], int incy )
//****************************************************************************80
//
//  Purpose:
//
//    DSWAP interchanges two vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of Y.
//
{
	int i, ix, iy, m;
	double temp;

	if ( n <= 0 ){	}
	else if ( incx == 1 && incy == 1 ){
		m = n % 3;

		for ( i = 0; i < m; i++ ){
			temp = x[i];			x[i] = y[i];			y[i] = temp;
		}

		for ( i = m; i < n; i = i + 3 ){
			temp = x[i];			x[i] = y[i];			y[i] = temp;
			temp = x[i+1];			x[i+1] = y[i+1];		y[i+1] = temp;
			temp = x[i+2];			x[i+2] = y[i+2];		y[i+2] = temp;
		}
	}else{
		if ( 0 <= incx ){	ix = 0;					}
		else			{	ix = ( - n + 1 ) * incx;}

		if ( 0 <= incy ){	iy = 0;					}
		else			{	iy = ( - n + 1 ) * incy;}

		for ( i = 0; i < n; i++ ){
			temp = x[ix];
			x[ix] = y[iy];
			y[iy] = temp;
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return;
}
//****************************************************************************80
double ddot ( int n, double dx[], int incx, double dy[], int incy )
//****************************************************************************80
//
//  Purpose:
//
//    DDOT forms the dot product of two vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries in DX.
//
//    Input, double DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries in DY.
//
//    Output, double DDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
	double dtemp = ZERO;
	int i, ix, iy, m;

	if ( n <= 0 ){	return dtemp;	}
	/*!
	 * Code for unequal increments or equal increments
	 * not equal to 1.
	 */
	if ( incx != 1 || incy != 1 ){
		if ( 0 <= incx ){	ix = 0;					}
		else			{	ix = ( - n + 1 ) * incx;}

		if ( 0 <= incy ){	iy = 0;					}
		else			{	iy = ( - n + 1 ) * incy;}

		for ( i = 0; i < n; i++ ){
			dtemp = dtemp + dx[ix] * dy[iy];
			ix = ix + incx;
			iy = iy + incy;
		}
	}else{	//  Code for both increments equal to 1.
		m = n % 5;

		for ( i = 0; i < m; i++ ){	dtemp = dtemp + dx[i] * dy[i];	}

		for ( i = m; i < n; i = i + 5 ){
			dtemp = dtemp + dx[i  ] * dy[i  ]  + dx[i+1] * dy[i+1]  + dx[i+2] * dy[i+2]
																				   + dx[i+3] * dy[i+3]  + dx[i+4] * dy[i+4];
		}
	}
	return dtemp;
}
//****************************************************************************80
double dnrm2 ( int n, double x[], int incx )
//****************************************************************************80
//
//  Purpose:
//
//    DNRM2 returns the euclidean norm of a vector.
//
//  Discussion:
//
//     DNRM2 ( X ) = sqrt ( X' * X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector whose norm is to be computed.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, double DNRM2, the Euclidean norm of X.
//
{
	double absxi, norm, scale, ssq;
	int i, ix;

	if ( n < 1 || incx < 1 ){	norm = ZERO;				}
	else if ( n == 1 )		{	norm = fabs ( x[0] );	}
	else{
		scale = ZERO;		ssq = ONE;		ix = 0;

		for ( i = 0; i < n; i++ ){
			if ( x[ix] != ZERO ){
				absxi = fabs ( x[ix] );
				if ( scale < absxi ){
					ssq = ONE + ssq * ( scale / absxi ) * ( scale / absxi );
					scale = absxi;
				}else{
					ssq = ssq + ( absxi / scale ) * ( absxi / scale );
				}
			}
			ix = ix + incx;
		}

		norm  = scale * sqrt ( ssq );
	}
	return norm;
}
//****************************************************************************80
double dasum ( int n, double x[], int incx )
//****************************************************************************80
//
//  Purpose:
//
//    DASUM takes the sum of the absolute values of a vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector to be examined.
//
//    Input, int INCX, the increment between successive entries of X.
//    INCX must not be negative.
//
//    Output, double DASUM, the sum of the absolute values of X.
//
{
	int i, j = 0;
	double value = ZERO;

	for ( i = 0; i < n; i++ ){
		value = value + fabs ( x[j] );
		j = j + incx;
	}
	return value;
}
//****************************************************************************80
int idamax ( int n, double dx[], int incx )
//****************************************************************************80
//
//  Purpose:
//
//    IDAMAX finds the index of the vector element of maximum absolute value.
//
//  Discussion:
//
//    WARNING: This index is a 1-based index, not a 0-based index!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector to be examined.
//
//    Input, int INCX, the increment between successive entries of SX.
//
//    Output, int IDAMAX, the index of the element of maximum
//    absolute value.
//
{
	double dmax;
	int i, ix, value;

	value = 0;

	if ( n < 1 || incx <= 0 ){	return value;	}

	value = 1;

	if ( n == 1 )			{	return value;	}

	if ( incx == 1 ){
		dmax = fabs ( dx[0] );
		for ( i = 1; i < n; i++ ){
			if ( dmax < fabs ( dx[i] ) ){
				value = i + 1;
				dmax = fabs ( dx[i] );
			}
		}
	}else{
		ix = 0;
		dmax = fabs ( dx[0] );
		ix = ix + incx;

		for ( i = 1; i < n; i++ ){
			if ( dmax < fabs ( dx[ix] ) ){
				value = i + 1;
				dmax = fabs ( dx[ix] );
			}
			ix = ix + incx;
		}
	}
	return value;
}
//****************************************************************************80
void dcopy ( int n, double dx[], int incx, double dy[], int incy )
//****************************************************************************80
//
//  Purpose:
//
//    DCOPY copies a vector X to a vector Y.
//
//  Discussion:
//
//    The routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Output, double DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
	int i, ix, iy, m;

	if ( n <= 0 ){	return;	}

	if ( incx == 1 && incy == 1 ){
		m = n % 7;

		if ( m != 0 ){	for ( i = 0; i < m; i++ ){	dy[i] = dx[i];	}	}

		for ( i = m; i < n; i = i + 7 )	{
			dy[i] = dx[i];
			dy[i + 1] = dx[i + 1];			dy[i + 2] = dx[i + 2];
			dy[i + 3] = dx[i + 3];			dy[i + 4] = dx[i + 4];
			dy[i + 5] = dx[i + 5];			dy[i + 6] = dx[i + 6];
		}
	}else{
		if ( 0 <= incx ){	ix = 0;					}
		else			{	ix = ( -n + 1 ) * incx;	}

		if ( 0 <= incy ){	iy = 0;					}
		else			{	iy = ( -n + 1 ) * incy;	}

		for ( i = 0; i < n; i++ ){
			dy[iy] = dx[ix];
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return;
}
//****************************************************************************80
int i4_min ( int i1, int i2 )
//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
	int value;

	if ( i1 < i2 )	{	value = i1;	}
	else			{	value = i2;	}
	return value;
}
//****************************************************************************80
int i4_max ( int i1, int i2 )
//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
	int value;

	if ( i2 < i1 )	{	value = i1;	}
	else			{	value = i2;	}
	return value;
}
//****************************************************************************80
double r8_abs ( double x )
//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Discussion:
//
//    The C++ math library provides the function fabs() which is preferred.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
	double value;

	if ( ZERO <= x )	{	value = + x;	}
	else			{	value = - x;	}
	return value;
}
//****************************************************************************80
double r8_max ( double x, double y )
//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Discussion:
//
//    The C++ math library provides the function fmax() which is preferred.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
	double value;

	if ( y < x ){	value = x;	}
	else		{	value = y;	}
	return value;
}
//****************************************************************************80
double r8_sign ( double x )
//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
	double value;

	if ( x < ZERO )	{	value = -ONE;	}
	else			{	value = ONE;	}
	return value;
}

bool lsame_(char *ca, char *cb)
/*  -- LAPACK auxiliary routine (version 2.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
       Courant Institute, Argonne National Lab, and Rice University
       January 31, 1994

    Purpose
    =======

    LSAME returns .TRUE. if CA is the same letter as CB regardless of
    case. */
{
	/* System generated locals */
	bool ret_val;

	/* Local variables */
	int inta, intb, zcode;

	ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
	if (ret_val) {
		return ret_val;
	}
	/*     Now test for equivalence if both characters are alphabetic. */
	zcode = 'Z';
	/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
       machines, on which ICHAR returns a value with bit 8 set.
       ICHAR('A') on Prime machines returns 193 which is the same as
       ICHAR('A') on an EBCDIC machine. */

	inta = *(unsigned char *)ca;
	intb = *(unsigned char *)cb;

	if ((zcode == 90) || (zcode == 122)) {

		/* ASCII is assumed - ZCODE is the ASCII code of either lower or upper case 'Z'. */

		if ((inta >= 97) && (inta <= 122)) {	inta += -32;	}
		if ((intb >= 97) && (intb <= 122)) {	intb += -32;	}

	} else if ((zcode == 233) || (zcode == 169)) {

		/* EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or upper case 'Z'. */

		if (((inta >= 129) && (inta <= 137)) || ((inta >= 145) && (inta <= 153)) || ((inta >= 162) && (inta <= 169))) {	inta += 64;	}
		if (((intb >= 129) && (intb <= 137)) || ((intb >= 145) && (intb <= 153)) || ((intb >= 162) && (intb <= 169))) {	intb += 64;	}

	} else if (zcode == 218 || zcode == 250) {
		/* ASCII is assumed, on Prime machines - ZCODE is the ASCII code
          plus 128 of either lower or upper case 'Z'. */

		if (inta >= 225 && inta <= 250) {	inta += -32;	}
		if (intb >= 225 && intb <= 250) {	intb += -32;	}
	}
	ret_val = inta == intb;

	return ret_val;
} /* lsame_ */
/* Subroutine */ int dscal_(int *n, double *da, double *dx, int *incx)
/*     scales a vector by a constant.
       uses unrolled loops for increment equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 3/93 to return if incx .le. 0.
       modified 12/3/93, array(1) declarations changed to array(*)*/
{
	/* System generated locals */
	int i__1, i__2;
	/* Local variables */
	static int i, m, nincx, mp1;
	/*  Function Body */
#define DX(I) dx[(I)-1]
	if (*n <= 0 || *incx <= 0) {	return 0;	}
	if (*incx == 1) {	goto L20;	}
	/*        code for increment not equal to 1 */
	nincx = *n * *incx;
	i__1 = nincx;	i__2 = *incx;
	for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {	DX(i) = *da * DX(i);	}
	return 0;
	/*        code for increment equal to 1. clean-up loop */
	L20:
	m = *n % 5;
	if (m == 0) {	goto L40;	}
	i__2 = m;
	for (i = 1; i <= m; ++i) {	DX(i) = *da * DX(i);	}
	if (*n < 5) {	return 0;	}
	L40:
	mp1 = m + 1;	i__2 = *n;
	for (i = mp1; i <= *n; i += 5) {
		DX(i) = *da * DX(i);
		DX(i + 1) = *da * DX(i + 1);
		DX(i + 2) = *da * DX(i + 2);
		DX(i + 3) = *da * DX(i + 3);
		DX(i + 4) = *da * DX(i + 4);
	}
	return 0;
} /* dscal_ */
/* Subroutine */ int dlabad_(double *small, double *large)
/*  -- LAPACK auxiliary routine (version 3.1) --
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
       November 2006

    Purpose
    =======

    DLABAD takes as input the values computed by DLAMCH for underflow and
    overflow, and returns the square root of each of these values if the
    log of LARGE is sufficiently large.  This subroutine is intended to
    identify machines with a large exponent range, such as the Crays, and
    redefine the underflow and overflow limits to be the square roots of
    the values computed by DLAMCH.  This subroutine is needed because
    DLAMCH does not compensate for poor arithmetic in the upper half of
    the exponent range, as is found on a Cray.*/
{
	/* Builtin functions */
#define log10e 0.43429448190325182765
	double tmp = log10e * log(*large);
	if (tmp > 2e3) {
		*small = sqrt(*small);
		*large = sqrt(*large);
	}

	return 0;
} /* dlabad_ */
/* Subroutine */ int drscl_(int *n, double *sa, double *sx, int *incx)
/*  -- LAPACK auxiliary routine (version 3.1) --
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
       November 2006


    Purpose
    =======

    DRSCL multiplies an n-element real vector x by the real scalar 1/a.
    This is done without overflow or underflow as long as
    the final result x/a does not overflow or underflow.*/
{
	double mul, cden;
	bool done;
	double cnum, cden1, cnum1;
	double bignum, smlnum;

	--sx;
	char S[] = "S";
	/* Function Body */
	if (*n <= 0) {		return 0;	}
	/*     Get machine parameters */
	smlnum = 1.17549E-38;
	bignum = 1. / smlnum;
	//	dlabad_(&smlnum, &bignum);
#define log10e 0.43429448190325182765
	double tmp = log10e * log(bignum);
	if (tmp > 2e3) {
		smlnum = sqrt(smlnum);
		bignum = sqrt(bignum);
	}
	/*     Initialize the denominator to SA and the numerator to 1. */
	cden = *sa;		cnum = 1.;
	L10:
	cden1 = cden * smlnum;		cnum1 = cnum / bignum;
	if (abs(cden1) > abs(cnum) && cnum != 0.) {
		/* Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. */
		mul = smlnum;		done = false;		cden = cden1;
	} else if (abs(cnum1) > abs(cden)) {

		/* Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. */

		mul = bignum;		done = false;		cnum = cnum1;
	} else {
		/* Multiply X by CNUM / CDEN and return. */
		mul = cnum / cden;		done = true;
	}
	/* Scale the vector X by MUL */
	dscal_(n, &mul, &sx[1], incx);
	if (! done) {	goto L10;	}

	return 0;
} /* drscl_ */
double dasum_(int *n, double *dx, int *incx)
/*     takes the sum of the absolute values.
   jack dongarra, linpack, 3/11/78.
   modified 3/93 to return if incx .le. 0.
   modified 12/3/93, array(1) declarations changed to array(*)*/
{
	/* System generated locals */
	int i__1, i__2;
	double ret_val, d__1, d__2, d__3, d__4, d__5, d__6;
	/* Local variables */
	int i, m;
	double dtemp;
	int nincx, mp1;
#define DX(I) dx[(I)-1]
	ret_val = 0.;	dtemp = 0.;
	if (*n <= 0 || *incx <= 0) {	return ret_val;	}
	if (*incx == 1) {		goto L20;	}
	/*        code for increment not equal to 1 */
	nincx = *n * *incx;		i__1 = nincx;		i__2 = *incx;
	for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
		dtemp += (d__1 = DX(i), abs(d__1));
	}
	ret_val = dtemp;
	return ret_val;
	/*  code for increment equal to 1. clean-up loop */
	L20:
	m = *n % 6;
	if (m == 0) {	goto L40;	}
	i__2 = m;
	for (i = 1; i <= m; ++i) {	dtemp += (d__1 = DX(i), abs(d__1));	}
	if (*n < 6) {	goto L60;	}
	L40:
	mp1 = m + 1;	i__2 = *n;
	for (i = mp1; i <= *n; i += 6) {
		dtemp = dtemp + (d__1 = DX(i), abs(d__1)) + (d__2 = DX(i + 1), abs(
				d__2)) + (d__3 = DX(i + 2), abs(d__3)) + (d__4 = DX(i + 3),
						abs(d__4)) + (d__5 = DX(i + 4), abs(d__5)) + (d__6 = DX(i + 5)
								, abs(d__6));
	}
	L60:
	ret_val = dtemp;
	return ret_val;
} /* dasum_ */
/* Subroutine */ int dcopy_(int *n, double *dx, int *incx, double *dy, int *incy)
/*     copies a vector, x, to a vector, y.
       uses unrolled loops for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*)*/
{
	/* System generated locals */
	int i__1;
	/* Local variables */
	static int i, m, ix, iy, mp1;
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]
	if (*n <= 0) {	return 0;	}
	if (*incx == 1 && *incy == 1) {	goto L20;	}
	/* code for unequal increments or equal increments not equal to 1 */
	ix = 1;		iy = 1;
	if (*incx < 0) {	ix = (-(*n) + 1) * *incx + 1;	}
	if (*incy < 0) {	iy = (-(*n) + 1) * *incy + 1;	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
		DY(iy) = DX(ix);
		ix += *incx;
		iy += *incy;
	}
	return 0;
	/* code for both increments equal to 1. clean-up loop */
	L20:
	m = *n % 7;
	if (m == 0) {	goto L40;	}
	i__1 = m;
	for (i = 1; i <= m; ++i) {	DY(i) = DX(i);	}
	if (*n < 7) {	return 0;	}
	L40:
	mp1 = m + 1;	i__1 = *n;
	for (i = mp1; i <= *n; i += 7) {
		DY(i) = DX(i);
		DY(i + 1) = DX(i + 1);
		DY(i + 2) = DX(i + 2);
		DY(i + 3) = DX(i + 3);
		DY(i + 4) = DX(i + 4);
		DY(i + 5) = DX(i + 5);
		DY(i + 6) = DX(i + 6);
	}
	return 0;
} /* dcopy_ */
int idamax_(int *n, double *dx, int *incx)
/*     finds the index of element having max. absolute value.
       jack dongarra, linpack, 3/11/78.
       modified 3/93 to return if incx .le. 0.
       modified 12/3/93, array(1) declarations changed to array(*) */
{
	/* System generated locals */
	int ret_val, i__1;
	double d__1;
	/* Local variables */
	double dmax__;
	int i, ix;
#define DX(I) dx[(I)-1]
	ret_val = 0;
	if (*n < 1 || *incx <= 0) {	return ret_val;	}
	ret_val = 1;
	if (*n == 1) {	return ret_val;	}
	if (*incx == 1) {	goto L20;	}
	/*        code for increment not equal to 1 */
	ix = 1;
	dmax__ = abs(DX(1));
	ix += *incx;		i__1 = *n;
	for (i = 2; i <= *n; ++i) {
		if ((d__1 = DX(ix), abs(d__1)) <= dmax__) {		goto L5;	}
		ret_val = i;
		dmax__ = (d__1 = DX(ix), abs(d__1));
		L5:
		ix += *incx;
		/* L10: */
	}
	return ret_val;
	/*        code for increment equal to 1 */
	L20:
	dmax__ = abs(DX(1));		i__1 = *n;
	for (i = 2; i <= *n; ++i) {
		if ((d__1 = DX(i), abs(d__1)) <= dmax__) {		goto L30;	}
		ret_val = i;
		dmax__ = (d__1 = DX(i), abs(d__1));
		L30:
		;
	}
	return ret_val;
} /* idamax_ */
/* Subroutine */ int dlacn2_(int *n, double *v, double *x,
		int *isgn, double *est, int *kase, int *isave)

/*  -- LAPACK auxiliary routine (version 3.1) --
   Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   November 2006
Purpose
=======

DLACN2 estimates the 1-norm of a square, real matrix A.
Reverse communication is used for evaluating matrix-vector products.*/
{
	/* Table of constant values */
	int c__1 = 1;
	double c_b11 = 1.;

	/* System generated locals */
	int i__1;
	double d__1;
	/* Builtin functions */
#define d_sign(a,b) ( *b >= 0 ? (*a >= 0 ? *a : - *a) : -(*a >= 0 ? *a : - *a))
#define i_dnnt(x) (int)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x))
	/* Local variables */
	int i__;
	double temp;
	int jlast;
	double altsgn, estold;
	--isave;	--isgn;		--x;	--v;
	/* Function Body */
	if (*kase == 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {	x[i__] = 1. / (double) (*n);	}
		*kase = 1;		isave[1] = 1;
		return 0;
	}
	switch (isave[1]) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L110;
	case 5:  goto L140;
	}
	/*     ................ ENTRY   (ISAVE( 1 ) = 1)
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */
	L20:
	if (*n == 1) {
		v[1] = x[1];
		*est = abs(v[1]);
		/*        ... QUIT */
		goto L150;
	}
	*est = dasum_(n, &x[1], &c__1);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__] = d_sign(&c_b11, &x[i__]);
		isgn[i__] = i_dnnt(&x[i__]);
	}
	*kase = 2;		isave[1] = 2;
	return 0;
	/*     ................ ENTRY   (ISAVE( 1 ) = 2)
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
	L40:
	isave[2] = idamax_(n, &x[1], &c__1);
	isave[3] = 2;
	/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */
	L50:
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {		x[i__] = 0.;	}
	x[isave[2]] = 1.;		*kase = 1;		isave[1] = 3;
	return 0;
	/*     ................ ENTRY   (ISAVE( 1 ) = 3)
       X HAS BEEN OVERWRITTEN BY A*X. */
	L70:
	dcopy_(n, &x[1], &c__1, &v[1], &c__1);
	estold = *est;
	*est = dasum_(n, &v[1], &c__1);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		d__1 = d_sign(&c_b11, &x[i__]);
		if (i_dnnt(&d__1) != isgn[i__]) {	goto L90;	}
	}
	/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
	goto L120;
	L90:
	/*     TEST FOR CYCLING. */
	if (*est <= estold) {	goto L120;	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__] = d_sign(&c_b11, &x[i__]);
		isgn[i__] = i_dnnt(&x[i__]);
	}
	*kase = 2;
	isave[1] = 4;
	return 0;
	/*     ................ ENTRY   (ISAVE( 1 ) = 4)
       X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
	L110:
	jlast = isave[2];
	isave[2] = idamax_(n, &x[1], &c__1);
	if (x[jlast] != (d__1 = x[isave[2]], abs(d__1)) && isave[3] < 5) {
		++isave[3];
		goto L50;
	}
	/*     ITERATION COMPLETE.  FINAL STAGE. */
	L120:
	altsgn = 1.;	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__] = altsgn * ((double) (i__ - 1) / (double) (*n - 1) +
				1.);
		altsgn = -altsgn;
		/* L130: */
	}
	*kase = 1;
	isave[1] = 5;
	return 0;
	/*     ................ ENTRY   (ISAVE( 1 ) = 5)
       X HAS BEEN OVERWRITTEN BY A*X. */
	L140:
	temp = dasum_(n, &x[1], &c__1) / (double) (*n * 3) * 2.;
	if (temp > *est) {
		dcopy_(n, &x[1], &c__1, &v[1], &c__1);
		*est = temp;
	}
	L150:
	*kase = 0;
	return 0;
} /* dlacn2_ */
/* Subroutine */ int daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy)
/*     constant times a vector plus a vector.
       uses unrolled loops for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*)  */
{
	/* System generated locals */
	int i__1;
	/* Local variables */
	int i, m, ix, iy, mp1;
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]
	if (*n <= 0) 	{	return 0;	}
	if (*da == 0.)  {	return 0;	}
	if (*incx == 1 && *incy == 1) {		goto L20;	}
	/* code for unequal increments or equal increments not equal to 1 */
	ix = 1;		iy = 1;
	if (*incx < 0) {	ix = (-(*n) + 1) * *incx + 1;	}
	if (*incy < 0) {	iy = (-(*n) + 1) * *incy + 1;	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
		DY(iy) += *da * DX(ix);
		ix += *incx;		iy += *incy;
	}
	return 0;
	/* code for both increments equal to 1. clean-up loop */
	L20:
	m = *n % 4;
	if (m == 0) {	goto L40;	}
	i__1 = m;
	for (i = 1; i <= m; ++i) {	DY(i) += *da * DX(i);	}
	if (*n < 4) {	return 0;	}
	L40:
	mp1 = m + 1;
	i__1 = *n;
	for (i = mp1; i <= *n; i += 4) {
		DY(i) += *da * DX(i);
		DY(i + 1) += *da * DX(i + 1);
		DY(i + 2) += *da * DX(i + 2);
		DY(i + 3) += *da * DX(i + 3);
	}
	return 0;
} /* daxpy_ */
double ddot_(int *n, double *dx, int *incx, double *dy, int *incy)
/*     forms the dot product of two vectors.
       uses unrolled loops for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*)*/
{
	/* System generated locals */
	int i__1;
	double ret_val;
	/* Local variables */
	int i, m;
	double dtemp;
	int ix, iy, mp1;
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]
	ret_val = 0.;		dtemp = 0.;
	if (*n <= 0) {	return ret_val;	}
	if (*incx == 1 && *incy == 1) {	goto L20;	}

	/* code for unequal increments or equal increments not equal to 1 */

	ix = 1;		iy = 1;
	if (*incx < 0) {	ix = (-(*n) + 1) * *incx + 1;	}
	if (*incy < 0) {	iy = (-(*n) + 1) * *incy + 1;	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
		dtemp += DX(ix) * DY(iy);
		ix += *incx;
		iy += *incy;
	}
	ret_val = dtemp;
	return ret_val;
	/* code for both increments equal to 1. clean-up loop */
	L20:
	m = *n % 5;
	if (m == 0) {	goto L40;	}
	i__1 = m;
	for (i = 1; i <= m; ++i) {	dtemp += DX(i) * DY(i);	}
	if (*n < 5) {	goto L60;	}
	L40:
	mp1 = m + 1;	i__1 = *n;
	for (i = mp1; i <= *n; i += 5) {
		dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) *
				DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
	}
	L60:
	ret_val = dtemp;
	return ret_val;
} /* ddot_ */
/* Subroutine */ int dtrsv_(char *uplo, char *trans, char *diag, int *n,
		double *a, int *lda, double *x, int *incx)
/*  Purpose
=======

DTRSV  solves one of the systems of equations

   A*x = b,   or   A'*x = b,

where b and x are n element vectors and A is an n by n unit, or
non-unit, upper or lower triangular matrix.

No test for singularity or near-singularity is included in this
routine. Such tests must be performed before calling this routine.  */
{
	/* System generated locals */
	int i__1, i__2;
	/* Local variables */
	int info;
	double temp;
	int i, j;
	int ix, jx, kx;
	bool nounit;
#define X(I) x[(I)-1]
#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
	info = 0;
	char U[] = "U", L[] = "L", N[] = "N", T[] = "T", C[] = "C";
	if (! lsame_(uplo, U) && ! lsame_(uplo, L)) {	info = 1;	}
	else if (! lsame_(trans, N) && ! lsame_(trans, T) &&! lsame_(trans, C)) {	info = 2;	}
	else if (! lsame_(diag, U) && ! lsame_(diag, N)) {	info = 3;	}
	else if (*n < 0) {	info = 4;	}
	else if (*lda < max(1,*n)) {	info = 6;	}
	else if (*incx == 0) {	info = 8;	}
	if (info != 0) {
		//xerbla_("DTRSV ", &info);
		std::cout << "** On entry to DTRSV, parameter number " << info <<  "had an illegal value\n";
		return 0;
	}
	/*     Quick return if possible. */
	if (*n == 0) {	return 0;	}
	nounit = lsame_(diag, N);
	/*     Set up the start point in X if the increment is not unity. This
       will be  ( N - 1 )*INCX  too small for descending loops. */
	if (*incx <= 0) {		kx = 1 - (*n - 1) * *incx;	}
	else if (*incx != 1) {	kx = 1;	}
	/*     Start the operations. In this version the elements of A are
       accessed sequentially with one pass through A. */
	if (lsame_(trans, N)) {
		/*        Form  x := inv( A )*x. */
		if (lsame_(uplo, U)) {
			if (*incx == 1) {
				for (j = *n; j >= 1; --j) {
					if (X(j) != 0.) {
						if (nounit) {	X(j) /= A(j,j);		}
						temp = X(j);
						for (i = j - 1; i >= 1; --i) {	X(i) -= temp * A(i,j);	}
					}
				}
			} else {
				jx = kx + (*n - 1) * *incx;
				for (j = *n; j >= 1; --j) {
					if (X(jx) != 0.) {
						if (nounit) {	X(jx) /= A(j,j);	}
						temp = X(jx);		ix = jx;
						for (i = j - 1; i >= 1; --i) {
							ix -= *incx;	X(ix) -= temp * A(i,j);
						}
					}
					jx -= *incx;
				}
			}
		} else {
			if (*incx == 1) {
				i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					if (X(j) != 0.) {
						if (nounit) {	X(j) /= A(j,j);		}
						temp = X(j);		i__2 = *n;
						for (i = j + 1; i <= *n; ++i) {	X(i) -= temp * A(i,j);	}
					}
				}
			} else {
				jx = kx;		i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					if (X(jx) != 0.) {
						if (nounit) {	X(jx) /= A(j,j);	}
						temp = X(jx);	ix = jx;	i__2 = *n;
						for (i = j + 1; i <= *n; ++i) {
							ix += *incx;		X(ix) -= temp * A(i,j);
						}
					}
					jx += *incx;
				}
			}
		}
	} else {
		/*        Form  x := inv( A' )*x. */
		if (lsame_(uplo, U)) {
			if (*incx == 1) {
				i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					temp = X(j);		i__2 = j - 1;
					for (i = 1; i <= j-1; ++i) {temp -= A(i,j) * X(i);	}
					if (nounit) {	temp /= A(j,j);	}
					X(j) = temp;
				}
			} else {
				jx = kx;		i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					temp = X(jx);		ix = kx;		i__2 = j - 1;
					for (i = 1; i <= j-1; ++i) {
						temp -= A(i,j) * X(ix);			ix += *incx;
					}
					if (nounit) {	temp /= A(j,j);		}
					X(jx) = temp;		jx += *incx;
				}
			}
		} else {
			if (*incx == 1) {
				for (j = *n; j >= 1; --j) {
					temp = X(j);		i__1 = j + 1;
					for (i = *n; i >= j+1; --i) {	temp -= A(i,j) * X(i);	}
					if (nounit) {	temp /= A(j,j);		}
					X(j) = temp;
				}
			} else {
				kx += (*n - 1) * *incx;		jx = kx;
				for (j = *n; j >= 1; --j) {
					temp = X(jx);		ix = kx;		i__1 = j + 1;
					for (i = *n; i >= j+1; --i) {
						temp -= A(i,j) * X(ix);			ix -= *incx;
					}
					if (nounit) {	temp /= A(j,j);	}
					X(jx) = temp;	jx -= *incx;
				}
			}
		}
	}
	return 0;
} /* dtrsv_ */
/* Subroutine */ int dlatrs_(char *uplo, char *trans, char *diag, char *
		normin, int *n, double *a, int *lda, double *x,
		double *scale, double *cnorm, int *info)
/*  -- LAPACK auxiliary routine (version 3.1) --
   Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   November 2006


Purpose
=======

DLATRS solves one of the triangular systems

   A *x = s*b  or  A'*x = s*b

with scaling to prevent overflow.  Here A is an upper or lower
triangular matrix, A' denotes the transpose of A, x and b are
n-element vectors, and s is a scaling factor, usually less than
or equal to 1, chosen so that the components of x will be less than
the overflow threshold.  If the unscaled problem will not cause
overflow, the Level 2 BLAS routine DTRSV is called.  If the matrix A
is singular (A(j,j) = 0 for some j), then s is set to 0 and a
non-trivial solution to A*x = 0 is returned. */
{
	/* Table of constant values */
	int c__1 = 1;
	double c_b36 = .5;

	/* System generated locals */
	int a_dim1, a_offset, i__1, i__2, i__3;
	double d__1, d__2, d__3;
	/* Local variables */
	int i__, j;
	double xj, rec, tjj;
	int jinc;
	double xbnd;
	int imax;
	double tmax, tjjs, xmax, grow, sumj;
	double tscal, uscal;
	int jlast;
	bool upper;
	double bignum;
	bool notran;
	int jfirst;
	double smlnum;
	bool nounit;
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--x;	--cnorm;

	/* Function Body */
	*info = 0;
	char U[] = "U", L[] = "L", N[] = "N", T[] = "T", C[] = "C", Y[] = "Y";
	upper = lsame_(uplo, U);
	notran = lsame_(trans, N);
	nounit = lsame_(diag, N);
	/*     Test the input parameters. */
	if (! upper && ! lsame_(uplo, L)) {	*info = -1;	}
	else if (! notran && ! lsame_(trans, T) && !	lsame_(trans, C)) {	*info = -2;	}
	else if (! nounit && ! lsame_(diag, U)) {	*info = -3;	}
	else if (! lsame_(normin, Y) && ! lsame_(normin, N)) {	*info = -4;	}
	else if (*n < 0) {	*info = -5;	}
	else if (*lda < max(1,*n)) {	*info = -7;	}
	if (*info != 0) {
		i__1 = -(*info);
		//xerbla_("DLATRS", &i__1);
		std::cout << "** On entry to DLATRS, parameter number " << i__1 <<  "had an illegal value\n";
		return 0;
	}
	/*     Quick return if possible */
	if (*n == 0) {	return 0;	}
	/*     Determine machine dependent parameters to control overflow. */
	smlnum = 1.17549E-38 / 1.19209E-07;
	bignum = 1. / smlnum;
	*scale = 1.;
	if (lsame_(normin, N)) {
		/*        Compute the 1-norm of each column, not including the diagonal. */
		if (upper) {
			/*           A is upper triangular. */
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = j - 1;
				cnorm[j] = dasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
				/* L10: */
			}
		} else {
			/*           A is lower triangular. */
			i__1 = *n - 1;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *n - j;
				cnorm[j] = dasum_(&i__2, &a[j + 1 + j * a_dim1], &c__1);
				/* L20: */
			}
			cnorm[*n] = 0.;
		}
	}
	/*     Scale the column norms by TSCAL if the maximum element in CNORM is
       greater than BIGNUM. */
	imax = idamax_(n, &cnorm[1], &c__1);
	tmax = cnorm[imax];
	if (tmax <= bignum) {	tscal = 1.;	}
	else {
		tscal = 1. / (smlnum * tmax);
		dscal_(n, &tscal, &cnorm[1], &c__1);
	}

	/*     Compute a bound on the computed solution vector to see if the
       Level 2 BLAS routine DTRSV can be used. */

	j = idamax_(n, &x[1], &c__1);
	xmax = (d__1 = x[j], abs(d__1));
	xbnd = xmax;
	if (notran) {
		/*        Compute the growth in A * x = b. */
		if (upper) {
			jfirst = *n;		jlast = 1;			jinc = -1;
		} else {
			jfirst = 1;			jlast = *n;			jinc = 1;
		}
		if (tscal != 1.) {
			grow = 0.;			goto L50;
		}
		if (nounit) {
			/*           A is non-unit triangular.
             Compute GROW = 1/G(j) and XBND = 1/M(j).
             Initially, G(0) = max{x(i), i=1,...,n}. */
			grow = 1. / max(xbnd,smlnum);
			xbnd = grow;			i__1 = jlast;			i__2 = jinc;
			for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
				/*              Exit the loop if the growth factor is too small. */
				if (grow <= smlnum) {	goto L50;	}
				/*              M(j) = G(j-1) / abs(A(j,j)) */
				tjj = (d__1 = a[j + j * a_dim1], abs(d__1));
				/* Computing MIN */
				d__1 = xbnd, d__2 = min(1.,tjj) * grow;
				xbnd = min(d__1,d__2);
				if (tjj + cnorm[j] >= smlnum) {
					/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */
					grow *= tjj / (tjj + cnorm[j]);
				} else {
					/*                 G(j) could overflow, set GROW to 0. */
					grow = 0.;
				}
			}
			grow = xbnd;
		} else {
			/*           A is unit triangular.
             Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
   Computing MIN */
			d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
			grow = min(d__1,d__2);
			i__2 = jlast;
			i__1 = jinc;
			for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
				/*              Exit the loop if the growth factor is too small. */
				if (grow <= smlnum) {	goto L50;	}

				/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

				grow *= 1. / (cnorm[j] + 1.);
				/* L40: */
			}
		}
		L50:

		;
	} else {

		/*        Compute the growth in A' * x = b. */

		if (upper) {
			jfirst = 1;
			jlast = *n;
			jinc = 1;
		} else {
			jfirst = *n;
			jlast = 1;
			jinc = -1;
		}

		if (tscal != 1.) {
			grow = 0.;
			goto L80;
		}

		if (nounit) {

			/*           A is non-unit triangular.

             Compute GROW = 1/G(j) and XBND = 1/M(j).
             Initially, M(0) = max{x(i), i=1,...,n}. */

			grow = 1. / max(xbnd,smlnum);
			xbnd = grow;
			i__1 = jlast;
			i__2 = jinc;
			for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

				/*              Exit the loop if the growth factor is too small. */

				if (grow <= smlnum) {
					goto L80;
				}

				/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

				xj = cnorm[j] + 1.;
				/* Computing MIN */
				d__1 = grow, d__2 = xbnd / xj;
				grow = min(d__1,d__2);

				/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

				tjj = (d__1 = a[j + j * a_dim1], abs(d__1));
				if (xj > tjj) {
					xbnd *= tjj / xj;
				}
				/* L60: */
			}
			grow = min(grow,xbnd);
		} else {

			/*           A is unit triangular.

             Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.

   Computing MIN */
			d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
			grow = min(d__1,d__2);
			i__2 = jlast;
			i__1 = jinc;
			for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

				/*              Exit the loop if the growth factor is too small. */

				if (grow <= smlnum) {
					goto L80;
				}

				/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

				xj = cnorm[j] + 1.;
				grow /= xj;
				/* L70: */
			}
		}
		L80:
		;
	}

	if (grow * tscal > smlnum) {

		/*        Use the Level 2 BLAS solve if the reciprocal of the bound on
          elements of X is not too small. */

		dtrsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1);
	} else {

		/*        Use a Level 1 BLAS solve, scaling intermediate results. */

		if (xmax > bignum) {

			/*           Scale X so that its components are less than or equal to
             BIGNUM in absolute value. */

			*scale = bignum / xmax;
			dscal_(n, scale, &x[1], &c__1);
			xmax = bignum;
		}

		if (notran) {

			/*           Solve A * x = b */

			i__1 = jlast;
			i__2 = jinc;
			for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

				/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

				xj = (d__1 = x[j], abs(d__1));
				if (nounit) {
					tjjs = a[j + j * a_dim1] * tscal;
				} else {
					tjjs = tscal;
					if (tscal == 1.) {
						goto L100;
					}
				}
				tjj = abs(tjjs);
				if (tjj > smlnum) {

					/*                    abs(A(j,j)) > SMLNUM: */

					if (tjj < 1.) {
						if (xj > tjj * bignum) {

							/*                          Scale x by 1/b(j). */

							rec = 1. / xj;
							dscal_(n, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}
					x[j] /= tjjs;
					xj = (d__1 = x[j], abs(d__1));
				} else if (tjj > 0.) {

					/*                    0 < abs(A(j,j)) <= SMLNUM: */

					if (xj > tjj * bignum) {

						/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                         to avoid overflow when dividing by A(j,j). */

						rec = tjj * bignum / xj;
						if (cnorm[j] > 1.) {

							/*                          Scale by 1/CNORM(j) to avoid overflow when
                            multiplying x(j) times column j. */

							rec /= cnorm[j];
						}
						dscal_(n, &rec, &x[1], &c__1);
						*scale *= rec;
						xmax *= rec;
					}
					x[j] /= tjjs;
					xj = (d__1 = x[j], abs(d__1));
				} else {

					/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                      scale = 0, and compute a solution to A*x = 0. */

					i__3 = *n;
					for (i__ = 1; i__ <= i__3; ++i__) {
						x[i__] = 0.;
						/* L90: */
					}
					x[j] = 1.;
					xj = 1.;
					*scale = 0.;
					xmax = 0.;
				}
				L100:

				/*              Scale x if necessary to avoid overflow when adding a
                multiple of column j of A. */

				if (xj > 1.) {
					rec = 1. / xj;
					if (cnorm[j] > (bignum - xmax) * rec) {

						/*                    Scale x by 1/(2*abs(x(j))). */

						rec *= .5;
						dscal_(n, &rec, &x[1], &c__1);
						*scale *= rec;
					}
				} else if (xj * cnorm[j] > bignum - xmax) {

					/*                 Scale x by 1/2. */

					dscal_(n, &c_b36, &x[1], &c__1);
					*scale *= .5;
				}

				if (upper) {
					if (j > 1) {

						/*                    Compute the update
                         x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

						i__3 = j - 1;
						d__1 = -x[j] * tscal;
						daxpy_(&i__3, &d__1, &a[j * a_dim1 + 1], &c__1, &x[1],
								&c__1);
						i__3 = j - 1;
						i__ = idamax_(&i__3, &x[1], &c__1);
						xmax = (d__1 = x[i__], abs(d__1));
					}
				} else {
					if (j < *n) {

						/*                    Compute the update
                         x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

						i__3 = *n - j;
						d__1 = -x[j] * tscal;
						daxpy_(&i__3, &d__1, &a[j + 1 + j * a_dim1], &c__1, &
								x[j + 1], &c__1);
						i__3 = *n - j;
						i__ = j + idamax_(&i__3, &x[j + 1], &c__1);
						xmax = (d__1 = x[i__], abs(d__1));
					}
				}
				/* L110: */
			}

		} else {

			/*           Solve A' * x = b */

			i__2 = jlast;
			i__1 = jinc;
			for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

				/*              Compute x(j) = b(j) - sum A(k,j)*x(k).
                                      k<>j */

				xj = (d__1 = x[j], abs(d__1));
				uscal = tscal;
				rec = 1. / max(xmax,1.);
				if (cnorm[j] > (bignum - xj) * rec) {

					/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

					rec *= .5;
					if (nounit) {
						tjjs = a[j + j * a_dim1] * tscal;
					} else {
						tjjs = tscal;
					}
					tjj = abs(tjjs);
					if (tjj > 1.) {

						/*                       Divide by A(j,j) when scaling x if A(j,j) > 1.

   Computing MIN */
						d__1 = 1., d__2 = rec * tjj;
						rec = min(d__1,d__2);
						uscal /= tjjs;
					}
					if (rec < 1.) {
						dscal_(n, &rec, &x[1], &c__1);
						*scale *= rec;
						xmax *= rec;
					}
				}

				sumj = 0.;
				if (uscal == 1.) {

					/*                 If the scaling needed for A in the dot product is 1,
                   call DDOT to perform the dot product. */

					if (upper) {
						i__3 = j - 1;
						sumj = ddot_(&i__3, &a[j * a_dim1 + 1], &c__1, &x[1],
								&c__1);
					} else if (j < *n) {
						i__3 = *n - j;
						sumj = ddot_(&i__3, &a[j + 1 + j * a_dim1], &c__1, &x[
																			  j + 1], &c__1);
					}
				} else {

					/*                 Otherwise, use in-line code for the dot product. */

					if (upper) {
						i__3 = j - 1;
						for (i__ = 1; i__ <= i__3; ++i__) {
							sumj += a[i__ + j * a_dim1] * uscal * x[i__];
							/* L120: */
						}
					} else if (j < *n) {
						i__3 = *n;
						for (i__ = j + 1; i__ <= i__3; ++i__) {
							sumj += a[i__ + j * a_dim1] * uscal * x[i__];
							/* L130: */
						}
					}
				}

				if (uscal == tscal) {

					/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                   was not used to scale the dotproduct. */

					x[j] -= sumj;
					xj = (d__1 = x[j], abs(d__1));
					if (nounit) {
						tjjs = a[j + j * a_dim1] * tscal;
					} else {
						tjjs = tscal;
						if (tscal == 1.) {
							goto L150;
						}
					}

					/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

					tjj = abs(tjjs);
					if (tjj > smlnum) {

						/*                       abs(A(j,j)) > SMLNUM: */

						if (tjj < 1.) {
							if (xj > tjj * bignum) {

								/*                             Scale X by 1/abs(x(j)). */

								rec = 1. / xj;
								dscal_(n, &rec, &x[1], &c__1);
								*scale *= rec;
								xmax *= rec;
							}
						}
						x[j] /= tjjs;
					} else if (tjj > 0.) {

						/*                       0 < abs(A(j,j)) <= SMLNUM: */

						if (xj > tjj * bignum) {

							/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

							rec = tjj * bignum / xj;
							dscal_(n, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
						x[j] /= tjjs;
					} else {

						/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                         scale = 0, and compute a solution to A'*x = 0. */

						i__3 = *n;
						for (i__ = 1; i__ <= i__3; ++i__) {
							x[i__] = 0.;
							/* L140: */
						}
						x[j] = 1.;
						*scale = 0.;
						xmax = 0.;
					}
					L150:
					;
				} else {

					/*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                   product has already been divided by 1/A(j,j). */

					x[j] = x[j] / tjjs - sumj;
				}
				/* Computing MAX */
				d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
				xmax = max(d__2,d__3);
				/* L160: */
			}
		}
		*scale /= tscal;
	}

	/*     Scale the column norms by 1/TSCAL for return. */

	if (tscal != 1.) {
		d__1 = 1. / tscal;
		dscal_(n, &d__1, &cnorm[1], &c__1);
	}

	return 0;

	/*     End of DLATRS */

} /* dlatrs_ */
/* Subroutine */ int dlassq_(int *n, double *x, int *incx, double *scale, double *sumsq)
/*  -- LAPACK auxiliary routine (version 3.1) --
   Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   November 2006


Purpose
=======

DLASSQ  returns the values  scl  and  smsq  such that

   ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,

where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
assumed to be non-negative and  scl  returns the value

   scl = max( scale, abs( x( i ) ) ).

scale and sumsq must be supplied in SCALE and SUMSQ and
scl and smsq are overwritten on SCALE and SUMSQ respectively.

The routine makes only one pass through the vector x.  */
{
	/* System generated locals */
	int i__1, i__2;
	double d__1;
	/* Local variables */
	int ix;
	double absxi;
	--x;
	/* Function Body */
	if (*n > 0) {
		i__1 = (*n - 1) * *incx + 1;
		i__2 = *incx;
		for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
			if (x[ix] != 0.) {
				absxi = (d__1 = x[ix], abs(d__1));
				if (*scale < absxi) {
					/* Computing 2nd power */
					d__1 = *scale / absxi;
					*sumsq = *sumsq * (d__1 * d__1) + 1;
					*scale = absxi;
				} else {
					/* Computing 2nd power */
					d__1 = absxi / *scale;
					*sumsq += d__1 * d__1;
				}
			}
		}
	}
	return 0;
} /* dlassq_ */
double dlange_(char *norm, int *m, int *n, double *a, int *lda, double *work)
/*  -- LAPACK auxiliary routine (version 3.1) --
   Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   November 2006


Purpose
=======

DLANGE  returns the value of the one norm,  or the Frobenius norm, or
the  infinity norm,  or the  element of  largest absolute value  of a
real matrix A.

Description
===========

DLANGE returns the value

   DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
            (
            ( norm1(A),         NORM = '1', 'O' or 'o'
            (
            ( normI(A),         NORM = 'I' or 'i'
            (
            ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

where  norm1  denotes the  one norm of a matrix (maximum column sum),
normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
normF  denotes the  Frobenius norm of a matrix (square root of sum of
squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm. */
{
	/* Table of constant values */
	int c__1 = 1;
	/* System generated locals */
	int a_dim1, a_offset, i__1, i__2;
	double ret_val, d__1, d__2, d__3;
	/* Local variables */
	int i__, j;
	double sum, scale;
	double value;
	char O[] = "O", M[] = "M", I[] = "I", F[] = "F", E[] = "E";
	a_dim1 = *lda;			a_offset = 1 + a_dim1;
	a -= a_offset;			--work;
	/* Function Body */
	if (min(*m,*n) == 0) {		value = 0.;	}
	else if (lsame_(norm, M)) {
		/*        Find max(abs(A(i,j))). */
		value = 0.;		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
				/* Computing MAX */
				d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
				value = max(d__2,d__3);
			}
		}
	} else if (lsame_(norm, O) || *(unsigned char *) norm == '1') {
		/*        Find norm1(A). */
		value = 0.;		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			sum = 0.;			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
				sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
			}
			value = max(value,sum);
		}
	} else if (lsame_(norm, I)) {
		/*        Find normI(A). */
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__) {		work[i__] = 0.;		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
				work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
			}
		}
		value = 0.;
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* Computing MAX */
			d__1 = value, d__2 = work[i__];
			value = max(d__1,d__2);
		}
	} else if (lsame_(norm, F) || lsame_(norm, E)) {
		/*        Find normF(A). */
		scale = 0.;		sum = 1.;		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			dlassq_(m, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
		}
		value = scale * sqrt(sum);
	}
	ret_val = value;
	return ret_val;
} /* dlange_ */
/* Subroutine */ int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
		double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc)
/*  Purpose
=======

DGEMM  performs one of the matrix-matrix operations

   C := alpha*op( A )*op( B ) + beta*C,

where  op( X ) is one of

   op( X ) = X   or   op( X ) = X',

alpha and beta are scalars, and A, B and C are matrices, with op( A )

an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */
{
	/* System generated locals */
	int i__1, i__2, i__3;

	/* Local variables */
	int info;
	bool nota, notb;
	double temp;
	int i, j, l, ncola;
	int nrowa, nrowb;
#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]
	char N[] = "N", C[] = "C", T[] = "T";
	nota = lsame_(transa, N);
	notb = lsame_(transb, N);
	if (nota) 	{		nrowa = *m;		ncola = *k;	}
	else 		{		nrowa = *k;		ncola = *m;	}
	if (notb) 	{		nrowb = *k;					}
	else 		{		nrowb = *n;					}
	/*     Test the input parameters. */
	info = 0;
	if (! nota && ! lsame_(transa, C) && ! lsame_(transa, T)) 		{	info = 1;	}
	else if (! notb && ! lsame_(transb, C) && ! lsame_(transb, T)) 	{	info = 2;	}
	else if (*m < 0) {	info = 3;	}
	else if (*n < 0) {	info = 4;	}
	else if (*k < 0) {	info = 5;	}
	else if (*lda < max(1,nrowa)) 	{	info = 8;	}
	else if (*ldb < max(1,nrowb)) 	{	info = 10;	}
	else if (*ldc < max(1,*m)) 		{	info = 13;	}
	if (info != 0) {
		//xerbla_("DGEMM ", &info);
		std::cout << "** On entry to DGEMM, parameter number " << info <<  "had an illegal value\n";
		return 0;
	}
	/*     Quick return if possible. */
	if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {	return 0;	}
	/*     And if  alpha.eq.zero. */
	if (*alpha == 0.) {
		if (*beta == 0.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
				i__2 = *m;
				for (i = 1; i <= *m; ++i) {	C(i,j) = 0.; }
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
				i__2 = *m;
				for (i = 1; i <= *m; ++i) {	C(i,j) = *beta * C(i,j);}
			}
		}
		return 0;
	}
	/*     Start the operations. */
	if (notb) {
		if (nota) {
			/*           Form  C := alpha*A*B + beta*C. */
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
				if (*beta == 0.) {
					i__2 = *m;
					for (i = 1; i <= *m; ++i) {	C(i,j) = 0.;}
				} else if (*beta != 1.) {
					i__2 = *m;
					for (i = 1; i <= *m; ++i) {	C(i,j) = *beta * C(i,j);}
				}
				i__2 = *k;
				for (l = 1; l <= *k; ++l) {
					if (B(l,j) != 0.) {
						temp = *alpha * B(l,j);
						i__3 = *m;
						for (i = 1; i <= *m; ++i) {	C(i,j) += temp * A(i,l);}
					}
				}
			}
		} else {
			/*           Form  C := alpha*A'*B + beta*C */
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
				i__2 = *m;
				for (i = 1; i <= *m; ++i) {
					temp = 0.;
					i__3 = *k;
					for (l = 1; l <= *k; ++l) {	temp += A(l,i) * B(l,j);	}
					if (*beta == 0.) {	C(i,j) = *alpha * temp;	}
					else {	C(i,j) = *alpha * temp + *beta * C(i,j);	}
				}
			}
		}
	} else {
		if (nota) {
			/*           Form  C := alpha*A*B' + beta*C */
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
				if (*beta == 0.) {
					i__2 = *m;
					for (i = 1; i <= *m; ++i) {	C(i,j) = 0.;}
				} else if (*beta != 1.) {
					i__2 = *m;
					for (i = 1; i <= *m; ++i) {	C(i,j) = *beta * C(i,j);	}
				}
				i__2 = *k;
				for (l = 1; l <= *k; ++l) {
					if (B(j,l) != 0.) {
						temp = *alpha * B(j,l);
						i__3 = *m;
						for (i = 1; i <= *m; ++i) {
							C(i,j) += temp * A(i,l);
						}
					}
				}
			}
		} else {
			/*           Form  C := alpha*A'*B' + beta*C */
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {
				i__2 = *m;
				for (i = 1; i <= *m; ++i) {
					temp = 0.;
					i__3 = *k;
					for (l = 1; l <= *k; ++l) {
						temp += A(l,i) * B(j,l);
					}
					if (*beta == 0.) {
						C(i,j) = *alpha * temp;
					} else {
						C(i,j) = *alpha * temp + *beta * C(i,j);
					}
				}
			}
		}
	}
	return 0;
} /* dgemm_ */
/* Subroutine */ int dtrsm_(char *side, char *uplo, char *transa, char *diag,
		int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb)
/*  Purpose
=======

DTRSM  solves one of the matrix equations

   op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,

where alpha is a scalar, X and B are m by n matrices, A is a unit, or

non-unit,  upper or lower triangular matrix  and  op( A )  is one  of


   op( A ) = A   or   op( A ) = A'.

The matrix X is overwritten on B. */
{
	/* System generated locals */
	int i__1, i__2, i__3;

	/* Local variables */
	int info;
	double temp;
	int i, j, k;
	bool lside;
	int nrowa;
	bool upper;
	bool nounit;

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
	char L[] = "L", N[] = "N", U[] = "U", R[] = "R", T[] = "T", C[] = "C";
	lside = lsame_(side, L);
	if (lside) {
		nrowa = *m;
	} else {
		nrowa = *n;
	}
	nounit = lsame_(diag, N);
	upper = lsame_(uplo, U);

	info = 0;
	if (! lside && ! lsame_(side, R)) {
		info = 1;
	} else if (! upper && ! lsame_(uplo, L)) {
		info = 2;
	} else if (! lsame_(transa, N) && ! lsame_(transa, T)
			&& ! lsame_(transa, C)) {
		info = 3;
	} else if (! lsame_(diag, U) && ! lsame_(diag, N)) {
		info = 4;
	} else if (*m < 0) {
		info = 5;
	} else if (*n < 0) {
		info = 6;
	} else if (*lda < max(1,nrowa)) {
		info = 9;
	} else if (*ldb < max(1,*m)) {
		info = 11;
	}
	if (info != 0) {
		//xerbla_("DTRSM ", &info);
		std::cout << "** On entry to DTRSM, parameter number " << info <<  "had an illegal value\n";
		return 0;
	}

	/*     Quick return if possible. */

	if (*n == 0) {
		return 0;
	}

	/*     And when  alpha.eq.zero. */

	if (*alpha == 0.) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
				B(i,j) = 0.;
				/* L10: */
			}
			/* L20: */
		}
		return 0;
	}

	/*     Start the operations. */

	if (lside) {
		if (lsame_(transa, N)) {

			/*           Form  B := alpha*inv( A )*B. */

			if (upper) {
				i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					if (*alpha != 1.) {
						i__2 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,j) = *alpha * B(i,j);
							/* L30: */
						}
					}
					for (k = *m; k >= 1; --k) {
						if (B(k,j) != 0.) {
							if (nounit) {
								B(k,j) /= A(k,k);
							}
							i__2 = k - 1;
							for (i = 1; i <= k-1; ++i) {
								B(i,j) -= B(k,j) * A(i,k);
								/* L40: */
							}
						}
						/* L50: */
					}
					/* L60: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					if (*alpha != 1.) {
						i__2 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,j) = *alpha * B(i,j);
							/* L70: */
						}
					}
					i__2 = *m;
					for (k = 1; k <= *m; ++k) {
						if (B(k,j) != 0.) {
							if (nounit) {
								B(k,j) /= A(k,k);
							}
							i__3 = *m;
							for (i = k + 1; i <= *m; ++i) {
								B(i,j) -= B(k,j) * A(i,k);
								/* L80: */
							}
						}
						/* L90: */
					}
					/* L100: */
				}
			}
		} else {

			/*           Form  B := alpha*inv( A' )*B. */

			if (upper) {
				i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					i__2 = *m;
					for (i = 1; i <= *m; ++i) {
						temp = *alpha * B(i,j);
						i__3 = i - 1;
						for (k = 1; k <= i-1; ++k) {
							temp -= A(k,i) * B(k,j);
							/* L110: */
						}
						if (nounit) {
							temp /= A(i,i);
						}
						B(i,j) = temp;
						/* L120: */
					}
					/* L130: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					for (i = *m; i >= 1; --i) {
						temp = *alpha * B(i,j);
						i__2 = *m;
						for (k = i + 1; k <= *m; ++k) {
							temp -= A(k,i) * B(k,j);
							/* L140: */
						}
						if (nounit) {
							temp /= A(i,i);
						}
						B(i,j) = temp;
						/* L150: */
					}
					/* L160: */
				}
			}
		}
	} else {
		if (lsame_(transa, N)) {

			/*           Form  B := alpha*B*inv( A ). */

			if (upper) {
				i__1 = *n;
				for (j = 1; j <= *n; ++j) {
					if (*alpha != 1.) {
						i__2 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,j) = *alpha * B(i,j);
							/* L170: */
						}
					}
					i__2 = j - 1;
					for (k = 1; k <= j-1; ++k) {
						if (A(k,j) != 0.) {
							i__3 = *m;
							for (i = 1; i <= *m; ++i) {
								B(i,j) -= A(k,j) * B(i,k);
								/* L180: */
							}
						}
						/* L190: */
					}
					if (nounit) {
						temp = 1. / A(j,j);
						i__2 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,j) = temp * B(i,j);
							/* L200: */
						}
					}
					/* L210: */
				}
			} else {
				for (j = *n; j >= 1; --j) {
					if (*alpha != 1.) {
						i__1 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,j) = *alpha * B(i,j);
							/* L220: */
						}
					}
					i__1 = *n;
					for (k = j + 1; k <= *n; ++k) {
						if (A(k,j) != 0.) {
							i__2 = *m;
							for (i = 1; i <= *m; ++i) {
								B(i,j) -= A(k,j) * B(i,k);
								/* L230: */
							}
						}
						/* L240: */
					}
					if (nounit) {
						temp = 1. / A(j,j);
						i__1 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,j) = temp * B(i,j);
							/* L250: */
						}
					}
					/* L260: */
				}
			}
		} else {

			/*           Form  B := alpha*B*inv( A' ). */

			if (upper) {
				for (k = *n; k >= 1; --k) {
					if (nounit) {
						temp = 1. / A(k,k);
						i__1 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,k) = temp * B(i,k);
							/* L270: */
						}
					}
					i__1 = k - 1;
					for (j = 1; j <= k-1; ++j) {
						if (A(j,k) != 0.) {
							temp = A(j,k);
							i__2 = *m;
							for (i = 1; i <= *m; ++i) {
								B(i,j) -= temp * B(i,k);
								/* L280: */
							}
						}
						/* L290: */
					}
					if (*alpha != 1.) {
						i__1 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,k) = *alpha * B(i,k);
							/* L300: */
						}
					}
					/* L310: */
				}
			} else {
				i__1 = *n;
				for (k = 1; k <= *n; ++k) {
					if (nounit) {
						temp = 1. / A(k,k);
						i__2 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,k) = temp * B(i,k);
							/* L320: */
						}
					}
					i__2 = *n;
					for (j = k + 1; j <= *n; ++j) {
						if (A(j,k) != 0.) {
							temp = A(j,k);
							i__3 = *m;
							for (i = 1; i <= *m; ++i) {
								B(i,j) -= temp * B(i,k);
								/* L330: */
							}
						}
						/* L340: */
					}
					if (*alpha != 1.) {
						i__2 = *m;
						for (i = 1; i <= *m; ++i) {
							B(i,k) = *alpha * B(i,k);
							/* L350: */
						}
					}
					/* L360: */
				}
			}
		}
	}

	return 0;

	/*     End of DTRSM . */

} /* dtrsm_ */
/* Subroutine */ int dger_(int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda)
/*  Purpose
    =======

    DGER   performs the rank 1 operation

       A := alpha*x*y' + A,

    where alpha is a scalar, x is an m element vector, y is an n element

    vector and A is an m by n matrix.*/
{
	/* System generated locals */
	int i__1, i__2;

	/* Local variables */
	int info;
	double temp;
	int i, j, ix, jy, kx;
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

	info = 0;
	if (*m < 0) {
		info = 1;
	} else if (*n < 0) {
		info = 2;
	} else if (*incx == 0) {
		info = 5;
	} else if (*incy == 0) {
		info = 7;
	} else if (*lda < max(1,*m)) {
		info = 9;
	}
	if (info != 0) {
		//xerbla_("DGER  ", &info);
		std::cout << "** On entry to DGER, parameter number " << info <<  "had an illegal value\n";
		return 0;
	}

	/*     Quick return if possible. */

	if (*m == 0 || *n == 0 || *alpha == 0.) {
		return 0;
	}

	/*     Start the operations. In this version the elements of A are
       accessed sequentially with one pass through A. */

	if (*incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (*n - 1) * *incy;
	}
	if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
			if (Y(jy) != 0.) {
				temp = *alpha * Y(jy);
				i__2 = *m;
				for (i = 1; i <= *m; ++i) {
					A(i,j) += X(i) * temp;
					/* L10: */
				}
			}
			jy += *incy;
			/* L20: */
		}
	} else {
		if (*incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (*m - 1) * *incx;
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
			if (Y(jy) != 0.) {
				temp = *alpha * Y(jy);
				ix = kx;
				i__2 = *m;
				for (i = 1; i <= *m; ++i) {
					A(i,j) += X(ix) * temp;
					ix += *incx;
					/* L30: */
				}
			}
			jy += *incy;
			/* L40: */
		}
	}

	return 0;

	/*     End of DGER  . */

} /* dger_ */
/* Subroutine */ int dswap_(int *n, double *dx, int *incx, double *dy, int *incy)
/*     interchanges two vectors.
       uses unrolled loops for increments equal one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*)*/
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i, m;
	double dtemp;
	int ix, iy, mp1;
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


	if (*n <= 0) {
		return 0;
	}
	if (*incx == 1 && *incy == 1) {
		goto L20;
	}

	/*       code for unequal increments or equal increments not equal
           to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
		ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
		iy = (-(*n) + 1) * *incy + 1;
	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
		dtemp = DX(ix);
		DX(ix) = DY(iy);
		DY(iy) = dtemp;
		ix += *incx;
		iy += *incy;
		/* L10: */
	}
	return 0;

	/*       code for both increments equal to 1


         clean-up loop */

	L20:
	m = *n % 3;
	if (m == 0) {
		goto L40;
	}
	i__1 = m;
	for (i = 1; i <= m; ++i) {
		dtemp = DX(i);
		DX(i) = DY(i);
		DY(i) = dtemp;
		/* L30: */
	}
	if (*n < 3) {
		return 0;
	}
	L40:
	mp1 = m + 1;
	i__1 = *n;
	for (i = mp1; i <= *n; i += 3) {
		dtemp = DX(i);
		DX(i) = DY(i);
		DY(i) = dtemp;
		dtemp = DX(i + 1);
		DX(i + 1) = DY(i + 1);
		DY(i + 1) = dtemp;
		dtemp = DX(i + 2);
		DX(i + 2) = DY(i + 2);
		DY(i + 2) = dtemp;
		/* L50: */
	}
	return 0;
} /* dswap_ */
/* Subroutine */ int dgetf2_(int *m, int *n, double *a, int *lda, int *ipiv, int *info)
/*  -- LAPACK routine (version 3.1) --
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
       November 2006


    Purpose
    =======

    DGETF2 computes an LU factorization of a general m-by-n matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 2 BLAS version of the algorithm.*/
{
	/* Table of constant values */
	int c__1 = 1;
	double c_b8 = -1.;

	/* System generated locals */
	int a_dim1, a_offset, i__1, i__2, i__3;
	double d__1;
	/* Local variables */
	int i__, j, jp;
	double sfmin;

	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipiv;

	/* Function Body */
	*info = 0;
	if (*m < 0) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	} else if (*lda < max(1,*m)) {
		*info = -4;
	}
	if (*info != 0) {
		i__1 = -(*info);
		//xerbla_("DGETF2", &i__1);
		std::cout << "** On entry to DGETF2, parameter number " << i__1 <<  "had an illegal value\n";
		return 0;
	}

	/*     Quick return if possible */

	if (*m == 0 || *n == 0) {
		return 0;
	}

	/*     Compute machine safe minimum */

	sfmin = 1.17549E-38;

	i__1 = min(*m,*n);
	for (j = 1; j <= i__1; ++j) {

		/*        Find pivot and test for singularity. */

		i__2 = *m - j + 1;
		jp = j - 1 + idamax_(&i__2, &a[j + j * a_dim1], &c__1);
		ipiv[j] = jp;
		if (a[jp + j * a_dim1] != 0.) {

			/*           Apply the interchange to columns 1:N. */

			if (jp != j) {
				dswap_(n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
			}

			/*           Compute elements J+1:M of J-th column. */

			if (j < *m) {
				if ((d__1 = a[j + j * a_dim1], abs(d__1)) >= sfmin) {
					i__2 = *m - j;
					d__1 = 1. / a[j + j * a_dim1];
					dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
				} else {
					i__2 = *m - j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						a[j + i__ + j * a_dim1] /= a[j + j * a_dim1];
						/* L20: */
					}
				}
			}

		} else if (*info == 0) {

			*info = j;
		}

		if (j < min(*m,*n)) {

			/*           Update trailing submatrix. */

			i__2 = *m - j;
			i__3 = *n - j;
			dger_(&i__2, &i__3, &c_b8, &a[j + 1 + j * a_dim1], &c__1, &a[j + (
					j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * a_dim1], lda);
		}
		/* L10: */
	}
	return 0;

	/*     End of DGETF2 */

} /* dgetf2_ */
void s_copy(register char *a, register char *b, int la, int lb)
{
	register char *aend, *bend;
	aend = a + la;
	if(la <= lb)
#ifndef NO_OVERWRITE
		if (a <= b || a >= b + la)
#endif
			while(a < aend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else
			for(b += la; a < aend; )
				*--aend = *--b;
#endif

	else {
		bend = b + lb;
#ifndef NO_OVERWRITE
		if (a <= b || a >= bend)
#endif
			while(b < bend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else {
			a += lb;
			while(b < bend)
				*--a = *--bend;
			a += lb;
		}
#endif
		while(a < aend)
			*a++ = ' ';
	}
}
int s_cmp(char *a0, char *b0, int la, int lb)
{
	register unsigned char *a, *aend, *b, *bend;
	a = (unsigned char *)a0;
	b = (unsigned char *)b0;
	aend = a + la;
	bend = b + lb;

	if(la <= lb)
	{
		while(a < aend)
			if(*a != *b)
				return( *a - *b );
			else
			{ ++a; ++b; }

		while(b < bend)
			if(*b != ' ')
				return( ' ' - *b );
			else	++b;
	}

	else
	{
		while(b < bend)
			if(*a == *b)
			{ ++a; ++b; }
			else
				return( *a - *b );
		while(a < aend)
			if(*a != ' ')
				return(*a - ' ');
			else	++a;
	}
	return(0);
}
int ieeeck_(int *ispec, float *zero, float *one)
/*  -- LAPACK auxiliary routine (version 3.1) --
   Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   November 2006


Purpose
=======

IEEECK is called from the ILAENV to verify that Infinity and
possibly NaN arithmetic is safe (i.e. will not trap). */
{
	int ret_val;
	/* Local variables */
	static float nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, negzro,
	newzro;


	ret_val = 1;

	posinf = *one / *zero;
	if (posinf <= *one) {
		ret_val = 0;
		return ret_val;
	}

	neginf = -(*one) / *zero;
	if (neginf >= *zero) {
		ret_val = 0;
		return ret_val;
	}

	negzro = *one / (neginf + *one);
	if (negzro != *zero) {
		ret_val = 0;
		return ret_val;
	}

	neginf = *one / negzro;
	if (neginf >= *zero) {
		ret_val = 0;
		return ret_val;
	}

	newzro = negzro + *zero;
	if (newzro != *zero) {
		ret_val = 0;
		return ret_val;
	}

	posinf = *one / newzro;
	if (posinf <= *one) {
		ret_val = 0;
		return ret_val;
	}

	neginf *= posinf;
	if (neginf >= *zero) {
		ret_val = 0;
		return ret_val;
	}

	posinf *= posinf;
	if (posinf <= *one) {
		ret_val = 0;
		return ret_val;
	}
	/*     Return if we were only asked to check infinity arithmetic */
	if (*ispec == 0) {
		return ret_val;
	}

	nan1 = posinf + neginf;

	nan2 = posinf / neginf;

	nan3 = posinf / posinf;

	nan4 = posinf * *zero;

	nan5 = neginf * negzro;

	nan6 = nan5 * 0.f;

	if (nan1 == nan1) {
		ret_val = 0;
		return ret_val;
	}

	if (nan2 == nan2) {
		ret_val = 0;
		return ret_val;
	}

	if (nan3 == nan3) {
		ret_val = 0;
		return ret_val;
	}

	if (nan4 == nan4) {
		ret_val = 0;
		return ret_val;
	}

	if (nan5 == nan5) {
		ret_val = 0;
		return ret_val;
	}

	if (nan6 == nan6) {
		ret_val = 0;
		return ret_val;
	}

	return ret_val;
} /* ieeeck_ */
int iparmq_(int *ispec, char *name__, char *opts, int *n, int *ilo, int *ihi, int *lwork){
	/* System generated locals */
	int ret_val, i__1, i__2;
	float r__1;
	/* Local variables */
	static int nh, ns;

	if (*ispec == 15 || *ispec == 13 || *ispec == 16) {

		/*        ==== Set the number simultaneous shifts ==== */

		nh = *ihi - *ilo + 1;
		ns = 2;
		if (nh >= 30) {
			ns = 4;
		}
		if (nh >= 60) {
			ns = 10;
		}
		if (nh >= 150) {
			/* Computing MAX */
			r__1 = log((float) nh) / log(2.f);
			i__1 = 10, i__2 = nh / (int)(r__1 >= 0 ? floor(r__1 + .5) : -floor(.5 - r__1));
			ns = max(i__1,i__2);
		}
		if (nh >= 590) {
			ns = 64;
		}
		if (nh >= 3000) {
			ns = 128;
		}
		if (nh >= 6000) {
			ns = 256;
		}
		/* Computing MAX */
		i__1 = 2, i__2 = ns - ns % 2;
		ns = max(i__1,i__2);
	}

	if (*ispec == 12) {


		/*        ===== Matrices of order smaller than NMIN get sent
          .     to xLAHQR, the classic double shift algorithm.
          .     This must be at least 11. ==== */

		ret_val = 75;

	} else if (*ispec == 14) {

		/*        ==== INIBL: skip a multi-shift qr iteration and
          .    whenever aggressive early deflation finds
          .    at least (NIBBLE*(window size)/100) deflations. ==== */

		ret_val = 14;

	} else if (*ispec == 15) {

		/*        ==== NSHFTS: The number of simultaneous shifts ===== */

		ret_val = ns;

	} else if (*ispec == 13) {

		/*        ==== NW: deflation window size.  ==== */

		if (nh <= 500) {
			ret_val = ns;
		} else {
			ret_val = ns * 3 / 2;
		}

	} else if (*ispec == 16) {

		/*        ==== IACC22: Whether to accumulate reflections
          .     before updating the far-from-diagonal elements
          .     and whether to use 2-by-2 block structure while
          .     doing it.  A small amount of work could be saved
          .     by making this choice dependent also upon the
          .     NH=IHI-ILO+1. */

		ret_val = 0;
		if (ns >= 14) {
			ret_val = 1;
		}
		if (ns >= 14) {
			ret_val = 2;
		}

	} else {
		/*        ===== invalid value of ispec ===== */
		ret_val = -1;

	}

	/*     ==== End of IPARMQ ==== */

	return ret_val;
} /* iparmq_ */
/* Table of constant values */

static int c__0 = 0;
static float c_b163 = 0.f;
static float c_b164 = 1.f;
static int c__1 = 1;

int ilaenv_(int *ispec, char *name__, char *opts, int *n1, int *n2, int *n3, int *n4, int name_len, int opts_len)
/*  -- LAPACK auxiliary routine (version 3.1.1) --
   Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   January 2007


Purpose
=======

ILAENV is called from the LAPACK routines to choose problem-dependent
parameters for the local environment.  See ISPEC for a description of
the parameters.

ILAENV returns an INTEGER
if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.

This version provides a set of parameters which should give good,
but not optimal, performance on many of the currently available
computers.  Users are encouraged to modify this subroutine to set
the tuning parameters for their particular machine using the option
and problem size information in the arguments.

This routine will not function correctly if it is converted to all
lower case.  Converting it to all upper case is allowed. */
{
	/* System generated locals */
	int ret_val;

	/* Builtin functions
       integer s_cmp(char *, char *, int, int);

       /* Local variables */
	static int i__;
	static char c1[1], c2[2], c3[3], c4[2];
	static int ic, nb, iz, nx;
	static bool cname;
	static int nbmin;
	static bool sname;
	static char subnam[6];

	char TRF[] = "TRF", QRF[] = "QRF", RQF[] = "RQF", LQF[] = "LQF", QLF[] = "QLF";
	char HRD[] = "HRD", BRD[] = "BRD", TRI[] = "TRI", TRD[] = "TRD", GST[] = "GST";
	char UUM[] = "UUM", EBZ[] = "EBZ";
	char OR[] = "OR", PO[] = "PO", SY[] = "SY", HE[] = "HE", GE[] = "GE", QR[] = "QR";
	char RQ[] = "RQ", LQ[] = "LQ", QL[] = "QL", HR[] = "HR", TR[] = "TR", BR[] = "BR";
	char UN[] = "UN", GB[] = "GB", PB[] = "PB", LA[] = "LA", ST[] = "ST";

	switch (*ispec) {
	case 1:  goto L10;
	case 2:  goto L10;
	case 3:  goto L10;
	case 4:  goto L80;
	case 5:  goto L90;
	case 6:  goto L100;
	case 7:  goto L110;
	case 8:  goto L120;
	case 9:  goto L130;
	case 10:  goto L140;
	case 11:  goto L150;
	case 12:  goto L160;
	case 13:  goto L160;
	case 14:  goto L160;
	case 15:  goto L160;
	case 16:  goto L160;
	}

	/*     Invalid value for ISPEC */

	ret_val = -1;
	return ret_val;

	L10:

	/*     Convert NAME to upper case if the first character is lower case. */

	ret_val = 1;
	s_copy(subnam, name__, (int)6, name_len);
	ic = *(unsigned char *)subnam;
	iz = 'Z';
	if (iz == 90 || iz == 122) {

		/*        ASCII character set */

		if (ic >= 97 && ic <= 122) {
			*(unsigned char *)subnam = (char) (ic - 32);
			for (i__ = 2; i__ <= 6; ++i__) {
				ic = *(unsigned char *)&subnam[i__ - 1];
				if (ic >= 97 && ic <= 122) {
					*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
				}
				/* L20: */
			}
		}

	} else if (iz == 233 || iz == 169) {

		/*        EBCDIC character set */

		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 &&
				ic <= 169) {
			*(unsigned char *)subnam = (char) (ic + 64);
			for (i__ = 2; i__ <= 6; ++i__) {
				ic = *(unsigned char *)&subnam[i__ - 1];
				if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >=
						162 && ic <= 169) {
					*(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
				}
				/* L30: */
			}
		}

	} else if (iz == 218 || iz == 250) {

		/*        Prime machines:  ASCII+128 */

		if (ic >= 225 && ic <= 250) {
			*(unsigned char *)subnam = (char) (ic - 32);
			for (i__ = 2; i__ <= 6; ++i__) {
				ic = *(unsigned char *)&subnam[i__ - 1];
				if (ic >= 225 && ic <= 250) {
					*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
				}
				/* L40: */
			}
		}
	}

	*(unsigned char *)c1 = *(unsigned char *)subnam;
	sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
	cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
	if (! (cname || sname)) {
		return ret_val;
	}
	s_copy(c2, subnam + 1, (int)2, (int)2);
	s_copy(c3, subnam + 3, (int)3, (int)3);
	s_copy(c4, c3 + 1, (int)2, (int)2);

	switch (*ispec) {
	case 1:  goto L50;
	case 2:  goto L60;
	case 3:  goto L70;
	}

	L50:

	/*     ISPEC = 1:  block size

       In these examples, separate code is provided for setting NB for
       real and complex.  We assume that NB will take the same value in
       single or double precision. */
	nb = 1;

	if (s_cmp(c2, GE, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRF, (int)3, (int)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		} else if (s_cmp(c3, QRF, (int)3, (int)3) == 0 || s_cmp(c3,
				RQF, (int)3, (int)3) == 0 || s_cmp(c3, LQF, (int)3, (int)3) == 0 || s_cmp(c3, QLF, (int)3, (int)3)
				== 0) {
			if (sname) {
				nb = 32;
			} else {
				nb = 32;
			}
		} else if (s_cmp(c3, HRD, (int)3, (int)3) == 0) {
			if (sname) {
				nb = 32;
			} else {
				nb = 32;
			}
		} else if (s_cmp(c3, BRD, (int)3, (int)3) == 0) {
			if (sname) {
				nb = 32;
			} else {
				nb = 32;
			}
		} else if (s_cmp(c3, TRI, (int)3, (int)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		}
	} else if (s_cmp(c2, PO, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRF, (int)3, (int)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		}
	} else if (s_cmp(c2, SY, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRF, (int)3, (int)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		} else if (sname && s_cmp(c3, TRD, (int)3, (int)3) == 0) {
			nb = 32;
		} else if (sname && s_cmp(c3, GST, (int)3, (int)3) == 0) {
			nb = 64;
		}
	} else if (cname && s_cmp(c2, HE, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRF, (int)3, (int)3) == 0) {
			nb = 64;
		} else if (s_cmp(c3, TRD, (int)3, (int)3) == 0) {
			nb = 32;
		} else if (s_cmp(c3, GST, (int)3, (int)3) == 0) {
			nb = 64;
		}
	} else if (sname && s_cmp(c2, OR, (int)2, (int)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nb = 32;
			}
		} else if (*(unsigned char *)c3 == 'M') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nb = 32;
			}
		}
	} else if (cname && s_cmp(c2, UN, (int)2, (int)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nb = 32;
			}
		} else if (*(unsigned char *)c3 == 'M') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nb = 32;
			}
		}
	} else if (s_cmp(c2, GB, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRF, (int)3, (int)3) == 0) {
			if (sname) {
				if (*n4 <= 64) {
					nb = 1;
				} else {
					nb = 32;
				}
			} else {
				if (*n4 <= 64) {
					nb = 1;
				} else {
					nb = 32;
				}
			}
		}
	} else if (s_cmp(c2, PB, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRF, (int)3, (int)3) == 0) {
			if (sname) {
				if (*n2 <= 64) {
					nb = 1;
				} else {
					nb = 32;
				}
			} else {
				if (*n2 <= 64) {
					nb = 1;
				} else {
					nb = 32;
				}
			}
		}
	} else if (s_cmp(c2, TR, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRI, (int)3, (int)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		}
	} else if (s_cmp(c2, LA, (int)2, (int)2) == 0) {
		if (s_cmp(c3, UUM, (int)3, (int)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		}
	} else if (sname && s_cmp(c2, ST, (int)2, (int)2) == 0) {
		if (s_cmp(c3, EBZ, (int)3, (int)3) == 0) {
			nb = 1;
		}
	}
	ret_val = nb;
	return ret_val;

	L60:

	/*     ISPEC = 2:  minimum block size */

	nbmin = 2;
	if (s_cmp(c2, GE, (int)2, (int)2) == 0) {
		if (s_cmp(c3, QRF, (int)3, (int)3) == 0 || s_cmp(c3, RQF, (
				int)3, (int)3) == 0 || s_cmp(c3, LQF, (int)3, (
						int)3) == 0 || s_cmp(c3, QLF, (int)3, (int)3) == 0)
		{
			if (sname) {
				nbmin = 2;
			} else {
				nbmin = 2;
			}
		} else if (s_cmp(c3, HRD, (int)3, (int)3) == 0) {
			if (sname) {
				nbmin = 2;
			} else {
				nbmin = 2;
			}
		} else if (s_cmp(c3, BRD, (int)3, (int)3) == 0) {
			if (sname) {
				nbmin = 2;
			} else {
				nbmin = 2;
			}
		} else if (s_cmp(c3, TRI, (int)3, (int)3) == 0) {
			if (sname) {
				nbmin = 2;
			} else {
				nbmin = 2;
			}
		}
	} else if (s_cmp(c2, SY, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRF, (int)3, (int)3) == 0) {
			if (sname) {
				nbmin = 8;
			} else {
				nbmin = 8;
			}
		} else if (sname && s_cmp(c3, TRD, (int)3, (int)3) == 0) {
			nbmin = 2;
		}
	} else if (cname && s_cmp(c2, HE, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRD, (int)3, (int)3) == 0) {
			nbmin = 2;
		}
	} else if (sname && s_cmp(c2, OR, (int)2, (int)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nbmin = 2;
			}
		} else if (*(unsigned char *)c3 == 'M') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nbmin = 2;
			}
		}
	} else if (cname && s_cmp(c2, UN, (int)2, (int)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nbmin = 2;
			}
		} else if (*(unsigned char *)c3 == 'M') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nbmin = 2;
			}
		}
	}
	ret_val = nbmin;
	return ret_val;

	L70:

	/*     ISPEC = 3:  crossover point */

	nx = 0;
	if (s_cmp(c2, GE, (int)2, (int)2) == 0) {
		if (s_cmp(c3, QRF, (int)3, (int)3) == 0 || s_cmp(c3, RQF, (
				int)3, (int)3) == 0 || s_cmp(c3, LQF, (int)3, (
						int)3) == 0 || s_cmp(c3, QLF, (int)3, (int)3) == 0)
		{
			if (sname) {
				nx = 128;
			} else {
				nx = 128;
			}
		} else if (s_cmp(c3, HRD, (int)3, (int)3) == 0) {
			if (sname) {
				nx = 128;
			} else {
				nx = 128;
			}
		} else if (s_cmp(c3, BRD, (int)3, (int)3) == 0) {
			if (sname) {
				nx = 128;
			} else {
				nx = 128;
			}
		}
	} else if (s_cmp(c2, SY, (int)2, (int)2) == 0) {
		if (sname && s_cmp(c3, TRD, (int)3, (int)3) == 0) {
			nx = 32;
		}
	} else if (cname && s_cmp(c2, HE, (int)2, (int)2) == 0) {
		if (s_cmp(c3, TRD, (int)3, (int)3) == 0) {
			nx = 32;
		}
	} else if (sname && s_cmp(c2, OR, (int)2, (int)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nx = 128;
			}
		}
	} else if (cname && s_cmp(c2, UN, (int)2, (int)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, QR, (int)2, (int)2) == 0 || s_cmp(c4, RQ,
					(int)2, (int)2) == 0 || s_cmp(c4, LQ, (int)2, (
							int)2) == 0 || s_cmp(c4, QL, (int)2, (int)2) ==
									0 || s_cmp(c4, HR, (int)2, (int)2) == 0 || s_cmp(
											c4, TR, (int)2, (int)2) == 0 || s_cmp(c4, BR, (int)2, (int)2) == 0) {
				nx = 128;
			}
		}
	}
	ret_val = nx;
	return ret_val;

	L80:

	/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

	ret_val = 6;
	return ret_val;

	L90:

	/*     ISPEC = 5:  minimum column dimension (not used) */

	ret_val = 2;
	return ret_val;

	L100:

	/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

	ret_val = (int) ((float) min(*n1,*n2) * 1.6f);
	return ret_val;

	L110:

	/*     ISPEC = 7:  number of processors (not used) */

	ret_val = 1;
	return ret_val;

	L120:

	/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

	ret_val = 50;
	return ret_val;

	L130:

	/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
                   computation tree in the divide-and-conquer algorithm
                   (used by xGELSD and xGESDD) */

	ret_val = 25;
	return ret_val;

	L140:

	/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap

       ILAENV = 0 */
	ret_val = 1;
	if (ret_val == 1) {
		ret_val = ieeeck_(&c__0, &c_b163, &c_b164);
	}
	return ret_val;

	L150:

	/*     ISPEC = 11: infinity arithmetic can be trusted not to trap

       ILAENV = 0 */
	ret_val = 1;
	if (ret_val == 1) {
		ret_val = ieeeck_(&c__1, &c_b163, &c_b164);
	}
	return ret_val;

	L160:

	/*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. */
	//int iparmq_(ispec, name__, opts, n, ilo, ihi, lwork)
	//ret_val = iparmq_(ispec, name__, opts, n1, n2, n3, n4, name_len, opts_len);
	ret_val = iparmq_(ispec, name__, opts, n1, n2, n3, n4/*, name_len, opts_len*/);
	return ret_val;

	/*     End of ILAENV */

} /* ilaenv_ */
/* Subroutine */ int dlaswp_(int *n, double *a, int *lda, int
		*k1, int *k2, int *ipiv, int *incx)
/*  -- LAPACK auxiliary routine (version 3.1) --
   Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   November 2006


Purpose
=======

DLASWP performs a series of row interchanges on the matrix A.
One row interchange is initiated for each of rows K1 through K2 of A. */
{
	/* System generated locals */
	int a_dim1, a_offset, i__1, i__2, i__3, i__4;
	/* Local variables */
	static int i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
	static double temp;

	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipiv;

	/* Function Body */
	if (*incx > 0) {
		ix0 = *k1;
		i1 = *k1;
		i2 = *k2;
		inc = 1;
	} else if (*incx < 0) {
		ix0 = (1 - *k2) * *incx + 1;
		i1 = *k2;
		i2 = *k1;
		inc = -1;
	} else {
		return 0;
	}

	n32 = *n / 32 << 5;
	if (n32 != 0) {
		i__1 = n32;
		for (j = 1; j <= i__1; j += 32) {
			ix = ix0;
			i__2 = i2;
			i__3 = inc;
			for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3)
			{
				ip = ipiv[ix];
				if (ip != i__) {
					i__4 = j + 31;
					for (k = j; k <= i__4; ++k) {
						temp = a[i__ + k * a_dim1];
						a[i__ + k * a_dim1] = a[ip + k * a_dim1];
						a[ip + k * a_dim1] = temp;
						/* L10: */
					}
				}
				ix += *incx;
				/* L20: */
			}
			/* L30: */
		}
	}
	if (n32 != *n) {
		++n32;
		ix = ix0;
		i__1 = i2;
		i__3 = inc;
		for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
			ip = ipiv[ix];
			if (ip != i__) {
				i__2 = *n;
				for (k = n32; k <= i__2; ++k) {
					temp = a[i__ + k * a_dim1];
					a[i__ + k * a_dim1] = a[ip + k * a_dim1];
					a[ip + k * a_dim1] = temp;
					/* L40: */
				}
			}
			ix += *incx;
			/* L50: */
		}
	}

	return 0;

	/*     End of DLASWP */

} /* dlaswp_ */
/* Subroutine */ int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info)
/*  -- LAPACK routine (version 3.1) --
   Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
   November 2006


Purpose
=======

DGETRF computes an LU factorization of a general M-by-N matrix A
using partial pivoting with row interchanges.

The factorization has the form
   A = P * L * U
where P is a permutation matrix, L is lower triangular with unit
diagonal elements (lower trapezoidal if m > n), and U is upper
triangular (upper trapezoidal if m < n).

This is the right-looking Level 3 BLAS version of the algorithm. */
{
	/* Table of constant values */
	int c__1 = 1;
	int c_n1 = -1;
	double c_b16 = 1.;
	double c_b19 = -1.;

	/* System generated locals */
	int a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
	/* Local variables */
	int i__, j, jb, nb;
	int iinfo;
	char DGETRF[] = "DGETRF", SPACE[] = " ", Unit[] = "Unit";
	char Left[] = "Left", Lower[] = "Lower", No_transpose[] = "No transpose";

	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipiv;

	/* Function Body */
	*info = 0;
	if (*m < 0) 				{	*info = -1;	}
	else if (*n < 0) 			{	*info = -2;	}
	else if (*lda < max(1,*m))	{	*info = -4;	}
	if (*info != 0) {
		i__1 = -(*info);
		//xerbla_("DGETRF", &i__1);
		std::cout << "** On entry to DGETRF, parameter number " << i__1 <<  "had an illegal value\n";
		return 0;
	}
	/*     Quick return if possible */
	if (*m == 0 || *n == 0) {	return 0;	}
	/*     Determine the block size for this environment. */
	nb = ilaenv_(&c__1, DGETRF, SPACE, m, n, &c_n1, &c_n1, (int)6, (int)1);
	if (nb <= 1 || nb >= min(*m,*n)) {
		/*        Use unblocked code. */
		dgetf2_(m, n, &a[a_offset], lda, &ipiv[1], info);
	} else {
		/*        Use blocked code. */
		i__1 = min(*m,*n);		i__2 = nb;
		for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
			/* Computing MIN */
			i__3 = min(*m,*n) - j + 1;
			jb = min(i__3,nb);
			/*           Factor diagonal and subdiagonal blocks and test for exact singularity. */
			i__3 = *m - j + 1;
			dgetf2_(&i__3, &jb, &a[j + j * a_dim1], lda, &ipiv[j], &iinfo);
			/*           Adjust INFO and the pivot indices. */
			if (*info == 0 && iinfo > 0) {	*info = iinfo + j - 1;	}
			/* Computing MIN */
			i__4 = *m, i__5 = j + jb - 1;
			i__3 = min(i__4,i__5);
			for (i__ = j; i__ <= i__3; ++i__) {	ipiv[i__] = j - 1 + ipiv[i__];	}
			/*           Apply interchanges to columns 1:J-1. */
			i__3 = j - 1;			i__4 = j + jb - 1;
			dlaswp_(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);
			if (j + jb <= *n) {
				/*              Apply interchanges to columns J+JB:N. */
				i__3 = *n - j - jb + 1;
				i__4 = j + jb - 1;
				dlaswp_(&i__3, &a[(j + jb) * a_dim1 + 1], lda, &j, &i__4, &ipiv[1], &c__1);
				/*              Compute block row of U. */
				i__3 = *n - j - jb + 1;
				dtrsm_(Left, Lower, No_transpose, Unit, &jb, &i__3, &
						c_b16, &a[j + j * a_dim1], lda, &a[j + (j + jb) * a_dim1], lda);
				if (j + jb <= *m) {
					/*                 Update trailing submatrix. */
					i__3 = *m - j - jb + 1;					i__4 = *n - j - jb + 1;
					dgemm_(No_transpose, No_transpose, &i__3, &i__4, &jb,
							&c_b19, &a[j + jb + j * a_dim1], lda, &a[j + (j +
									jb) * a_dim1], lda, &c_b16, &a[j + jb + (j + jb) * a_dim1], lda);
				}
			}
		}
	}
	return 0;
} /* dgetrf_ */
/* Subroutine */ int dgecon_(char *norm, int *n, double *a, int *
		lda, double *anorm, double *rcond, double *work, int *
		iwork, int *info)
{
	/*  -- LAPACK routine (version 3.1) --
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
       November 2006

       Modified to call DLACN2 in place of DLACON, 5 Feb 03, SJH.


    Purpose
    =======

    DGECON estimates the reciprocal of the condition number of a general
    real matrix A, in either the 1-norm or the infinity-norm, using
    the LU factorization computed by DGETRF.

    An estimate is obtained for norm(inv(A)), and the reciprocal of the
    condition number is computed as
       RCOND = 1 / ( norm(A) * norm(inv(A)) ). */

	/* Table of constant values */
	int c__1 = 1;

	/* System generated locals */
	int a_dim1, a_offset, i__1;
	double d__1;
	/* Local variables */
	double sl;
	int ix;
	double su;
	int kase, kase1;
	double scale;
	int isave[3];
	double ainvnm;
	bool onenrm;
	char normin[1];
	double smlnum;


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--work;
	--iwork;

	/* Function Body */
	*info = 0;
	char O[] = "O", I[] = "I";
	char Lower[] = "Lower", No_transpose[] = "No transpose", Unit[] = "Unit";
	char Upper[] = "Upper", Non_unit[] = "Non-unit", Transpose[] = "Transpose";
	onenrm = *(unsigned char *)norm == '1' || lsame_(norm, O);
	if (! onenrm && ! lsame_(norm, I)) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	} else if (*lda < max(1,*n)) {
		*info = -4;
	} else if (*anorm < 0.) {
		*info = -5;
	}
	if (*info != 0) {
		i__1 = -(*info);
		//xerbla_("DGECON", &i__1);
		std::cout << "** On entry to DGECON, parameter number " << i__1 <<  "had an illegal value\n";
		printf("** On entry to %6s, parameter number %2i had an illegal value\n", "DGECON", *info);
		return 0;
	}

	/*     Quick return if possible */

	*rcond = 0.;
	if (*n == 0) {
		*rcond = 1.;
		return 0;
	} else if (*anorm == 0.) {
		return 0;
	}

	smlnum = 1.17549E-38;

	/*     Estimate the norm of inv(A). */

	ainvnm = 0.;
	*(unsigned char *)normin = 'N';
	if (onenrm) {
		kase1 = 1;
	} else {
		kase1 = 2;
	}
	kase = 0;
	L10:
	dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
	if (kase != 0) {
		if (kase == kase1) {

			/*           Multiply by inv(L). */

			dlatrs_(Lower, No_transpose, Unit, normin, n, &a[a_offset],
					lda, &work[1], &sl, &work[(*n << 1) + 1], info);

			/*           Multiply by inv(U). */

			dlatrs_(Upper, No_transpose, Non_unit, normin, n, &a[
																 a_offset], lda, &work[1], &su, &work[*n * 3 + 1], info);
		} else {

			/*           Multiply by inv(U'). */

			dlatrs_(Upper, Transpose, Non_unit, normin, n, &a[a_offset],
					lda, &work[1], &su, &work[*n * 3 + 1], info);

			/*           Multiply by inv(L'). */

			dlatrs_(Lower, Transpose, Unit, normin, n, &a[a_offset],
					lda, &work[1], &sl, &work[(*n << 1) + 1], info);
		}

		/*        Divide X by 1/(SL*SU) if doing so will not cause overflow. */

		scale = sl * su;
		*(unsigned char *)normin = 'Y';
		if (scale != 1.) {
			ix = idamax_(n, &work[1], &c__1);
			if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.)
			{
				goto L20;
			}
			drscl_(n, &scale, &work[1], &c__1);
		}
		goto L10;
	}

	/*     Compute the estimate of the reciprocal condition number. */

	if (ainvnm != 0.) {
		*rcond = 1. / ainvnm / *anorm;
	}

	L20:
	return 0;

	/*     End of DGECON */

} /* dgecon_ */
