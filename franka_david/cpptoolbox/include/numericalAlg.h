/*
 * numericalAlg.h
 *
 *  Created on: Nov 10, 2015
 *      Author: fares
 */

#ifndef INCLUDE_NUMERICALALG_H_
#define INCLUDE_NUMERICALALG_H_

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include "mathUtility.hpp"

//!< LINPACK functions
int dchdc ( double a[], int lda, int p, double work[], int ipvt[], int job );
int dchdd ( double r[], int ldr, int p, double x[], double z[], int ldz,
		int nz, double y[], double rho[], double c[], double s[] );
void dchex ( double r[], int ldr, int p, int k, int l, double z[], int ldz,
		int nz, double c[], double s[], int job );
void dchud ( double r[], int ldr, int p, double x[], double z[], int ldz,
		int nz, double y[], double rho[], double c[], double s[] );
double dgbco ( double abd[], int lda, int n, int ml, int mu, int ipvt[],
		double z[] );
void dgbdi ( double abd[], int lda, int n, int ml, int mu, int ipvt[],
		double det[2] );
int dgbfa ( double abd[], int lda, int n, int ml, int mu, int ipvt[] );
void dgbsl ( double abd[], int lda, int n, int ml, int mu, int ipvt[],
		double b[], int job );
double dgeco ( double a[], int lda, int n, int ipvt[], double z[] );
void dgedi ( double a[], int lda, int n, int ipvt[], double det[],
		double work[], int job );
int dgefa ( double a[], int lda, int n, int ipvt[] );
void dgesl ( double a[], int lda, int n, int ipvt[], double b[], int job );
int dgtsl ( int n, double c[], double d[], double e[], double b[] );
double dpbco ( double abd[], int lda, int n, int m, double z[] );
void dpbdi ( double abd[], int lda, int n, int m, double det[] );
int dpbfa ( double abd[], int lda, int n, int m );
void dpbsl ( double abd[], int lda, int n, int m, double b[] );
double dpoco ( double a[], int lda, int n, double z[] );
void dpodi ( double a[], int lda, int n, double det[], int job );
int dpofa ( double a[], int lda, int n );
void dposl ( double a[], int lda, int n, double b[] );
double dppco ( double ap[], int n, double z[] );
void dppdi ( double ap[], int n, double det[2], int job );
int dppfa ( double ap[], int n );
void dppsl ( double ap[], int n, double b[] );
void dptsl ( int n, double d[], double e[], double b[] );
void dqrdc ( double a[], int lda, int n, int p, double qraux[], int jpvt[],
		double work[], int job );
int dqrsl ( double a[], int lda, int n, int k, double qraux[], double y[],
		double qy[], double qty[], double b[], double rsd[], double ab[], int job );
double dsico ( double a[], int lda, int n, int kpvt[], double z[] );
void dsidi ( double a[], int lda, int n, int kpvt[], double det[2],
		int inert[3], double work[], int job );
int dsifa ( double a[], int lda, int n, int kpvt[] );
void dsisl ( double a[], int lda, int n, int kpvt[], double b[] );
double dspco ( double ap[], int n, int kpvt[], double z[] );
void dspdi ( double ap[], int n, int kpvt[], double det[2], int inert[3],
		double work[], int job );
int dspfa ( double ap[], int n, int kpvt[] );
void dspsl ( double ap[], int n, int kpvt[], double b[] );
int dsvdc ( double a[], int lda, int m, int n, double s[], double e[],
		double u[], int ldu, double v[], int ldv, double work[], int job );
double dtrco ( double t[], int ldt, int n, double z[], int job );
int dtrdi ( double t[], int ldt, int n, double det[], int job );
int dtrsl ( double t[], int ldt, int n, double b[], int job );
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
void dscal ( int n, double sa, double x[], int incx );
void drot ( int n, double x[], int incx, double y[], int incy, double c,
  double s );
void drotg ( double *sa, double *sb, double *c, double *s );
void dswap ( int n, double x[], int incx, double y[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
double dnrm2 ( int n, double x[], int incx );
double dasum ( int n, double x[], int incx );
int idamax ( int n, double dx[], int incx );
void dcopy ( int n, double dx[], int incx, double dy[], int incy );

int i4_min ( int i1, int i2 );
int i4_max ( int i1, int i2 );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_sign ( double x );

//!< LANPACK functions
int dgecon_(char *norm, int *n, double *a, int *lda, double *anorm, double *rcond, double *work, int *iwork, int *info);
double dlange_(char *norm, int *m, int *n, double *a, int *lda, double *work);
int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

#endif /* INCLUDE_NUMERICALALG_H_ */
