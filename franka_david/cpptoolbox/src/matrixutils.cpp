//=============================================================================
//==============================================================================

//	file:	matrixutils.cpp

//	author:	Fares J. Abu-Dakka
//	date:	2006 - 2015

//==============================================================================
#include "matrix.hpp"

namespace math{


//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
ColumnVector cross(const ColumnVector & v1, const ColumnVector & v2){
	if (v1.size() != 3 || v2.size() != 3)
		throw Exception("error: CrossProduct(v1,v2): vectors should be of size 3.");

	ColumnVector v3(3);
	v3(1) = v1(2)*v2(3) - v2(2)*v1(3);
	v3(2) = v2(1)*v1(3) - v1(1)*v2(3);
	v3(3) = v1(1)*v2(2) - v2(1)*v1(2);
	return v3;
}
//-------------------------------------------------------------------------------------------
RowVector cross(const RowVector & v1, const RowVector & v2){
	if (v1.size() != 3 || v2.size() != 3)
		throw Exception("error: CrossProduct(v1,v2): vectors should be of size 3.");

	RowVector v3(3);
	const Real *a = v1.getData(), *b = v2.getData();	Real *c = v3.getData();
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
	return v3;
}
//-------------------------------------------------------------------------------------------
Real dot(const ColumnVector & v1, const ColumnVector & v2){
	if (v1.size() != 3 || v2.size() != 3)
		throw Exception("error: DotProduct(v1,v2): vectors should be of size 3.");

	Real sum = 0;
	int n = 3;
	const Real *a = v1.getData(), *b = v2.getData();
	while (n--) sum += *a++ * *b++;
	return sum;
}
//-------------------------------------------------------------------------------------------
Real dot(const RowVector & v1, const RowVector & v2){
	if (v1.size() != 3 || v2.size() != 3)
		throw Exception("error: DotProduct(v1,v2): vectors should be of size 3.");

	Real sum = 0;
	int n = 3;
	const Real *a = v1.getData(), *b = v2.getData();
	while (n--) sum += *a++ * *b++;
	return sum;
}
//-------------------------------------------------------------------------------------------
//! Distance Between Two Points in 3D
Real distance(const ColumnVector & v1, const ColumnVector & v2){
	if (v1.getnRows() != v2.getnRows())
		throw InvalidIndex();
	Real ss = ZERO;
	for (int i = 1; i <= v1.getnRows(); i++)
		ss += (v1(i) - v2(i))*(v1(i) - v2(i));

	return sqrt(ss);
}
//-------------------------------------------------------------------------------------------
ColumnVector UnitVector(const ColumnVector & v){ return v*(1 / sqrt(dot(v, v))); }
//-------------------------------------------------------------------------------------------
Matrix x_prod_matrix(const ColumnVector & v)
//!  @brief Cross product matrix
{
	if (v.getnRows() != 3)
		throw Exception("error: x_prod_matrix(v): vector must be of size 3.");
	Matrix S(3, 3); S = 0.0;
	S(1, 2) = -v(3); S(1, 3) = v(2);
	S(2, 1) = v(3);                 S(2, 3) = -v(1);
	S(3, 1) = -v(2); S(3, 2) = v(1);

	return (S);
}
//-------------------------------------------------------------------------------------------
Matrix x_prod_matrix(const RowVector & v)
//!  @brief Cross product matrix
{
	if (v.getnColumns() != 3)
		throw Exception("error: x_prod_matrix(v): vector must be of size 3.");

	return x_prod_matrix(ColumnVector(v(1), v(2), v(3)));
}
//-------------------------------------------------------------------------------------------
Matrix cholesky(const Matrix &A, ColumnVector &p){
	/*!
	 * Cholesky decomposition of  a  symmetric  positive-definite matrix.
	 * The result of an algorithm is a representation  of  A  as
	 * A=U^T*U  or A=L*L^T
	 */
	if (!A.isSquare())
		throw Exception("Error: Matrix A is not a square matrix in Cholesky decomposition");
	int i,j,k, n=A.getnRows();
	Matrix a(A);
	p = ColumnVector(n);
	Real sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			sum = a(i+1,j+1);
			for (k = i - 1; k >= 0; k--)	sum -= a(i+1,k+1) * a(j+1,k+1);
			if (i == j) {
				if (sum <= 0)	printf(" a is not positive definite!\n");
				p(i+1) = sqrt(sum);
			}else	a(j+1,i+1) = sum / p(i+1);
		}
	}
	for (i = 0; i < n; i++) {
		a(i+1,i+1) = p(i+1);
		for (j = i + 1; j < n; j++)		a(i+1,j+1) = 0;
	}
	return a;
}
//-------------------------------------------------------------------------------------------
Matrix icholesky(const Matrix &A){
	/*!
	 * Inverse of Cholesky decomposition.
	 */
	if (!A.isSquare())
		throw Exception("Error: Matrix A is not a square matrix in Inverse of Cholesky decomposition");
	int i,j,k, n=A.getnRows();
	Real sum;
	ColumnVector p;
	Matrix a = cholesky(A,p);
	for (i = 0; i < n; i++) {
		a(i+1,i+1) = 1 / p(i+1);
		for (j = i + 1; j < n; j++) {
			sum = 0;
			for (k = i; k < j; k++)		sum -= a(j+1,k+1) * a(k+1,i+1);
			a(j+1,i+1) = sum / p(j+1);
		}
	}
	return a;
}
//-------------------------------------------------------------------------------------------
Matrix SVD(const math::Matrix &A, math::Matrix &U, math::Matrix &V)
/*!
 * @brief computes the singular value decomposition of a real rectangular matrix.
 * 		  based on DSVDC: Original FORTRAN77 version by Jack Dongarra, Cleve Moler,
 * 		  Jim Bunch, Pete Stewart.
 * 		  C++ version by John Burkardt.
 *
 * @Description: The form of the singular value decomposition
 * 				 A(MxN) = U(MxM) * S(MxN) * V(NxN)'
 */
{
	int m = A.getnRows(), n = A.getnColumns();
	int job = 11, lda = m, ldu = m, ldv = n;
	double *a = new double[m*n], *s = new double[m+n], *sigma = new double[m*n];
	//  E must be dimensioned at least maximum(M+1,N).
	double *e = new double[m+n];
	//  S must be dimensioned at least maximum(M+1,N).
	double *u = new double[m*m]; double *uu = new double[m*m];
	double *v = new double[n*n]; double *vv = new double[n*n];
	double *work = new double[m];
	for (int i = 1; i <= m; i++ )
		for (int j = 1; j <= n; j++ )
			a[(i-1)+(j-1)*m] = A(i,j);

	int info = dsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job );

	for (int i = 1; i <= m; i++ )
		for (int j = 1; j <= m; j++ )
			uu[(i-1)*m+(j-1)] =u[(i-1)+(j-1)*m];
	for (int i = 1; i <= n; i++ )
		for (int j = 1; j <= n; j++ )
			vv[(i-1)*n+(j-1)] =v[(i-1)+(j-1)*n];
	Matrix S(m,n);
	U = Matrix(m,m,uu);
	V = Matrix(n,n,vv);
	for ( int j = 1; j <= m; j++ )		S(j,j) = s[j-1];
	if ( info != 0 ){		cout << "  value of the error flag, INFO = " << info << "\n";		return Matrix();	}

	delete [] a;		delete [] s;		delete [] sigma;		delete [] e;
	delete [] u;		delete [] v;		delete [] work;			delete [] uu;
	delete [] vv;
	return S;
}
//-------------------------------------------------------------------------------------------
Matrix chol(const Matrix A, const string &tri)
/*!
 * @brief computes the Cholesky decomposition of a positive definite matrix.
 * 		  based on DCHDC: Original FORTRAN77 version by Jack Dongarra, Cleve Moler,
 * 		  Jim Bunch, Pete Stewart.
 * 		  C++ version by John Burkardt.
 * @Description: uses only the diagonal and upper triangle of A.
 * 				 If A is positive definite, then U = chol(A) produces an upper triangular R
 * 				 so that U'*U = A.
 * 				 If A is not positive definite, an error message is printed.
 */
{
	if (!A.isSquare())
		throw Exception("Error using chol: Input must be a square matrix.");
	int N = A.getnRows();
	int LDA = N;
	double *a = new double[LDA*N], *b = new double[LDA*N], *work = new double[N];
	int i, info, *ipvt = new int[N], j, job = 0;
	for ( i = 1; i <= N; i++ ){
		ipvt[i] = 0;
		for ( j = 1; j <= N; j++ ){	a[i-1+(j-1)*LDA] = A(i,j);}
	}
	info = dchdc ( a, LDA, N, work, ipvt, job );
	if ( info != N ){
		cout << "\n";
		cout << "  DCHDC returned INFO = " << info << "\n";
		cout << "  This means the matrix is not positive definite.\n";
		return Matrix();
	}
	//  Zero out the lower diagonal.
	for ( i = 2; i <= N; i++ )
		for ( j = 1; j <= i-1; j++ ){a[i-1+(j-1)*LDA] = 0.0;}
	double *aa = new double[LDA*N];
	for (int i=1; i<=N; i++){for (int j=1; j<=N; j++){aa[(i-1)*N+(j-1)] =a[(i-1)+(j-1)*N];}}
	delete[] a;
	delete[] b;
	delete[] work;
	delete[] ipvt;
	if (tri == "lower") {
		Matrix res(N, aa);
		delete[] aa;
		return ~res;
	}
	Matrix res(N, aa);
	delete[] aa;
	return res;
}
//-------------------------------------------------------------------------------------------
ColumnVector linsolve(const Matrix &A, const ColumnVector &B)
/*!
 * @breif solves a real general linear system A * X = B.
 * 		  based on DGESL which can solve either of the systems
 * 		  A * X = B or A' * X = B.
 * @param A is m x n matrix.
 * @param B is n x 1 vector.
 * @return n x 1 vector (the solution X)
 */
{
	if (A.getnColumns() != B.getnRows())
		throw Exception("Error using linsolve: Matrix dimensions must agree.");
	int N = A.getnRows();
	int LDA = N;
	int info, *ipvt=new int[N], job;
	double *a = new double[LDA*N], *x = new double[N];
	for (int i = 1; i <= LDA; i++ ){
		x[i-1] = B(i);
		for (int j = 1; j <= LDA; j++ )	{
			a[(i-1)+(j-1)*LDA] = A.getData()[(i-1)*LDA+(j-1)];
		}
	}
	//  Factor the matrix.
	info = dgefa ( a, LDA, N, ipvt );
	if ( info != 0 ){
		cout << "  DGEFA returned an error flag INFO = " << info << "\n";
		return ColumnVector();
	}
	//  Solve the system.
	job = 0;
	dgesl ( a, LDA, N, ipvt, x, job );
	ColumnVector res(N, x);
	delete[] ipvt;
	delete[] a;
	delete[] x;
	return res;
}
//-------------------------------------------------------------------------------------------
double rcond(const Matrix &A){
	if (!A.isSquare())
		throw Exception("Error using rcond: Input must be a square matrix.");
	char norm[] = "O";
	int N = A.getnRows();
	int LDA = N, *ipvt = new int[N], info = 0, *iwork = new int[N];
	double *a = new double[N*N], *z = new double[N], *work = new double[4*N], rcond;
	for (int i = 1; i <= N; i++ ){
		ipvt[i] = 0;
		for (int j = 1; j <= N; j++ ){
			a[(i-1)+(j-1)*LDA] = A(i,j);
		}
	}
	double anorm = dlange_(norm, &N, &N, a, &LDA, work);
	dgetrf_(&N, &N, a, &LDA, ipvt, &info);
	if ( info != 0 ){
		if (info<0)
			cout << "\nLAPACKE_dgetrf failed because the " << info << " argument had an illegal value\n";
		if (info>0)
			cout << "\nLAPACKE_dgetrf failed because U(i,i); where i = " << info << ", is exactly zero.\n"
			<< "The factorization has been completed, but the factor U is exactly singular, and \n"
			<< "division by zero will occur if it is used to solve a system of equations.\n";
		delete[] ipvt;
		delete[] work;
		delete[] iwork;
		delete[] a;
		delete[] z;
		return 0;
	}
	dgecon_(norm, &N, a, &LDA, &anorm, &rcond, work, iwork, &info);
	delete[] ipvt;
	delete[] work;
	delete[] iwork;
	delete[] a;
	delete[] z;
	if ( info != 0 ){
		cout << "\nLAPACKE_dgecon failed because the " << info << " argument had an illegal value\n";
		return 0;
	}
	return rcond;
}
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
Real my_pythag(Real a, Real b){
	//Computes .a2 Cb2/1 = 2 without destructive underflow or overflow.
	//#define SQR(a) ((a)*(a))
	Real absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb)
		return absa * ::sqrt(ONE + SQR(absb / absa));
	else
		return (absb == ZERO ? ZERO : absb * ::sqrt(ONE + SQR(absa / absb)));
}
math::Matrix SVDcmpN(const math::Matrix &a, math::Matrix &u, math::Matrix &v){
	//Given the matrix A stored in u(0..m-1, 0..n-1), this routine computes its singular value
	//decomposition, A = U* W *V^Real and stores the results in the matrices u and v, and the vector w.
	/* Give a matrix A, with logical dimensions M by N and physical
		dimensions MP by NP, this routine computes its singular value
		decomposition, A = U * W * transpose V. The matrix U replaces
		A on output. The diagonal matrix of singular values, W, is output
		as a vector W. The matrix V (not the transpose of V) is output as
		V. M must be greater or equal to N. If it is smaller then A should
		be filled up to square with zero rows.

		adapted from:
		Numerical Recipes in C, The Art of Scientific Computing,
		Press, William H. and Flannery, Brian P. and Teukolsky, Saul A.
		and Vetterling, William T., Cambridge University Press, 1988.
	 */
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

	bool flag;
	int m = a.getnRows();
	int n = a.getnColumns();
	int i, its, j, jj, k, l = 0, nm = 0;
	Real anorm, c, f, g, h, s, scale, x, y, z;
	Real aux;
	math::ColumnVector rv1(n);
	math::Matrix w(n, n);
	// 	v.reConstruct(math::Matrix(n,n));
	// 	u.reConstruct(math::Matrix(a));
	v = math::Matrix(n, n);
	u = math::Matrix(a);

	//added to fix pb of convergence with modern pentiums
	Real eps = ::pow(static_cast<Real>(2.0), static_cast<Real>(-22.0));

	g = scale = anorm = ZERO;	//Householder reduction to bidiagonal form.
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1(i) = scale * g;
		g = s = scale = ZERO;
		if (i <= m) {
			for (k = i; k <= m; k++) scale += fabs(u(k, i));
			if (scale) {
				for (k = i; k <= m; k++) {
					u(k, i) /= scale;
					s += u(k, i)*u(k, i);
				}
				f = u(i, i);
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				u(i, i) = f - g;
				for (j = l; j <= n; j++) {
					for (s = ZERO, k = i; k <= m; k++) s += u(k, i)*u(k, j);
					f = s / h;
					for (k = i; k <= m; k++) u(k, j) += f * u(k, i);
				}
				for (k = i; k <= m; k++) u(k, i) *= scale;
			}
		}
		w(i, i) = scale * g;
		g = s = scale = ZERO;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++) scale += fabs(u(i, k));
			if (scale) {
				for (k = l; k <= n; k++) {
					u(i, k) /= scale;
					s += u(i, k)*u(i, k);
				}
				f = u(i, l);
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				u(i, l) = f - g;
				for (k = l; k <= n; k++) rv1(k) = u(i, k) / h;
				for (j = l; j <= m; j++) {
					for (s = ZERO, k = l; k <= n; k++) s += u(j, k) * u(i, k);
					for (k = l; k <= n; k++) u(j, k) += s * rv1(k);
				}
				for (k = l; k <= n; k++) u(i, k) *= scale;
			}
		}
		anorm = max(anorm, (fabs(w(i, i)) + fabs(rv1(i))));
	}
	for (i = n; i >= 1; i--) {	//Accumulation of right-hand transformations.
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++)	//Double division to avoid possible under ow.
					v(j, i) = (u(i, j) / u(i, l)) / g;
				for (j = l; j <= n; j++) {
					for (s = ZERO, k = l; k <= n; k++) s += u(i, k) * v(k, j);
					for (k = l; k <= n; k++) v(k, j) += s * v(k, i);
				}
			}
			for (j = l; j <= n; j++) v(i, j) = v(j, i) = ZERO;
		}
		v(i, i) = ONE;
		g = rv1(i);
		l = i;
	}
	for (i = min(m, n); i >= 1; i--) {	//Accumulation of left-hand transformations.
		l = i + 1;
		g = w(i, i);
		for (j = l; j <= n; j++) u(i, j) = ZERO;
		if (g) {
			g = ONE / g;
			for (j = l; j <= n; j++) {
				for (s = ZERO, k = l; k <= m; k++) s += u(k, i) * u(k, j);
				f = (s / u(i, i)) * g;
				for (k = i; k <= m; k++) u(k, j) += f * u(k, i);
			}
			for (j = i; j <= m; j++) u(j, i) *= g;
		}
		else for (j = i; j <= m; j++) u(j, i) = ZERO;
		++u(i, i);
	}
	for (k = n; k >= 1; k--) {	//Diagonalization of the bidiagonal form: Loop over
		for (its = 1; its <= 30; its++) {	//singular values, and over allowed iterations.
			flag = 1;
			for (l = k; l >= 1; l--) {	//Test for splitting.
				nm = l - 1;
				if ((Real)(fabs(rv1(l))) <= anorm*eps) {
					flag = 0;
					break;
				}
				if ((Real)(fabs(w(nm, nm))) <= anorm*eps) break;
			}
			if (flag) {
				c = ZERO;	//Cancellation of rv1(l), if l > 0.
				s = ONE;
				for (i = l; i <= k; i++) {
					f = s * rv1(i);
					rv1(i) = c * rv1(i);
					if ((Real)(fabs(f)) <= anorm*eps) break;
					g = w(i, i);
					h = my_pythag(f, g);
					w(i, i) = h;
					h = ONE / h;
					c = g * h;
					s = -f * h;
					for (j = 1; j <= m; j++) {
						y = u(j, nm);
						z = u(j, i);
						u(j, nm) = y * c + z * s;
						u(j, i) = z * c - y * s;
					}
				}
			}
			z = w(k, k);
			if (l == k) {	//Convergence.
				if (z < ZERO) {	//Singular value is made nonnegative.
					w(k, k) = -z;
					for (j = 1; j <= n; j++) v(j, k) = -v(j, k);
				}
				break;
			}
			if (its == 30) ErrorMsg("no convergence in 30 svdcmp iterations");
			x = w(l, l);	//Shift from bottom 2-by-2 minor.
			nm = k - 1;
			y = w(nm, nm);
			g = rv1(nm);
			h = rv1(k);
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (TWO * h * y);
			g = my_pythag(f, ONE);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;	//Next QR transformation:
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1(i);
				y = w(i, i);
				h = s * g;
				g = c * g;
				z = my_pythag(f, h);
				rv1(j) = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1; jj <= n; jj++) {
					x = v(jj, j);
					z = v(jj, i);
					v(jj, j) = x * c + z * s;
					v(jj, i) = z * c - x * s;
				}
				z = my_pythag(f, h);
				w(j, j) = z;	//Rotation can be arbitrary if z == 0.
				if (z) {
					z = ONE / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1; jj <= m; jj++) {
					y = u(jj, j);
					z = u(jj, i);
					u(jj, j) = y * c + z * s;
					u(jj, i) = z * c - y * s;
				}
			}
			rv1(l) = ZERO;
			rv1(k) = f;
			w(k, k) = x;
		}
	}

	//Given the output of decompose, this routine sorts the singular values, and corresponding columns
	//of u and v, by decreasing magnitude. Also, signs of corresponding columns are
	//flipped so as to maximize the number of positive elements.
	// sort the columns of u and v according to w coefficients
	for (j = 1; j <= n; ++j)
		for (i = 1; i <= n - j; ++i)
			if (w(i, i) < w(i + 1, i + 1)) {
				for (jj = 1; jj <= n; ++jj) {
					aux = v(jj, i);
					v(jj, i) = v(jj, i + 1);
					v(jj, i + 1) = aux;
				}
				for (jj = 1; jj <= m; ++jj) {
					aux = u(jj, i);
					u(jj, i) = u(jj, i + 1);
					u(jj, i + 1) = aux;
				}
				aux = w(i, i);
				w(i, i) = w(i + 1, i + 1);
				w(i + 1, i + 1) = aux;
			}
	/*	for (k = 1; k<=n; k++){ //Flip signs.
		int s = 0;
		for (i = 1; i<=m; i++) if (u(i, k) < 0.) s++;
		for (j = 1; j<=n; j++) if (v(j, k) < 0.) s++;
		if (s > (m+n)/2){
		for (i = 1; i<=m; i++) u(i, k) = -u(i, k);
		for (j = 1; j<=n; j++) v(j, k) = -v(j, k);
		}
		}*/
	return w;
}
//-------------------------------------------------------------------------------------------
math::ColumnVector SVDcmp(const math::Matrix &a, math::Matrix &u, math::Matrix &v, int maxitr){
	//Given the matrix A stored in u(0..m-1, 0..n-1), this routine computes its singular value
	//decomposition, A = U* W *V^Real and stores the results in the matrices u and v, and the vector w.
	/* Give a matrix A, with logical dimensions M by N and physical
		dimensions MP by NP, this routine computes its singular value
		decomposition, A = U * W * transpose V. The matrix U replaces
		A on output. The diagonal matrix of singular values, W, is output
		as a vector W. The matrix V (not the transpose of V) is output as
		V. M must be greater or equal to N. If it is smaller then A should
		be filled up to square with zero rows.

		adapted from:
		Numerical Recipes in C, The Art of Scientific Computing,
		Press, William H. and Flannery, Brian P. and Teukolsky, Saul A.
		and Vetterling, William T., Cambridge University Press, 1988.
	 */
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

	bool flag;
	int m = a.getnRows();
	int n = a.getnColumns();
	int i, its, j, jj, k, l = 0, nm = 0;
	Real anorm, c, f, g, h, s, scale, x, y, z;
	Real aux;
	math::ColumnVector rv1(n), w(n);
	// 	v.reConstruct(math::Matrix(n,n));
	// 	u.reConstruct(math::Matrix(a));
	v = math::Matrix(n, n);
	u = math::Matrix(a);

	//added to fix pb of convergence with modern pentiums
	Real eps = ::pow(static_cast<Real>(2.0), static_cast<Real>(-22.0));

	g = scale = anorm = ZERO;	//Householder reduction to bidiagonal form.
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1(i) = scale * g;
		g = s = scale = ZERO;
		if (i <= m) {
			for (k = i; k <= m; k++) scale += fabs(u(k, i));
			if (scale) {
				for (k = i; k <= m; k++) {
					u(k, i) /= scale;
					s += u(k, i)*u(k, i);
				}
				f = u(i, i);
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				u(i, i) = f - g;
				for (j = l; j <= n; j++) {
					for (s = ZERO, k = i; k <= m; k++) s += u(k, i)*u(k, j);
					f = s / h;
					for (k = i; k <= m; k++) u(k, j) += f * u(k, i);
				}
				for (k = i; k <= m; k++) u(k, i) *= scale;
			}
		}
		w(i) = scale * g;
		g = s = scale = ZERO;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++) scale += fabs(u(i, k));
			if (scale) {
				for (k = l; k <= n; k++) {
					u(i, k) /= scale;
					s += u(i, k)*u(i, k);
				}
				f = u(i, l);
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				u(i, l) = f - g;
				for (k = l; k <= n; k++) rv1(k) = u(i, k) / h;
				for (j = l; j <= m; j++) {
					for (s = ZERO, k = l; k <= n; k++) s += u(j, k) * u(i, k);
					for (k = l; k <= n; k++) u(j, k) += s * rv1(k);
				}
				for (k = l; k <= n; k++) u(i, k) *= scale;
			}
		}
		anorm = max(anorm, (fabs(w(i)) + fabs(rv1(i))));
	}
	for (i = n; i >= 1; i--) {	//Accumulation of right-hand transformations.
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++)	//Double division to avoid possible under ow.
					v(j, i) = (u(i, j) / u(i, l)) / g;
				for (j = l; j <= n; j++) {
					for (s = ZERO, k = l; k <= n; k++) s += u(i, k) * v(k, j);
					for (k = l; k <= n; k++) v(k, j) += s * v(k, i);
				}
			}
			for (j = l; j <= n; j++) v(i, j) = v(j, i) = ZERO;
		}
		v(i, i) = ONE;
		g = rv1(i);
		l = i;
	}
	for (i = min(m, n); i >= 1; i--) {	//Accumulation of left-hand transformations.
		l = i + 1;
		g = w(i);
		for (j = l; j <= n; j++) u(i, j) = ZERO;
		if (g) {
			g = ONE / g;
			for (j = l; j <= n; j++) {
				for (s = ZERO, k = l; k <= m; k++) s += u(k, i) * u(k, j);
				f = (s / u(i, i)) * g;
				for (k = i; k <= m; k++) u(k, j) += f * u(k, i);
			}
			for (j = i; j <= m; j++) u(j, i) *= g;
		}
		else for (j = i; j <= m; j++) u(j, i) = ZERO;
		++u(i, i);
	}
	for (k = n; k >= 1; k--) {	//Diagonalization of the bidiagonal form: Loop over
		for (its = 1; its <= maxitr; its++) {	//singular values, and over allowed iterations.
			flag = 1;
			for (l = k; l >= 1; l--) {	//Test for splitting.
				nm = l - 1;
				if ((Real)(fabs(rv1(l))) <= anorm*eps) {
					flag = 0;
					break;
				}
				if ((Real)(fabs(w(nm))) <= anorm*eps) break;
			}
			if (flag) {
				c = ZERO;	//Cancellation of rv1(l), if l > 0.
				s = ONE;
				for (i = l; i <= k; i++) {
					f = s * rv1(i);
					rv1(i) = c * rv1(i);
					if ((Real)(fabs(f)) <= anorm*eps) break;
					g = w(i);
					h = my_pythag(f, g);
					w(i) = h;
					h = ONE / h;
					c = g * h;
					s = -f * h;
					for (j = 1; j <= m; j++) {
						y = u(j, nm);
						z = u(j, i);
						u(j, nm) = y * c + z * s;
						u(j, i) = z * c - y * s;
					}
				}
			}
			z = w(k);
			if (l == k) {	//Convergence.
				if (z < ZERO) {	//Singular value is made nonnegative.
					w(k) = -z;
					for (j = 1; j <= n; j++) v(j, k) = -v(j, k);
				}
				break;
			}
			if (its == maxitr) ErrorMsg("no convergence in 30 svdcmp iterations");
			x = w(l);	//Shift from bottom 2-by-2 minor.
			nm = k - 1;
			y = w(nm);
			g = rv1(nm);
			h = rv1(k);
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (TWO * h * y);
			g = my_pythag(f, ONE);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;	//Next QR transformation:
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1(i);
				y = w(i);
				h = s * g;
				g = c * g;
				z = my_pythag(f, h);
				rv1(j) = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1; jj <= n; jj++) {
					x = v(jj, j);
					z = v(jj, i);
					v(jj, j) = x * c + z * s;
					v(jj, i) = z * c - x * s;
				}
				z = my_pythag(f, h);
				w(j) = z;	//Rotation can be arbitrary if z == 0.
				if (z) {
					z = ONE / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1; jj <= m; jj++) {
					y = u(jj, j);
					z = u(jj, i);
					u(jj, j) = y * c + z * s;
					u(jj, i) = z * c - y * s;
				}
			}
			rv1(l) = ZERO;
			rv1(k) = f;
			w(k) = x;
		}
	}

	//Given the output of decompose, this routine sorts the singular values, and corresponding columns
	//of u and v, by decreasing magnitude. Also, signs of corresponding columns are
	//flipped so as to maximize the number of positive elements.
	// sort the columns of u and v according to w coefficients
	for (j = 1; j <= n; ++j)
		for (i = 1; i <= n - j; ++i)
			if (w(i) < w(i + 1)) {
				for (jj = 1; jj <= n; ++jj) {
					aux = v(jj, i);
					v(jj, i) = v(jj, i + 1);
					v(jj, i + 1) = aux;
				}
				for (jj = 1; jj <= m; ++jj) {
					aux = u(jj, i);
					u(jj, i) = u(jj, i + 1);
					u(jj, i + 1) = aux;
				}
				aux = w(i);
				w(i) = w(i + 1);
				w(i + 1) = aux;
			}
	/*	for (k = 1; k<=n; k++){ //Flip signs.
		int s = 0;
		for (i = 1; i<=m; i++) if (u(i, k) < 0.) s++;
		for (j = 1; j<=n; j++) if (v(j, k) < 0.) s++;
		if (s > (m+n)/2){
		for (i = 1; i<=m; i++) u(i, k) = -u(i, k);
		for (j = 1; j<=n; j++) v(j, k) = -v(j, k);
		}
		}*/
	return w;
}
//-------------------------------------------------------------------------------------------
math::ColumnVector SVDcmp(const math::Matrix &a){
	math::Matrix u, v;
	return SVDcmp(a, u, v);
}
//-------------------------------------------------------------------------------------------

void svbksb(Matrix &u, ColumnVector &w, Matrix &v, ColumnVector &b, ColumnVector &x){
	int jj, j, i;
	Real s;

	int m = u.getnRows();
	int n = u.getnColumns();
	ColumnVector tmp(n);
	for (j = 1; j <= n; j++) {
		s = 0.0;
		if (w(j) != 0.0) {
			for (i = 1; i <= m; i++)	
				s += u(i, j) * b(i);
			s /= w(j);
		}
		tmp(j) = s;
	}
	for (j = 1; j<n; j++) {
		s = 0.0;
		for (jj = 1; jj<n; jj++)	s += v(j, jj) * tmp(jj);
		x(j) = s;
	}
}
//-------------------------------------------------------------------------------------------
void ErrorMsg(string error_text){
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text.c_str());
	fprintf(stderr, "...now exiting to system...\n");
	fprintf(stderr, "\n        hit return to continue...\n");
	getchar();
	exit(-1);
}
//-------------------------------------------------------------------------------------------
void QR(const Matrix &M, Matrix &R, Matrix &Q)
/*!
 * @Brief computes the QR factorization of a real rectangular matrix.
 * 		  based on DQRDC: Original FORTRAN77 version by Jack Dongarra, Cleve Moler,
 * 		  Jim Bunch, Pete Stewart.
 * 		  C++ version by John Burkardt.
 */
{
	int N = M.getnRows(), P = M.getnColumns();
	int LDA = N;
	double 	 *b = new double[LDA*P],    *q = new double[N*N], *qraux = new double[P];
	double *qty = new double[N]	   ,   *qy = new double[N]	, 	  *r = new double[N*P];
	double *rsd = new double[N]	   , *work = new double[P]	, 	 *xb = new double[N];
	double 	 *y = new double[N]	   , 	*a = new double[LDA*P], *rr = new double[N*P];
	double *qq = new double[N*N];
	int i, info, *ipvt = new int[P], j, job;
	for (int i = 1; i <= N; i++ ){
		for (int j = 1; j <= P; j++ ){
			a[(i-1)+(j-1)*LDA] = M.getData()[(i-1)*P+(j-1)];
			if (i==1)	ipvt[j-1] = 0;
		}
	}
	//  Decompose the matrix.
	job = 0;
	dqrdc ( a, LDA, N, P, qraux, ipvt, work, job );
	//  The resulting R factor.
	for ( i = 1; i <= N; i++ ){
		for ( j = 1; j <= P; j++ ){
			if ( j < i ){	r[i-1+(j-1)*N] = 0.0;				}
			else{			r[i-1+(j-1)*N] = a[i-1+(j-1)*LDA];	}
		}
	}
	//double *rr = new double[N*P];
	for ( i = 1; i <= N; i++ ){
		for ( j = 1; j <= P; j++ ){
			rr[(i-1)*P+(j-1)] = r[(i-1)+(j-1)*LDA];
		}
	}
	R = Matrix(N,P,rr);
	/*!
	 * Call DQRSL to extract the information about the Q matrix.
	 * We do this, essentially, by asking DQRSL to tell us the
	 * value of Q*Y, where Y is a column of the identity matrix.
	 */
	job = 10000;
	for ( i = 1; i <= N; i++ ){
		//  Set the vector Y.
		for ( j = 1; j <= N; j++ ){		y[j-1] = 0.0;	}
		y[i-1] = 1.0;
		//  Ask DQRSL to tell us what Q*Y is.
		info = dqrsl ( a, LDA, N, P, qraux, y, qy, qty, b, rsd, xb, job );

		if ( info != 0 ){
			cout << "  Error!  DQRSL returns INFO = " << info << "\n";
			return;
		}
		//  Copy QY into the appropriate column of Q.
		for ( j = 1; j <= N; j++ ){		q[j-1+(i-1)*N] = qy[j-1];	}
	}
	//  The Q matrix we have extracted.
	for ( i = 1; i <= N; i++ ){
		for ( j = 1; j <= N; j++ ){
			qq[(i-1)*N+(j-1)] = q[(i-1)+(j-1)*N];
		}
	}
	Q = Matrix(N,N,qq);

	delete [] qraux;
	delete [] ipvt;
	delete [] work;
	delete [] qty;
	delete [] rsd;
	delete [] qy;
	delete [] xb;
	delete [] rr;
	delete [] qq;
	delete [] b;
	delete [] q;
	delete [] r;
	delete [] y;
	delete [] a;
}
//-------------------------------------------------------------------------------------------
void QRZ(math::Matrix &X, math::Matrix &U){//U is UpperTriangularMatrix
	//REPORT
	//Tracer et("QRZ(1)");
	int n = X.getnRows(), s = X.getnColumns();
	U.resize(s, s, ZERO);
	if (n == 0 || s == 0)	return;

	Real *xi0 = X.getData(), *u0 = U.getData(), *u;
	int j, k, J = s, i = s;
	while (i--){
		Real *xj0 = xi0, *xi = xi0;
		k = n;
		if (k) for (;;){
			u = u0;
			Real Xi = *xi, *xj = xj0;
			j = J;
			while (j--)		*u++ += Xi * *xj++;

			if (!(--k))		break;

			xi += s;
			xj0 += s;
		}

		Real sum = sqrt(*u0); *u0 = sum; u = u0 + 1;
		if (sum == 0.0){
			//REPORT
			j = J - 1;
			while (j--)	*u++ = 0.0;

			xj0 = xi0++;
			k = n;
			if (k)
				for (;;){
					*xj0 = 0.0;
					if (!(--k)) break;

					xj0 += s;
				}
			u0 += J--;
		}
		else{
			int J1 = J - 1; j = J1;
			while (j--)	*u++ /= sum;

			xj0 = xi0;
			xi = xi0++;
			k = n;
			if (k)
				for (;;){
					u = u0 + 1;
					Real Xi = *xi, *xj = xj0;
					Xi /= sum;
					*xj++ = Xi;
					j = J1;

					while (j--)		*xj++ -= *u++ * Xi;

					if (!(--k))		break;

					xi += s;
					xj0 += s;
				}
			u0 += J--;
		}
	}
	// 		cout << "\nU = \n" << setw(7) << setprecision(3) << U;
	int num = 0;
	if (s % 2 == 0)
		num = (s / 2) * s + (s / 2);
	else
		num = ((s + 1) / 2) * s;
	Matrix temp = U;
	U.zeros();
	int kk = 1, mm = 1;
	for (int i = 1; i <= s; i++){
		for (int j = i; j <= s; j++){
			U(i, j) = temp(mm, kk);
			kk++;
			if (kk > s){ kk = 1; mm++; }
		}
		if (kk>num) break;
	}
}


void QRZ(math::Matrix &X, math::Matrix &Y, math::ColumnVector &M){
	//REPORT
	//Tracer et("QRZ(2)");
	int n = X.getnRows(), s = X.getnColumns(), t = Y.getnColumns();
	if (Y.getnRows() != n)
		throw Exception("Unequal column lengths");
	M.resize(s, ZERO);
	Real *m0 = M.getData(), *m, *xi0 = X.getData();
	int j, k; int i = s;
	while (i--){
		Real *xj0 = Y.getData(), *xi = xi0; k = n;
		if (k)
			for (;;){
				m = m0;
				Real Xi = *xi, *xj = xj0;
				j = t;
				while (j--)	*m++ += Xi * *xj++;

				if (!(--k)) break;

				xi += s;
				xj0 += t;
			}

		xj0 = Y.getData();
		xi = xi0++;
		k = n;
		if (k) for (;;){
			m = m0;
			Real Xi = *xi, *xj = xj0;
			j = t;
			while (j--)	*xj++ -= *m++ * Xi;

			if (!(--k)) break;

			xi += s;
			xj0 += t;
		}
		m0 += t;
	}
}

void odeint(Matrix(*xdot)(Real time, const Matrix & xin),
		Matrix & xo, Real to, Real tf, Real eps, Real h1, Real hmin,
		int & nok, int & nbad,
		RowVector & tout, Matrix & xout, Real dtsav)
/*!
		@brief Integrate the ordinary differential equation xdot from time to
		to time tf using an adaptive step size strategy

		adapted from:
		Numerical Recipes in C, The Art of Scientific Computing,
		Press, William H. and Flannery, Brian P. and Teukolsky, Saul A.
		and Vetterling, William T., Cambridge University Press, 1988.
 */
{
#define MAXSTP 10000
#define TINY 1.0e-30
	Real tsav, t, hnext, hdid, h;
	RowVector tv(1);

	Matrix xscal, x, dxdt;

	tv = to;
	tout = tv;
	xout = xo;
	xscal = xo;
	t = to;
	h = (tf > to) ? fabs(h1) : -fabs(h1);
	nok = (nbad) = 0;
	x = xo;
	tsav = t;
	for (int nstp = 1; nstp <= MAXSTP; nstp++){
		dxdt = (*xdot)(t, x);
		for (int i = 1; i <= x.getnRows(); i++)
			xscal(i, 1) = fabs(x(i)) + fabs(dxdt(i)*h) + TINY;
		if ((t + h - tf)*(t + h - to) > 0.0) h = tf - t;
		rkqc(x, dxdt, t, h, eps, xscal, hdid, hnext, xdot);
		if (hdid == h) ++(nok); else ++(nbad);
		if ((t - tf)*(tf - to) >= 0.0) {
			xo = x;
			tv = t;
			tout = tout | tv;
			xout = xout | x;
			return;
		}
		if (fabs(t - tsav) > fabs(dtsav)) {
			tv = t;
			tout = tout | tv;
			xout = xout | x;
			tsav = t;
		}
		if (fabs(hnext) <= hmin) {
			cerr << "Step size too small in ODEINT\n";
			cerr << setw(7) << setprecision(3) << (tout & xout).Transpose();
			exit(1);
		}
		h = hnext;
	}
	cerr << "Too many step in routine ODEINT\n";
	exit(1);
}

void rkqc(Matrix & x, Matrix & dxdt, Real & t, Real htry,
		Real eps, Matrix & xscal, Real & hdid, Real & hnext,
		Matrix(*xdot)(Real time, const Matrix & xin))
/*!
		@brief Compute one adaptive step based on two rk4.

		adapted from:
		Numerical Recipes in C, The Art of Scientific Computing,
		Press, William H. and Flannery, Brian P. and Teukolsky, Saul A.
		and Vetterling, William T., Cambridge University Press, 1988.
 */
{
#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666
#define SAFETY 0.9
#define ERRCON 6.0E-4
	Real tsav, hh, h, temp, errmax;
	Matrix dxsav, xsav, xtemp;

	tsav = t;
	xsav = x;
	dxsav = dxdt;
	h = htry;
	for (;;) {
		hh = 0.5*h;
		xtemp = rk4(xsav, dxsav, tsav, hh, xdot);
		t = tsav + hh;
		dxdt = (*xdot)(t, xtemp);
		x = rk4(xtemp, dxdt, t, hh, xdot);
		t = tsav + h;
		if (t == tsav) {
			cerr << "Step size too small in routine RKQC\n";
			exit(1);
		}
		xtemp = rk4(xsav, dxsav, tsav, h, xdot);
		errmax = 0.0;
		xtemp = x - xtemp;
		for (int i = 1; i <= x.getnRows(); i++) {
			temp = fabs(xtemp(i, 1) / xscal(i, 1));
			if (errmax < temp) errmax = temp;
		}
		errmax /= eps;
		if (errmax <= 1.0) {
			hdid = h;
			hnext = (errmax > ERRCON ?
					SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
			break;
		}
		h = SAFETY*h*exp(PSHRNK*log(errmax));
	}
	x += xtemp*FCOR;
}

Matrix rk4(const Matrix & x, const Matrix & dxdt, Real t, Real h,
		Matrix(*xdot)(Real time, const Matrix & xin))
/*!
 * @brief Compute one Runge-Kutta fourth order step
 *
 * adapted from:
 * Numerical Recipes in C, The Art of Scientific Computing,
 * Press, William H. and Flannery, Brian P. and Teukolsky, Saul A.
 * and Vetterling, William T., Cambridge University Press, 1988.
 */
{
	Matrix xt, xout, dxm, dxt;
	Real th, hh, h6;

	hh = h*static_cast<Real>(0.5);
	h6 = h / static_cast<Real>(6.0);
	th = t + hh;
	xt = x + hh*dxdt;
	dxt = (*xdot)(th, xt);
	xt = x + hh*dxt;
	dxm = (*xdot)(th, xt);
	xt = x + h*dxm;
	dxm += dxt;
	dxt = (*xdot)(t + h, xt);
	xout = x + h6*(dxdt + dxt + static_cast<Real>(2.0)*dxm);

	return xout;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
ColumnVector Integ_Trap(const ColumnVector & present, ColumnVector & past, const Real dt){
	ColumnVector integration(present.getnRows());
	integration = (past+present)*0.5*dt;
	past = present;
	return integration;
}

};
