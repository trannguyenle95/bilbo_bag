//=============================================================================
//==============================================================================

//	file:	matrix.cpp

//	author:	Fares J. Abu-Dakka
//	date:	2006 - 2015

//==============================================================================
#include <float.h>
#include "matrix.hpp"

#define TRUE 1
#define FALSE 0

namespace math{
    //-------------------------------------------------------------------------------------------
    Real getTime()	{	return Real(clock())/CLOCKS_PER_SEC;		}
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
	bool BaseMatrix::in_range(int r, int c) const
	/*!
	  @brief Check whether element (r,c) is within the matrix
	*/
	{
		return ((r >= 1) && (r <= _nRows) && (c >= 1) && (c <= _nColumns));
	}
	//-------------------------------------------------------------------------------------------
	bool BaseMatrix::row_in_range(int r) const
	/*!
	  @brief Check whether row \c r is in the allowed range
	*/
	{
		return ((r >= 1) && (r <= _nRows));
	}
	//-------------------------------------------------------------------------------------------
	bool BaseMatrix::col_in_range(int c) const
	/*!
	  @brief Check whether column \c c is in the allowed range
	*/
	{
		return ((c >= 1) && (c <= _nColumns));
	}
	//-------------------------------------------------------------------------------------------
	bool BaseMatrix::in_range(int i) const
	/*!
	  @brief Check whether element \c i is in the allowed range
	*/
	{
		return ((i >= 0) && (i <= size()));
	}
	//-------------------------------------------------------------------------------------------
	const BaseMatrix::assignMatix BaseMatrix::operator << (const Real& val)
	/*!
	  @brief Assign the given value to every spot in this matrix
	*/
	{
		for (long r = 1; r <= _nRows; ++r)
			for (long c = 1; c <= _nColumns; ++c)
				(*this)(r,c) = val;
		return BaseMatrix::assignMatix(this);
	}
	//-------------------------------------------------------------------------------------------
    BaseMatrix::BaseMatrix(int r, int c, Real x)
	/*!
	  @brief Constructor: Creates a matrix with all values equal to x.
	*/
	{
		if ( r<1 || c<1 )
            throw InvalidIndex();

		_nRows = r;			_r=1;	has_been_used = false;
		_nColumns = c;		_c=1;
		_data = new Real[_nRows*_nColumns];    //! Allocate memory
		assert(_data != 0);          //! Check that memory was allocated
		for ( int i=0; i<_nRows*_nColumns; i++)
			_data[i] = x;
    }
	//-------------------------------------------------------------------------------------------
    BaseMatrix::BaseMatrix(int r, int c, const Real *x)
	/*!
		@brief Constructor: builds a matrix from a given array.
	*/
	{
		if ( r<1 || c<1 )	throw InvalidIndex();
		if (!x)				throw Exception("Invalid array");
		_nRows = r;		_nColumns = c;
		_r = 1;			_c = 1;			has_been_used = false;

		_data = new Real[_nRows*_nColumns];
		assert(_data != 0);          //! Check that memory was allocated
		for ( int i=0; i<_nRows*_nColumns; i++)		// I think I should add an exception to check the array size
			_data[i] = x[i];
    }
	//-------------------------------------------------------------------------------------------
    BaseMatrix::BaseMatrix(int r, int c, const string &str)
	/*!
		@brief Constructor: builds a matrix from a given array.
	*/
	{
		if ( r<1 || c<1 )	throw InvalidIndex();
		if (!(str.compare("NULL")) || !(str.compare("null")) || !(str.compare("0"))){
			_nRows = r;		_nColumns = c;
			_r = 1;			_c = 1;			has_been_used = false;
			_data = new Real[_nRows*_nColumns];
			assert(_data != 0);          //! Check that memory was allocated
		}else
			throw Exception("Invalid string");
    }
    //-------------------------------------------------------------------------------------------
    void BaseMatrix::copy(const BaseMatrix & m)
	/*!
		@brief Copy Matrix m to this
	*/
	{
		(*this)._nRows = m._nRows;
		(*this)._nColumns = m._nColumns;
		(*this)._data = new Real[_nColumns * _nRows];
		assert(_data != 0);          //! Check that memory was allocated
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++)
			(*this)._data[i] = m._data[i];
    }
    //-------------------------------------------------------------------------------------------
    BaseMatrix::~BaseMatrix(){
	if (_data)
            delete [] _data;
    }
   //-------------------------------------------------------------------------------------------
    BaseMatrix & BaseMatrix::operator =(const BaseMatrix & m){
		if ((*this)._data)	delete [] (*this)._data;
		copy(m);
		return *this;
    }
	//-------------------------------------------------------------------------------------------
	BaseMatrix & BaseMatrix::operator >> (const Real *val ){
		int i = 0;
		for (long r = 1; r <= _nRows; ++r)
			for (long c = 1; c <= _nColumns; ++c)
				(*this)(r,c) = val[i++];
		return *this;
	}
    //-------------------------------------------------------------------------------------------
    BaseMatrix & BaseMatrix::operator =(Real x){
	if (!_data)
            throw Exception("BaseMatrix::operator =(Real x): Invalid Matrix");
	for (int i=0; i<_nRows*_nColumns ; i++)
            _data[i] = x;
	return *this;
    }
    //-------------------------------------------------------------------------------------------
	const Real & BaseMatrix::operator () (int i, int j) const{
		if ( !in_range(i,j) )
			throw Exception("BaseMatrix::operator () (int i, int j): Invalid Index");
		return _data[(i-1) * _nColumns + (j-1)];
	}
    //-------------------------------------------------------------------------------------------
	Real & BaseMatrix::operator () (int i, int j){
		if ( !in_range(i,j) )
			throw Exception("BaseMatrix::operator () (int i, int j): Invalid Index");
		return _data[(i-1) * _nColumns + (j-1)];
	}
	//-------------------------------------------------------------------------------------------
	Real & BaseMatrix::operator () (int i, int j, Real e) {
		if ( !in_range(i,j) )
			throw Exception("BaseMatrix::operator () (int i, int j, Real e) : Invalid Index");
		return _data[(i-1) * _nColumns + (j-1)] = e;
	}
	//-------------------------------------------------------------------------------------------
	BaseMatrix & BaseMatrix::operator () (Real *array, int n){
		for (int i=0; i<n && i<_nRows*_nColumns; i++)// maybe I should add warning if the array size is less or more the matrix size
			_data[i] = array[i];
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	void BaseMatrix::setElement( int i, int j, Real e){
		// Set an Element in the Matrix
		if ( !in_range(i,j) )
			_data[(i-1) * _nColumns + (j-1)] = e;
	}
	//-------------------------------------------------------------------------------------------
	ostream & operator << (ostream & o, const BaseMatrix & m){
		for (int i=1; i<=m.getnRows(); i++){
			for (int j=1; j<=m.getnColumns(); j++)
				o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << m(i,j) << " ";
			o << endl;
		}
		//o << endl;
		return o;
	}
	//-------------------------------------------------------------------------------------------
	void BaseMatrix::print (char const* ch, int prec, int w) const{
		int n = 0;
		if (ch){
			string s(ch);
			std::cout << ch << "= ";
			n = s.size();
		}
		//cout << (*this);
		for (int i=1; i<=_nRows; i++){
			for (int j=1; j<=_nColumns; j++)
				cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << (*this)(i,j) << " ";
			if (i<_nRows)
				cout << endl << std::setw(n+2) << " ";
			else
				cout << endl;
		}
		//cout << endl;
	}
    //-------------------------------------------------------------------------------------------
	void BaseMatrix::zeros(){
		// Replaces all values in the Matrix with zero.
		for ( int i=0; i<_nRows * _nColumns; i++ )
			_data[i] = 0;
	}
	//-------------------------------------------------------------------------------------------
	void BaseMatrix::ones (){
		// Replaces all values in the Matrix with one.
		for ( int i=0; i<_nRows * _nColumns; i++ )
			_data[i] = 1;
	}
	//-------------------------------------------------------------------------------------------
	void BaseMatrix::fill(Real val){
		for ( int i=0; i<_nRows * _nColumns; i++ )
			_data[i] = val;
	}
	//-------------------------------------------------------------------------------------------
	bool BaseMatrix::operator ==(const BaseMatrix & m)const{
		if (!( _nRows == m._nRows && _nColumns == m._nColumns )){
			std::cout << "warning: operator ==: The matrices have different sizes...";
			return false;
		}
		for ( int i=0; i<_nRows*_nColumns; i++ )
			if ( _data[i] != m._data[i] )
				return false;
		return true;
	}
    //-------------------------------------------------------------------------------------------
	Real BaseMatrix::sum(const BaseMatrix &m){
		Real s = 0;
		for ( int i=0; i<m._nRows * m._nColumns; i++ )
			s += m._data[i];
		return s;
	}
    //-------------------------------------------------------------------------------------------
	Real BaseMatrix::maxV(const BaseMatrix &m){
		Real maxv = m._data[0];
		for ( int i=1; i<m._nRows * m._nColumns; i++ )
			if (maxv < m._data[i])		maxv = m._data[i];
		return maxv;
	}
    //-------------------------------------------------------------------------------------------
	Real BaseMatrix::maxV(const BaseMatrix &m,int &ind1,int &ind2){
		Real maxv = m(1,1);
		for ( int i=1; i<=m.getnRows(); i++ )
			for ( int j=1; j<=m.getnRows(); j++ )
				if (maxv < m(i,j)){
					maxv = m(i,j);
					ind1 = i;
					ind2 = j;
				}
		return maxv;
	}
    //-------------------------------------------------------------------------------------------
	Real BaseMatrix::minV(const BaseMatrix &m){
		Real minv = m._data[0];
		for ( int i=1; i<m._nRows * m._nColumns; i++ )
			if (minv > m._data[i])		minv = m._data[i];
		return minv;
	}
    //-------------------------------------------------------------------------------------------
	Real BaseMatrix::minV(const BaseMatrix &m,int &ind1,int &ind2){
		Real minv = m(1,1);
		for ( int i=1; i<=m.getnRows(); i++ )
			for ( int j=1; j<=m.getnRows(); j++ )
				if (minv > m(i,j)){
					minv = m(i,j);
					ind1 = i;
					ind2 = j;
				}
		return minv;
	}
    //-------------------------------------------------------------------------------------------
	Real BaseMatrix::MaximumAbsoluteValue(const BaseMatrix &m){
		Real maxv = abs(m._data[0]);
		for ( int i=1; i<m._nRows * m._nColumns; i++ )
			if (maxv > abs(m._data[i]))		maxv = abs(m._data[i]);
		return maxv;
	}
    //-------------------------------------------------------------------------------------------
	Real BaseMatrix::MinimumAbsoluteValue(const BaseMatrix &m){
		Real minv = abs(m._data[0]);
		for ( int i=1; i<m._nRows * m._nColumns; i++ )
			if (minv < abs(m._data[i]))		minv = abs(m._data[i]);
		return minv;
	}
	Real BaseMatrix::SumSquare()const{
		Real sum = ZERO;
		for ( int i=0; i<_nRows * _nColumns; i++ )
                        sum += SQR(*(_data+i));
                return sum;
	}
   //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    //------------------------------------------Matrix-------------------------------------------
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    /************************** LU transformation ****************************/
	void Matrix::ludcmp(Matrix &a, int *indx, Real &d)const{
		int i, imax, j, k;
		Real big, dum, sum, temp;

		int n = a.getnRows();
		ColumnVector vv(n);
		d = ONE;
		for (i=0; i<n; i++) {
			big = ZERO;
			for (j=0; j<n; j++)
				if ((temp = fabs(a(i+1,j+1))) > big) big = temp;
			if (big < EPSilon) throw Exception("Singular matrix in routine ludcmp");
			vv(i+1) = ONE/big;
		}
		for (j=0; j<n; j++) {
			for (i=0; i<j; i++) {
				sum = a(i+1,j+1);
				for (k=0; k<i; k++)		sum -= a(i+1,k+1) * a(k+1,j+1);
				a(i+1,j+1) = sum;
			}
			big = ZERO;
			for (i=j; i<n; i++) {
				sum = a(i+1,j+1);
				for (k=0; k<j; k++)		sum -= a(i+1,k+1) * a(k+1,j+1);
				a(i+1,j+1) = sum;
				if ((dum = vv(i+1) * fabs(sum)) >= big) {
					big = dum;
					imax = i;
				}
			}
			if (j != imax) {
				for (k=0; k<n; k++) {
					dum = a(imax+1,k+1);
					a(imax+1,k+1) = a(j+1,k+1);
					a(j+1,k+1) = dum;
				}
				d = -d;
				vv(imax+1) = vv(j+1);
			}
			indx[j] = imax;
			if (fabs(a(j+1,j+1)) < EPSilon)		a(j+1,j+1) = EPSilon;
			if (j != n-1) {
				dum = ONE/(a(j+1,j+1));
				for (i=j+1; i<n; i++)	a(i+1,j+1) *= dum;
			}
		}
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::lubksb(Matrix &a, int *indx, ColumnVector &b)const{
		int i, ii=0, ip, j;
		Real sum = ZERO;

		int n = a.getnRows();
		for (i=0; i<n; i++) {
			ip = indx[i];
			sum = b(ip+1);
			b(ip+1) = b(i+1);
			if (ii != 0)
				for (j=ii-1; j<i; j++) sum -= a(i+1,j+1) * b(j+1);
			else if (sum != ZERO)
				ii = i+1;
			b(i+1) = sum;
		}
		for (i=n-1; i>=0; i--) {
			sum = b(i+1);
			for (j=i+1; j<n; j++) sum -= a(i+1,j+1) * b(j+1);
			b(i+1) = sum/a(i+1,i+1);
		}
	}
    //-------------------------------------------------------------------------------------------
	Matrix Matrix::InverseLU()const{
		if (!isSquare())
			throw Exception("InverseLU()::The matrix is NOT square matrix");
		Matrix Ainv = *this, A = *this;
		int *indx, n = Ainv.getnRows();	indx = new int[Ainv.getnRows()];
		Real d;
		ludcmp(A,indx,d);
		for (int i=1; i<=n; i++){
			ColumnVector v(n);
			v(i) = 1;
			lubksb(A,indx,v);
			Ainv.setCol(v,i);
		}
		return Ainv;
	}
    //-------------------------------------------------------------------------------------------
	Matrix Matrix::Inverse()const{
		if (!isSquare())
			throw Exception("Inverse()::The matrix is NOT square matrix");
		int N = (*this).getnRows();
		//! Copy matrix to ensure Ainv is same size
		Matrix Ainv(*this), A(*this);

		int i, j, k;
		ColumnVector scale(N);	//! Scale factor
		Matrix b(N);		//! Work array
		int *index;  index = new int [N+1];

		//!* Matrix b is initialized to the identity matrix
		b = eye(N);

		//!* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
		for( i=1; i<=N; i++ ) {
			index[i] = i;			  //! Initialize row index list
			Real scalemax = 0.;
			for( j=1; j<=N; j++ )
				scalemax = (scalemax > fabs(A(i,j))) ? scalemax : fabs(A(i,j));
			scale(i) = scalemax;
		}

		//!* Loop over rows k = 1, ..., (N-1)
		int signDet = 1;
		for( k=1; k<=N-1; k++ ) {
			//!* Select pivot row from max( |a(j,k)/s(j)| )
			Real ratiomax = ZERO;
			int jPivot = k;
			for( i=k; i<=N; i++ ) {
				Real ratio = fabs(A(index[i],k))/scale(index[i]);
				if( ratio > ratiomax ) {
					jPivot=i;
					ratiomax = ratio;
				}
			}
			//!* Perform pivoting using row index list
			int indexJ = index[k];
			if( jPivot != k ) {	          //! Pivot
				indexJ = index[jPivot];
				index[jPivot] = index[k];   //! Swap index jPivot and k
				index[k] = indexJ;
				signDet *= -1;			  //! Flip sign of determinant
			}
			//!* Perform forward elimination
			for( i=k+1; i<=N; i++ ) {
				Real coeff = A(index[i],k)/A(indexJ,k);
				for( j=k+1; j<=N; j++ )
					A(index[i],j) -= coeff*A(indexJ,j);
				A(index[i],k) = coeff;
				for( j=1; j<=N; j++ )
					b(index[i],j) -= A(index[i],k)*b(indexJ,j);
			}
		}
		//!* Compute determinant as product of diagonal elements
		Real determ = static_cast<Real>(signDet);	   //! Sign of determinant
		for( i=1; i<=N; i++ )
			determ *= A(index[i],i);

		//!* Perform back-substitution
		for( k=1; k<=N; k++ ) {
			Ainv(N,k) = b(index[N],k)/A(index[N],N);
			for( i=N-1; i>=1; i--) {
				Real sum = b(index[i],k);
				for( j=i+1; j<=N; j++ )
					sum -= A(index[i],j)*Ainv(j,k);
				Ainv(i,k) = sum/A(index[i],i);
			}
		}

		delete [] index;	//! Release allocated memory
		return Ainv;
	}
    //-------------------------------------------------------------------------------------------
	Matrix Matrix::InverseCholesky()const{
		if (!isSquare())
			throw Exception("InverseCholesky()::The matrix is NOT square matrix");
		int i,j,k;
		Matrix a = icholesky(*this);
	/*	for (i = 0; i < _nColumns; i++) {
			for (j = i + 1; j < _nColumns; j++) {
				a(i+1,j+1) = ZERO;
			}
		}*/
		for (i = 0; i < _nColumns; i++) {
			a(i+1,i+1) *= a(i+1,i+1);
			for (k = i + 1; k < _nColumns; k++)
				a(i+1,i+1) += a(k+1,i+1) * a(k+1,i+1);
			for (j = i + 1; j < _nColumns; j++)
				for (k = j; k < _nColumns; k++)
					a(i+1,j+1) += a(k+1,i+1) * a(k+1,j+1);
		}
		for (i = 0; i < _nColumns; i++)
			for (j = 0; j < i; j++)
				a(i+1,j+1) = a(j+1,i+1);
		return a;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::InverseGJ()const{
		if (!isSquare())
			throw Exception("Inverse()::The matrix is NOT square matrix");
		Matrix a(*this), b(*this);
		int n = a.getnRows();
		int m = a.getnColumns();
		int *indxc, *indxr, *ipiv;
		int i, icol, irow, j, k, l, ll;
		Real big, dum, pivinv;

		indxc = new int[n];
		indxr = new int[n];
		ipiv = new int[n];
		for (j=0; j<n; j++) ipiv[j] = 0;
		for (i=0; i<n; i++) {
			big = ZERO;
			for (j=0; j<n; j++)
				if (ipiv[j] != 1)
					for (k=0; k<n; k++) {
						if (ipiv[k] == 0) {
							if (fabs(a(j+1,k+1)) >= big) {
								big = fabs(a(j+1,k+1));
								irow = j;
								icol = k;
							}
						} else if (ipiv[k] > 1) {
							throw Exception("error: Matrix::InverseGJ(): Gauss-Jordan: Singular Matrix-1");
							exit(1);
						}
					}
			++(ipiv[icol]);
			if (irow != icol) {
				for (l=0; l<n; l++) Swap(a(irow+1,l+1),a(icol+1,l+1));
				for (l=0; l<m; l++) Swap(b(irow+1,l+1),b(icol+1,l+1));
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if (a(icol+1,icol+1) == ZERO) cout << "gaussj: Singular Matrix-2" << endl;
			pivinv = ONE/a(icol+1,icol+1);
			a(icol+1,icol+1) = ONE;
			for (l=0;  l<n; l++) a(icol+1,l+1) *= pivinv;
			for (l=0;  l<m; l++) b(icol+1,l+1) *= pivinv;
			for (ll=0;ll<n;ll++)
				if (ll != icol) {
					dum = a(ll+1,icol+1);
					a(ll+1,icol+1) = ZERO;
					for (l=0; l<n; l++) a(ll+1,l+1) -= a(icol+1,l+1) * dum;

					for (l=0; l<m; l++) b(ll+1,l+1) -= b(icol+1,l+1) * dum;
				}
		}
		for (l=n-1; l>=0; l--) {
			if (indxr[l] != indxc[l])
				for (k=0; k<n; k++)
					Swap(a(k+1,indxr[l]+1),a(k+1,indxc[l]+1));
		}
		delete [] ipiv;
		delete [] indxr;
		delete [] indxc;
		return a;
	}
    //-------------------------------------------------------------------------------------------
    void Matrix::set(const std::string &str){
		// actual row counter
		int rowsCount = 0;

		Real xx[1000] = {0};
		int v_size = 0; int cc = 0, dataSize = 0;
		// variable to store the start of a current vector
		std::string::size_type beg = 0, end = 0;
		while (end != std::string::npos) {
			// find next occurrence of a semicolon in string str
			end = str.find(';', beg);
			// parse first row into a vector v
			std::string data = str.substr(beg, end - beg);
			istringstream ss(data);
			while ((ss >> xx[dataSize])){
				v_size++;
				dataSize++;
			}
			if (cc == 0){
				(*this)._nColumns = v_size;
			}
			cc++;
			for (int i=dataSize; i<cc*(*this)._nColumns; i++)
				xx[i] = 0;
			dataSize = cc * (*this)._nColumns;
			// this check is necessary to parse the following two strings as the
			// same matrix: "1 0 1; ; 1 1; " and "1 0 1; 0 0 0; 1 1 0"
			if ((end != std::string::npos) || (v_size > 0))		rowsCount++;
			// update the starting position of the next vector in the parsed string
			beg = end + 1;
			v_size = 0;
		} // if ((end != std::string::npos) || (v.size > 0))
		(*this)._nRows = rowsCount;

		(*this)._data = new Real[(*this)._nRows * (*this)._nColumns];
		for (int i=0; i<(*this)._nRows * (*this)._nColumns; i++)
			(*this)._data[i] = xx[i];
	}
    //-------------------------------------------------------------------------------------------
	Real Matrix::operator () (int i) const{
		if ( !(*this).in_range(i) )		throw Exception("Matrix::operator () (int i): Invalid Index");
		return (*this)._data[i-1];
	}
	//-------------------------------------------------------------------------------------------
	Real & Matrix::operator () (int i){
		if ( !(*this).in_range(i) )		throw Exception("Matrix::operator () (int i): Invalid Index");
		return (*this)._data[i-1];
	}
    //-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator +(const Matrix &m)const{
		if(	(*this)._nColumns != m._nColumns || (*this)._nRows != m._nRows )
			throw Exception("Matrix::operator +(const Matrix &m): Invalid Matrix size");

		Matrix mat((*this)._nRows,(*this)._nColumns);
		const Real *pA = (*this).getData(), *pB = m.getData();
		Real *pC = mat.getData();
		int size = (*this)._nRows * (*this)._nColumns;
		for ( int i=0; i<size; i++ ){
			*pC = *pA + *pB;
			pA++;	pB++;	pC++;
		}
		return mat;
	}
	const Matrix Matrix::divideScalarbyMatrix(Real d)const{
		Matrix mat((*this)._nRows,(*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			if ((*this)._data[i] != 0)
				mat._data[i] = d / (*this)._data[i];
			else
				Exception("division by zero");
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator -(const Matrix &m)const{
		// Matrix Subtraction
		if(	(*this)._nColumns != m._nColumns || (*this)._nRows != m._nRows )
			throw Exception("Matrix::operator -(const Matrix &m): Invalid Matrix size");

		Matrix mat((*this)._nRows,(*this)._nColumns);
		const Real *pA = (*this).getData(), *pB = m.getData();
		Real *pC = mat.getData();
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ ){
			*pC = *pA - *pB;
			pA++;	pB++;	pC++;
		}
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator *(const Matrix & m)const{
		// Matrix Multiplication
		if ( (*this)._nColumns != m._nRows )
			throw Exception("Matrix::operator *(const Matrix &m): Invalid Matrix size");

		Matrix mat(_nRows, m._nColumns);
		Real *pC = mat._data;
		for ( int i = 0; i < _nRows; i++ )
			for ( int j = 0; j < m._nColumns; j++ ){
				Real sigma = 0;
				Real *pA = _data + i * _nColumns,
						*pB = m._data + j;

				for ( int k = 0; k < _nColumns; k++ ){
					sigma += *pA * *pB;
					pA++;
					pB += m._nColumns;
				}
				*pC++ = sigma;
			}
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	const ColumnVector Matrix::operator * (const ColumnVector &v)const{
		if ( (*this)._nColumns != v.getnRows() )
			throw Exception("Matrix::operator *(const ColumnVector &m): Invalid Matrix size");

		ColumnVector vec((*this)._nRows);
		Real *pC = vec.getData();
		for ( int i = 0; i < (*this)._nRows; i++ )
			for ( int j = 0; j < v.getnColumns(); j++ ){
				Real sigma = 0;
				Real *pA = (*this)._data + i * (*this)._nColumns;
				const Real *pB = v.getData() + j;

				for ( int k = 0; k < (*this)._nColumns; k++ ){
					sigma += *pA * *pB;
					pA++; pB += v.getnColumns();
				}
				*pC++ = sigma;
			}
		return vec;
	}
    //-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator * (const RowVector &v)const{
		if ( (*this)._nColumns != v.getnRows() )
			throw Exception("Matrix::operator *(const RowVector &m): Invalid Matrix size");

		Matrix vec((*this)._nRows, v.getnColumns());
		Real *pC = vec.getData();
		for ( int i = 0; i < (*this)._nRows; i++ )
			for ( int j = 0; j < v.getnColumns(); j++ ){
				Real sigma = 0;
				Real *pA = (*this)._data + i * (*this)._nColumns;
				const Real *pB = v.getData() + j;

				for ( int k = 0; k < (*this)._nColumns; k++ ){
					sigma += *pA * *pB;
					pA++; pB += v.getnColumns();
				}
				*pC++ = sigma;
			}
		return vec;
	}
    //-------------------------------------------------------------------------------------------
	const Matrix & Matrix::operator *=(const Matrix &m){
		if ( (*this)._nColumns != m.getnColumns() )
			throw Exception("Matrix::operator *=(const Matrix &m): Invalid Matrix size");
		return *this = *this * m;
	}
	//-------------------------------------------------------------------------------------------
	Matrix & Matrix::operator ++(){
		*this = *this + ONE;
		return *this;
	}
	Matrix   Matrix::operator ++(int){
		*this = *this + ONE;
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	Matrix & Matrix::operator --(){
		*this = *this - ONE;
		return *this;
	}

	Matrix   Matrix::operator --(int){
		*this = *this - ONE;
		return *this;
	}
    //-------------------------------------------------------------------------------------------
	const Matrix Matrix::ElementWiseMul(const Matrix &m)const{
		if ((*this)._nRows != m._nRows && (*this)._nColumns != m._nColumns)
			throw Exception("Matrix::ElementWiseMul(): The two Matrices have different sizes...");
		Matrix mat((*this)._nRows,(*this)._nColumns);
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++)
			mat._data[i] = (*this)._data[i] * m._data[i];
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const Matrix Matrix::ElementWiseDiv(const Matrix &m)const{
		if ((*this)._nRows != m._nRows && (*this)._nColumns != m._nColumns)
			throw Exception("Matrix::ElementWiseDiv(): The two Matrices have different sizes...");
		Matrix mat((*this)._nRows,(*this)._nColumns);
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++)
			if (::abs(m._data[i]) != 0)
				mat._data[i] = (*this)._data[i] / m._data[i];
			else
				throw Exception("Matrix::ElementWiseDiv(): division by zero...");
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator + (const Real & n)const{
		//! Matrix Scalar Multiplication mat * n
		Matrix mat((*this)._nRows,(*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			mat._data[i] = (*this)._data[i] + n;
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator - (const Real & n)const{
		//! Matrix Scalar Multiplication mat * n
		Matrix mat((*this)._nRows,(*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			mat._data[i] = (*this)._data[i] - n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator * (const Real & n)const{
		//! Matrix Scalar Multiplication mat * n
		Matrix mat((*this)._nRows,(*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			mat._data[i] = (*this)._data[i] * n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator / (const Real & n)const{
		//! Matrix Scalar divide mat / n
		if ( n == 0 )
			throw Exception("Matrix::operator / (const Real & n): n must not equal 0");
		Matrix mat((*this)._nRows,(*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			mat._data[i] = (*this)._data[i] / n;
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	const Matrix Matrix::operator^(const Real & n)const{
		// rise each element to the power n, mat(i)^n
		if ((*this)._nRows < 1 || (*this)._nColumns < 1)
			throw InvalidMatrixSize();
		Matrix mat((*this)._nRows,(*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			mat._data[i] = ::pow((*this)._data[i], n);
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	void Matrix::fastAdd(const Matrix &a, const Matrix &b, Matrix &c){
		if(	a.getnColumns() != b.getnColumns() || a.getnColumns() != c.getnColumns() ||
			a.getnRows() != b.getnRows() || a.getnRows() != c.getnRows() )
			throw InvalidMatrixSize();
		int nr = a.getnRows(), nc = a.getnColumns(), n = nr * nc;
		const Real *Pa = a.getData(), *Pb = b.getData();
		Real *Pc = c.getData();
		for (int i=0; i<n; i++){
			*Pc = *Pa + *Pb;
			Pc++;	Pa++;	Pb++;
		}
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::fastSub(const Matrix &a, const Matrix &b, Matrix &c){
		if(	a.getnColumns() != b.getnColumns() || a.getnColumns() != c.getnColumns() ||
			a.getnRows() != b.getnRows() || a.getnRows() != c.getnRows() )
			throw InvalidMatrixSize();
		int nr = a.getnRows(), nc = a.getnColumns(), n = nr * nc;
		const Real *Pa = a.getData(), *Pb = b.getData();
		Real *Pc = c.getData();
		for (int i=0; i<n; i++){
			*Pc = *Pa - *Pb;
			Pc++;	Pa++;	Pb++;
		}
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::operator ~()const{
		// Matrix Transpose
		if ((*this)._nRows && (*this)._nColumns){
			Matrix mat((*this)._nColumns,(*this)._nRows);
			for (int i=1; i<=(*this)._nColumns; i++)
				for (int j=1; j<=(*this)._nRows; j++)
					mat(i,j, (*this)(j,i));
			return mat;
		}else
			return Matrix();
	}
    //-------------------------------------------------------------------------------------------
	RowVector Matrix::getRow(int r) const{
		if (!(*this).row_in_range(r))
			throw Exception("Matrix::getRow(): Indexing out of range");
		RowVector row((*this)._nColumns);
		for ( int i=1; i<=(*this)._nColumns; i++)
			row(i) = getElement(r,i);	//_data[(i-1)*_nRows+(n-1)];  or   Row(i) = (*this) (i,n);
		return row;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::getRows(int r1, int r2) const{
		if (!(*this).row_in_range(r1) || !(*this).row_in_range(r2) || r1 > r2)
			throw Exception("Matrix::getRow(): Indexing out of range");
		Matrix mat(r2-r1+1, (*this)._nColumns);
		for (int i=(r1-1)*(*this)._nColumns, j=0; i<(r2-1)*(*this)._nColumns+(*this)._nColumns; i++, j++)
			mat._data[j] = (*this)._data[i];
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	ColumnVector Matrix::getColumn(int n) const{
		if (!(*this).col_in_range(n))
			throw Exception("Matrix::getColumn(): Indexing out of range");
		ColumnVector col((*this)._nRows);
		for ( int i=1; i<=(*this)._nRows; i++)
			col(i) = (*this).getElement(i,n);	//_data[(i-1)*_nColumns+(n-1)];  or   col(i) = (*this) (i,n);
		return col;
	}
    //-------------------------------------------------------------------------------------------
	Matrix Matrix::getColumns(int c1, int c2) const{
		if (!(*this).col_in_range(c1) || !(*this).col_in_range(c2) || c1 > c2)
			throw Exception("Matrix::getColumn(): Indexing out of range");
		int nc = c2 - c1, xx = c2-1;
		Matrix mat((*this)._nRows, nc+1);
		for (int m=c1-1, k=0; k<(*this)._nRows*(nc+1); k++, m++){
			mat._data[k] = (*this)._data[m];
			if (m == xx){
				m += (*this)._nColumns-nc-1;
				xx += (*this)._nColumns;
			}
		}
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::subMatrix(int r1, int r2, int c1, int c2) const{
		if (!(*this).row_in_range(r1) || !(*this).row_in_range(r2) || r1 > r2)
			throw Exception("Matrix::subMatrix(): Indexing of rows out of range");
		if (!(*this).col_in_range(c1) || !(*this).col_in_range(c2) || c1 > c2)
			throw Exception("Matrix::subMatrix(): Indexing of columns out of range");
		Matrix mat(r2-r1+1, (*this)._nColumns);
		mat = (*this).getRows(r1,r2);
		mat = mat.getColumns(c1,c2);
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	ColumnVector Matrix::getRowAsCol(int m) const{
		// Return a specific Row from a Matrix as a Column
		ColumnVector col((*this)._nColumns);
		for ( int i=1; i<=(*this)._nColumns; i++)
			col(i) = (*this).getElement(m,i);	//_data[(i-1)*_nRows+(n-1)];  or   Row(i) = (*this) (i,n);
		return col;
	}
	//-------------------------------------------------------------------------------------------
	RowVector Matrix::getColAsRow(int m) const{
		// Return a specific Row from a Matrix as a Column
		RowVector row((*this)._nRows);
		for ( int i=1; i<=(*this)._nRows; i++)
			row(i) = (*this).getElement(i,m);	//_data[(i-1)*_nRows+(n-1)];  or   Row(i) = (*this) (i,n);
		return row;
	}
	ColumnVector Matrix::getMatrixAsColumn() const{
		int n = (*this)._nRows*(*this)._nColumns;
		ColumnVector col(n);
		for(int i=0; i<n; i++)
			col(i+1) = (*this)._data[i];
		return col;
	}
	//-------------------------------------------------------------------------------------------
	void Matrix::setColAsRow(const ColumnVector &col, int i) {
		// Set a Column as a Row in a Matrix
		for ( int j=1; j<=(*this)._nColumns; j++)
			if ( i>0 && i <= (*this)._nRows && j>0 && j <= (*this)._nColumns)
				(*this)._data[(i-1) * (*this)._nColumns + (j-1)] = col(j);
	}
	//-------------------------------------------------------------------------------------------
	void Matrix::setCol(const ColumnVector &col, int c){
		for (int i=1; i<=(*this).getnRows(); i++){
			(*this)._data[(i-1) * (*this)._nColumns + (c-1) ] = col(i);
		}
	}
	//-------------------------------------------------------------------------------------------
	void Matrix::setRow(const ColumnVector &col, int r){
		for (int i=1; i<=(*this).getnColumns(); i++){
			(*this)._data[ (r-1) * (*this)._nColumns + (i-1) ] = col(i);
		}
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::swap_rows(int r1, int r2){
		if ( !(*this).row_in_range(r1) && !(*this).row_in_range(r2) )
			throw Exception("Matrix::swap_rows: Index out of range");
		Real tmp = 0;
		for (int i=1; i<=(*this)._nColumns; i++){
			tmp = (*this)(r1,i);
			(*this)(r1,i) = (*this)(r2,i);
			(*this)(r2,i) = tmp;
		}
	}
	//-------------------------------------------------------------------------------------------
	void Matrix::swap_cols(int c1, int c2){
		if ( !(*this).col_in_range(c1) && !(*this).col_in_range(c2) )
			throw Exception("Matrix::swap_cols: Indexes are out of range");
		Real tmp = 0;
		for (int i=1; i<=(*this)._nRows; i++){
			tmp = (*this)(i,c1);
			(*this)(i,c1) = (*this)(i,c2);
			(*this)(i,c2) = tmp;
		}
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::repmat(const ColumnVector &v, int n){
		if (n < 1)
			throw Exception("Matrix::Matrix::repmat: n should be >= 1.");
		Matrix mat(v.getnRows(), n);
		for (int i=1; i<=v.getnRows(); i++)
			for (int j=1; j<=n; j++)
				mat(i,j) = v(i);
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	Matrix Matrix::repmat(const RowVector &v, int m){
		if (m < 1)
			throw Exception("Matrix::Matrix::repmat: n should be >= 1.");
		Matrix mat(m, v.getnColumns());
		for (int i=1; i<=m; i++)
			for (int j=1; j<=v.getnColumns(); j++)
				mat(i,j) = v(j);
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::asDiagonal(){
		if ( !(_nColumns != 1 || _nRows != 1) )
			throw Exception("Matrix::Matrix::asDiagonal: the matrix must be a row or a column.");
		int n;
		Matrix m;

		if (_nColumns > 1){
			n = _nColumns;
			m = Matrix(n);
			for (int i = 1; i <= n; i++)	m(i, i) = (*this)(1, i);
		}

		if (_nRows > 1){
			n = _nRows;
			m = Matrix(n);
			for (int i = 1; i <= n; i++)	m(i, i) = (*this)(i, 1);
		}
		return m;
	}
	//-------------------------------------------------------------------------------------------
	void Matrix::pushBackROW(const RowVector &rv){
		if ((*this).getnColumns() != rv.getnColumns())
			throw Exception("Matrix::pushBackROW(const RowVector &rv): Invalid size");
		int oldSize = (*this)._nRows * (*this)._nColumns;
		Real *oldData;	oldData = new Real[oldSize];
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++)		oldData[i] = (*this)._data[i];
		(*this)._nRows += 1;
		delete [] (*this)._data;
		(*this)._data = new Real[(*this)._nRows * (*this)._nColumns];
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++){
			if (i<oldSize)
				(*this)._data[i] = oldData[i];
			else
				(*this)._data[i] = rv(i-oldSize+1);
		}
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::pushBackROW(){
		int oldSize = (*this)._nRows * (*this)._nColumns;
		Real *oldData;	oldData = new Real[oldSize];
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++)		oldData[i] = (*this)._data[i];
		(*this)._nRows += 1;
		delete [] (*this)._data;
		(*this)._data = new Real[(*this)._nRows*(*this)._nColumns];
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++){
			if (i<oldSize)
				(*this)._data[i] = oldData[i];
			else
				(*this)._data[i] = 0;
		}
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::pushBackCOL(const ColumnVector &cv){
		if ((*this).getnRows() != cv.getnRows())
			throw Exception("error: pushBackCOL: Invalid vector size...");
		int oldSize = (*this)._nRows * (*this)._nColumns;
		Real *oldData;	oldData = new Real[oldSize];
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++)		oldData[i] = (*this)._data[i];

		(*this)._nColumns += 1;			int cc = 0;
		delete [] (*this)._data;
		(*this)._data = new Real[(*this)._nRows * (*this)._nColumns];
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++){
			if (((i+1)%(*this)._nColumns)!=0)
				(*this)._data[i] = oldData[i-cc];
			if (i>0)
				if (((i+1)%(*this)._nColumns)==0){
					(*this)._data[i] =  cv(cc+1);
					cc++;
				}
		}
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::pushBackCOL(){
		int oldSize = (*this)._nRows * (*this)._nColumns;
		Real *oldData;	oldData = new Real[oldSize];
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++)		oldData[i] = (*this)._data[i];

		(*this)._nColumns += 1;			int cc = 0;
		delete [] (*this)._data;
		(*this)._data = new Real[(*this)._nRows*(*this)._nColumns];
		for (int i=0; i<(*this)._nRows*(*this)._nColumns; i++){
			if (((i+1)%(*this)._nColumns)!=0)
				(*this)._data[i] = oldData[i-cc];
			if (i>0)
				if (((i+1)%(*this)._nColumns)==0){
					(*this)._data[i] =  0;
					cc++;
				}
		}
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::popBackROW(){
		if ((*this)._nRows <= 1)
			throw Exception("error: popBackROW: can't popBack the rows of this matrix; nRows = 1");
		(*this)._nRows -= 1;
		int sz = (*this)._nRows * (*this)._nColumns;
		Real *oldData;	oldData = new Real[sz];
		for (int i=0; i<sz; i++)		oldData[i] = (*this)._data[i];
		delete [] (*this)._data;
		(*this)._data = new Real[sz];
		for (int i=0; i<sz; i++)
			(*this)._data[i] = oldData[i];
	}
	//-------------------------------------------------------------------------------------------
	void Matrix::popBackCOL(){
		if ((*this)._nColumns <= 1)
			throw Exception("error: popBackCOL: can't popBack the columns of this matrix; nColumns = 1");
		int ncol = (*this)._nColumns - 1;
		int sz = (*this)._nRows * ncol;
		Real *oldData;	oldData = new Real[sz];
		for (int i=0, j=0; i<sz; i++,j++){
			if (((j+1)%(*this)._nColumns)!=0)		oldData[i] = (*this)._data[j];
			if (i>0)
				if (((j+1)%(*this)._nColumns)==0)		i--;
		}
		(*this)._nColumns--;
		delete [] (*this)._data;
		(*this)._data = new Real[sz];
		for (int i=0; i<sz; i++)	(*this)._data[i] = oldData[i];
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::eraseROW(int ro){
		if (ro < 1 || ro > (*this)._nRows || (*this)._nRows == 1)
			throw Exception("error: eraseROW: Index out of range");
		int nro = (*this)._nRows - 1;
		int sz = (*this)._nColumns*nro;
		Real *oldData;	oldData = new Real[sz];
		int i=0;
		for (i=0; i<(ro-1)*(*this)._nColumns; i++)
			oldData[i] = (*this)._data[i];
		for (int k=i,j=ro*(*this)._nColumns; k<sz; j++,k++)
			oldData[k] = (*this)._data[j];

		(*this)._nRows--;
		delete [] (*this)._data;
		(*this)._data = new Real[sz];
		for (int i=0; i<sz; i++)	(*this)._data[i] = oldData[i];
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::eraseCOL(int col){
		if (col < 1 || col > (*this)._nColumns || (*this)._nColumns == 1)
			throw Exception("error: eraseCOL: Index out of range");
		int ncol = (*this)._nColumns - 1;
		int sz = (*this)._nRows * ncol, cc = 0;
		Real *oldData;	oldData = new Real[sz];
		for (int i=0, j=0; i<sz; i++,j++){
			if (((j+1)%(col+cc))!=0)
				oldData[i] = (*this)._data[j];
			if (i>0 || col == 1)
				if (((j+1)%(col+cc))==0){
					i--;
					cc += (*this)._nColumns;
				}
		}
		(*this)._nColumns--;
		delete [] (*this)._data;
		(*this)._data = new Real[sz];
		for (int i=0; i<sz; i++)	(*this)._data[i] = oldData[i];
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::insertROW(int r, const RowVector &rv){
		if ((*this)._nColumns != rv.getnColumns() || r < 1 || r > (*this)._nRows+1)
			throw Exception("error: insertROW: Index out of range");
		int nro = (*this)._nRows + 1;
		int sz = (*this)._nColumns * nro;
		Real *newData;	newData = new Real[sz];
		int i=0, k=0, j=0;
		for (i=0; i<(r-1)*(*this)._nColumns; i++)
			newData[i] = (*this)._data[i];
		for ( k=i, j=0; k<i+(*this)._nColumns; j++,k++)
			newData[k] = rv(j+1);
		for (int m=k,j=(r-1)*(*this)._nColumns; m<sz; j++,m++)
			newData[m] = (*this)._data[j];

		(*this)._nRows++;
		delete [] (*this)._data;
		(*this)._data = new Real[sz];
		for (int i=0; i<sz; i++)	(*this)._data[i] = newData[i];
	}
    //-------------------------------------------------------------------------------------------
	void Matrix::insertCOL(int col, const ColumnVector &cv){
		if ((*this)._nRows != cv.getnRows() || col < 1 || col > (*this)._nColumns+1)
			throw Exception("error: insertCOL: Invalid vector size...");
		int ncol = (*this)._nColumns + 1;
		int sz = (*this)._nRows * ncol, cc = 0, count = 0;
		Real *newData;	newData = new Real[sz];
		for (int i=0, j=0; i<sz; i++,j++){
			if (((j+1)%(col+cc))!=0)
				newData[i] = (*this)._data[j];
			if (i>0 || col == 1)
				if (((j+1)%(col+cc))==0){
					newData[i] = cv(count+1);
					cc += (*this)._nColumns;
					count++;
					j--;
				}
		}
		(*this)._nColumns++;
		delete [] (*this)._data;
		(*this)._data = new Real[sz];
		for (int i=0; i<sz; i++)	(*this)._data[i] = newData[i];
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::HorizontallyConcatenated(const Matrix& mat){
		int i, j;
		if (_nRows != mat.getnRows())
			throw Exception("Matrix::HorizontallyConcatenated: Number of rows must be the same for horizontal concatination");
		int c = _nColumns + mat.getnColumns();
		Matrix res(_nRows, c);
		for (i = 1; i <= _nRows; i++)
			for(j=1; j<=c; j++)
				if (j <= _nColumns)	res(i, j) = (*this)(i, j);
				else res(i, j) = mat(i, j-_nColumns);
		return res;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::HorizontallyConcatenated(const ColumnVector& mat){
		int i, j;
		if (_nRows != mat.getnRows())
			throw Exception("Matrix::HorizontallyConcatenated: Number of rows must be the same for horizontal concatination");
		int ct = (*this)._nColumns;
		int c = ct + mat.getnColumns();
		Matrix res(_nRows, c);
		for (i = 1; i <= _nRows; i++){
			res(i,c) = mat(i);
			for(j=1; j<=(*this)._nColumns; j++)
				res(i,j) = (*this)(i,j);
		}
		return res;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::VerticallyConcatenated(const Matrix& mat){
		if (_nColumns != mat.getnColumns())
			throw Exception("Matrix::VerticallyConcatenated: Number of Col. must be the same for horizontal concatination");
		int i, j;
		int r = _nRows + mat.getnRows();
		Matrix res(r, _nColumns);
		for (i = 1; i <= _nColumns; i++)
			for (j = 1; j <= r; j++)
				if (j <= _nRows)	res(j, i) = (*this)(j, i);
				else res(j, i) = mat(j-_nRows, i);
		return res;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::VerticallyConcatenated(const RowVector& mat){
		int c = (*this)._nColumns;
		if (c != mat.getnColumns())
			throw Exception("Matrix::VerticallyConcatenated: Number of Col. must be the same for horizontal concatination");
		int rt = (*this)._nRows, i, j;
		int r = rt + mat.getnRows();
		Matrix res(r,c);
		for(i=1; i<=c; i++){
			for(j=1; j<=(*this)._nRows; j++){
				res(j,i) = (*this)(j,i);
				if(j==1)
					res(j+rt,i) = mat(i);
			}
		}
		return res;
	}
    //************************************
    // Method:      reshape
    // FullName:    myMat::Matrix::reshape
    // Access:      public
    // Returns:     void
    // Parameter:   int r, c
    // Description: RESHAPE(M,N) returns the M-by-N matrix whose elements
    //				are taken column-wise from this matrix object.  An error
    //				results if X does not have M*N elements.
    //************************************
	Matrix Matrix::reshape (int r, int c){
		// Resizes the matrix, initializes new values to zero.
		if ( (*this)._nRows * (*this)._nColumns != r * c)
			throw InvalidIndex();

		Matrix temp(r, c);
		int i, j, ii = 1, jj = 1;
		for (j = 1; j <= (*this)._nColumns; j++) {
			for (i = 1; i < (*this)._nRows; i++) {
				temp(ii++, jj) = (*this)(i, j);
				if (ii == r) {
					jj++;
					ii = 1;
				}
			}
		}
		return temp;
	}
    //-------------------------------------------------------------------------------------------
	void  Matrix::resize (int r, int c){
		// Resizes the matrix, initializes new values to zero.
		// 	if ( r<1 || c<1 || r<_nRows || c<_nColumns )
		if ( r<1 || c<1 )
			throw InvalidIndex();
		int old_nRows = (*this)._nRows, old_nColumns = (*this)._nColumns;
		Real *nData;		nData = new Real[old_nRows*old_nColumns];
		for ( int i=0; i<old_nRows * old_nColumns; i++ )
			nData[i] = (*this)._data[i];

		(*this)._nRows = r;		(*this)._nColumns = c;
		delete [] (*this)._data;	(*this)._data = new Real[r*c];
		for ( int i=0; i<(*this)._nRows * (*this)._nColumns; i++ )
			if (i < old_nRows * old_nColumns)
				(*this)._data[i] = nData[i];
			else
				(*this)._data[i] = 0;
		delete [] nData;
	}
	//-------------------------------------------------------------------------------------------
	void  Matrix::resize (int r, int c, Real val){
		// Resizes the matrix, initializes new values to zero.
		// 	if ( r<1 || c<1 || r<_nRows || c<_nColumns )
		if ( r<1 || c<1 )
			throw InvalidIndex();

		(*this)._nRows = r;		(*this)._nColumns = c;
		delete [] (*this)._data;
		(*this)._data = new Real[r*c];
		for ( int i=0; i<(*this)._nRows * (*this)._nColumns; i++ )
			(*this)._data[i] = val;
	}
	//-------------------------------------------------------------------------------------------
	void  Matrix::reConstruct (const Matrix &m){
		// reConstruct the Matrix base to a given one.
		if ((*this)._data)
			delete [] (*this)._data;
		(*this)._nRows = m._nRows;
		(*this)._nColumns = m._nColumns;
		(*this)._data = new Real[(*this)._nRows * (*this)._nColumns];
		for ( int i=0; i<m._nRows * m._nColumns; i++ )
			(*this)._data[i] = m._data[i];
	}
	//-------------------------------------------------------------------------------------------
	void Matrix::ReverseElements()
	/*!
		@brief Reversing in place
	*/
	{
		int n = size(); Real* x = getData(); Real* rx = x + n;
		n /= 2;
		while (n--) { Real t = *(--rx); *rx = *x; *(x++) = t; }
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::ReverseElements(const Matrix & mat)
	/*!
		@brief Reversing into a new matrix
	*/
	{
		Matrix n_mat(mat.getnRows(),mat.getnColumns());
		int n = mat.size(); Real* rx = n_mat.getData() + n; const Real* x = mat.getData();
		while (n--) *(--rx) = *(x++);
		return n_mat;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::abs()const{
		Matrix mat((*this)._nRows, (*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			mat._data[i] = ::abs((*this)._data[i]);
		return mat;
	}
	//-------------------------------------------------------------------------------------------

	Matrix Matrix::Sqr()const{
		Matrix mat((*this)._nRows, (*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			mat._data[i] = (*this)._data[i] * (*this)._data[i];
		return mat;
	}
	//-------------------------------------------------------------------------------------------

	Matrix Matrix::Sqrt()const{
		Matrix mat((*this)._nRows, (*this)._nColumns);
		for ( int i=0; i<(*this)._nRows*(*this)._nColumns; i++ )
			mat._data[i] = ::sqrt((*this)._data[i]);
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::exp()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ )
			mat._data[i] = ::exp((*this)._data[i]);
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	Matrix Matrix::ceil()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ )
			mat._data[i] = ::ceil((*this)._data[i]);
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::floor()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ )
			mat._data[i] = ::floor((*this)._data[i]);
		return mat;
	}
	//-------------------------------------------------------------------------------------------

	Matrix Matrix::round()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ ){
			mat._data[i] = ((*this)._data[i] > ZERO) ? ::floor((*this)._data[i] + static_cast<Real>(0.5)) : ::ceil((*this)._data[i] - static_cast<Real>(0.5));
		}
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	Matrix Matrix::log()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ )
			mat._data[i] = ::log((*this)._data[i]);
		return mat;
	}
	//-------------------------------------------------------------------------------------------

	Matrix Matrix::log10()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ )
			mat._data[i] = ::log10((*this)._data[i]);
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::log1p()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ )
			mat._data[i] = ::log(ONE + (*this)._data[i]);
		return mat;
	}
    //-------------------------------------------------------------------------------------------
	Matrix Matrix::log2()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ )
			mat._data[i] = ::log((*this)._data[i])/::log(static_cast<Real>(2));
		return mat;
	}
	//-------------------------------------------------------------------------------------------
/*
	Matrix Matrix::logb()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		Matrix mat(r,c);
		for ( int i=0; i<r*c; i++ )
			mat._data[i] = ::logb((*this)._data[i]);
		return mat;
	}
*/
	//-------------------------------------------------------------------------------------------
	Real Matrix::trace()const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		if (r != c)		throw Exception("error: trace(): matrix must be square sized...");
		Real s = 0;
		for (int i=1; i<=r; i++)	s += (*this)(i,i);
		return s;
	}
    //-------------------------------------------------------------------------------------------
	ColumnVector Matrix::diag(int k)const{
		int r = (*this).getnRows(), c = (*this).getnColumns();
		if (r != c)		throw Exception("error: diag(): matrix must be square sized...");
		if (::abs(k) >= r)	throw Exception("error: diag(): requested diagonal out of bounds");
		ColumnVector v;
		if (k == 0){
			v = ColumnVector(r);
			for (int i=1; i<=r; i++)	v(i) = (*this)(i,i);
		}else if (k > 0){
			v = ColumnVector(r-k);
			for (int i=1; i<=r-k; i++)
				v(i) = (*this)(i,i+k);
		}else{
			v = ColumnVector(r-::abs(k));
			for (int i=1; i<=r-::abs(k); i++)	v(i) = (*this)(i+::abs(k),i);
		}
		return v;
	}
	ColumnVector Matrix::diag(const Matrix &m, int k){return m.diag(k);	}
    //-------------------------------------------------------------------------------------------
	bool Matrix::isSquare()const{
		if ((*this)._nRows == (*this)._nColumns)	return true;
		return false;
	}
	//-------------------------------------------------------------------------------------------
	Real Matrix::norm()const{
		ColumnVector s((*this).getnColumns());
		s = SVDcmp(*this);
		Real max = s(1);
		for (int i=2; i<=(*this).getnColumns(); i++)
			if (s(i) > max)
				max = s(i);
		return max;
	}
 	Real Matrix::normINF()const{
		ColumnVector s = sum(2);
		Real max = s(1);
		for (int i=2; i<=(*this).getnColumns(); i++)
			if (s(i) > max)
				max = s(i);
		return max;
	}
   //-------------------------------------------------------------------------------------------
	ColumnVector Matrix::sum(int dim)const{
		if ((dim != 1) || (dim != 2))
			throw Exception("error: Matrix::sum(dim): dimension need to be 1 or 2");
		ColumnVector v;
		if (dim == 1){
			v = ColumnVector((*this)._nColumns);
			for (int i = 1; i <= (*this)._nColumns; i++)
				v(i) = BaseMatrix::sum((*this).getColumn(i));
		}else{
			v = ColumnVector((*this)._nRows);
			for (int i = 1; i <= (*this)._nRows; i++)
				v(i) = BaseMatrix::sum((*this).getRow(i));
		}
		return v;
	}
	ColumnVector Matrix::sum(const Matrix &m, int dim)	{return m.sum(dim);	}
    //-------------------------------------------------------------------------------------------
    Matrix Matrix::PseudoInverse(math::Matrix &A){
		int m = A.getnRows();
		int n = A.getnColumns();
		math::Matrix X;
		Real s = 0;

		if (n>m){
			Matrix Aa = ~A;
			return ~PseudoInverse(Aa);
		}

		math::ColumnVector S(n);
		math::Matrix U(A), V(n,n);
		S = math::SVDcmp(A, U, V);
		if (m > 1){
			math::ColumnVector s(n);
			s = S;
			Real tol = MAX(m,n) * ::pow(static_cast<Real>(2), static_cast<Real>(-22));
			Real r = 0;
			for (int i=1; i<=n; i++)
				r += s(i) > tol;
			if (r == 0)
				return math::Matrix(n,m);
			else{
				int rr = (int)r;
				math::Matrix s2(rr,rr), v(n,rr), u(m,rr);
				for (int i=1; i<=rr; i++)
					s2(i,i) = ONE/s(i);
				for (int i=1; i<=rr; i++){
					for (int j=1; j<=n; j++)
						v(j,i) = V(j,i);
					for (int j=1; j<=m; j++)
						u(j,i) = U(j,i);
				}
				X = math::Matrix(n,m);
				X = v * s2 * ~u;
			}
			return X;
		}else if (m == 1)
			s = S(1);
		else
			s = 0;

		Real tol = MAX(m,n) * ::pow(static_cast<Real>(2), static_cast<Real>(-22));
		Real r = 0;
		r = s > tol;
		if (r == 0)
			return math::Matrix(n,m);
		else{
			int rr = (int)r;
			math::Matrix v(n,rr), u(m,rr);
			s = ONE/S(1);
			for (int i=1; i<=rr; i++){
				for (int j=1; j<=n; j++)
					v(j,i) = V(j,i);
				for (int j=1; j<=m; j++)
					u(j,i) = U(j,i);
			}
			X = math::Matrix(n,m);
			X = s * v * ~u;
		}
		return X;
    }
    //-------------------------------------------------------------------------------------------
	int Matrix::rank(Real tol){
		int rows = (*this)._nRows;
		int cols = (*this)._nColumns;
		if ((rows == 0) || (cols == 0))
			return 0;

		ColumnVector sing_val = SVDcmp(*this);
		// Calculate default tolerance
		if (tol < ZERO)
			tol = EPSilon * sing_val(1) * (rows > cols ? rows : cols);
		// Count number of nonzero singular values
		int r = 1;
		while ((r <= sing_val.size()) && (sing_val(r) > tol))
			r++;
		return r;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::eye(int n){
		if (n <= 0)
			throw Exception("error:Matrix::eye(n): the value of n must be n>=1.");
		if (n==1)
			return Matrix(n) = 1;
		Matrix m(n);
		for (int i=1; i<=n; i++)
			m(i,i) = 1;
		return m;
	}
    //-------------------------------------------------------------------------------------------
	Real Matrix::determinant()const{
		if (!isSquare())
			throw Exception("determinant()::The matrix is NOT square matrix");
		int n = (*this).getnRows(), i, j, i_count, j_count, count=0;
		Matrix array(n-1);
		Real det = 0;

		if(n<1){
			puts("Error");
			exit(1);
		}
		if(n == 1) return (*this)(1,1);
		if(n == 2) return ((*this)(1,1)*(*this)(2,2) - (*this)(1,2)*(*this)(2,1));

		for(count=1; count<=n; count++){
			//Creating array of Minors
			i_count = 1;
			for(i=2; i<=n; i++){
				j_count = 1;
				for(j=1; j<=n; j++){
					if(j == count) continue;
					array(i_count,j_count) = (*this)(i,j);
					j_count++;
				}
				i_count++;
			}
			det += ::pow(-ONE, count+1) * (*this)(1,count) * determinant(array);	//Recursive call
		}
		return det;
	}
    //-------------------------------------------------------------------------------------------
	Real Matrix::determinantLU()const{
		if (!isSquare())
			throw Exception("determinantLU()::The matrix is NOT square matrix");
		Matrix temp = (*this);
		Real d = ONE;
		int j, *indx;

		indx = new int[_nRows];

		ludcmp(temp, indx, d);
		for (j = 0; j < _nRows; j++) d *= temp(j+1,j+1);

		delete [] indx;

		return d;
	}
    //-------------------------------------------------------------------------------------------
	Real Matrix::detCholesky()const{
		if (!isSquare())
			throw Exception("detCholesky()::The matrix is NOT square matrix");
		ColumnVector p;
		Matrix A = cholesky(*this,p);
		Real d=1; int i;

		for (i = 0; i < _nColumns; i++)  d *= A(i+1,i+1);
		return d * d;
	}
    //-------------------------------------------------------------------------------------------
	ColumnVector Matrix::solve(const ColumnVector &b)const{
		if (!isSquare())
			throw Exception("solve()::The matrix is NOT square matrix");
		Matrix temp(*this);
		int *indx, n = temp.getnRows();
		Real d;
		ColumnVector x = b;
		indx = new int[n];

		ludcmp (temp, indx, d);
		lubksb (temp, indx, x);

		delete [] indx;
		return x;
	}
	//-------------------------------------------------------------------------------------------
	Matrix Matrix::LU()const{
		if (!isSquare())
			throw Exception("LU()::The matrix is NOT square matrix");
		Matrix temp(*this);
		int *indx, dim = temp.getnRows();
		Real d = ONE;
		indx = new int[dim];
		ludcmp (temp, indx, d);
		// 	for (i = 1; i <= dim; i++) d *= temp(i,i);
		// 	determinant = d;
		delete [] indx;
		return temp;
	}
    //-------------------------------------------------------------------------------------------
	Matrix random(int m, int n){
		if (m<1 || n<1)
			throw Exception("random(m,n): m and n must be >= 1");
		Matrix mat(m,n);
		for (int i=1; i<=m*n; i++)
			random(mat(i));
		return mat;
	}
	Matrix magic(int m, int n) { return random(m,n);	}
    //-------------------------------------------------------------------------------------------
	DiagonalMatrix::DiagonalMatrix(int n, const ColumnVector &v) : Matrix(v.getnRows(), n){
		for (int i = 1; i <= (v.getnRows()<n ? v.getnRows() : n); i++) (*this)(i) = v(i);
	}
    //-------------------------------------------------------------------------------------------
	Matrix inv(const Matrix & mat){
		if (!mat.isSquare())
			throw Exception("Error using inv: Input must be a square matrix.");
		int N = mat.getnRows();
		int LDA = N;
		double det[2], *work = new double[N];
		int info, *ipvt = new int[N], job;
		double *a = new double[LDA*N];
		for (int i = 1; i <= LDA; i++ ){
			for (int j = 1; j <= LDA; j++ )	{
				a[(i-1)+(j-1)*LDA] = mat.getData()[(i-1)*LDA+(j-1)];
			}
		}
		info = dgefa ( a, LDA, N, ipvt );
		if ( info != 0 ){
			cout << "  Error!  The matrix is nearly singular!\n";
			return Matrix();
		}
		job = 11;
		dgedi ( a, LDA, N, ipvt, det, work, job );
		double *b = new double[LDA*N];
		for (int i = 1; i <= LDA; i++ ){
			for (int j = 1; j <= LDA; j++ )	{
				b[(i-1)*LDA+(j-1)] = a[(i-1)+(j-1)*LDA];
			}
		}
		Matrix res(N, b);
		delete[] ipvt;
		delete[] work;
		delete[] a;
		delete[] b;
		return res;
	}
    //-------------------------------------------------------------------------------------------
	Real det(const Matrix & mat){
		return Matrix::detCholesky(mat);
	}

};
