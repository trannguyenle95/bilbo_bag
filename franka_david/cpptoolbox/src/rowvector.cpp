
#include "matrix.hpp"

namespace math{


	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//-----------------------------------------RowVector-----------------------------------------
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	RowVector::RowVector(Real x, Real y, Real z){
		*this = RowVector(3);
		(*this)._data[0] = x;	(*this)._data[1] = y;	(*this)._data[2] = z;
	}
	//-------------------------------------------------------------------------------------------
	void RowVector::set(const std::string &str){
		// actual row counter
		int rowsCount = 0;

		Real xx[1000] = { 0 };
		int v_size = 0; int cc = 0, dataSize = 0;
		RowVector v;
		// variable to store the start of a current vector
		std::string::size_type beg = 0;
		std::string::size_type end = 0;
		if (end != std::string::npos) {
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
			for (int i = dataSize; i<cc*(*this)._nColumns; i++)
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

		if ((*this)._nRows != 1)
			throw Exception("error: RowVector::set: Invalid String; no. of rows != 1");

		(*this)._data = new Real[(*this)._nRows * (*this)._nColumns];
		for (int i = 0; i<(*this)._nRows * (*this)._nColumns; i++)		(*this)._data[i] = xx[i];
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::operator +(const RowVector &m)const{
		if ((*this)._nColumns != m._nColumns)	throw InvalidMatrixSize();

		RowVector v((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)	v._data[i] = (*this)._data[i] + m._data[i];
		return v;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::operator -(const RowVector &m)const{
		if ((*this)._nColumns != m._nColumns)	throw InvalidMatrixSize();

		RowVector v((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)	v._data[i] = (*this)._data[i] - m._data[i];
		return v;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::operator *(const Matrix	  &m)const{
		if ((*this)._nColumns != m.getnRows())
			throw InvalidMatrixSize();

		RowVector rv(m.getnColumns());
		Real *pC = rv.getData();
		for (int i = 0; i < (*this)._nRows; i++)
			for (int j = 0; j < m.getnColumns(); j++){
			Real sigma = 0;
			Real *pA = (*this)._data + i * (*this)._nColumns;
			const Real *pB = m.getData() + j;

			for (int k = 0; k < (*this)._nColumns; k++){
				sigma += *pA * *pB;
				pA++; pB += m.getnColumns();
			}
			*pC++ = sigma;
			}
		return rv;
	}
	//-------------------------------------------------------------------------------------------
	const Real RowVector::operator * (const ColumnVector & rv)const{
		if ((*this)._nColumns != rv.getnRows())
			throw InvalidMatrixSize();

		Real m = 0;
		for (int i = 1; i <= (*this)._nColumns; i++)
			m += (*this)(i)* rv(i);
		return m;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector & RowVector::operator *=(const Matrix &v){
		if ((*this)._nColumns != v.getnColumns())
			throw InvalidMatrixSize();
		return *this = *this * v;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::ElementWiseMul(const RowVector &m)const{
		if ((*this)._nColumns != m._nColumns)
			throw Exception("error: ElementWiseMul: The two Matrices have different sizes...");
		RowVector mat((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)
			mat._data[i] = (*this)._data[i] * m._data[i];
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::ElementWiseDiv(const RowVector &m)const{
		if ((*this)._nColumns != m._nColumns)
			throw Exception("error: ElementWiseDiv: The two Matrices have different sizes...");
		RowVector mat((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)
			if (::abs(m._data[i]) != 0)
				mat._data[i] = (*this)._data[i] / m._data[i];
			else
				throw Exception("error: ElementWiseDiv: division by zero...");
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::operator + (const Real & n)const{
		// Matrix Scalar Multiplication mat * n
		RowVector mat((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)
			mat._data[i] = (*this)._data[i] + n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::operator - (const Real & n)const{
		// Matrix Scalar Multiplication mat * n
		RowVector mat((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)
			mat._data[i] = (*this)._data[i] - n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::operator * (const Real & n)const{
		// Matrix Scalar Multiplication mat * n
		RowVector mat((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)
			mat._data[i] = (*this)._data[i] * n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::operator / (const Real & n)const{
		// Matrix Scalar divide mat / n
		if (n == 0)
			throw InvalidIndex();
		RowVector mat((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)
			mat._data[i] = (*this)._data[i] / n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const RowVector RowVector::operator^(const Real & n)const{
		// rise each element to the power n, mat(i)^n
		if ((*this)._nRows < 1 || (*this)._nColumns < 1)
			throw InvalidMatrixSize();
		RowVector mat((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)
			mat._data[i] = ::pow((*this)._data[i], n);
		return mat;
	}
	const RowVector RowVector::divideScalarbyRow(Real d)const{
		RowVector mat((*this)._nColumns);
		for (int i = 0; i<(*this)._nColumns; i++)
			if ((*this)._data[i] != 0)
				mat._data[i] = d / (*this)._data[i];
			else
				Exception("division by zero");
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	ColumnVector RowVector::operator ~()const{
		// Transpose Row Vector to Column Vector
		ColumnVector mat((*this).getnColumns());
		for (int i = 1; i <= (*this).getnColumns(); i++)
			mat(i) = (*this)(i);
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	ColumnVector RowVector::Transpose()const{ return ~(*this); }
	//-------------------------------------------------------------------------------------------
	void RowVector::ReverseElements()
		/*!
		@brief Reversing in place
		*/
	{
		int n = size(); Real* x = getData(); Real* rx = x + n;
		n /= 2;
		while (n--) { Real t = *(--rx); *rx = *x; *(x++) = t; }
	}
	//-------------------------------------------------------------------------------------------
	RowVector RowVector::ReverseElements(const RowVector & mat)
		/*!
		@brief Reversing into a new matrix
		*/
	{
		RowVector n_mat(mat.getnColumns());
		int n = mat.size(); Real* rx = n_mat.getData() + n; const Real* x = mat.getData();
		while (n--) *(--rx) = *(x++);
		return n_mat;
	}
	//-------------------------------------------------------------------------------------------
	// Calculates an absolute value of each matrix element.
	RowVector RowVector::abs()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = ::abs((*this)._data[i]);
		return v;
	}
	// Computes a square value of each matrix element.
	RowVector RowVector::Sqr()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = (*this)._data[i] * (*this)._data[i];
		return v;
	}
	// Computes a square root of each matrix element.
	RowVector RowVector::Sqrt()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = ::sqrt((*this)._data[i]);
		return v;
	}
	// Computes an exponent of each matrix element.
	RowVector RowVector::exp()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = ::exp((*this)._data[i]);
		return v;
	}
	// Changes each element to the smallest integer greater than or equal to that element.
	RowVector RowVector::ceil()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = ::ceil((*this)._data[i]);
		return v;
	}
	// Changes each element to the largest integer less than or equal to that element.
	RowVector RowVector::floor()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = ::floor((*this)._data[i]);
		return v;
	}
	// Rounds each element to the nearest integral value.
	RowVector RowVector::round()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++){
			v._data[i] = ((*this)._data[i] > ZERO) ? ::floor((*this)._data[i] + static_cast<Real>(0.5)) : ::ceil((*this)._data[i] - static_cast<Real>(0.5));
		}
		return v;
	}
	// Computes the natural logarithm of each element in the matrix.
	RowVector RowVector::log()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = ::log((*this)._data[i]);
		return v;
	}
	// Computes the logarithm of each element in the matrix to base 10.
	RowVector RowVector::log10()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = ::log10((*this)._data[i]);
		return v;
	}
	// Computes the natural logarithm of each element in the matrix plus 1.
	RowVector RowVector::log1p()const{
		int c = (*this).getnColumns();
		RowVector v(c);
		for (int i = 0; i<c; i++)
			v._data[i] = ::log(ONE + (*this)._data[i]);
		return v;
	}
	// Computes the logarithm of each element in the matrix to base 2. log(n)/log(2) is log2.
	RowVector RowVector::log2()const{
		int r = (*this).getnColumns();
		RowVector v(r);
		for (int i = 0; i<r; i++)
			v._data[i] = ::log((*this)._data[i]) / ::log(static_cast<Real>(2));
		return v;
	}
	// Computes the exponent of each element in the matrix.
	/*
	RowVector RowVector::logb()const{
	int r = (*this).getnColumns();
	RowVector v(r);
	for ( int i=0; i<r; i++ )
	v._data[i] = ::logb((*this)._data[i]);
	return v;
	}
	*/
	//-------------------------------------------------------------------------------------------
	Real RowVector::norm(int p)const{
		Real sum = 0;
		for (int i = 0; i<(*this)._nRows; i++)
			sum += ::pow(fabs((*this)._data[i]), static_cast<Real>(p));
		return std::pow(sum, ONE / static_cast<Real>(p));
	}
	//-------------------------------------------------------------------------------------------
	Real RowVector::norm(const std::string &s)const{
		if (s != "fro")
			throw Exception("error:  ColumnVector::norm: Unrecognised norm...");
		return norm(2);
	}
	//-------------------------------------------------------------------------------------------
	void RowVector::pushBack(const Real e){
		Real *x;		x = new Real[(*this)._nColumns];
		for (int i = 0; i<(*this)._nColumns; i++)	x[i] = (*this)._data[i];
		(*this)._nColumns += 1;
		delete[](*this)._data;	(*this)._data = new Real[(*this)._nColumns];
		for (int i = 0; i<(*this)._nColumns - 1; i++)	(*this)._data[i] = x[i];
		(*this)._data[(*this)._nColumns - 1] = e;
	}
	//-------------------------------------------------------------------------------------------
	void RowVector::popBack(){
		(*this)._nColumns -= 1;
		Real *x;		x = new Real[(*this)._nColumns];
		for (int i = 0; i<(*this)._nColumns; i++)	x[i] = (*this)._data[i];
		delete[](*this)._data;	(*this)._data = new Real[(*this)._nColumns];
		for (int i = 0; i<(*this)._nColumns; i++)	(*this)._data[i] = x[i];
	}
	//-------------------------------------------------------------------------------------------
	RowVector RowVector::HorizontallyConcatenated(const RowVector& mat){
		int ct = (*this)._nColumns, i;
		int cm = mat.getnColumns();
		RowVector res(ct + cm);
		for (i = 1; i <= ct; i++)
			res(i) = (*this)(i);
		for (i = 1; i <= cm; i++)
			res(i + ct) = mat(i);
		return res;
	}
	//-------------------------------------------------------------------------------------------
	Matrix RowVector::VerticallyConcatenated(const RowVector& mat){
		int c = (*this)._nColumns;
		if (c != mat.getnColumns())
			throw Exception("error: insertCOL: Invalid vector size...");
		int i;
		Matrix res(2, c);
		for (i = 1; i <= c; i++){
			res(1, i) = (*this)(i);
			res(2, i) = mat(i);
		}
		return res;
	}
	//-------------------------------------------------------------------------------------------
	Matrix RowVector::VerticallyConcatenated(const Matrix& mat){
		int ct = (*this)._nColumns;
		if (ct != mat.getnColumns())
			throw Exception("error: insertCOL: Invalid vector size...");
		int i, j;
		int rm = mat.getnRows();
		Matrix res(rm + 1, ct);
		for (i = 1; i <= rm; i++)
			for (j = 1; j <= ct; j++){
			res(i + 1, j) = mat(i, j);
			if (i == 1) res(i, j) = (*this)(j);
			}
		return res;
	}
	//-------------------------------------------------------------------------------------------
	void RowVector::insert(const Real & n, int k){
		if (k<0 || k>(*this)._nColumns + 1)
			throw Exception("error:  ColumnVector::insert: Invalid Index...");

		Real *x;		x = new Real[(*this)._nColumns];
		for (int i = 0; i<(*this)._nColumns; i++)	x[i] = (*this)._data[i];
		(*this)._nColumns += 1;
		delete[](*this)._data;	(*this)._data = new Real[(*this)._nColumns];
		for (int i = 0, j = 0; i<(*this)._nColumns; i++, j++)
			if ((i + 1) == k){
			(*this)._data[i] = n;
			j--;
			}
			else
				(*this)._data[i] = x[j];
	}
	//-------------------------------------------------------------------------------------------
	void RowVector::swap(int i, int j){
		Real temp = (*this)(i);
		(*this)(i) = (*this)(j);
		(*this)(j) = temp;
	}
	//-------------------------------------------------------------------------------------------
	RowVector RowVector::FromTo(int i, int j)const{
		if (i < 1 || i > j || j > (*this)._nColumns)
			throw InvalidIndex();
		RowVector v(j - i + 1);
		for (int k = 1; k <= (j - i + 1); k++)
			v(k) = (*this)(i + k - 1);
		return v;
	}


};