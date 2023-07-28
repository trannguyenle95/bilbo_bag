
#include "matrix.hpp"

namespace math{


	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//----------------------------------------ColumnVector---------------------------------------
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	ColumnVector::ColumnVector(Real x, Real y, Real z){
		*this = ColumnVector(3);
		(*this)._data[0] = x;	(*this)._data[1] = y;	(*this)._data[2] = z;
	}
	//-------------------------------------------------------------------------------------------
	void ColumnVector::set(const std::string &str){
		// actual row counter
		int colsCount = 0;

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
				(*this)._nRows = v_size;
			}
			cc++;
			for (int i = dataSize; i<cc*(*this)._nColumns; i++)
				xx[i] = 0;
			dataSize = cc * (*this)._nColumns;
			// this check is necessary to parse the following two strings as the
			// same matrix: "1 0 1; ; 1 1; " and "1 0 1; 0 0 0; 1 1 0"
			if ((end != std::string::npos) || (v_size > 0))		colsCount++;
			// update the starting position of the next vector in the parsed string
			beg = end + 1;
			v_size = 0;
		} // if ((end != std::string::npos) || (v.size > 0))
		(*this)._nColumns = colsCount;

		if ((*this)._nColumns != 1)
			throw Exception("error: ColumnVector::set: Invalid String; no. of rows != no. of columns");

		(*this)._data = new Real[(*this)._nRows * (*this)._nColumns];
		for (int i = 0; i<(*this)._nRows * (*this)._nColumns; i++)
			(*this)._data[i] = xx[i];
	}
	//-------------------------------------------------------------------------------------------
	const ColumnVector ColumnVector::operator +(const ColumnVector &m)const{
		if ((*this)._nRows != m._nRows)	throw InvalidMatrixSize();

		ColumnVector v((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)	v._data[i] = (*this)._data[i] + m._data[i];
		return v;
	}
	//-------------------------------------------------------------------------------------------
	const ColumnVector ColumnVector::operator -(const ColumnVector &m)const{
		if ((*this)._nRows != m._nRows)	throw InvalidMatrixSize();

		ColumnVector v((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)	v._data[i] = (*this)._data[i] - m._data[i];
		return v;
	}
	//-------------------------------------------------------------------------------------------
	const Matrix ColumnVector::operator * (const RowVector & rv)const{

		Matrix mat((*this)._nRows, rv.getnColumns());
		Real *pC = mat.getData();
		for (int i = 0; i < (*this)._nRows; i++)
			for (int j = 0; j < rv.getnColumns(); j++)
			{
			Real sigma = 0;
			Real *pA = (*this)._data + i * (*this)._nColumns;
			const Real *pB = rv.getData() + j;

			for (int k = 0; k < (*this)._nColumns; k++)
			{
				sigma += *pA * *pB;
				pA++; pB += rv.getnColumns();
			}
			*pC++ = sigma;
			}
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const ColumnVector ColumnVector::ElementWiseMul(const ColumnVector &v)const{
		if ((*this)._nRows != v._nRows)
			throw Exception("error: ElementWiseMul: The two Matrices have different sizes...");
		ColumnVector vec((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			vec._data[i] = (*this)._data[i] * v._data[i];
		return vec;
	}
	//-------------------------------------------------------------------------------------------
	const ColumnVector ColumnVector::ElementWiseDiv(const ColumnVector &v)const{
		if ((*this)._nRows != v._nRows)
			throw Exception("error: ElementWiseDiv: The two Matrices have different sizes...");
		ColumnVector vec((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			if (std::abs(v._data[i]) != 0)
				vec._data[i] = (*this)._data[i] / v._data[i];
			else
				throw Exception("error: ElementWiseDiv: division by zero...");
		return vec;
	}
	//-------------------------------------------------------------------------------------------
	const ColumnVector ColumnVector::operator + (const Real & n)const{
		// Matrix Scalar Multiplication mat * n
		ColumnVector mat((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			mat._data[i] = (*this)._data[i] + n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const ColumnVector ColumnVector::operator - (const Real & n)const{
		// Matrix Scalar Multiplication mat * n
		ColumnVector mat((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			mat._data[i] = (*this)._data[i] - n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------

	const ColumnVector ColumnVector::operator * (const Real & n)const{
		// Matrix Scalar Multiplication mat * n
		ColumnVector mat((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			mat._data[i] = (*this)._data[i] * n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const ColumnVector ColumnVector::operator / (const Real & n)const{
		// Matrix Scalar divide mat / n
		if (n == 0)
			throw InvalidIndex();
		ColumnVector mat((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			mat._data[i] = (*this)._data[i] / n;
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	const ColumnVector ColumnVector::operator^(const Real & n)const{
		// rise each element to the power n, mat(i)^n
		if ((*this)._nRows < 1)
			throw InvalidMatrixSize();
		ColumnVector mat((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			mat._data[i] = ::pow((*this)._data[i], n);
		return mat;
	}
	const ColumnVector ColumnVector::divideScalarbyCol(Real d)const{
		ColumnVector mat((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			if ((*this)._data[i] != 0)
				mat._data[i] = d / (*this)._data[i];
			else
				Exception("division by zero");
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	RowVector ColumnVector::operator ~()const{
		// Transpose Column Vector to Row Vector
		RowVector mat((*this).getnRows());
		for (int i = 1; i <= (*this).getnRows(); i++)
			mat(i) = (*this)(i);
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	RowVector ColumnVector::Transpose()const{ return ~(*this); }
	//-------------------------------------------------------------------------------------------
	void ColumnVector::ReverseElements()
		/*!
		@brief Reversing in place
		*/
	{
		int n = size(); Real* x = getData(); Real* rx = x + n;
		n /= 2;
		while (n--) { Real t = *(--rx); *rx = *x; *(x++) = t; }
	}
	//-------------------------------------------------------------------------------------------
	ColumnVector ColumnVector::ReverseElements(const ColumnVector & mat)
		/*!
		@brief Reversing into a new matrix
		*/
	{
		ColumnVector n_mat(mat.getnRows());
		int n = mat.size(); Real* rx = n_mat.getData() + n; const Real* x = mat.getData();
		while (n--) *(--rx) = *(x++);
		return n_mat;
	}
	//-------------------------------------------------------------------------------------------
	// Calculates an absolute value of each matrix element.
	ColumnVector ColumnVector::abs()const{
		ColumnVector v((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			v._data[i] = std::abs((*this)._data[i]);
		return v;
	}
	// Computes a square value of each matrix element.
	ColumnVector ColumnVector::Sqr()const{
		ColumnVector v((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			v._data[i] = (*this)._data[i] * (*this)._data[i];
		return v;
	}
	// Computes a square root of each matrix element.
	ColumnVector ColumnVector::Sqrt()const{
		ColumnVector v((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			v._data[i] = ::sqrt((*this)._data[i]);
		return v;
	}
	// Computes an exponent of each matrix element.
	ColumnVector ColumnVector::exp()const{
		ColumnVector v((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			v._data[i] = ::exp((*this)._data[i]);
		return v;
	}
	// Changes each element to the smallest integer greater than or equal to that element.
	ColumnVector ColumnVector::ceil()const{
		ColumnVector v((*this)._nRows);
		for (int i = 0; i<(*this)._nRows; i++)
			v._data[i] = ::ceil((*this)._data[i]);
		return v;
	}
	// Changes each element to the largest integer less than or equal to that element.
	ColumnVector ColumnVector::floor()const{
		int r = (*this).getnRows();
		ColumnVector v(r);
		for (int i = 0; i<r; i++)
			v._data[i] = ::floor((*this)._data[i]);
		return v;
	}
	// Rounds each element to the nearest integral value.
	ColumnVector ColumnVector::round()const{
		int r = (*this).getnRows();
		ColumnVector v(r);
		for (int i = 0; i<r; i++){
			v._data[i] = ((*this)._data[i] > ZERO) ? ::floor((*this)._data[i] + static_cast<Real>(0.5)) : ::ceil((*this)._data[i] - static_cast<Real>(0.5));
		}
		return v;
	}
	// Computes the natural logarithm of each element in the matrix.
	ColumnVector ColumnVector::log()const{
		int r = (*this).getnRows();
		ColumnVector v(r);
		for (int i = 0; i<r; i++)
			v._data[i] = ::log((*this)._data[i]);
		return v;
	}
	// Computes the logarithm of each element in the matrix to base 10.
	ColumnVector ColumnVector::log10()const{
		int r = (*this).getnRows();
		ColumnVector v(r);
		for (int i = 0; i<r; i++)
			v._data[i] = ::log10((*this)._data[i]);
		return v;
	}
	// Computes the natural logarithm of each element in the matrix plus 1.
	ColumnVector ColumnVector::log1p()const{
		int r = (*this).getnRows();
		ColumnVector v(r);
		for (int i = 0; i<r; i++)
			v._data[i] = ::log(ONE + (*this)._data[i]);
		return v;
	}
	// Computes the logarithm of each element in the matrix to base 2. log(n)/log(2) is log2.
	ColumnVector ColumnVector::log2()const{
		int r = (*this).getnRows();
		ColumnVector v(r);
		for (int i = 0; i<r; i++)
			v._data[i] = ::log((*this)._data[i]) / ::log(static_cast<Real>(2));
		return v;
	}
	// Computes the exponent of each element in the matrix.
	/*
	ColumnVector ColumnVector::logb()const{
	int r = (*this).getnRows();
	ColumnVector v(r);
	for ( int i=0; i<r; i++ )
	v._data[i] = ::logb((*this)._data[i]);
	return v;
	}
	*/
	//-------------------------------------------------------------------------------------------
	Real ColumnVector::norm(int p)const{
		Real sum = 0;
		for (int i = 0; i<(*this)._nRows; i++)
			sum += ::pow(fabs((*this)._data[i]), static_cast<Real>(p));
		return ::pow(sum, ONE / static_cast<Real>(p));
	}
	//-------------------------------------------------------------------------------------------
	Real ColumnVector::norm(const std::string &s)const{
		if (s == "fro")
			return norm(2);
		if (s == "inf"){
			Real max = 0;
			for (int i=0; i<_nRows; i++)
				max = MAX(max,std::abs(_data[i]));
			return max;
		}
		throw Exception("error:  ColumnVector::norm: Unrecognised norm...");
	}
	//-------------------------------------------------------------------------------------------
	void ColumnVector::pushBack(const Real e){
		Real *x;		x = new Real[(*this)._nRows];
		for (int i = 0; i<(*this)._nRows; i++)	x[i] = (*this)._data[i];
		(*this)._nRows += 1;
		delete[](*this)._data;	(*this)._data = new Real[(*this)._nRows];
		for (int i = 0; i<(*this)._nRows - 1; i++)	(*this)._data[i] = x[i];
		(*this)._data[(*this)._nRows - 1] = e;
	}
	//-------------------------------------------------------------------------------------------
	void ColumnVector::popBack(){
		(*this)._nRows -= 1;
		Real *x;		x = new Real[(*this)._nRows];
		for (int i = 0; i<(*this)._nRows; i++)	x[i] = (*this)._data[i];
		delete[](*this)._data;	(*this)._data = new Real[(*this)._nRows];
		for (int i = 0; i<(*this)._nRows; i++)	(*this)._data[i] = x[i];
	}
	//-------------------------------------------------------------------------------------------
	Matrix ColumnVector::HorizontallyConcatenated(const ColumnVector& mat){
		int r = (*this)._nRows, i;
		if (r != mat.getnRows())
			throw Exception("error: insertCOL: Invalid vector size...");
		Matrix res(r, 2);
		for (i = 1; i <= r; i++){
			res(i, 1) = (*this)(i);
			res(i, 2) = mat(i);
		}
		return res;
	}
	Matrix ColumnVector::HorizontallyConcatenated(const Matrix& mat){
		int r = (*this)._nRows, i, j;
		if (r != mat.getnRows())
			throw Exception("error: insertCOL: Invalid vector size...");
		int ct = (*this)._nColumns;
		int c = ct + mat.getnColumns();
		Matrix res(r, c);
		for (i = 1; i <= r; i++){
			for (j = 1; j <= mat.getnColumns(); j++){
				if (j == 1) res(i, j) = (*this)(i);
				res(i, j + ct) = mat(i, j);
			}
		}
		return res;
	}
	ColumnVector ColumnVector::VerticallyConcatenated(const ColumnVector& mat){
		int ct = (*this)._nRows, i;
		int cm = mat.getnRows();
		ColumnVector res(ct + cm);
		for (i = 1; i <= ct; i++)
			res(i) = (*this)(i);
		for (i = 1; i <= cm; i++)
			res(i + ct) = mat(i);
		return res;
	}

	//-------------------------------------------------------------------------------------------
	void ColumnVector::insert(const Real & n, int k){
		if (k<1 || k>(*this)._nRows + 1)
			throw Exception("error:  ColumnVector::insert: Invalid Index...");

		Real *x;		x = new Real[(*this)._nRows];
		for (int i = 0; i<(*this)._nRows; i++)	x[i] = (*this)._data[i];
		(*this)._nRows += 1;
		delete[](*this)._data;	(*this)._data = new Real[(*this)._nRows];
		for (int i = 0, j = 0; i<(*this)._nRows; i++, j++)
			if ((i + 1) == k){
			(*this)._data[i] = n;
			j--;
			}
			else
				(*this)._data[i] = x[j];
	}
	//-------------------------------------------------------------------------------------------
	void ColumnVector::swap(int i, int j){
		Real temp = (*this)(i);
		(*this)(i) = (*this)(j);
		(*this)(j) = temp;
	}
	//-------------------------------------------------------------------------------------------
	ColumnVector ColumnVector::FromTo(int i, int j)const{
		if (i < 1 || i > j || j > (*this)._nRows)
			throw InvalidIndex();
		ColumnVector v(j - i + 1);
		for (int k = 1; k <= (j - i + 1); k++)
			v(k) = (*this)(i + k - 1);
		return v;
	}
	//-------------------------------------------------------------------------------------------
	Matrix ColumnVector::skewMatrix()const{
		if ((*this)._nRows != 3)
			throw Exception("error: ColumnVector::skewMatrix(): Vector size must be = 3.");
		Matrix s(3, 3);
		s(1, 2) = -(*this)._data[2];	s(2, 1) = (*this)._data[2];
		s(1, 3) = (*this)._data[1];	s(3, 1) = -(*this)._data[1];
		s(2, 3) = -(*this)._data[0];	s(3, 2) = (*this)._data[0];
		return s;
	}
	//-------------------------------------------------------------------------------------------

};
