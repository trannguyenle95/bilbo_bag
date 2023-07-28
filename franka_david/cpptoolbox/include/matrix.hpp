#ifndef ___MATRIX_H
#define ___MATRIX_H

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <assert.h>
#include "include.h"
#include "mathUtility.hpp"
#include "Exceptions.h"
#include "numericalAlg.h"

#ifdef _MSC_VER
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif // !_CRT_SECURE_NO_WARNINGS
#endif

using namespace std;
namespace math{

	Real getTime();

	class BaseMatrix;
	class Matrix;
	class RowVector;
	class ColumnVector;

	Matrix inv(const Matrix & mat);

	class InvalidIndex : public Exception {
	public:
		InvalidIndex() : Exception( "Index out of range" ) { }
	};
	class InvalidMatrixSize : public Exception {
	public:
		InvalidMatrixSize() : Exception( "Invalid Matrix Size" ) { }
	};
	class InvalidMatrix : public Exception {
	public:
		InvalidMatrix() : Exception( "Invalid Matrix" ) { }
	};
	class InvalidSquareMatrix : public Exception {
	public:
		InvalidSquareMatrix() : Exception( "Invalid Square Matrix" ) { }
	};

	class BaseMatrix{
	protected:
		Real *_data;
		int _nRows;
		int _nColumns;
		bool in_range(int r, int c) const;
		bool row_in_range(int r) const;
		bool col_in_range(int c) const;
		bool in_range(int i) const;
	private:
		mutable long _r;
		mutable long _c;
		mutable bool has_been_used;
		struct assignMatix{
			assignMatix(const assignMatix& mat) : _m(mat._m), _r(mat._r), _c(mat._c), _used(false) {}
			assignMatix(BaseMatrix* m_): _m(m_), _r(1), _c(1),_used(false) {next();}
			~assignMatix(){ }
			const assignMatix& operator, (const Real& val) const {
				(*_m)(_r,_c) = val;
				next();
				_used = true;
				return *this;
			}
		private:
			void next () const{
				++_c;
				if (_c > _m->_nColumns){_c = 1;		++_r;}
			}
			BaseMatrix * _m;
			mutable long _r;
			mutable long _c;
			mutable bool _used;
        };
	public:
		//! assign the given value to every spot in this matrix
		const assignMatix operator << (const Real& val );
		//! Default constructor. Creates a 1 by 1 matrix; sets value to zero.
		BaseMatrix() : _r(1), _c(1), has_been_used(false), _data(new Real[1]), _nRows(1), _nColumns(1) {_data[0] = 0;}
		//! Create a matrix of size (rows, cols). sets all values to val.
		BaseMatrix(int r, int c, Real val=0);
		//! Constructor. Create a r by c matrix. Fill it from a given array x.
		BaseMatrix(int r, int c, const Real *x);
		//! Constructor. Create a r by c matrix. Fill it from a given array x.
		BaseMatrix(int r, int c, const string &str);
		//! Copy Constructor
		BaseMatrix(const BaseMatrix & m) {	copy(m);    }
		//! Copy Constructor
		virtual void copy(const BaseMatrix & m);
		//! Destructor
		virtual ~BaseMatrix();

		//! Assignment operator function..
		BaseMatrix & operator =(const BaseMatrix & m);
		BaseMatrix & operator >> (const Real *val );
		//! Set all elements of the matrix equal to x.
		BaseMatrix & operator =(Real x);
		//! Set all elements of the matrix equal to x.
		BaseMatrix & operator =(const std::string &str)		{	set(str);	return *this;	}
		//! Set matrix equal to values in the string str
		virtual void set(const std::string &str) = 0;

		//! Get element (i,j) from the matrix
		const Real & operator () (int i, int j) const;
		//! Get or/and assign element (i,j) in the matrix
		Real & operator () (int i, int j);
		//! Set an Element in the Matrix
		Real & operator () (int i, int j, Real e);
		//! Set an Array in the Matrix
		BaseMatrix & operator () (Real *array, int n);
		//! Get Element (i,j) from the Matrix
		virtual const Real & getElement(int i, int j) const = 0;
		void setElement( int i, int j, Real e);	//! Set an Element in the Matrix
		//! Set an Array in /(static_cast<Real>(1.0) + temp * xx)the Matrix
		void setElements(Real *array, int n){(*this)(array,n);}

		friend ostream & operator << (ostream & o, const BaseMatrix & m);
		void print (char const* ch = 0, int prec=4, int w=8) const;

		//! The number of rows
		inline int getnRows() const			{ return _nRows;			}
		//! The number of columns
		inline int getnColumns() const		{ return _nColumns;			}
		//! The number of elements
		inline int size() const				{ return _nRows * _nColumns;}
		//! The data array
		inline Real * getData()				{ return _data;				}
		//! The data array
		inline const Real * getData()const	{ return _data;				}
		//! Set matrix equal to the all zero matrix.
		void zeros ();
		//! Set matrix equal to the all zero matrix.
		void clear()						{ zeros();					}
		//! Set matrix equal to the all one  matrix.
		void ones ();
		//! Set matrix equal to the all val  matrix
		void fill(Real val);

		//! Compare two matrices. False if wrong sizes or different values
		bool operator ==(const BaseMatrix & m)const;
		//! Compare two matrices. True if different
		bool operator !=(const BaseMatrix & m)const	{ return !((*this) == m);	}

		void add(const BaseMatrix &m1, const BaseMatrix &m2, BaseMatrix &m3);

		//! Sum of all elements in the matrix
		static Real sum(const BaseMatrix &m);
		Real sum()const							{ return sum(*this);		}
		//! The average value
		Real average()const						{return sum(*this) / size();}
		static Real average(const BaseMatrix &v){return sum(v) / v.size();	}
		//! Maximum value of the matrix
		static Real maxV(const BaseMatrix &m);
		Real maxV()const							{ return maxV(*this);		}
		//! Maximum value of matrix, also returns the index position of max value
		static Real maxV(const BaseMatrix &m,int &i,int &j);
		Real maxV(int &i,int &j)const				{ return maxV(*this,i,j);	}
		//! Minimum value of the matrix
		static Real minV(const BaseMatrix &m);
		Real minV()const							{ return minV(*this);		}
		//! Minimum value of matrix, also returns the index position of max value
		static Real minV(const BaseMatrix &m,int &i,int &j);
		Real minV(int &i,int &j)const				{ return minV(*this,i,j);	}
		//! Maximum absolute value of the matrix
		static Real MaximumAbsoluteValue(const BaseMatrix &m);
		Real MaximumAbsoluteValue()const			{ return MaximumAbsoluteValue(*this);	}
		//! Minimum absolute value of the matrix
		static Real MinimumAbsoluteValue(const BaseMatrix &m);
		Real MinimumAbsoluteValue()const			{ return MinimumAbsoluteValue(*this);	}

		static Real SumSquare(const BaseMatrix &m){ return m.SumSquare();				}
		Real SumSquare()const;

		//! Norm of the matrix. Calculate the 2-norm: is defined in each matrix type
		virtual Real norm()const = 0;
    };
    class Matrix : public BaseMatrix{
    private:
        //! LU-Decomposition. (cf. Numerical Recipes)
        void ludcmp(Matrix & a, int * index, Real &d)const;
        //! LU-Back-substitution. (cf. Numerical Recipes)
        void lubksb(Matrix & a, int * index, ColumnVector & b)const;
    public:
		Matrix() : BaseMatrix(){}
		Matrix(int r, Real val=0) : BaseMatrix(r,r,val)	{}
		Matrix(int r, int c, Real val=0) : BaseMatrix(r,c,val)	{}
		Matrix(int r, int c, const Real *x) : BaseMatrix(r,c,x)	{}
		Matrix(int r, int c, const string &str) : BaseMatrix(r,c,str)	{}
		Matrix(int r, const Real *x) : BaseMatrix(r,r,x)	{}
		//! Set matrix equal to values in string
		Matrix(const std::string &str)	{set(str);}

		Matrix(const Matrix & m) : BaseMatrix(m)	{}
		void copy(const Matrix & m)					{BaseMatrix::copy(m);							}

		Matrix & operator =(const Matrix & m)		{return (Matrix&)((*(BaseMatrix*)this) = m);	}
		Matrix & operator =(Real x)					{return (Matrix&)((*(BaseMatrix*)this) = x);	}
		Matrix & operator =(const std::string &str)	{return (Matrix&)((*(BaseMatrix*)this) = str);	}
		bool operator ==(const Matrix & m)const		{return (*(BaseMatrix*)this) == m;				}
		bool operator !=(const Matrix & m)const		{ return !((*this) == m);						}

		void set(const std::string &str);

		//! Get element i using linear addressing (by rows)
		Real operator () (int i) const;
		Real & operator () (int i);
		//! Get & Set element i, j
		Real operator () (int i, int j) const{return (*(BaseMatrix *)this)(i,j);				}
		Real & operator () (int i, int j)			 {return (*(BaseMatrix *)this)(i,j);				}
		Real & operator () (int i, int j, Real e)	 {return (*(BaseMatrix *)this)(i,j,e);				}
		const Real & getElement(int i, int j) const	 {return (*(BaseMatrix *)this)(i,j);					}
		void setElement( int i, int j, Real e)		 {(*(BaseMatrix *)this).setElement(i,j,e);			}
		//! Fill a matrix by a giving array
		Matrix & operator () (Real *array, int n)	 {return (Matrix&)((*(BaseMatrix*)this)(array,n));	}

		//! Matrix Addition m1 + m2
		const Matrix operator +(const Matrix &m)const;
		//! Matrix Subtraction m1 - m2
		const Matrix operator -(const Matrix &m)const;
		//! Matrix Multiplication m1 * m2
		const Matrix operator *(const Matrix &m)const;
		//! Matrix and Column Vector Multiplication m(n,m) * cv(m)
		const ColumnVector operator *(const ColumnVector &cv)const;
		//! Matrix and Row Vector Multiplication m(n,m) * rv(m)
		const Matrix operator *(const RowVector &rv)const;

		//! Element-wise Division m1 ./ m2
		const Matrix   operator /  (const Matrix &m)const	{return ElementWiseDiv(m);}
		const Matrix & operator /= (const Matrix &m)		{return *this = *this / m;}
		const Matrix ElementWiseDiv(const Matrix &m)const;
		//! Element-wise Multiplication m1 .* m2
		const Matrix   operator %  (const Matrix &m)const	{return ElementWiseMul(m);}
		const Matrix & operator %= (const Matrix &m)		{return *this = *this % m;}
		const Matrix ElementWiseMul(const Matrix &m)const;

		//! m1 += m2	==>	m1 = m1 + m2
		const Matrix & operator +=(const Matrix &m)	{ return *this = *this + m;		}
		//! m1 -= m2	==>	m1 = m1 - m2
		const Matrix & operator -=(const Matrix &m)	{ return *this = *this - m;		}
		//! m1 *= m2	==>	m1 = m1 * m2
		const Matrix & operator *=(const Matrix &m);
		//! m++	==>	m = m + 1
		Matrix & operator ++();
		Matrix   operator ++(int);
		//! m++	==>	m = m - 1
		Matrix & operator --();
		Matrix   operator --(int);

		//! Matrix Scalar Addition mat + n
		const Matrix operator + (const Real & n)const;
		//! Matrix Scalar Subtraction mat - n
		const Matrix operator - (const Real & n)const;
		//! Matrix Scalar Multiplication mat * n
		const Matrix operator * (const Real & n)const;
		//! Matrix Scalar divide mat / n
		const Matrix operator / (const Real & n)const;
		//! Raises every matrix element to a power.
		const Matrix operator ^ (const Real & n)const;
		//! Unary operator
		const Matrix operator - ( )const	{ return *this * -ONE;}

		//! m1 += m2	==>	m1 = m1 + m2
		const Matrix & operator +=(const Real & n)	{ return *this = *this + n;	}
		//! m1 -= m2	==>	m1 = m1 - m2
		const Matrix & operator -=(const Real & n)	{ return *this = *this - n;	}
		//! m1 *= m2	==>	m1 = m1 * m2
		const Matrix & operator *=(const Real & n)	{ return *this = *this * n;	}
		//! m1 *= m2	==>	m1 = m1 / m2
		const Matrix & operator /=(const Real & n)	{ return *this = *this / n;	}
		//! m1 *= m2	==>	m1 = m1 ^ m2
		const Matrix & operator ^=(const Real & n)	{ return *this = *this ^ n;	}
		static Matrix divideScalarbyMatrix(Real d, const Matrix &m){return m.divideScalarbyMatrix(d);}
		const Matrix divideScalarbyMatrix(Real d)const;

		//! Matrix Scalar Addition n + mat
		friend Matrix operator + (Real d, const Matrix & m){return m + d;}
		//! Matrix Scalar Subtraction n - mat
		friend Matrix operator - (Real d, const Matrix & m){return (m * -ONE) + d;}
		//! Matrix Scalar Multiplication n * mat
		friend Matrix operator * (Real d, const Matrix & m){return m * d;}
		//! Matrix Scalar division n ./ mat
		friend Matrix operator / (Real d, const Matrix & m) {return m.divideScalarbyMatrix(d);}

		static void fastAdd(const Matrix &a, const Matrix &b, Matrix &c);
		static void fastSub(const Matrix &a, const Matrix &b, Matrix &c);

		//! Matrix Transpose
		Matrix operator ~ ()const;
		Matrix Transpose()const			{return ~(*this);				}

		static Matrix ones(int r, int c)	{return Matrix(r,c,1);		}
		void ones()							{((BaseMatrix&)*this).ones();}

		//! Get row r from a matrix
		RowVector getRow(int r) const;
		//! Get rows r1 through r2
		Matrix getRows(int r1, int r2) const;
		//! Get column c from a matrix
		ColumnVector getColumn(int c) const;
		//! Get column c1 through c2
		Matrix getColumns(int c1, int c2) const;
		//! Get Matrix(r1:r2, c1:c2);
		Matrix subMatrix(int r1, int r2, int c1, int c2) const;
		//! Get row r and return it as a column
		ColumnVector getRowAsCol(int m) const;
		//! Get column c and return it as a row
		RowVector getColAsRow(int m) const;
		//! Get Matrix as ColumnVector
		ColumnVector getMatrixAsColumn() const;
		//! Set ColumnVector as a Row in a matrix
		void setColAsRow(const ColumnVector &col, int i);
		//! Set column c to the matrix
		void setCol(const ColumnVector &c, int i);
		//! Set row r to to the matrix
		void setRow(const ColumnVector &r, int i);
		//! Swap the rows r1 and r2
		void swap_rows(int r1, int r2);
		//! Swap the columns c1 and c2
		void swap_cols(int c1, int c2);
		//! Create a Matrix from n copies of column vector v
		static Matrix repmat(const ColumnVector &v, int n);
		//! Create a Matrix from m copies of row vector v
		static Matrix repmat(const RowVector &v, int m);
		Matrix asDiagonal();

		//! add a new row at the end of the matrix
		void pushBackROW();
		void pushBackROW(const RowVector &rv);
		//! add a new column at the end of the matrix (from the right)
		void pushBackCOL();
		void pushBackCOL(const ColumnVector &cv);
		//! Delete the last row from the matrix
		void popBackROW();
		//! Delete the last col from the matrix (from the right)
		void popBackCOL();
		//! Delete a specific row from the matrix
		void eraseROW(int i);
		//! Delete a specific column from the matrix
		void eraseCOL(int i);
		//! Insert new row to the matrix, nRows++.
		void insertROW(int i, const RowVector &rv);
		//! Insert new column to the matrix, nColumn++.
		void insertCOL(int i, const ColumnVector &cv);

		Matrix HorizontallyConcatenated(const Matrix&);
		Matrix operator |(const Matrix& mat)		{return HorizontallyConcatenated(mat);	}
		Matrix HorizontallyConcatenated(const ColumnVector&);
		Matrix operator |(const ColumnVector& mat)	{return HorizontallyConcatenated(mat);	}

		Matrix VerticallyConcatenated(const Matrix& mat);
		Matrix operator &(const Matrix& mat)		{return VerticallyConcatenated(mat);	}
		Matrix VerticallyConcatenated(const RowVector& mat);
		Matrix operator &(const RowVector& mat)		{return VerticallyConcatenated(mat);	}

		//! Returns the r-by-c matrix whose elements	are taken column-wise from this object
		Matrix reshape(int r, int c);
		static Matrix reshape(Matrix &m, int r, int c){return m.reshape(r,c);}
		//! Resizes the Matrix, initializes new values to zero.
		void  resize (int r, int c);
		void  resize (int r, int c, Real val);
		//! reConstruct the Matrix base to a given one.
		void  reConstruct (const Matrix &m);

		// reversing in place
		void ReverseElements();
		// reversing into a new matrix
		static Matrix ReverseElements(const Matrix&);

		//! Calculates an absolute value of each matrix element.
		Matrix abs()const;
		static Matrix abs(const Matrix &m)	{return m.abs();	}
		//! Computes a square value of each matrix element.
		Matrix Sqr()const;
		static Matrix Sqr(const Matrix &m)	{return m.Sqr();	}
		//! Computes a square root of each matrix element.
		Matrix Sqrt()const;
		static Matrix Sqrt(const Matrix &m)	{return m.Sqrt();	}
		//! Raises every matrix element to a power.
		Matrix pow(const Real & n)const		{return (*this)^n;	}
		//! Computes an exponent of each matrix element.
		Matrix exp()const;
		static Matrix exp(const Matrix &m)	{return m.exp();	}
		//! Changes each element to the smallest integer greater than or equal to that element.
		Matrix ceil()const;
		static Matrix ceil(const Matrix &m)	{return m.ceil();	}
		//! Changes each element to the largest integer less than or equal to that element.
		Matrix floor()const;
		static Matrix floor(const Matrix &m){return m.floor();	}
		//! Rounds each element to the nearest integral value.
		Matrix round()const;
		static Matrix round(const Matrix &m){return m.round();	}
		//! Computes the natural logarithm of each element in the matrix.
		Matrix log()const;
		static Matrix log(const Matrix &m)	{return m.log();	}
		//! Computes the logarithm of each element in the matrix to base 10.
		Matrix log10()const;
		static Matrix log10(const Matrix &m){return m.log10();	}
		//! Computes the natural logarithm of each element in the matrix plus 1.
		Matrix log1p()const;
		static Matrix log1p(const Matrix &m){return m.log1p();	}
		//! Computes the logarithm of each element in the matrix to base 2. log(n)/log(2) is log2.
		Matrix log2()const;
		static Matrix log2(const Matrix &m)	{return m.log2();	}
		//! Computes the exponent of each element in the matrix.
/*
		Matrix logb()const;
		static Matrix logb(const Matrix &m)	{return m.logb();	}
*/

		//! Computes sum of diagonal elements.
		Real trace()const;
		static Real trace(const Matrix &m)	{return m.trace();	}
		//! Returns the k-th diagonal in a matrix.
		ColumnVector diag(int k = 0)const;
		static ColumnVector diag(const Matrix &m, int k = 0);
		//! Checks if the Matrix is Square Matrix or not.
		bool isSquare()const;
		//! Compute the 2-norm of matrix
		Real norm()const;
		static Real norm(const Matrix &m)	{return m.norm();	}
		Real normINF()const;
		//! Sum of elements in the matrix m, either along columns or rows
		//! sum(1) returns a vector where the elements are sum over each column,
		//! whereas (2) returns a vector where the elements are sum over each row.
		ColumnVector sum(int dim)const;
		static ColumnVector sum(const Matrix &m, int dim);

		static Matrix PseudoInverse(Matrix &A);
		//! Calculate the rank of the matrix. tol Tolerance used for comparing the singular values with zero.
		int rank(Real tol = -1.0);

		//! Matrix-Inversion using determinant
		Matrix Inverse()const;
		static Matrix Inverse(const Matrix &mat)	{return mat.Inverse();		}
		//! Matrix-Inversion using Cholesky decomposition
		Matrix InverseCholesky()const;
		Matrix i()const										{return math::inv(*this);		}
		static Matrix InverseCholesky(const Matrix &mat)	{return mat.InverseCholesky();	}
		//! Matrix-Inversion using LU-Decomposition.
		Matrix InverseLU()const;
		static Matrix InverseLU(const Matrix &mat)	{return mat.InverseLU();	}
		//! Matrix-Inversion using Gauss-Jordan.
		Matrix InverseGJ()const;
		static Matrix InverseGJ(const Matrix &mat)	{return mat.InverseGJ();	}
		//! Matrix determinant using recursive method
		Real determinant()const;
		static Real determinant(const Matrix &mat)	{return mat.determinant();	}
		//! Matrix determinant using LU-Decomposition.
		Real determinantLU()const;
		static Real determinantLU(const Matrix &mat){return mat.determinantLU();}
		//! Determinant of the matrix using Cholesky decomposition
		Real detCholesky()const;
		Real det()const								{return detCholesky();		}
		static Real detCholesky(const Matrix &mat)	{return mat.detCholesky();	}
		//! Solve the inhomogeneous linear system A*x = b with A quadratic using LU-decomposition
		ColumnVector solve(const ColumnVector &b)const;
		//! LU factorization.
		Matrix LU()const;

		static Matrix eye(int n);
		void eye()	{ *this = eye((*this).getnRows());}
		static Matrix Zeros(int m, int n = -1){ return (n>0) ? Matrix(m,n) : Matrix(m); }
    };

    class SquareMatrix : public Matrix{
    public:
    	SquareMatrix() : Matrix(){}
		SquareMatrix(int r, Real val=0) : Matrix(r,r,val)	{}
		SquareMatrix(int r, const Real *x) : Matrix(r,r,x)	{}
		SquareMatrix(int r, const string &str) : Matrix(r,r,str)	{}
		//! Set matrix equal to values in string
		SquareMatrix(const std::string &str)	{Matrix::set(str);}
		SquareMatrix(const Matrix & m) : Matrix(m)	{}
    };
    class DiagonalMatrix : public Matrix{
    public:
    	DiagonalMatrix() : Matrix(){}
    	DiagonalMatrix(int r, Real val=0) : Matrix(r)	{for(int i=1; i<=r; i++) (*this)(i) = val;}
		DiagonalMatrix(int n, const ColumnVector &v);
		//DiagonalMatrix(int r, const Real *x) : Matrix(r,r,x)	{}
		//! Set matrix equal to values in string
		DiagonalMatrix(const std::string &str)	{Matrix::set(str);}
		DiagonalMatrix(const Matrix & m) : Matrix(m)	{}
		//! Get element i using linear addressing (by rows)
		Real operator () (int i) const{ return Matrix::operator ()(i,i);}
		Real & operator () (int i){ return Matrix::operator ()(i,i);}
    };

    class ColumnVector : public BaseMatrix{
    public:
		ColumnVector (int n=1, Real val=0) : BaseMatrix(n,1,val)	{}
		ColumnVector( int n, const Real *x ) : BaseMatrix (n,1,x)	{}
		ColumnVector( int n, const string &str ) : BaseMatrix (n,1,str)	{}
		ColumnVector(const ColumnVector & cv) : BaseMatrix (cv._nRows,1){copy(cv);				}
		ColumnVector(const std::string &str)	{set(str);	}
		ColumnVector( Real x, Real y, Real z );

		void copy(const ColumnVector & cv)		{BaseMatrix::copy(cv);	}

		ColumnVector & operator =(const ColumnVector & cv)	{return (ColumnVector&)((*(BaseMatrix*)this) = cv);	}
		ColumnVector & operator =(Real x)					{return (ColumnVector&)((*(BaseMatrix*)this) = x);	}
		ColumnVector & operator =(const std::string &str)	{return (ColumnVector&)((*(BaseMatrix*)this) = str);}
		void set(const std::string &str);

		bool operator ==(const ColumnVector & m)const    {return (*(BaseMatrix*)this) == m;}
		bool operator !=(const ColumnVector & m)const    { return !((*this) == m);	}

		const Real & operator () (int i) const			{return (*(BaseMatrix *)this)(i,1);						}
		Real & operator () (int i)						{return (*(BaseMatrix *)this)(i,1);						}
		ColumnVector & operator () (Real *array, int n)	{return (ColumnVector&)((*(BaseMatrix*)this)(array,n));	}
		const Real & getElement(int i, int j) const		{return (*(BaseMatrix *)this)(i,1);						}
		const Real & getElement(int i) const			{return (*(BaseMatrix *)this)(i,1);						}
		void setElement( int i, Real e)					{(*(BaseMatrix *)this).setElement(i,1,e);				}

		//! Column Vector Addition m1 + m2
		const ColumnVector operator +(const ColumnVector &v)const;
		//! Column Vector Subtraction m1 - m2
		const ColumnVector operator -(const ColumnVector &v)const;
		//! Column Vector * Row Vector
		const Matrix operator *(const RowVector &v)const;

		//! Element-wise Multiplication of vectors m1 .* m2
		const ColumnVector   operator % (const ColumnVector &v)const{return ElementWiseMul(v);	}
		const ColumnVector & operator %=(const ColumnVector &v)		{return *this = *this % v;	}
		const ColumnVector ElementWiseMul(const ColumnVector &v)const;
		static ColumnVector ElementWiseMul(const ColumnVector &v1, const ColumnVector &v2){return v1.ElementWiseMul(v2);}
		//! Element-wise Division of vectors m1 ./ m2
		const ColumnVector   operator / (const ColumnVector &v)const{return ElementWiseDiv(v);	}
		const ColumnVector & operator /=(const ColumnVector &v)		{return *this = *this / v;	}
		const ColumnVector ElementWiseDiv(const ColumnVector &v)const;
		static ColumnVector ElementWiseDiv(const ColumnVector &v1, const ColumnVector &v2){return v1.ElementWiseDiv(v2);}

		//! m1 += m2	==>	m1 = m1 + m2
		const ColumnVector & operator +=(const ColumnVector &v)		{ return *this = *this + v; }
		//! m1 -= m2	==>	m1 = m1 - m2
		const ColumnVector & operator -=(const ColumnVector &v)		{ return *this = *this - v; }
		//! m1 ./= m2	==>	m1 = m1 ./ m2

		//! Column Vector Scalar Addition v + n
		const ColumnVector operator + (const Real & n)const;
		//! Column Vector Scalar Subtraction v - n
		const ColumnVector operator - (const Real & n)const;
		//! Column Vector Scalar Multiplication mat * n
		const ColumnVector operator * (const Real & n)const;
		//! Column Vector Scalar divide v / n
		const ColumnVector operator / (const Real & n)const;
		//! Raises every matrix element to a power.
		const ColumnVector operator ^ (const Real & n)const;
		//! Unary operator
		const ColumnVector operator - ()const	{ return *this * -ONE;		}
		static ColumnVector divideScalarbyCol(Real d, const ColumnVector &m){return m.divideScalarbyCol(d);}
		const ColumnVector divideScalarbyCol(Real d)const;

		//! m1 += m2	==>	m1 = m1 + m2
		const ColumnVector & operator +=(const Real & n)	{ return *this = *this + n;	}
		//! m1 -= m2	==>	m1 = m1 - m2
		const ColumnVector & operator -=(const Real & n)	{ return *this = *this - n;	}
		//! m1 *= m2	==>	m1 = m1 * m2
		const ColumnVector & operator *=(const Real & n)	{ return *this = *this * n;	}
		//! m1 /= m2	==>	m1 = m1 / m2
		const ColumnVector & operator /=(const Real & n)	{ return *this = *this / n;	}
		//! m1 ^= m2	==>	m1 = m1 ^ m2
		const ColumnVector & operator ^=(const Real & n)	{ return *this = *this ^ n;	}
		//! v++	==>	v = v + 1
		const ColumnVector & operator ++() { return *this = *this + ONE;	}
		//! v++	==>	v = v - 1
		const ColumnVector & operator --() { return *this = *this - ONE;	}

		//! Column Vector Scalar Addition n + v
		friend ColumnVector operator + (Real d, const ColumnVector & m){return m + d;}
		//! Column Vector Scalar Subtraction n - v
		friend ColumnVector operator - (Real d, const ColumnVector & m){return (m * -ONE) + d;}
		//! Column Vector Scalar Multiplication n * v
		friend ColumnVector operator * (Real d, const ColumnVector & m){return m * d;}
		friend ColumnVector operator / (Real d, const ColumnVector & m){return m.divideScalarbyCol(d);}

		//! Column Vector Transpose
		RowVector operator ~ ()const;
		RowVector Transpose()const;

		// reversing in place
		void ReverseElements();
		// reversing into a new matrix
		static ColumnVector ReverseElements(const ColumnVector&);

		static ColumnVector ones(int r){return ColumnVector(r,1);}
		void ones(){return ((BaseMatrix&)*this).ones();}

		//! Calculates an absolute value of each matrix element.
		ColumnVector abs()const;
		static ColumnVector abs(const ColumnVector &v)	{return v.abs();	}
		//! Computes a square value of each matrix element.
		ColumnVector Sqr()const;
		static ColumnVector Sqr(const ColumnVector &v)	{return v.Sqr();	}
		//! Computes a square root of each matrix element.
		ColumnVector Sqrt()const;
		static ColumnVector Sqrt(const ColumnVector &v)	{return v.Sqrt();	}
		//! Raises every matrix element to a power.
		ColumnVector pow(const Real & n)const			{return (*this)^n;	}
		//! Computes an exponent of each matrix element.
		ColumnVector exp()const;
		static ColumnVector exp(const ColumnVector &v)	{return v.exp();	}
		//! Changes each element to the smallest integer greater than or equal to that element.
		ColumnVector ceil()const;
		static ColumnVector ceil(const ColumnVector &v)	{return v.ceil();	}
		//! Changes each element to the largest integer less than or equal to that element.
		ColumnVector floor()const;
		static ColumnVector floor(const ColumnVector &v){return v.floor();	}
		//! Rounds each element to the nearest integral value.
		ColumnVector round()const;
		static ColumnVector round(const ColumnVector &v){return v.round();	}
		//! Computes the natural logarithm of each element in the matrix.
		ColumnVector log()const;
		static ColumnVector log(const ColumnVector &v)	{return v.log();	}
		//! Computes the logarithm of each element in the matrix to base 10.
		ColumnVector log10()const;
		static ColumnVector log10(const ColumnVector &v){return v.log10();	}
		//! Computes the natural logarithm of each element in the matrix plus 1.
		ColumnVector log1p()const;
		static ColumnVector log1p(const ColumnVector &v){return v.log1p();	}
		//! Computes the logarithm of each element in the matrix to base 2. log(n)/log(2) is log2.
		ColumnVector log2()const;
		static ColumnVector log2(const ColumnVector &v)	{return v.log2();	}
		//! Computes the exponent of each element in the matrix.
/*
		ColumnVector logb()const;
		static ColumnVector logb(const ColumnVector &v)	{return v.logb();	}
*/

		//! Compute the 2-norm of vector
		Real norm()const								{return norm(2);	}
		static Real norm(const ColumnVector &v)			{return v.norm();	}
		//! Compute the p-norm. v.norm(p)=sum(abs(v).^2)^(1/p)
		Real norm(int p)const;
		//! Compute the Frobenius norm for s = "fro" (equal to 2-norm)
		Real norm(const std::string &s)const;
		//! Maximum value of the vector
		static Real maxV(const ColumnVector &v,int &i)	{int j = 1; return BaseMatrix::maxV(v,i,j);	}
		Real maxV(int &i)const							{return maxV(*this,i);						}
		Real maxV()const								{return BaseMatrix::maxV(*this);			}
		//! Minimum value of the vector
		static Real minV(const ColumnVector &v,int &i)	{int j = 1; return BaseMatrix::minV(v,i,j);	}
		Real minV(int &i)const							{return minV(*this,i);						}
		Real minV()const								{return BaseMatrix::minV(*this);			}

		//! Add element at the end
		void pushBack(const Real e);
		//! Delete last element
		void popBack();

		Matrix HorizontallyConcatenated(const ColumnVector&);
		Matrix operator |(const ColumnVector& mat)		{return HorizontallyConcatenated(mat);}
		Matrix HorizontallyConcatenated(const Matrix&);
		Matrix operator |(const Matrix& mat)			{return HorizontallyConcatenated(mat);}
		ColumnVector VerticallyConcatenated(const ColumnVector& mat);
		ColumnVector operator &(const ColumnVector& mat){return VerticallyConcatenated(mat);}

		//! Insert an element at position k; ==> k = 1, ... , V.size()
		void insert(const Real & n, int k = 1);
		void resize(int n)					{(*(Matrix*)this).resize(n,1);		}
		void resize(int n, Real val)		{(*(Matrix*)this).resize(n,1,val);	}
		void swap(int i, int j);
		ColumnVector FromTo(int i, int j)const;

		//!Creates a skew symmetric matrix from a PositionVector3D.
		Matrix skewMatrix()const;
    };

    class RowVector : public BaseMatrix{
    public:
		RowVector( int n=1, Real val=0 ) : BaseMatrix (1,n,val)	{}
		RowVector( int n, const Real *x ) : BaseMatrix (1,n,x)	{}
		RowVector( int n, const string &str ) : BaseMatrix (1,n,str)	{}
		RowVector(const RowVector & rv) : BaseMatrix (1,rv._nColumns) {copy(rv);}
		RowVector(const std::string &str) {set(str);}
		RowVector( Real x, Real y, Real z );

		void copy(const RowVector & rv)					{BaseMatrix::copy(rv);								}

		RowVector & operator =(const RowVector & rv)	{return (RowVector&)((*(BaseMatrix*)this) = rv);}
		RowVector & operator =(Real x)					{return (RowVector&)((*(BaseMatrix*)this) = x);	}
		RowVector & operator =(const std::string &str)	{return (RowVector&)((*(BaseMatrix*)this) = str);}
		void set(const std::string &str);

		bool operator ==(const RowVector & m)const	{return (*(BaseMatrix*)this) == m;	}
		bool operator !=(const RowVector & m)const	{ return !((*this) == m);			}

		const Real & operator () (int i) const		{return (*(BaseMatrix *)this)(1,i);	}
		Real & operator () (int i)					{return (*(BaseMatrix *)this)(1,i);	}
		RowVector& operator () (Real *array, int n)	{return (RowVector&)((*(BaseMatrix*)this)(array,n));}
		const Real & getElement(int i, int j) const	{return (*(BaseMatrix *)this)(1,j);	}
		const Real & getElement(int i) const		{return (*(BaseMatrix *)this)(1,i);	}
		void setElement( int i, Real e)				{(*(BaseMatrix *)this).setElement(1,i,e);}

		const RowVector operator +(const RowVector	&v)const;
		const RowVector operator -(const RowVector	&v)const;
		const RowVector operator *(const Matrix		&v)const;
		const Real operator *(const ColumnVector &v)const;

		//! Element-wise Multiplication m1 .* m2
		const RowVector   operator %  (const RowVector &v)const	{return ElementWiseMul(v);	}
		const RowVector & operator %= (const RowVector &v)		{return *this = *this % v;	}
		const RowVector ElementWiseMul(const RowVector &v)const;
		//! Element-wise Division m1 ./ m2
		const RowVector   operator /  (const RowVector &v)const	{return ElementWiseDiv(v);	}
		const RowVector & operator /= (const RowVector &v)		{return *this = *this / v;	}
		const RowVector ElementWiseDiv(const RowVector &v)const;

		const RowVector & operator +=(const RowVector    &v)	{ return *this = *this + v;	}
		const RowVector & operator -=(const RowVector    &v)	{ return *this = *this - v;	}
		const RowVector & operator *=(const Matrix	   &v);
		const RowVector & operator ++()	{ return *this = *this + ONE;	}
		const RowVector & operator --()	{ return *this = *this - ONE;	}

		const RowVector operator + (const Real & n)const;	//! Matrix Scalar Addition mat + n
		const RowVector operator - (const Real & n)const;	//! Matrix Scalar Subtraction mat - n
		const RowVector operator * (const Real & n)const;	//! Matrix Scalar Multiplication mat * n
		const RowVector operator / (const Real & n)const;	//! Matrix Scalar divide mat / n
		//! Raises every matrix element to a power.
		const RowVector operator ^ (const Real & n)const;
		const RowVector operator - ( )const	{ return (*this) * -1;	}	//! Unary operator
		static RowVector divideScalarbyRow(Real d, const RowVector &m){return m.divideScalarbyRow(d);}
		const RowVector divideScalarbyRow(Real d)const;

		const RowVector & operator +=(const Real & n) { return (*this) = (*this) + n;	}	//! m1 += n
		const RowVector & operator -=(const Real & n) { return (*this) = (*this) - n;	}	//! m1 -= n
		const RowVector & operator *=(const Real & n) { return (*this) = (*this) * n;	}	//! m1 *= n
		const RowVector & operator /=(const Real & n) { return (*this) = (*this) / n;	}	//! m1 /= n
		const RowVector & operator ^=(const Real & n) { return (*this) = (*this) ^ n;	}	//! m1 /= n

		friend RowVector operator + (Real d, const RowVector & m){return m + d;}
		friend RowVector operator - (Real d, const RowVector & m){return (m * -ONE) + d;}
		friend RowVector operator * (Real d, const RowVector & m){return m * d;}
		friend RowVector operator / (Real d, const RowVector & m){return m.divideScalarbyRow(d);}

		ColumnVector operator ~ ()const;		//! Matrix Transpose
		ColumnVector Transpose()const;

		// reversing in place
		void ReverseElements();
		// reversing into a new matrix
		static RowVector ReverseElements(const RowVector&);

		static RowVector ones(int r){return RowVector(r,1);}
		void ones(){return ((BaseMatrix&)*this).ones();}

		//! Calculates an absolute value of each matrix element.
		RowVector abs()const;
		static RowVector abs(const RowVector &v)	{return v.abs();	}
		//! Computes a square value of each matrix element.
		RowVector Sqr()const;
		static RowVector Sqr(const RowVector &v)	{return v.Sqr();	}
		//! Computes a square root of each matrix element.
		RowVector Sqrt()const;
		static RowVector Sqrt(const RowVector &v)	{return v.Sqrt();	}
		//! Raises every matrix element to a power.
		RowVector pow(const Real & n)const			{return (*this)^n;	}
		//! Computes an exponent of each matrix element.
		RowVector exp()const;
		static RowVector exp(const RowVector &v)	{return v.exp();	}
		//! Changes each element to the smallest integer greater than or equal to that element.
		RowVector ceil()const;
		static RowVector ceil(const RowVector &v)	{return v.ceil();	}
		//! Changes each element to the largest integer less than or equal to that element.
		RowVector floor()const;
		static RowVector floor(const RowVector &v)	{return v.floor();	}
		//! Rounds each element to the nearest integral value.
		RowVector round()const;
		static RowVector round(const RowVector &v)	{return v.round();	}
		//! Computes the natural logarithm of each element in the matrix.
		RowVector log()const;
		static RowVector log(const RowVector &v)	{return v.log();	}
		//! Computes the logarithm of each element in the matrix to base 10.
		RowVector log10()const;
		static RowVector log10(const RowVector &v)	{return v.log10();	}
		//! Computes the natural logarithm of each element in the matrix plus 1.
		RowVector log1p()const;
		static RowVector log1p(const RowVector &v)	{return v.log1p();	}
		//! Computes the logarithm of each element in the matrix to base 2. log(n)/log(2) is log2.
		RowVector log2()const;
		static RowVector log2(const RowVector &v)	{return v.log2();	}
		//! Computes the exponent of each element in the matrix.
/*
		RowVector logb()const;
		static RowVector logb(const RowVector &v)	{return v.logb();	}
*/

		//! Compute the 2-norm of vector
		Real norm()const							{return norm(2);	}
		static Real norm(const RowVector &v)		{return v.norm();	}
		//! Compute the p-norm. v.norm(p)=sum(abs(v).^2)^(1/p)
		Real norm(int p)const;
		//! Compute the Frobenius norm for s = "fro" (equal to 2-norm)
		Real norm(const std::string &s)const;
		//! Maximum value of the vector
		static Real maxV(const RowVector &m,int &i)	{int j = 1; return BaseMatrix::maxV(m,j,i);	}
		Real maxV(int &i)const						{return maxV(*this,i);						}
		Real maxV()const							{return BaseMatrix::maxV(*this);			}
		static Real minV(const RowVector &m,int &i)	{int j = 1; return BaseMatrix::minV(m,j,i);	}
		Real minV(int &i)const						{return minV(*this,i);						}
		Real minV()const							{return BaseMatrix::minV(*this);			}

		void pushBack(const Real e);	//!Add element at the end
		void popBack();		//!Delete last element

		RowVector HorizontallyConcatenated(const RowVector&);
		RowVector operator |(const RowVector& mat)	{return HorizontallyConcatenated(mat);	}
		Matrix VerticallyConcatenated(const RowVector& mat);
		Matrix operator &(const RowVector& mat)		{return VerticallyConcatenated(mat);	}
		Matrix VerticallyConcatenated(const Matrix& mat);
		Matrix operator &(const Matrix& mat)		{return VerticallyConcatenated(mat);	}

		void insert(const Real & n, int k = 0);
		void resize(int n)			{(*(Matrix*)this).resize(n,1);		}
		void resize(int n, Real val){(*(Matrix*)this).resize(n,1,val);	}
		void swap(int i, int j);
		RowVector FromTo(int i, int j)const;
	};

	/*!
	  @brief Array of matrices.
	*/
	class Matrix3D{
	private:
		Matrix *_m3d;
		int _nArr;
		void copy(const Matrix3D & orig);
	public:
		Matrix3D();
		Matrix3D(int nArr);
		Matrix3D(int nArr, int nRow, Real val = ZERO);
		Matrix3D(int nArr, int nRow, int nCol, Real val = ZERO);
		Matrix3D(int nArr, int nRow, int nCol, const string &str);
		Matrix3D(const Matrix3D & orig);
		~Matrix3D();

		Matrix3D & operator=(const Matrix3D & orig);
		Matrix3D & operator=(Real val);
		const Matrix & operator()(int i)const;
		Matrix & operator()(int i);
		const Real & operator()(int i, int j, int k)const;
		Real & operator()(int i, int j, int k);

		int getnArr()const{ return _nArr; }
	};

	//-------------------------------------------------------------------------------------------
	ColumnVector cross(const ColumnVector & v1, const ColumnVector & v2);
	//-------------------------------------------------------------------------------------------
	RowVector cross(const RowVector & v1, const RowVector & v2);
	//-------------------------------------------------------------------------------------------
	Real dot(const ColumnVector & v1, const ColumnVector & v2);
	//-------------------------------------------------------------------------------------------
	Real dot(const RowVector & v1, const RowVector & v2);
	//-------------------------------------------------------------------------------------------
	//! Distance Between Two ColumnVectos
	Real distance(const ColumnVector & v1, const ColumnVector & v2);
	Matrix x_prod_matrix(const ColumnVector & v);
	Matrix x_prod_matrix(const RowVector & v);
	//! Return the Unit Vector for a given Vector v.
	ColumnVector UnitVector(const ColumnVector & v);
	Real my_pythag(Real a, Real b);
	//! Solve the Singular Value Decomposition. Adopted from numerical recipes
	ColumnVector SVDcmp(const math::Matrix &a, math::Matrix &u, math::Matrix &v, int maxitr = 30);
	math::Matrix SVDcmpN(const math::Matrix &a, math::Matrix &u, math::Matrix &v);
	ColumnVector SVDcmp(const math::Matrix &a);
	//! The correspondent back substitution for the singular value decomposition. Adopted from numerical recipes
	void svbksb(Matrix &u, ColumnVector &w, Matrix &v, ColumnVector &b, ColumnVector &x);
	void ErrorMsg(string error_text);

	void QR(const Matrix &M, Matrix &R, Matrix &Q);
	void QRZ(math::Matrix& X, math::Matrix& U);
	void QRZ(math::Matrix& X, math::Matrix& Y, math::ColumnVector& M);

	Matrix rk4(const Matrix & x, const Matrix & dxdt, Real t, Real h,
		Matrix (*xdot)(Real time, const Matrix & xin));
	void rkqc(Matrix & x, Matrix & dxdt, Real & t, Real htry,
		Real eps, Matrix & xscal, Real & hdid, Real & hnext,
		Matrix (*xdot)(Real time, const Matrix & xin));
	void odeint(Matrix (*xdot)(Real time, const Matrix & xin), Matrix & xo, 
		Real to, Real tf, Real eps, Real h1, Real hmin, int & nok, int & nbad,
		RowVector & tout, Matrix & xout, Real dtsav);
	ColumnVector Integ_Trap(const ColumnVector & present, ColumnVector & past, const Real dt);
	void SingularValueDecomposition(int m, int n, int lda, double *aa);
	Matrix SVD(const math::Matrix &a, math::Matrix &u, math::Matrix &v);
	Matrix chol(const Matrix A, const string &tri="upper");
	Matrix cholesky(const Matrix &A, ColumnVector &p);
	Matrix icholesky(const Matrix &A);

	ColumnVector linsolve(const Matrix &A, const ColumnVector &B);
	Real det(const Matrix & mat);
	Matrix random(int m, int n);
	Matrix magic(int m, int n);
	double rcond(const Matrix &A);
};

#endif
