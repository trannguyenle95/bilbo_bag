/**============================================================================
//==============================================================================

//	file:	PositionVector3D.cpp

//	author:	Fares J. Abu-Dakka
//	date:	May. 2012

//	Description: Implements a 3x1 Position Vector
//==============================================================================*/

#include <cmath>
#include <stdlib.h>
#include "Frame3D.h"

namespace math{


	//-------------------------------------------------------------------------------------------
	PositionVector3D::PositionVector3D(Real x, Real y, Real z)
	/*!
	  @brief Creates a 3D position vector
	  @param $ x $ with default value = 0.0
	  @param $ y $ with default value = 0.0
	  @param $ z $ with default value = 0.0
	*/
	{
		_v[0] = x;		_v[1] = y;		_v[2] = z;
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D::PositionVector3D(const math::ColumnVector &v)
	/*!
	  @brief Creates a 3D position vector
	  @param ColumnVector$ v $
	*/
	{
		if (v.getnRows() != 3 || v.getnColumns() != 1)
			throw Exception("Invalid Vector");
		_v[0] = v(1);		_v[1] = v(2);		_v[2] = v(3);
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D::PositionVector3D(const math::RowVector &v)
	/*!
	  @brief Creates a 3D position vector
	  @param RowVector $ v $
	*/
	{
		if (v.getnRows() != 1 || v.getnColumns() != 3)
			throw Exception("Invalid Vector");
		_v[0] = v(1);		_v[1] = v(2);		_v[2] = v(3);
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D::PositionVector3D(const PositionVector3D &p)
	/*!
	  @brief Copy constructor.
	  @param $ x $ with default value = 0.0
	  @param $ y $ with default value = 0.0
	  @param $ z $ with default value = 0.0
	*/
	{
		_v[0] = p(1);		_v[1] = p(2);		_v[2] = p(3);
	}
	//-------------------------------------------------------------------------------------------
	void PositionVector3D::set(const std::string &str)
	/*!
	  @brief Set the vector values to a values given through a string.
	  @param $ str $ string: its elements separated by spaces.
	*/
	{
		istringstream ss(str);
		ss >> _v[0];		ss >> _v[1];		ss >> _v[2];
	}
	//-------------------------------------------------------------------------------------------
	Real PositionVector3D::operator () (int i) const
	/*!
	  @brief Returns reference to vector element
	  @param $ i $ index in the vector $i\in \{1,2,3\} $
	  @return const reference to element
	*/
	{
		switch (i){
		case 1: return _v[i-1];
		case 2: return _v[i-1];
		case 3: return _v[i-1];
		default:
			throw Exception("Value of i must be 1, 2, or 3 in operator () of PositionVector3D");
		}
	}
	//-------------------------------------------------------------------------------------------
	Real & PositionVector3D::operator () (int i)
	/*!
	  @brief Returns reference to vector element
	  @param $ i $ index in the vector $i\in \{1,2,3\} $
	  @return const reference to element
	*/
	{
		switch (i){
		case 1: return _v[i-1];
		case 2: return _v[i-1];
		case 3: return _v[i-1];
		default:
			throw Exception("Value of i must be 1, 2, or 3 in operator () of PositionVector3D");
		}
	}
	//-------------------------------------------------------------------------------------------
	Real PositionVector3D::operator [] (int i) const
	/*!
	  @brief Returns reference to vector element
	  @param $ i $ index in the vector $i\in \{0,1,2\} $
	  @return const reference to element
	*/
	{
		switch (i){
		case 0: return _v[i];
		case 1: return _v[i];
		case 2: return _v[i];
		default:
			throw Exception("Value of i must be 0, 1, or 2 in operator () of PositionVector3D");
		}
	}
	//-------------------------------------------------------------------------------------------
	Real & PositionVector3D::operator [] (int i)
	/*!
	  @brief Returns reference to vector element
	  @param $ i $ index in the vector $i\in \{0,1,2\} $
	  @return const reference to element
	*/
	{
		switch (i){
		case 0: return _v[i];
		case 1: return _v[i];
		case 2: return _v[i];
		default:
			throw Exception("Value of i must be 0, 1, or 2 in operator () of PositionVector3D");
		}
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D & PositionVector3D::operator =(const PositionVector3D & P)
	/*!
	  @brief Assignment operator between two position vectors
	*/
	{
		_v[0] = P(1);		_v[1] = P(2);		_v[2] = P(3);
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D & PositionVector3D::operator =(const math::ColumnVector & v)
	/*!
	  @brief Assignment operator between Position vectors and Column vector
	  @param ColumnVector $ \mathbf{v} $.
	*/
	{
		if (v.getnRows() != 3 || v.getnColumns() != 1)
			throw Exception("Invalid Vector");
		_v[0] = v(1);		_v[1] = v(2);		_v[2] = v(3);
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D & PositionVector3D::operator =(const math::RowVector & v)
	/*!
	  @brief Assignment operator between Position vectors and Row vector
	  @param RowVector $ \mathbf{v} $.
	*/
	{
		if (v.getnRows() != 1 || v.getnColumns() != 3)
			throw Exception("Invalid Vector");
		_v[0] = v(1);		_v[1] = v(2);		_v[2] = v(3);
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D & PositionVector3D::operator =(Real x)
	/*!
	  @brief Assignment operator. All vector elements assigned to a given value
	  @param Real $ x $.
	*/
	{
		_v[0] = x;		_v[1] = x;		_v[2] = x;
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	bool PositionVector3D::operator ==(const PositionVector3D &m)const
	/*!
	  @brief Compares \b a and \b b for equality.
	  @param a
	  @param b
	  @return True if a equals b, false otherwise.
	*/
	{
		return _v[0] == m(1) && _v[1] == m(2) && _v[2] == m(3);
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D PositionVector3D::operator / ( Real s) const
	/*!
	  @brief Scalar division.
	  @param Real $ x $.
	*/
	{
		if (s == 0)
			throw Exception("error: PositionVector3D::operator / ( Real s): Division by zero...");
		return PositionVector3D(_v[0] / s, _v[1] / s, _v[2] / s);
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D PositionVector3D::operator * ( Real s) const
	/*!
	  @brief Scalar multiplication. P = P1 * s
	  @param Real $ x $.
	*/
	{
		return PositionVector3D(_v[0] * s, _v[1] * s, _v[2] * s);
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D operator * (Real s, const PositionVector3D& v)
	/*!
	  @brief Scalar multiplication. P = s * P1
	  @param Real $ x $.
	*/
	{
		return v*s;
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D PositionVector3D::operator -(const PositionVector3D &v)const
	/*!
	  @brief Position vectors subtraction. P = P1 - P2
	  @param PositionVector3D $ \mathbf{v} $.
	*/
	{
		return PositionVector3D(_v[0] - v(1), _v[1] - v(2), _v[2] - v(3));
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D PositionVector3D::operator +(const PositionVector3D &v)const
	/*!
	  @brief Position vectors addition. P = P1 + P2
	  @param PositionVector3D $ \mathbf{v} $.
	*/
	{
		return PositionVector3D(_v[0] + v(1), _v[1] + v(2), _v[2] + v(3));
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D PositionVector3D::operator-() const
	/*!
	  @brief Unary minus
	*/
	{
		return PositionVector3D(-_v[0], -_v[1], -_v[2]);
	}
	void PositionVector3D::reverseSign()
	/*!
	  @brief reverse the sign for the this object
	*/
	{
		_v[0] = -_v[0];
		_v[1] = -_v[1];
		_v[2] = -_v[2];
	}
	//-------------------------------------------------------------------------------------------
	////// this function needs to be checked
	const Real PositionVector3D::operator * ( const PositionVector3D &s) const{
		return _v[0]*s(1) + _v[1]*s(2) + _v[2]*s(3);
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D operator * (const Matrix &s, const PositionVector3D& v)
	/*!
	  @brief Matrix multiplication with Position vector. P = mat * v
	*/
	{
		if (s.getnRows() != 3 || s.getnColumns() != 3)
			throw Exception("Matrix is not 3x3");
		return PositionVector3D(s(1,1)*v(1) + s(1,2)*v(2) + s(1,3)*v(3),
			s(2,1)*v(1) + s(2,2)*v(2) + s(2,3)*v(3),
			s(3,1)*v(1) + s(3,2)*v(2) + s(3,3)*v(3));
	}
	//-------------------------------------------------------------------------------------------
	const Real operator * (const RowVector &s, const PositionVector3D& v)
	/*!
	  @brief Row vector multiplication with Position vector. P = row * v
	*/
	{
		if (s.getnColumns() != 3)
			throw Exception("RowVector is not 1x3");
		return s(1)*v(1) + s(2)*v(2) + s(3)*v(3);
	}
 	//-------------------------------------------------------------------------------------------
	double PositionVector3D::Normalize(double eps)
	/*!
	  @brief Normalize this object with respect to its unitvector
	  @return the norm of this vector.
	*/
	{
		double v = this->norm2();
		if (v < eps) {
			*this = PositionVector3D(1,0,0);
			return v;
		} else {
			*this = (*this)/v;
			return v;
		}
	}
	//-------------------------------------------------------------------------------------------
	Real PositionVector3D::norm2() const
	/*!
	  @brief Returns the Euclidean norm (2-norm) of the vector
	  @return the norm
	*/
	{
		return sqrt(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]);
	}
	//-------------------------------------------------------------------------------------------
	Real PositionVector3D::norm1() const
	/*!
	  @brief Returns the Manhattan norm (1-norm) of the vector
	  @return the norm
	*/
	{
		return fabs(_v[0])+fabs(_v[1])+fabs(_v[2]);
	}
	//-------------------------------------------------------------------------------------------
	Real PositionVector3D::normInf() const
	/*!
	  @brief Returns the infinite norm ($\inf$-norm) of the vector
	  @return the norm
	*/
	{
		Real max = fabs(_v[0]);
		if (fabs(_v[1]) > max)
			max = fabs(_v[1]);
		if (fabs(_v[2]) > max)
			max = fabs(_v[2]);
		return max;
	}
	//-------------------------------------------------------------------------------------------
	bool PositionVector3D::equal(const PositionVector3D& p, Real precision)const
	/*!
	  @brief Compares rotations with a given precision
	  @param \b Rot: Rotation to compare with
	  @param \b precision: The precision to use for testing
	  @return True if all elements are less than \b precision apart.
	  
	  Performs an element wise comparison. Two elements are considered equal if the difference
	  are less than \b precision.
	*/
	{
		return (EQUAL(_v[0],p[0],precision)&&
                EQUAL(_v[1],p[1],precision)&&
                EQUAL(_v[2],p[2],precision)   );
	}
	//-------------------------------------------------------------------------------------------
	bool equal(const PositionVector3D& p1, const PositionVector3D& p2, Real precision)
	/*!
	  @brief Compares rotations with a given precision
	  @param \b Rot: Rotation to compare with
	  @param \b precision: The precision to use for testing
	  @return True if all elements are less than \b precision apart.
	  
	  Performs an element wise comparison. Two elements are considered equal if the difference
	  are less than \b precision.
	*/
	{
		return p1.equal(p2,precision);
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D PositionVector3D::ElementWiseDiv(const PositionVector3D &v)const
	/*!
	  @brief Element wise division
	*/
	{
		PositionVector3D vec;
		for (int i=1; i<=3; i++)
			if (::abs(v(i)) != 0)
				vec(i) = (*this)(i) / v(i);
			else
				throw Exception("error: ElementWiseDiv: division by zero...");
		return vec;
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D PositionVector3D::ElementWiseMul(const PositionVector3D &p1, const PositionVector3D &p2)
	/*!
	  @brief Element wise multiplication
	*/
	{
		return PositionVector3D(p1(1)*p2(1),p1(2)*p2(2),p1(3)*p2(3));
	}
	//-------------------------------------------------------------------------------------------
	Matrix PositionVector3D::skewMatrix()const
	/*!
	  @brief Creates a skew symmetric matrix from a PositionVector3D. Also
	  known as the cross product matrix of v.
	*/
	{
		Matrix s(3,3);
		s(1,2) = - _v[2];	s(2,1) =   _v[2];
		s(1,3) =   _v[1];	s(3,1) = - _v[1];
		s(2,3) = - _v[0];	s(3,2) =   _v[0];
		return s;
	}
	//-------------------------------------------------------------------------------------------
	ostream & operator << (ostream & o, const PositionVector3D & m){
		o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << m(1) << " ";
		o << endl;
		o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << m(2) << " ";
		o << endl;
		o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << m(3) << " ";
		o << endl;
		o << endl;
		return o;
	}
	//-------------------------------------------------------------------------------------------
	void PositionVector3D::print (char const* ch, int prec, int w) const{
		int n = 0;
		if (ch){
			string s(ch);
			std::cout << ch << ": ";
			n = s.size();
		}
		cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << _v[0] << " ";
		cout << endl << std::setw(n+2) << " ";
		cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << _v[1] << " ";
		cout << endl << std::setw(n+2) << " ";
		cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << _v[2] << " ";
		cout << endl << std::setw(n+2) << " ";
		cout << endl;
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D & PositionVector3D::operator >> (const Real *val ){
		_v[0] = val[0];		_v[1] = val[1];		_v[2] = val[2];
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	Real PositionVector3D::SumSquare()const{
		Real sum = ZERO;
		sum = SQR(_v[0]) + SQR(_v[1]) + SQR(_v[2]);
		return sum;
	}
	//-------------------------------------------------------------------------------------------
	/*std::istream& operator >> (std::istream& is, PositionVector3D& v)
	{
		IOTrace("Stream input Vector (vector or ZERO)");
		char storage[10];
		EatWord(is, "[]", storage, 10);
		if (strlen(storage) == 0) {
			Eat(is, '[');
			is >> v(1);
			Eat(is, ',');
			is >> v(2);
			Eat(is, ',');
			is >> v(3);
			EatEnd(is, ']');
			IOTracePop();
			return is;
		}
		if (strcmp(storage, "ZERO") == 0) {
			v = PositionVector3D();
			IOTracePop();
			return is;
		}
		throw Error_Frame_Vector_Unexpected_id();
	}*/
	//-------------------------------------------------------------------------------------------
	Real distance(const PositionVector3D & v1, const PositionVector3D & v2)
	/*!
	  @brief Calculates the distance between two position vectors
	  @param $ \mathbf{v1} $
	  @param $ \mathbf{v2} $
	  @return Real $ distance $
	*/
	{
		Real ss = ZERO;
		for (int i=1; i<=3; i++)
			ss += (v1(i)-v2(i))*(v1(i)-v2(i));

		return sqrt(ss);
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D cross(const PositionVector3D& v1, const PositionVector3D& v2)
	/*!
	  @brief Calculates the 3D vector cross product $ \mathbf{v1} \times \mathbf{v2} $
	  @param $ \mathbf{v1} $
	  @param $ \mathbf{v2} $
	  @return the 3D vector cross product $ \mathbf{v1} \times \mathbf{v2} $
	  
	  The 3D vector cross product is defined as:
	  $
	  \mathbf{v1} \times \mathbf{v2} = \left[\begin{array}{c}
		v1_y * v2_z - v1_z * v2_y \\
		v1_z * v2_x - v1_x * v2_z \\
		v1_x * v2_y - v1_y * v2_x
	  \end{array}\right]
	  $
	*/
	{
		PositionVector3D v3(3);
		v3(1) = v1(2)*v2(3) - v2(2)*v1(3);
		v3(2) = v2(1)*v1(3) - v1(1)*v2(3);
		v3(3) = v1(1)*v2(2) - v2(1)*v1(2);
		return v3;
	}
	//-------------------------------------------------------------------------------------------
	void cross(const PositionVector3D& v1, const PositionVector3D& v2, PositionVector3D& v3){
		v3 = cross(v1,v2);
	}
	//-------------------------------------------------------------------------------------------
	Real dot(const PositionVector3D& v1, const PositionVector3D& v2)
	/*!
	  @brief Calculates the dot product $ \mathbf{v1} . \mathbf{v2} $
	  @param $ \mathbf{v1} $
	  @param $ \mathbf{v2} $
	  @return the dot product $ \mathbf{v1} . \mathbf{v2} $
	*/
	{
		return v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3);
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D UnitVector(const PositionVector3D & v)
	/*!
	  @brief Calculates the unit vector $ \hat{\mathbf{v}} $
	  @param $ \mathbf{v} $
	  @return $ \hat{\mathbf{v}} $
	*/
	{
		return v*(1/sqrt(dot(v,v)));
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D normalize(const PositionVector3D& v)
	/*!
	  @brief Returns the normalized vector $\mathbf{n}=\frac{\mathbf{v}}{\|\mathbf{v}\|} $.
	  In case $ \|\mathbf{v}\| = 0$ the zero vector is returned.
	  @param $ \mathbf{v} $ which should be normalized
	  @return the normalized vector $ \mathbf{n} $
	*/
	{
		Real length = v.norm2();
		if (length != 0)
			return PositionVector3D(v(1)/length, v(2)/length, v(3)/length);
		else
			return PositionVector3D(0,0,0);
	}
	//-------------------------------------------------------------------------------------------
	Real angle(const PositionVector3D& v1, const PositionVector3D& v2, const PositionVector3D& n)
	/*!
	  @brief Calculates the angle from $ \mathbf{v1}$ to $ \mathbf{v2} $
	  around the axis defined by $ \mathbf{v1} \times \mathbf{v2} $ with n
	  determining the sign.
	  @param $ \mathbf{v1} $
	  @param $ \mathbf{v2} $
	  @param $ \mathbf{n} $
	  @return the angle
	*/
	{
		const PositionVector3D nv1 = normalize(v1);
		const PositionVector3D nv2 = normalize(v2);
		const PositionVector3D nn  = normalize(n);
		return atan2(dot(nn, cross(nv1, nv2)), dot(nv1, nv2));
	}
	//-------------------------------------------------------------------------------------------
	Real angle(const PositionVector3D& v1, const PositionVector3D& v2){
		PositionVector3D n = cross(v1, v2);
		return angle(v1,v2,n);
	}
	void random(PositionVector3D &v){
		math::random(v[1]);
		math::random(v[2]);
		math::random(v[3]);
	}
};
