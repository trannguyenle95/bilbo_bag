//***************************************************************************************************
//***************************************************************************************************

//	file:	Transform3D.cpp

//	author:	Fares J. Abu-Dakka
//	date:	May. 2012

//	Description: Transform3D represents a 3D transformation object
//				 It represents the transformation internally in a 4x4 Square Matrix.
//***************************************************************************************************

#include <math.h>
#include "Frame3D.h"
#include "quaternion.hpp"

namespace math{

	//***************************************************************************************************
	Transform3D::Transform3D() : _P(), _R()
	/*!
	  @brief Default Constructor.
	  Initializes with 0 translation and Identity matrix as rotation
	*/
	{ }
	//***************************************************************************************************
	Transform3D::Transform3D(const PositionVector3D& P, const RotationMatrix& R) : _P(P), _R(R)
	/*!
	  @brief Constructs a homogeneous transform.
	  @param $\mathbf{P}$ A 3x1 translation vector
	  @param $\mathbf{R}$ A 3x3 rotation matrix
	*/
	{ }
	//***************************************************************************************************
	Transform3D::Transform3D(const Real * CartPos) : _P(), _R()
	/*!
	  @brief Constructs a homogeneous transform.
	  @param CartPos is a 12 elements array.
	*/
	{
		_R(1,1) = CartPos[0];	_R(1,2) = CartPos[1];	_R(1,3) = CartPos[2];	_P = CartPos[3];
		_R(2,1) = CartPos[4];	_R(2,2) = CartPos[5];	_R(2,3) = CartPos[6];	_P = CartPos[7];
		_R(3,1) = CartPos[8];	_R(3,2) = CartPos[9];	_R(3,3) = CartPos[10];	_P = CartPos[11];
	}
	//***************************************************************************************************
	Transform3D::Transform3D(const RotationMatrix& R) : _P(PositionVector3D()), _R(R)
	/*!
	  @brief Constructs a homogeneous transform.
	  @param $\mathbf{R}$ A 3x3 rotation matrix
	  A homogeneous transform with a rotation of \b R and a translation of zero.
	*/
	{ }
	//***************************************************************************************************
	Transform3D::Transform3D(const PositionVector3D& P) : _P(P), _R(RotationMatrix())
	/*!
	  @brief Constructs a homogeneous transform.
	  @param $\mathbf{P}$ A 3x1 translation vector
	  A homogeneous transform with a rotation of zero and a translation of \b d.
	*/
	{ }
	//***************************************************************************************************
	const Transform3D Transform3D::DH(Real alpha, Real a, Real d, Real theta)
	/*!
	  @brief Constructs a homogeneous transform using the original Denavit-Hartenberg notation
	  @param $ \alpha_i $ is angle between $Z_{i}$ to $Z_{i+1}$ along $X_{i}$
	  @param $ a_i $ is the distance from $Z_{i}$ to $Z_{i+1}$ along $X_{i}$
	  @param $ d_i $ is the distance from $X_{i}$ to $X_{i+1}$ along $Z_{i+1}$
	  @param $ \theta_i $ is the angle between $X_{i}$ to $X_{i+1}$ along $X_{i+1}$
	  @return $ ^{i-1}\mathbf{T}_i $
	  
	  $
	  ^{i-1}\mathbf{T}_i=
	  \left[
		\begin{array}{cccc}
			c\theta_i & -s\theta_i c\alpha_i &  s\theta_i s\alpha_i & a_i c\theta_i \\
			s\theta_i &  c\theta_i c\alpha_i & -c\theta_i s\alpha_i & a_i s\theta_i \\
			0         &  s\alpha_i           &  c\alpha_i           & d_i \\
			0         &  0                   & 0                    & 1
	  \end{array}
	  \right]
	  $	
	  Reference: Denavit, J. and Hartenberg, R. S., A kinematic
      notation for lower-pair mechanisms based on matrices, ASME
      Journal of Applied Mechanics, 23:215-221, 1955.
	*/
	{
		double ct = cos(theta), st = sin(theta), sa = sin(alpha), ca = cos(alpha);
		return Transform3D(	
			PositionVector3D( a*ct,      a*st,  d),
			RotationMatrix(	
				ct,    -st*ca,   st*sa,
				st,     ct*ca,  -ct*sa,
				0 ,        sa,      ca		) );
	}
	//***************************************************************************************************
	const Transform3D Transform3D::modifiedDH(Real alpha, Real a, Real d, Real theta)
	/*!
	  @brief Constructs a homogeneous transform using the modified Denavit-Hartenberg notation
	  @param $ \alpha_i $ is angle between $Z_{i-1}$ to $Z_i$ along $X_{i-1}$
	  @param $ a_i $ is the distance from $Z_{i-1}$ to $Z_i$ along $X_{i-1}$
	  @param $ d_i $ is the distance from $X_{i-1}$ to $X_i$ along $Z_i$
	  @param $ \theta_i $ is the angle between $X_{i-1}$ to $X_i$ along $X_i$
	  @return $ ^{i-1}\mathbf{T}_i$
	  
	  @note The modified (Craig) Denavit-Hartenberg notation differs from
	  the original Denavit-Hartenberg notation and is given as
	  
	  $
	  ^{i-1}\mathbf{T}_i =
	  \left[
		\begin{array}{cccc}
			c\theta_i & -s\theta_i & 0 & a_{i-1} \\
			s\theta_i c\alpha_{i-1} & c\theta_i c\alpha_{i-1} & -s\alpha_{i-1} & -s\alpha_{i-1}d_i \\
			s\theta_i s\alpha_{i-1} & c\theta_i s\alpha_{i-1} &  c\alpha_{i-1} &  c\alpha_{i-1}d_i \\
			0 & 0 & 0 & 1
		\end{array}
	  \right]
	  $
	  Reference: Craig, J. J.,Introduction to Robotics: Mechanics and
	  Control, Addison-Wesley, isbn:0-201-10326-5, 1986.
	*/
	{
		double ct = cos(theta), st = sin(theta), sa = sin(alpha), ca = cos(alpha);
		return Transform3D( 
			PositionVector3D( a, -sa * d, ca * d), RotationMatrix(
														ct   ,    -st,    0 ,
														st*ca,  ct*ca,   -sa,
														st*sa,  ct*sa,    ca 	) );
	}
	//***************************************************************************************************
	Transform3D Transform3D::identity()
	/*!
	  @brief Constructs the identity transform
	  @return the identity transform
	  
	  $
	  \mathbf{T} =
	  \left[
		\begin{array}{cccc}
			1 & 0 & 0 & 0\\
			0 & 1 & 0 & 0\\
			0 & 0 & 1 & 0\\
			0 & 0 & 0 & 1
		\end{array}
	  \right]
	  $	
	*/
	{
		return Transform3D(PositionVector3D::zero(), RotationMatrix::identity());
	}
	//***************************************************************************************************
	bool Transform3D::equal(const Transform3D& T, Real precision)
	/*!
	  @brief Compares the transformations with a given precision
	  @param \b Rot: Rotation to compare with
	  @param \b precision: The precision to use for testing
	  @return True if all elements are less than \b precision apart.
	  
	  Performs an element wise comparison. Two elements are considered equal if the difference
	  are less than \b precision.
	*/
	{
		if(!(*this).R().equal(T.R(),precision) )
			return false;
		for(int i=1;i<=3;i++)
			if( ::abs((*this).P()(i)-T.P()(i)>precision) )
				return false;
		return true;
	}
	//***************************************************************************************************
	bool equal(const Transform3D& a,const Transform3D& b,double precision)
	/*!
	  @brief Compares the transformations with a given precision
	  @param \b Rot: Rotation to compare with
	  @param \b precision: The precision to use for testing
	  @return True if all elements are less than \b precision apart.

	  Performs an element wise comparison. Two elements are considered equal if the difference
	  are less than \b precision.
	*/
	{
		if(!a.R().equal(b.R(),precision) )
			return false;
		for(int i=1;i<=3;i++)
			if( ::abs(a.P()(i)-b.P()(i)>precision) )
				return false;
		return true;
	}
	//***************************************************************************************************
	const RotationMatrix & Transform3D::R()const
	/*!
	  @brief Gets the rotation part $ \mathbf{R} $ from $ \mathbf{T} $
	  @return $ \mathbf{R} $
	*/
	{
		return _R;
	}
	//***************************************************************************************************
	RotationMatrix & Transform3D::R()
	/*!
	  @brief Gets the rotation part $ \mathbf{R} $ from $ \mathbf{T} $
	  @return $ \mathbf{R} $
	*/
	{
		return _R;
	}
	//***************************************************************************************************
	const PositionVector3D & Transform3D::P()const
	/*!
	  @brief Gets the position part $ \mathbf{P} $ from $ \mathbf{T} $
	  @return $ \mathbf{P} $
	*/
	{
		return _P;
	}
	//***************************************************************************************************
	PositionVector3D & Transform3D::P()
	/*!
	  @brief Gets the position part $ \mathbf{P} $ from $ \mathbf{T} $
	  @return $ \mathbf{P} $
	*/
	{
		return _P;
	}
	//***************************************************************************************************
	Real & Transform3D::operator()(int i, int j)
	/*!
	  @brief Returns matrix element reference
	  @param row [in] row, row must be @f$ < 3 @f$
	  @param col [in] col, col must be @f$ < 4 @f$
	  @return reference to matrix element
	*/
	{
		if ( i<1 || i>3 || j<1 || j>4 )
			throw Exception("RotationMatrix::operator () (int i, int j): Invalid Index");
		if(i <= 3 && j <= 3)
			return _R(i,j);
		else
			return _P(i);
	}
	//***************************************************************************************************
	Real Transform3D::operator()(int i, int j) const
	/*!
	  @brief Returns const matrix element reference
	  @param row [in] row, row must be @f$ < 3 @f$
	  @param col [in] col, col must be @f$ < 4 @f$
	  @return const reference to matrix element
	*/
	{
		if ( i<1 || i>3 || j<1 || j>4 )
			throw Exception("RotationMatrix::operator () (int i, int j): Invalid Index");
	    if(i <= 3 && j <= 3)
	        return _R(i,j);
	    else
	        return _P(i);
	}
	//***************************************************************************************************
	Transform3D Transform3D::invMult(const Transform3D& t1, const Transform3D& t2)
	/*!
	  @brief computes the inverse of t1 and multiplies it with t2. t0 = inv(t1) * t2
	*/
	{
		Transform3D T0;
		const Real p0 = t1.P()(1), p1 = t1.P()(2), p2 = t1.P()(3);
		const Real r01 = t1.R()(1,2);
		const Real r12 = t1.R()(2,3);
		const Real r02 = t1.R()(1,3);

		T0(1,4) = (-p0 + t2.P()(1))*t1.R()(1,1) +
				 (-p1 + t2.P()(2))*t1.R()(2,1) +
				 (-p2 + t2.P()(3))*t1.R()(3,1);

		T0(2,4) = (-p0 + t2.P()(1))*r01 +
				 (-p1 + t2.P()(2))*t1.R()(2,2) +
				 (-p2 + t2.P()(3))*t1.R()(3,2);

		T0(3,4) = (-p0 + t2.P()(1))*r02 +
				 (-p1 + t2.P()(2))*r12 +
				 (-p2 + t2.P()(3))*t1.R()(3,3);

		T0(1,2) = t1.R()(1,1)*t2.R()(1,2) + t1.R()(2,1)*t2.R()(2,2) + t1.R()(3,1)*t2.R()(3,2);
		T0(1,3) = t1.R()(1,1)*t2.R()(1,3) + t1.R()(2,1)*t2.R()(2,3) + t1.R()(3,1)*t2.R()(3,3);
		T0(1,1) = t1.R()(1,1)*t2.R()(1,1) + t1.R()(2,1)*t2.R()(2,1) + t1.R()(3,1)*t2.R()(3,1);

		T0(2,1) = r01*t2.R()(1,1) + t1.R()(2,2)*t2.R()(2,1) + t1.R()(3,2)*t2.R()(3,1);
		T0(2,3) = r01*t2.R()(1,3) + t1.R()(2,2)*t2.R()(2,3) + t1.R()(3,2)*t2.R()(3,3);
		T0(2,2) = r01*t2.R()(1,2) + t1.R()(2,2)*t2.R()(2,2) + t1.R()(3,2)*t2.R()(3,2);

		T0(3,1) = r02*t2.R()(1,1) + r12*t2.R()(2,1) + t1.R()(3,3)*t2.R()(3,1);
		T0(3,2) = r02*t2.R()(1,2) + r12*t2.R()(2,2) + t1.R()(3,3)*t2.R()(3,2);
		T0(3,3) = r02*t2.R()(1,3) + r12*t2.R()(2,3) + t1.R()(3,3)*t2.R()(3,3);
		return T0;
	}
	//***************************************************************************************************
	Transform3D Transform3D::Inverse() const
	/*!
	  @brief Calculates $ ^{b}\mathbf{T}_a = ^{a}\mathbf{T}_b^{-1} $
	  @return $ ^{b}\mathbf{T}_a = ^{a}\mathbf{T}_b^{-1} $
	  $
	  ^{a}\mathbf{T}_b^{-1} =
	  \left[
		\begin{array}{cc}
			^{a}\mathbf{R}_b^{T} & - ^{a}\mathbf{R}_b^{T}\,  ^{a}\mathbf{d}_b \\
			\begin{array}{ccc}0 & 0 & 0 
			\end{array}& 1
		\end{array}
	  \right]
	  $
	*/
	{
		RotationMatrix Rt = ~((*this)._R);
		PositionVector3D P = -Rt * (*this)._P;

		return Transform3D(P,Rt);
	}
	//***************************************************************************************************
	Wrench Transform3D::Inverse(const Wrench& arg) const
	{
		Wrench tmp;
		tmp.Force() =  _R.Inverse(arg.Force());
		tmp.Torque() = _R.Inverse(arg.Torque()-_P*arg.Force());
		return tmp;
	}
	//***************************************************************************************************
	Twist Transform3D::Inverse(const Twist& arg) const
	{
		Twist tmp;
		tmp.Rot() = _R.Inverse(arg.Rot());
		tmp.Vel() = _R.Inverse(arg.Vel()-_P*arg.Rot());
		return tmp;
	}
	//***************************************************************************************************
	Transform3D & Transform3D::operator =(const Transform3D & T0)
	/*!
	  @brief Sets Transform3D equal to another Transform3D.
	*/
	{
		(*this)._R = T0.R();
		(*this)._P = T0.P();
		return (*this);
	}
	//***************************************************************************************************
	Transform3D & Transform3D::operator =(const RotationMatrix & Rot)
	/*!
	  @brief Sets the Rotation part equal to a given one.
	  @param $\mathbf{Rot}$ A 3x3 rotation matrix
	*/
	{
		(*this)._R = Rot;
		return (*this);
	}
	//***************************************************************************************************
	Transform3D & Transform3D::operator =(const PositionVector3D & P)
	/*!
	  @brief Sets the translation part equal to a given one.
	  @param $\mathbf{P}$ A 3x1 translation vector
	*/
	{
		_P = P;
		return *this;
	}
	//***************************************************************************************************
	const Transform3D Transform3D::operator*(const Transform3D& Ti) const
	/*!
	  @brief Calculates $^{a}\mathbf{T}_c =\,  ^{a}\mathbf{T}_b \: ^{b}\mathbf{T}_c $
	  @param $ ^{a}\mathbf{T}_b $
	  @param $ ^{b}\mathbf{T}_c $
	  @return $ ^{a}\mathbf{T}_c $
	  
	  $
	  ^{a}\mathbf{T}_c =
	  \left[
		\begin{array}{cc}
		^{a}\mathbf{R}_b \:  ^{b}\mathbf{R}_c & ^{a}\mathbf{d}_b + ^{a}\mathbf{R}_b \: ^{b}\mathbf{d}_c \\
		\begin{array}{ccc}0 & 0 & 0\end{array} & 1
		\end{array}
	  \right]
	  $	
	*/
	{
		return Transform3D(_P + _R * Ti.P(), _R * Ti.R());
	}
	//***************************************************************************************************
	const PositionVector3D Transform3D::operator*( const PositionVector3D& Pos) const
	/*!
	  @brief Calculates $ ^{a}\mathbf{P} = \, ^{a}\mathbf{T}_b \,  ^{b}\mathbf{p} $ 
			 thus transforming point $ \mathbf{p} $ from frame $ b $ to frame $ a $
	  @param $ ^{b}\mathbf{p} $
	  @return $ ^{a}\mathbf{p} $	
	*/
	{
		return _R * Pos + _P;
	}
	//***************************************************************************************************
	const Transform3D Transform3D::operator / ( Real s) const
	/*!
	  @brief Calculates $\mathbf{Tnew} =\,  \mathbf{T} / s $
	  @param $ s $
	  @return $ \mathbf{Tnew} $
	*/
	{
		return Transform3D(	_P/s, _R/s );
	}
	//***************************************************************************************************
	const Transform3D Transform3D::operator * ( Real s) const
	/*!
	  @brief Calculates $\mathbf{Tnew} =\,  \mathbf{T} * s $
	  @param $ s $
	  @return $ \mathbf{Tnew} $
	*/
	{
		return Transform3D(	_P*s, _R*s );
	}
	//***************************************************************************************************
	const Transform3D operator * (Real s, const Transform3D& v)
	/*!
	  @brief Calculates $\mathbf{Tnew} =\,  s * \mathbf{T} $
	  @param $ s $
	  @return $ \mathbf{Tnew} $
	*/
	{
		return v * s;
	}
	//***************************************************************************************************
	const Transform3D Transform3D::operator -(const Transform3D &Ti)const
	/*!
	  @brief Calculates $\mathbf{Tnew} =\,  \mathbf{T} - \mathbf{Ti} $
	  @param $ \mathbf{Ti} $
	  @return $ \mathbf{Tnew} $
	*/
	{
		return Transform3D(	P() - Ti.P(), R() - Ti.R() );
	}
	//***************************************************************************************************
	const Transform3D Transform3D::operator +(const Transform3D &Ti)const
	/*!
	  @brief Calculates $\mathbf{Tnew} =\,  \mathbf{T} + \mathbf{Ti} $
	  @param $ \mathbf{Ti} $
	  @return $ \mathbf{Tnew} $
	*/
	{
		return Transform3D(	P() + Ti.P(), R() + Ti.R() );
	}
	//***************************************************************************************************
	const Transform3D & Transform3D::operator *=(Real s)
	/*!
	  @brief Calculates $\mathbf{T} *= s $
	  @param $ s $
	  @return $ \mathbf{T} $
	*/
	{ 
		return (*this) = (*this) * s;
	}
	//***************************************************************************************************
	const Transform3D & Transform3D::operator /=(Real s)
	/*!
	  @brief Calculates $\mathbf{T} /= s $
	  @param $ s $
	  @return $ \mathbf{T} $
	*/
	{
		return (*this) = (*this) / s;
	}
	//***************************************************************************************************
	const Transform3D & Transform3D::operator +=(const Transform3D &Ti)
	/*!
	  @brief Calculates $\mathbf{T} += \mathbf{Ti} $
	  @param $ \mathbf{Ti} $
	  @return $ \mathbf{T} $
	*/
	{
		return (*this) = (*this) + Ti;
	}
	//***************************************************************************************************
	const Transform3D & Transform3D::operator -=(const Transform3D &Ti)
	/*!
	  @brief Calculates $\mathbf{T} -= \mathbf{Ti} $
	  @param $ \mathbf{Ti} $
	  @return $ \mathbf{T} $
	*/
	{
		return (*this) = (*this) - Ti;
	}
	//***************************************************************************************************
	ostream & operator << (ostream & o, const Transform3D & m)
	/*
      @brief Writes Transformation matrix to stream
      
      @param os output stream to use
      @param r rotation matrix to print
      @return the updated output stream
	*/
	{
		for (int i=1; i<=3; i++){
			for (int j=1; j<=4; j++)
				o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << m(i,j) << " ";
			o << endl;
		}
		o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << ZERO << " ";
		o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << ZERO << " ";
		o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << ZERO << " ";
		o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << ONE  << " ";
		o << endl;
		o << endl;
		return o;
	}
	//***************************************************************************************************
	void Transform3D::print (char const* ch, int prec, int w) const
	/*
      @brief Prints rotation matrix to screen
      
      @param ch variable name
	  @param prec output precision
	  @param w width of the element
	*/
	{
		int n = 0;
		if (ch){
			string s(ch);
			std::cout << ch << ": ";
			n = s.size();
		}
		for (int i=1; i<=3; i++){
			for (int j=1; j<=4; j++)
				cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << (*this)(i,j) << " ";
			cout << endl << std::setw(n+2) << " ";
		}
		cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << ZERO << " ";
		cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << ZERO << " ";
		cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << ZERO << " ";
		cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << ONE  << " ";
		cout << endl;
		cout << endl;
	}
	//***************************************************************************************************
	Wrench Transform3D::operator * (const Wrench& arg) const
	//! Transformation of both the force reference point
	//! and of the base to which the wrench is expressed.
	//! look at Rotation*Wrench operator for a transformation
	//! of only the base to which the twist is expressed.
	//!
	//! Complexity : 24M+18A
	{
		Wrench tmp;
		tmp.Force()  = _R*arg.Force();
		tmp.Torque() = _R*arg.Torque() + _P*tmp.Force();
		return tmp;
	}
	//***************************************************************************************************
	Twist Transform3D::operator * (const Twist& arg) const
	//! Transformation of both the velocity reference point
	//! and of the base to which the twist is expressed.
	//! look at Rotation*Twist for a transformation of only the
	//! base to which the twist is expressed.
	//!
	//! Complexity : 24M+18A
	{
		Twist tmp;
		tmp.Rot() = _R*arg.Rot();
		tmp.Vel() = _R*arg.Vel()+_P*tmp.Rot();
		return tmp;
	}
	//***************************************************************************************************
	Real * Transform3D::getCartPos()const
	/*!
	  @brief Fill and returns CartPos pinter from this Transformation3D object.
	  @return $ \mathbf{CartPos} $ pointer of 12 elements
	*/
	{
		Real * CartPos;
		CartPos = new Real[12];
		CartPos[0]=(*this)(1,1) ;	CartPos[1]=(*this)(1,2);	CartPos[2] =(*this)(1,3);	CartPos[3] =(*this)(1,4);
		CartPos[4]=(*this)(2,1) ;	CartPos[5]=(*this)(2,2);	CartPos[6] =(*this)(2,3);	CartPos[7] =(*this)(2,4);
		CartPos[8]=(*this)(3,1) ;	CartPos[9]=(*this)(3,2);	CartPos[10]=(*this)(3,3);	CartPos[11]=(*this)(3,4);
		return CartPos;
	}
	//***************************************************************************************************
	void Transform3D::getCartPos(Real CartPos[12])
	/*!
	  @brief Fill the CartPos array from this Transformation3D object.
	*/
	{
		CartPos[0]=(*this)(1,1) ;	CartPos[1]=(*this)(1,2);	CartPos[2] =(*this)(1,3);	CartPos[3] =(*this)(1,4);
		CartPos[4]=(*this)(2,1) ;	CartPos[5]=(*this)(2,2);	CartPos[6] =(*this)(2,3);	CartPos[7] =(*this)(2,4);
		CartPos[8]=(*this)(3,1) ;	CartPos[9]=(*this)(3,2);	CartPos[10]=(*this)(3,3);	CartPos[11]=(*this)(3,4);
	}
	//***************************************************************************************************
	Transform3D Transform3D::roty(Real t)
	/*!
	  @brief  Rotates about Y axis. T = roty(theta)
	  @Returns a homogeneous transformation representing a rotation of theta about the Y axis.
	*/
	{
		Real ct = cos(t), st = sin(t);
		Transform3D T0;
		T0(1,1) = ct;	T0(1,2) = 0;		T0(1,3) = st;	T0(1,4) = 0;
		T0(2,1) = 0;	T0(2,2) = 1;		T0(2,3) = 0;	T0(2,4) = 0;
		T0(3,1) = -st;	T0(3,2) = 0;		T0(3,3) = ct;	T0(3,4) = 0;
		return T0;
	}
	//***************************************************************************************************
	Transform3D Transform3D::rotx(Real t)
	/*!
	  @brief  Rotates about X axis. T = rotx(theta)
	  @Returns a homogeneous transformation representing a rotation of theta about the X axis.
	*/
	{
		Real ct = cos(t), st = sin(t);
		Transform3D T0;
		T0(1,1) = 1;		T0(1,2) = 0;	T0(1,3) = 0;	T0(1,4) = 0;
		T0(2,1) = 0;		T0(2,2) = ct;	T0(2,3) = -st;	T0(2,4) = 0;
		T0(3,1) = 0;		T0(3,2) = st;	T0(3,3) = ct;	T0(3,4) = 0;
		return T0;
	}
	//***************************************************************************************************
	Transform3D Transform3D::rotz(Real t)
	/*!
	  @brief  Rotates about Z axis. T = rotz(theta)
	  @Returns a homogeneous transformation representing a rotation of theta about the Z axis.
	*/
	{
		Real ct = cos(t) , st = sin(t);
		Transform3D T0;
		T0(1,1) = ct;	T0(1,2) = -st;	T0(1,3) = 0;		T0(1,4) = 0;
		T0(2,1) = st;	T0(2,2) = ct;	T0(2,3) = 0;		T0(2,4) = 0;
		T0(3,1) = 0;	T0(3,2) = 0;	T0(3,3) = 1;		T0(3,4) = 0;
		return T0;
	}
	//***************************************************************************************************
	Transform3D Transform3D::eul2tr(Real phi, Real theta, Real psi)
	/*!
	  @brief  Converts Euler angles to homogeneous transformation
	  			TR = eul2tr(PHI, THETA, PSI)
	  @Returns a homogeneous transformation for the specified Euler angles. 
	  These correspond to rotations about the Z, Y, Z axes respectively.
	*/
	{
		return rotz(phi) * roty(theta) * rotz(psi);
	}
	//***************************************************************************************************
	Transform3D Transform3D::rpy2tr(Real roll, Real pitch, Real yaw)
	/*!
	  @brief  Converts Roll/pitch/yaw to homogeneous transform
	  			TR = rpy2tr(R,P,Y)
	  @Returns a homogeneous transformation for the specified roll/pitch/yaw angles. 
	  These correspond to rotations about the Z, Y, X axes respectively.
	*/
	{
		return rotz(roll) * roty(pitch) * rotx(yaw);
	}
	//***************************************************************************************************
	PositionVector3D Transform3D::tr2eul(Transform3D T0)
	/*!
	  @brief  Converts a homogeneous transform matrix to Euler angle form
	  			<PHI, THETA, PSI> = tr2eul(TR)
	  @Returns a vector of Euler angles corresponding to the rotational part of
	  the homogeneous transform TR.  The 3 angles correspond to rotations about
	  the Z, Y and Z axes respectively.
	*/
	{
		PositionVector3D euler;
		euler(1) = atan2(T0(2,3), T0(1,3));
		Real sp = sin(euler(1));
		Real cp = cos(euler(1));
		euler(2) = atan2(cp*T0(1,3) + sp*T0(2,3), T0(3,3));
		euler(3) = atan2(-sp * T0(1,1) + cp * T0(2,1), -sp*T0(1,2) + cp*T0(2,2));
		return euler;
	}
	//***************************************************************************************************
	PositionVector3D Transform3D::tr2rpy(Transform3D Ti)
	/*!
	  @brief  Converts a homogeneous transform matrix to roll/pitch/yaw angles
	  			RPY = tr2rpy(TR)
	  @Returns a vector of roll/pitch/yaw angles corresponding to the rotational
	  part of the homogeneous transform TR.
	*/
	{
		return RotationMatrix::r2rpy(Ti.R());
	}
	//***************************************************************************************************
	Transform3D Transform3D::normalizeT(Transform3D Ti)
	/*!
	  @brief  Normalizes a homogeneous transformation.
	  			TN = normalizeT(Real)
	  @Returns a normalized homogeneous transformation matrix in which the rotation
	  sub-matrix is a proper orthogonal matrix. The O and V vectors are 
	  normalized and the normal vector is formed from O x A.
	*/
	{
		PositionVector3D O, A, i, j, k;
		O = Ti.R().getCol(2);
		A = Ti.R().getCol(3);
		PositionVector3D n = math::cross(O, A);	// N = O x A
		i = math::UnitVector(n);
		j = math::UnitVector(O);
		k = math::UnitVector(A);
		return Transform3D(Ti.P(), RotationMatrix(i,j,k));
	}
	//***************************************************************************************************
	Transform3D Transform3D::MakeT(RotationMatrix &R, PositionVector3D X)
	/*!
		@brief Makes transformation matrix from rotation matrix (R) and position (x).
	*/
	{
		return Transform3D(X,R);
	}
	//***************************************************************************************************
	Transform3D Transform3D::MakeT(Transform3D &T0, PositionVector3D X)
	/*!
		@brief Makes transformation matrix from rotation matrix (R) and position (x).
	*/
	{
		return Transform3D(X, T0.R());
	}
	//***************************************************************************************************
	Transform3D Transform3D::rotk(const Real theta, const PositionVector3D & k)
	/*!
		@brief Rotation around arbitrary axis.
	*/
	{
		Transform3D rot;
		Real c, s, vers, kx, ky, kz;

		vers = PositionVector3D::SumSquare(k);
		if (vers != ZERO) { // compute the rotation if the vector norm is not 0.0
			vers = sqrt(1/vers);
			kx = k(1)*vers;
			ky = k(2)*vers;
			kz = k(3)*vers;
			s = sin(theta);
			c = cos(theta);
			vers = 1-c;

			rot(1,1) = kx*kx*vers+c;
			rot(1,2) = kx*ky*vers-kz*s;
			rot(1,3) = kx*kz*vers+ky*s;
			rot(2,1) = kx*ky*vers+kz*s;
			rot(2,2) = ky*ky*vers+c;
			rot(2,3) = ky*kz*vers-kx*s;
			rot(3,1) = kx*kz*vers-ky*s;
			rot(3,2) = ky*kz*vers+kx*s;
			rot(3,3) = kz*kz*vers+c;
		}
		return rot;
	}
	//***************************************************************************************************
	Transform3D Transform3D::rotd(Real theta, const PositionVector3D & k1,
		const PositionVector3D & k2)
	/*!
		@brief Rotation around an arbitrary line.
	*/
	{
		Transform3D rot;
		rot = Transform3D(k1)*rotk(theta,k2-k1)*Transform3D(-k1);
		//rot.Release();
		return rot;
	}
	//***************************************************************************************************
	Transform3D Transform3D::rpy(const PositionVector3D & a)
	/*!
		@brief Roll Pitch Yaw rotation.
	*/
	{
		Transform3D rot;
		Real ca, sa, cb, sb, cc, sc;

		ca = cos(a(1));
		sa = sin(a(1));
		cb = cos(a(2));
		sb = sin(a(2));
		cc = cos(a(3));
		sc = sin(a(3));

		rot(1,1) = cb*cc;
		rot(1,2) = sa*sb*cc-ca*sc;
		rot(1,3) = ca*sb*cc+sa*sc;
		rot(2,1) = cb*sc;
		rot(2,2) = sa*sb*sc+ca*cc;
		rot(2,3) = ca*sb*sc-sa*cc;
		rot(3,1) = -sb;
		rot(3,2) = sa*cb;
		rot(3,3) = ca*cb;

		//rot.Release();
		return rot;
	}
	//***************************************************************************************************
	Transform3D Transform3D::eulzxz(const PositionVector3D & a)
	/*!
		@brief Euler ZXZ rotation.
	*/
	{
		Transform3D rot;
		Real c1, s1, ca, sa, c2, s2;

		c1 = cos(a(1));
		s1 = sin(a(1));
		ca = cos(a(2));
		sa = sin(a(2));
		c2 = cos(a(3));
		s2 = sin(a(3));

		rot(1,1) = c1*c2-s1*ca*s2;
		rot(1,2) = -c1*s2-s1*ca*c2;
		rot(1,3) = sa*s1;
		rot(2,1) = s1*c2+c1*ca*s2;
		rot(2,2) = -s1*s2+c1*ca*c2;
		rot(2,3) = -sa*c1;
		rot(3,1) = sa*s2;
		rot(3,2) = sa*c2;
		rot(3,3) = ca;

		//rot.Release();
		return rot;
	}
	//***************************************************************************************************
	ColumnVector Transform3D::irotk(const Transform3D & R)
	/*!
		@brief Obtain axis from a rotation matrix.
	*/
	{
		ColumnVector k(4);
		Real a, b, c;

		a = (R(3,2)-R(2,3));
		b = (R(1,3)-R(3,1));
		c = (R(2,1)-R(1,2));
		k(4) = atan2(sqrt(a*a + b*b + c*c),(R(1,1) + R(2,2) + R(3,3)-1));
		k(1) = (R(3,2)-R(2,3))/(2*sin(k(4)));
		k(2) = (R(1,3)-R(3,1))/(2*sin(k(4)));
		k(3) = (R(2,1)-R(1,2))/(2*sin(k(4)));

		//k.Release();
		return k;
	}
	//***************************************************************************************************
	PositionVector3D Transform3D::irpy(const Transform3D & R)
	/*!
		@brief Obtain Roll, Pitch and Yaw from a rotation matrix.
	*/
	{
		PositionVector3D k;

		if (R(3,1)==1) {
			k(1) = atan2(-R(1,2),-R(1,3));
			k(2) = -M_PI/2;
			k(3) = 0.0;
		} else if (R(3,1)==-1) {
			k(1) = atan2(R(1,2),R(1,3));
			k(2) = M_PI/2;
			k(3) = 0.0;
		} else {
			k(1) = atan2(R(3,2), R(3,3));
			k(2) = atan2(-R(3,1), sqrt(R(1,1)*R(1,1) + R(2,1)*R(2,1)));
			k(3) = atan2(R(2,1), R(1,1));
		}

		//k.Release();
		return k;
	}
	//***************************************************************************************************
	PositionVector3D Transform3D::ieulzxz(const Transform3D & R)
	/*!
		@brief Obtain Roll, Pitch and Yaw from a rotation matrix.
	*/
	{
		PositionVector3D  a;

		if ((R(3,3)==1)  || (R(3,3)==-1)) {
			a(1) = 0.0;
			a(2) = ((R(3,3) == 1) ? 0.0 : M_PI);
			a(3) = atan2(R(2,1),R(1,1));
		} else {
			a(1) = atan2(R(1,3), -R(2,3));
			a(2) = atan2(sqrt(R(1,3)*R(1,3) + R(2,3)*R(2,3)), R(3,3));
			a(3) = atan2(R(3,1), R(3,2));
		}

		return a;
	}
	//***************************************************************************************************
	Transform3D Transform3D::CV2T(const math::ColumnVector & v)
	/*!
		@brief Obtain Fill Transformation matrix from a Col. Vector.
	*/
	{
		Transform3D Ti;
		Ti(1,1) = v(1);		Ti(1,2) = v(2) ;		Ti(1,3) = v(3) ;		Ti(1,4) = v(4) ;
		Ti(2,1) = v(5);		Ti(2,2) = v(6) ;		Ti(2,3) = v(7) ;		Ti(2,4) = v(8) ;
		Ti(3,1) = v(9);		Ti(3,2) = v(10);		Ti(3,3) = v(11);		Ti(3,4) = v(12);
		return Ti;
	}
	//***************************************************************************************************
	ColumnVector Transform3D::T2CV(const Transform3D &Ti)
	{
	/*!
		@brief Return a Col Vector from a Transformation matrix.
	*/
		math::ColumnVector v(12);
		v(1) = Ti(1,1);		v(2) = Ti(1,2);		v(3) = Ti(1,3);		v(4) = Ti(1,4);
		v(5) = Ti(2,1);		v(6) = Ti(2,2);		v(7) = Ti(2,3);		v(8) = Ti(2,4);
		v(9) = Ti(3,1);		v(10)= Ti(3,2);		v(11)= Ti(3,3);		v(12)= Ti(3,4);
		return v;
	}
	//***************************************************************************************************
	Matrix Transform3D::subMatrix() const{
		Matrix mat(3, 4);
		for (int i=1; i<=3; i++)
			for (int j=1; j<=4; j++)
				if(j <= 3)
					mat(i,j) = _R(i,j);
				else
					mat(i,j) = _P(i);
		return mat;
	}
	//***************************************************************************************************
	void Transform3D::Integrate(const Twist& t_this,double samplefrequency)
	//! The twist <t_this> is expressed wrt the current
	//! frame.  This frame is integrated into an updated frame with
	//! <samplefrequency>.  Very simple first order integration rule.
	{
		double n = t_this.Rot().norm()/samplefrequency;
		if (n<EPSilon) {
			_P += _R*(t_this.Vel()/samplefrequency);
		} else {
			(*this) = (*this) * 
				Transform3D ( t_this.Vel()/samplefrequency, RotationMatrix::Rot( t_this.Rot(), n ) );
		}
	}
	//***************************************************************************************************
	/*std::istream& operator >> (std::istream& is, Transform3D& T)
	{
		IOTrace("Stream input Frame (Rotation,Vector) or DH[...]");
		char storage[10];
		EatWord(is, "[", storage, 10);
		if (strlen(storage) == 0) {
			Eat(is, '[');
			is >> T.R();
			is >> T.P();
			EatEnd(is, ']');
			IOTracePop();
			return is;
		}
		if (strcmp(storage, "DH") == 0) {
			double a, alpha, d, theta;
			const double deg2rad = 0.01745329251994329576923690768488;
			const double rad2deg = 57.2957795130823208767981548141052;
			Eat(is, '[');
			is >> a;
			Eat(is, ',');
			is >> alpha;
			Eat(is, ',');
			is >> d;
			Eat(is, ',');
			is >> theta;
			EatEnd(is, ']');
			T = Transform3D::DH(a, alpha*deg2rad, d, theta*deg2rad);
			IOTracePop();
			return is;
		}
		throw Error_Frame_Frame_Unexpected_id();
		return is;
	}*/
	//***************************************************************************************************
	/*
	ColumnVector Transform3D::getMatrixAsColumn() const{
		int n = (*this)._nRows*(*this)._nColumns;
		ColumnVector col(n);
		for(int i=0; i<n; i++)
			col(i+1) = (*this)._data[i];
		return col;
	}
*/
	//***************************************************************************************************
	Transform3D trinterp(const Transform3D & T0, const Transform3D & T1, Real r)
	/*!
		@brief  Interpolates homogeneous transformations
					TR = trinterp(T0, T1, R)
				Returns a homogeneous transform interpolation between T0 and T1 as
				R varies from 0 to 1.  Rotation is interpolated using quaternion
				spherical linear interpolation.
	*/
	{
		PositionVector3D p0, p1, pr;
		Quaternion q0, q1, qr;
		RotationMatrix Rr;

		q0 = T2Q(T0);
		q1 = T2Q(T1);

		p0 = T0.P();
		p1 = T1.P();

		qr = Quaternion::qinterp(q0, q1, r);
		pr = p0*(1-r) + r*p1;

		Rr = Q2R(qr);
		return Transform3D(pr, Rr);
	}
	//***************************************************************************************************
	Transform3D * ctraj(const Transform3D & T0, const Transform3D & T1, Real n)
	/*!
	 * @brief  Computes a Cartesian trajectory from pose T0 to T1.
	 * 			TC = ctraj(T0, T1, N)
	 *		   with N points that follow a trapezoidal velocity profile along the path.
	 * Copyright (C) 1993-2015, by Peter I. Corke
	*/
	{
		if ( (int)n <= 0)
			cout<<"\nInvalid n for ctraj\n";
		Transform3D *Ti = new Transform3D[(int)n+1];
		math::ColumnVector i((int)n), r((int)n);
		for (int k=1; k<=(int)n; k++){
			Ti[k-1] = Transform3D();
			i(k) = static_cast<Real>(k);
			if (k>1)
				r(k) = Real(i(k-1)/(n-1));
		}
		for (int k=0; k<(int)n; k++){
			Ti[k] = trinterp(T0, T1, r(k+1));
		}
		Ti[(int)n] = T1;
		return Ti;
	}
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	void jtraj(const math::ColumnVector &q0, const ColumnVector &q1, const ColumnVector &tv,
		Matrix &qt, Matrix &qdt,Matrix &qddt){
			int n = tv.getnRows();
			int m = q0.getnRows();
			Real tscal = tv.maxV();
			ColumnVector t = tv/tscal, qd0(m), qd1(m);
			//compute the polynomial coefficients
			ColumnVector A(m), B(m), C(m), H(m), F(m);
			A =  6 * (q1 - q0) - 3 * (qd1 + qd0) * tscal;
			B = -15* (q1 - q0) + (8 * qd0 + 7 * qd1) * tscal;
			C =  10* (q1 - q0) - (6 * qd0 + 4 * qd1) * tscal;
			H = qd0 * tscal; // as the t vector has been normalized
			F = q0;
			Matrix Real(n,m-1), c(m-1,m);
			ColumnVector t5(n), t4(n), t3(n), t2(n), ones(n);
			t5 = t.pow(5);
			t4 = t.pow(4);
			t3 = t.pow(3);
			t2 = t.pow(2);
			ones.ones();
			for (int i=1; i<=n;i++){
				Real(i,1) = t5(i);
				Real(i,2) = t4(i);
				Real(i,3) = t3(i);
				Real(i,4) = t2(i);
				Real(i,5) = tv(i);
				Real(i,6) = ones(i);
			}
			for (int i=1; i<=m; i++){

				c(1,i) = A(i);
				c(2,i) = B(i);
				c(3,i) = C(i);
				c(4,i) = 0.0;
				c(5,i) = H(i);
				c(6,i) = F(i);
			}
			// compute positions
			Matrix _qt = Real * c;
			for (int i=1; i<=m; i++){
				c(1,i) = 0.0;
				c(2,i) = 5 * A(i);
				c(3,i) = 4 * B(i);
				c(4,i) = 3 * C(i);
				c(5,i) = 0.0;
				c(6,i) = H(i);
			}
			// compute velocity
			Matrix _qdt = Real * c/tscal;
			for (int i=1; i<=m; i++){
				c(1,i) = 0.0;
				c(2,i) = 0.0;
				c(3,i) = 20 * A(i);
				c(4,i) = 12 * B(i);
				c(5,i) = 6  * C(i);
				c(6,i) = 0.0;
			}
			// compute acceleration
			Matrix _qddt = Real * c/(tscal*tscal);

			qt = Matrix(_qt);
			qdt = Matrix(_qdt);
			qddt = Matrix(_qddt);
			// 	qt.reConstruct(_qt);
			// 	qdt.reConstruct(_qdt);
			// 	qddt.reConstruct(_qddt);
	}
	void random(Transform3D& F) {
		random(F.R());
		random(F.P());
	}
	void distance(const Transform3D &T0, const Transform3D &T1, Real &dist, Real &angle){
		dist = distance(T0.P(), T1.P());
		angle = distance(T2Q(T0),T2Q(T1));
	}
};

