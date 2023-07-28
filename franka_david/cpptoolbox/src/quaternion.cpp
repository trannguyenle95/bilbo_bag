//=============================================================================
//==============================================================================

//	file:	quaternion.cpp

//	author:	Fares J. Abu-Dakka
//	date:	April. 2012

//	Description: Quaternion representation of the rotation in 3D
//==============================================================================

#include <math.h>
#include "quaternion.hpp"

namespace math{
	
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion() : _s(1), _v(math::RowVector(3))
	/*!
	  @brief Constructs default quaternion 1<0,0,0>.
	*/
	{ }
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const math::ColumnVector & v)
	/*!
	  @brief Constructs quaternion from 4D ColumnVector
	  @param 4D vector $ \mathbf{v} $.
	*/
	{
		if (v.getnRows() != 4 || v.getnColumns() != 1)
			throw Exception("Invalid Vector for Quaternion");
		_v = math::RowVector(3);
		s(v(1));	x(v(2));	y(v(3));	z(v(4));
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const math::PositionVector3D & v)
	/*!
	  @brief Constructs quaternion from 3D PositionVector. Sets scalar value to 1
	  @param 3D vector $ \mathbf{v} $.
	*/
	{
		_v = math::RowVector(3);
		s(ONE);	x(v(1));	y(v(2));	z(v(3));
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const math::RowVector & r)
	/*!
	  @brief Constructs quaternion from 4D RowVector
	  @param 4D vector $ \mathbf{v} $.
	*/
	{
		if (r.getnRows() != 1 || r.getnColumns() != 4)
			throw Exception("Invalid Vector for Quaternion");
		_v = math::RowVector(3);
		s(r(1));	x(r(2));	y(r(3));	z(r(4));
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const Quaternion & q) : _s(q.s()), _v(q.v())
	/*!
	  "brief Copy constructor
	*/
	{ }
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const math::RotationMatrix & R)
	/*
	  @brief Extracts a Quaternion from Rotation matrix using R2Q(const RotationMatrix& rot)
	  @param A 3x3 rotation matrix $ \mathbf{R} $
	*/
	{
		*this = R2Q(R);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const math::Transform3D & Ti)
	/*
	  @brief Extracts a Quaternion from Rotation matrix using T2Q(const Transform3D& Ti)
	  @param A 4x4 homogeneous transformation matrix $ \mathbf{Ti} $
	*/
	{
		*this = T2Q(Ti);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const Real s, const math::RowVector r)
	/*!
	  @brief Constructs quaternion from 3D RowVector
	  @param A scalar s;
	  @param 3D vector $ \mathbf{v} $.
	*/
	{
		if (r.getnRows() != 1 || r.getnColumns() != 3)
			throw Exception("Invalid Vector for Quaternion");
		_v = r;
		_s = s;
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const math::RowVector axis, const Real angle)
	/*!
	  @brief Constructs quaternion from a axis and angle.
	*/
	{
		if(axis.getnColumns() != 3){
			cerr << "Quaternion::Quaternion, size of axis != 3" << endl
				<< "Press enter to exit...." << endl;
			getchar();
			exit(1);
		}
		// make sure axis is a unit vector
		Real norm_axis = ::sqrt(math::dot(axis,axis));
		if(norm_axis != 1){
			cerr << "warning: Quaternion::Quaternion(angle, axis), axis is not unit" << endl;
			cerr << "Normalizing!!! Making the axis unit." << endl;
			_v = ::sin(angle/2) * axis/norm_axis;
		}
		else
			_v = ::sin(angle/2) * axis;
		_s = ::cos(angle/2);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const math::ColumnVector axis, const Real angle)
	/*!
	  @brief Constructs quaternion from a axis and angle.
	*/
	{
		if(axis.getnRows() != 3){
			cerr << "Quaternion::Quaternion, size of axis != 3" << endl
				<< "Press enter to exit...." << endl;
			getchar();
			exit(1);
		}
		_v = math::RowVector(3);
		// make sure axis is a unit vector
		Real norm_axis = ::sqrt(math::dot(axis,axis));
		if(norm_axis != 1){
			cerr << "warning: Quaternion::Quaternion(angle, axis), axis is not unit" << endl;
			cerr << "Normalizing!!! Making the axis unit." << endl;
			ColumnVector v = ::sin(angle/2) * axis/norm_axis;
			x(v(1));	y(v(2));	z(v(3));
		}else{
			ColumnVector v = ::sin(angle/2) * axis;
			x(v(1));	y(v(2));	z(v(3));
		}
		_s = ::cos(angle/2);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const math::PositionVector3D axis, const Real angle)
	/*!
	  @brief Constructs quaternion from a axis and angle.
	*/
	{
		_v = math::RowVector(3);
		// make sure axis is a unit vector
		Real norm_axis = ::sqrt(math::dot(axis,axis));
		if(norm_axis != 1){
			cerr << "warning: Quaternion::Quaternion(angle, axis), axis is not unit" << endl;
			cerr << "Normalizing!!! Making the axis unit." << endl;
			PositionVector3D v = ::sin(angle/2) * axis/norm_axis;
			x(v(1));	y(v(2));	z(v(3));
		}else{
			PositionVector3D v = ::sin(angle/2) * axis;
			x(v(1));	y(v(2));	z(v(3));
		}
		_s = ::cos(angle/2);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion(const Real s, const Real x, const Real y, const Real z) 
		: _s(s), _v(math::RowVector(x,y,z))
	/*!
	  @brief Constructs a Quaternion from 4 parameteres.
	  @param A scalar $ q_s $
	  @param x-value of quaternion vector $ q_x $
	  @param y-value of quaternion vector $ q_y $
	  @param z-value of quaternion vector $ q_z $
	*/
	{ }
	//-------------------------------------------------------------------------------------------
	Quaternion::Quaternion (const Real rzPsi, const Real ryTheta, const Real rxPhi)
	/*!
	  @brief Constructs quaternion by conversion rotation angles.
	*/
	{
		_v = math::RowVector(3);
		_s    = cos(rxPhi/TWO)*cos(ryTheta/TWO)*cos(rzPsi/TWO) + sin(rxPhi/TWO)*sin(ryTheta/TWO)*sin(rzPsi/TWO);
		_v(1) = sin(rxPhi/TWO)*cos(ryTheta/TWO)*cos(rzPsi/TWO) - cos(rxPhi/TWO)*sin(ryTheta/TWO)*sin(rzPsi/TWO);
		_v(2) = cos(rxPhi/TWO)*sin(ryTheta/TWO)*cos(rzPsi/TWO) + sin(rxPhi/TWO)*cos(ryTheta/TWO)*sin(rzPsi/TWO);
		_v(3) = cos(rxPhi/TWO)*cos(ryTheta/TWO)*sin(rzPsi/TWO) - sin(rxPhi/TWO)*sin(ryTheta/TWO)*cos(rzPsi/TWO);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion & Quaternion::operator =(const Quaternion & Q)
	/*!
	  @brief Assignment operator.
	*/
	{
		_s = Q.s(); 
		_v = Q.v(); 
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	Quaternion & Quaternion::operator =(const math::RotationMatrix & R)
	/*!
	  @brief Assignment operator.
	*/
	{	
		Quaternion q = R2Q(R);
		_s = q.s();
		_v = q.v();
		return *this;	
	}
	//-------------------------------------------------------------------------------------------
	const Quaternion Quaternion::operator  /(Real n)const
	/*!
	  @brief Scalar division.
	*/
	{
		if (n == 0)
			throw Exception("Quaternion, scalar division by zero.");
		return Quaternion(_s/n, _v/n);
	}
	//-------------------------------------------------------------------------------------------
	const Quaternion Quaternion::operator +(const Quaternion &Q)const
	/*!
	  @brief Quaternion addition. $q_1 + q_2 = <s_1, v_1> + <s_2, v_2> = <s_1+s_2, v_1+v_2>$
	  The result is not necessarily a unit quaternion even if q1 and
	  q2 are unit quaternions. So the result is normalized before it returns
	*/
	{
		Quaternion q(_s+Q.s(),_v+Q.v());
		q.Normalize();
		return q;
	}
	//-------------------------------------------------------------------------------------------
	const Quaternion Quaternion::operator -(const Quaternion &Q)const
	/*!
	  @brief Quaternion addition. $q_1 - q_2 = <s_1, v_1> - <s_2, v_2> = <s_1-s_2, v_1-v_2>$
	  The result is not necessarily a unit quaternion even if q1 and
	  q2 are unit quaternions. So the result is normalized before it returns
	*/
	{
		Quaternion q(_s-Q.s(),_v-Q.v());
		q.Normalize();
		return q;
	}
	//-------------------------------------------------------------------------------------------
	const Quaternion Quaternion::operator *(const Quaternion &r)const
	/*!
	  @brief Quaternion multiplication. 
	  q = q1*q2 = <s1*s2 - dotProduct(v1,v_2), crossProduct(v1,v2) + s1*v2 + s2*v1>
	  If q1 and q2 are unit quaternions, then q will also be a unit quaternion.
	*/
	{
		Real r0 = r.s(), r1 = r.x(), r2 = r.y(), r3 = r.z();
		Real q0 =   s(), q1 =   x(), q2 =   y(), q3 =   z();
		Real t0 = r0 * q0 - r1 * q1 - r2 * q2 - r3 * q3;
		Real t1 = r0 * q1 + r1 * q0 - r2 * q3 + r3 * q2;
		Real t2 = r0 * q2 + r1 * q3 + r2 * q0 - r3 * q1;
		Real t3 = r0 * q3 - r1 * q2 + r2 * q1 + r3 * q0;
		return Quaternion(t0,t1,t2,t3);
	}
	//-------------------------------------------------------------------------------------------
	const Quaternion Quaternion::operator /(const Quaternion &r)const
	/*!
	  @brief Quaternion division.
	*/
	{
		Real r0 = r.s(), r1 = r.x(), r2 = r.y(), r3 = r.z();
		Real q0 = s(), q1 = x(), q2 = y(), q3 = z();
		Real t0 = (r0 * q0 + r1 * q1 + r2 * q2 + r3 * q3)/(r0*r0 + r1*r1 + r2*r2 + r3*r3);
		Real t1 = (r0 * q1 - r1 * q0 - r2 * q3 + r3 * q2)/(r0*r0 + r1*r1 + r2*r2 + r3*r3);
		Real t2 = (r0 * q2 + r1 * q3 - r2 * q0 - r3 * q1)/(r0*r0 + r1*r1 + r2*r2 + r3*r3);
		Real t3 = (r0 * q3 - r1 * q2 + r2 * q1 - r3 * q0)/(r0*r0 + r1*r1 + r2*r2 + r3*r3);
		return Quaternion(t0,t1,t2,t3);
	}
	//-------------------------------------------------------------------------------------------
	void Quaternion::setQuaternion(const Real scalar, const math::RowVector &vector)
	/*!
	  @brief Sets quaternion from scalar and 3D RowVector
	  @param A scalar s;
	  @param 3D vector $ \mathbf{v} $.
	*/
	{
		if(vector.getnColumns() != 3){
			cerr << "Quaternion::setQuaternion, size of vector != 3" << endl
				<< "Press enter to exit...." << endl;
			getchar();
			exit(1);
		}
		_s = scalar;
		_v = vector;
	}
	//-------------------------------------------------------------------------------------------
	void Quaternion::setVector(const math::RowVector & vector)
	/*!
	  @brief Sets vector part of a quaternion from 3D RowVector
	  @param 3D vector $ \mathbf{v} $.
	*/
	{
		if(vector.getnColumns() != 3){
			cerr << "Quaternion::setVector, size of vector != 3" << endl
				<< "Press enter to exit...." << endl;
			getchar();
			exit(1);
		}
		_v = vector;
	}
	//-------------------------------------------------------------------------------------------
	ostream & operator << (ostream & o, const Quaternion & Q)
	/*
      @brief Writes a quaternion to stream
      
      @param os output stream to use
      @param Q quaternion to print
      @return the updated output stream
	*/
	{
		o << std::right << std::fixed<< std::setprecision(4) << std::setw(8)
		  << Q.s() << " < " << Q.v()(1) << " "
		  << std::right << std::fixed<< std::setprecision(4) << std::setw(8)
		  << Q.v()(2) << " "
		  << std::right << std::fixed<< std::setprecision(4) << std::setw(8)
		  << Q.v()(3) << " >\n\n";
		return o;
	}
	//-------------------------------------------------------------------------------------------
	void Quaternion::print (char const* ch, int prec, int w) const
	/*
      @brief Prints a quaternion to screen
      
      @param ch variable name
	  @param prec output precision
	  @param w width of the element
	*/
	{
		if (ch)
			std::cout << ch << ": ";

		cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8)
			 << _s << " < " << _v(1) << " "
			 << std::right << std::fixed<< std::setprecision(4) << std::setw(8)
			 << _v(2) << " "
			 << std::right << std::fixed<< std::setprecision(4) << std::setw(8)
			 << _v(3) << " > norm = " << norm() << "\n";
	}
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	bool Quaternion::is_inf() const {
	  	return
			( _s ==  std::numeric_limits<Real>::infinity() ) ||
			( _s == -std::numeric_limits<Real>::infinity() ) ||
			( _v ==  std::numeric_limits<Real>::infinity() ) ||
			( _v == -std::numeric_limits<Real>::infinity() );
	}
	//-------------------------------------------------------------------------------------------
	bool Quaternion::is_nan() const {
	  	return ( _s != _s) || ( _v != _v );
	}
	//-------------------------------------------------------------------------------------------
	bool Quaternion::is_neg_inf() const {
	  	return ( _s == -std::numeric_limits<Real>::infinity() ) && ( _v == ZERO );
	}
	//-------------------------------------------------------------------------------------------
	bool Quaternion::is_pos_inf() const {
	  	return ( _s ==  std::numeric_limits<Real>::infinity() ) && ( _v == ZERO );
	}
	//-------------------------------------------------------------------------------------------
	bool Quaternion::is_zero() const {
		return (_s == ZERO) && (_v == ZERO);
	}
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	Real Quaternion::normSquared()const
	/*!
	  @brief Return the squared of quaternion norm.
	
	  $
		normSquared(\mathbf{q}) = s^2 + \mathbf{v} \cdot \mathbf{v}
	  $
	*/
	{
		return _s*_s + _v(1)*_v(1) + _v(2)*_v(2) + _v(3)*_v(3);
	}
	//-------------------------------------------------------------------------------------------
	Real Quaternion::norm()const
	/*!
	  @brief Return the quaternion norm. A unit quaternion has norm of one.
	
	  $
		norm(\mathbf{q}) = \sqrt{s^2 + \mathbf{v} \cdot \mathbf{v}}
	  $
	*/
	{
		return sqrt(normSquared());
	}
	//-------------------------------------------------------------------------------------------
	Quaternion & Quaternion::Normalize()
	/*!
	  @brief Normalize a quaternion
	*/
	{
		Real n = norm();
		if (n > EPSilon)
			(*this) /= norm();
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::conjugate()const
	/*!
	  @brief Conjugate of a quaternion.
	  if $\mathbf{q} = <s, \mathbf{v}>$ then $\mathbf{q}^{*} = <s, -\mathbf{v}>$
	*/
	{
		return Quaternion(_s,_v*-1);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::qinv()const
	/*!
	  @brief Quaternion inverse.
	  $
	    \mathbf{q}^{-1} = \frac{\mathbf{q}^{*}}{norm(\mathbf{q})}
	  $

	  where $\mathbf{q}^{*}$ and $norm(\mathbf{q})$ are 
	  the quaternion conjugate and the quaternion norm respectively.

	  @note for unit quaternions, the inverse is equal to the conjugate.
	*/
	{
		return conjugate()/normSquared();
	}
	//-------------------------------------------------------------------------------------------
	void Quaternion::Quat2Angle(const Quaternion &q, Real &rzPsi, Real &ryTheta, Real &rxPhi)
	/*!
	  @brief Converts quaternion to rotation angles.
	*/
	{
		Real q0 = q.s(), q1 = q.x(), q2 = q.y(), q3 = q.z();
		rxPhi   = atan2( 2*(q0*q1 + q2*q3), 1 - 2*(q1*q1 + q2*q2) );
		ryTheta =  asin( 2*(q0*q2 - q3*q1) );
		rzPsi   = atan2( 2*(q0*q3 + q1*q2), 1 - 2*(q2*q2 + q3*q3) );
	}
	//-------------------------------------------------------------------------------------------
	void Quaternion::Angle2Quat(Real rzPsi, Real ryTheta, Real rxPhi, char * type)
	/*!
	  @brief Converts rotation angles to Quaternion.
	*/
	{
		math::PositionVector3D cang, sang;
		cang(1) = cos(rzPsi/TWO);		cang(2) = cos(ryTheta/TWO);		cang(3) = cos(rxPhi/TWO);
		sang(1) = sin(rzPsi/TWO);		sang(2) = sin(ryTheta/TWO);		sang(3) = sin(rxPhi/TWO);
		if (!type || strcmp(type, "zyx"))
			(*this) = Quaternion(rzPsi, ryTheta, rxPhi);
		else if (strcmp(type, "zyz")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*cang(2)*sang(3),
									cang(1)*sang(2)*sang(3) - sang(1)*sang(2)*cang(3),
									cang(1)*sang(2)*cang(3) + sang(1)*sang(2)*sang(3),
									sang(1)*cang(2)*cang(3) + cang(1)*cang(2)*sang(3)	);
		}
		else if (strcmp(type, "zxy")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*sang(2)*sang(3),
									cang(1)*sang(2)*cang(3) - sang(1)*cang(2)*sang(3),
									cang(1)*cang(2)*sang(3) + sang(1)*sang(2)*cang(3),
									cang(1)*sang(2)*sang(3) + sang(1)*cang(2)*cang(3)	);
		}
		else if (strcmp(type, "zxz")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*cang(2)*sang(3),
									cang(1)*sang(2)*cang(3) + sang(1)*sang(2)*sang(3),
									sang(1)*sang(2)*cang(3) - cang(1)*sang(2)*sang(3),
									cang(1)*cang(2)*sang(3) + sang(1)*cang(2)*cang(3)	);
		}
		else if (strcmp(type, "yxz")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) + sang(1)*sang(2)*sang(3),
									cang(1)*sang(2)*cang(3) + sang(1)*cang(2)*sang(3),
									sang(1)*cang(2)*cang(3) - cang(1)*sang(2)*sang(3),
									cang(1)*cang(2)*sang(3) - sang(1)*sang(2)*cang(3)	);
		}
		else if (strcmp(type, "yxy")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*cang(2)*sang(3),
									cang(1)*sang(2)*cang(3) + sang(1)*sang(2)*sang(3),
									sang(1)*cang(2)*cang(3) + cang(1)*cang(2)*sang(3),
									cang(1)*sang(2)*sang(3) - sang(1)*sang(2)*cang(3)	);
		}
		else if (strcmp(type, "yzx")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*sang(2)*sang(3),
									cang(1)*cang(2)*sang(3) + sang(1)*sang(2)*cang(3),
									cang(1)*sang(2)*sang(3) + sang(1)*cang(2)*cang(3),
									cang(1)*sang(2)*cang(3) - sang(1)*cang(2)*sang(3)	);
		}
		else if (strcmp(type, "yzy")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*cang(2)*sang(3),
									sang(1)*sang(2)*cang(3) - cang(1)*sang(2)*sang(3),
									cang(1)*cang(2)*sang(3) + sang(1)*cang(2)*cang(3),
									cang(1)*sang(2)*cang(3) + sang(1)*sang(2)*sang(3)	);
		}
		else if (strcmp(type, "xyz")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*sang(2)*sang(3),
									cang(1)*sang(2)*sang(3) + sang(1)*cang(2)*cang(3),
									cang(1)*sang(2)*cang(3) - sang(1)*cang(2)*sang(3),
									cang(1)*cang(2)*sang(3) + sang(1)*sang(2)*cang(3)	);
		}
		else if (strcmp(type, "xyx")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*cang(2)*sang(3),
									cang(1)*cang(2)*sang(3) + sang(1)*cang(2)*cang(3),
									cang(1)*sang(2)*cang(3) + sang(1)*sang(2)*sang(3),
									sang(1)*sang(2)*cang(3) - cang(1)*sang(2)*sang(3)	);
		}
		else if (strcmp(type, "xzy")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) + sang(1)*sang(2)*sang(3),
									sang(1)*cang(2)*cang(3) - cang(1)*sang(2)*sang(3),
									cang(1)*cang(2)*sang(3) - sang(1)*sang(2)*cang(3),
									cang(1)*sang(2)*cang(3) + sang(1)*cang(2)*sang(3)	);
		}
		else if (strcmp(type, "xzx")){
			(*this) = Quaternion(	cang(1)*cang(2)*cang(3) - sang(1)*cang(2)*sang(3),
									cang(1)*cang(2)*sang(3) + sang(1)*cang(2)*cang(3),
									cang(1)*sang(2)*sang(3) - sang(1)*sang(2)*cang(3),
									cang(1)*sang(2)*cang(3) + sang(1)*sang(2)*sang(3)	);
		}
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::qinterp(Quaternion Q1, Quaternion Q2, Real r)
	/*!
	  @brief Calculates quaternion interpolation using SLERP (Spherical Linear Interpolation) 
	  between \b Q1 and \b Q2.
	  
	  The SLERP interpolation ensures a constant velocity across the interpolation.
	  For $r=0$ the result is \b Q1 and for $r=1$ it is \b Q2.
	  
	  The quaternion $Q(r)$ interpolate the quaternions $Q_1$
	  and $Q_2$ given the parameter $r$ along the quaternion sphere.

	  $ Q(t) = c_1(r)Q_1 + c_2(t)Q_2 $

	  where $c_1$ and $c_2$ are float functions with $0\leq r \leq 1$.
	  As $r$ varies between 0 and 1. the values $Q(r)$ varies uniformly
	  along the circular arc from $Q_1$ and $Q_2$. The angle between
	  $Q(r)$ and $Q_1$ is $\cos(r\theta)$ and the angle between
	  $Q(r)$ and $Q_2$ is $\cos((1-r)\theta)$. Taking the dot product
	  of $Q(r)$ and $Q_1$ yields

	  $	  \cos(r\theta) = c_1(r) + \cos(\theta)c_2(r)  $

	  and taking the dot product of $Q(r)$ and $Q_2$ yields
	  
	  $   \cos((1-r)\theta) = \cos(\theta)c_1(r) + c_2(r)  $

	  These are two equations with $c_1$ and $c_2$. The solution is
	  
	  $
		c_1 = \frac{\sin((1-r)\theta)}{\sin(\theta)}
		c_2 = \frac{\sin(r\theta)}{sin(\theta)}
	  $

	  The interpolation is then
	  
	  $	  Slerp(Q_1, Q_2, r) = \frac{Q_1\sin((1-r)\theta)+Q_2\sin(r\theta)}{\sin(\theta)}  $

	  If $Q_1$ and $Q_2$ are unit quaternions the $Q(r)$ is also a unit
	  quaternions. For unit quaternions we have
	  
	  $	  Slerp(Q_1, Q_2, r) = Q_1(Q_1^{-1}Q_2)^r  $

	  For r = 0 and r = 1 we have
	  
	  $
		Q_1 = Slerp(Q_1, Q_2, 0)
		Q_2 = Slerp(Q_1, Q_2, 1)
	  $

	  It is customary to choose the sign G on Q2 so that Q1.GQ2 >= 0 (the angle
	  between q0 ang Gq1 is acute). This choice avoids extra spinning caused
	  by the interpolated rotations.

	  @note Algorithm and implementation is thanks to euclideanspace.com
	*/
	{
		if ((r < 0) || (r > 1))
			cerr << "qinterp(Q1, Q2, r): r < 0 or r > 1. r is set to 0." << endl;
		Real theta = Q1.dot(Q2);
		theta = acos( (theta > 1 ? 1 : (theta < -1 ? -1 : theta)) );
		Quaternion q;
		if (abs(theta) < EPSilon)
			q = Q1;
		else
			q = (sin((1-r)*theta) * Q1 + sin(r*theta) * Q2) / sin(theta) ;
		return q;
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::Slerp_prime(const Quaternion & q0, const Quaternion & q1, const Real t)
	/*!
	  @brief Spherical Linear Interpolation derivative.
	  
	  The derivative of the function \f$q^t\f$ where \f$q\f$ is a constant
	  unit quaternion is

	  $	  \frac{d}{dt}q^t = q^t log(q)	  $

	  Using the preceding equation the Slerp derivative is then
	  
	  $	  Slerp'(q_0, q_1, t) = q_0(q_0^{-1}q_1)^t log(q_0^{-1}q_1)	  $
	  
	  It is customary to choose the sign G on q1 so that q0.Gq1 >=0 (the angle
	  between q0 ang Gq1 is acute). This choice avoids extra spinning caused
	  by the interpolated rotations.
	  The result is not necessary a unit quaternion.
	*/

	{
		if ((t < 0) || (t > 1))
			cerr << "Slerp_prime(q0, q1, t): t < 0 or t > 1. t is set to 0." << endl;

		if (q0.dot(q1) >= 0)
			return qinterp(q0, q1, t)*(q0.qinv()*q1).ln();
		else
			return qinterp(q0, q1, t)*(q0.qinv()*-1 * q1).ln();
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::Squad(const Quaternion & p, const Quaternion & a, const Quaternion & b,
		const Quaternion & q, const Real t)
	/*!
	  @brief Spherical Cubic Interpolation.
	  
	  Let four quaternions be $q_i$ (p), $s_i$ (a), $s_{i+1}$ (b) and $q_{i+1}$
	  (q) be the ordered vertices of a quadrilateral. Obtain c from $q_i$ to $q_{i+1}$
	  interpolation. Obtain d from $s_i$ to $s_{i+1}$ interpolation. Obtain e,
	  the final result, from c to d interpolation.
	  
	  $
	  	Squad(q_i, s_i, s_{i+1}, q_{i+1}, t) = Slerp(Slerp(q_i,q_{i+1},t),Slerp(s_i,s_{i+1},t), 2t(1-t))
	  $
	  
	  The intermediate quaternion $s_i$ and $s_{i+1}$ are given by
	  
	  $	s_i = q_i exp\Big ( - \frac{log(q_i^{-1}q_{i+1}) + log(q_i^{-1}q_{i-1})}{4}\Big )	$
	*/
	{
		if ((t < 0) || (t > 1))
			cerr << "Squad(p,a,b,q, t): t < 0 or t > 1. t is set to 0." << endl;

		return qinterp(qinterp(p, q, t), qinterp(a, b, t), 2 * t*(1 - t));
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::Squad_prime(const Quaternion & p, const Quaternion & a, const Quaternion & b,
		const Quaternion & q, const Real t)
	/*!
	  @brief Spherical Cubic Interpolation derivative.
	  
	  The derivative of the function $q^t$ where $q$ is a constant unit quaternion is
	  
	  $	\frac{d}{dt}q^t = q^t log(q) $
	  
	  Recalling that $log(q) = [0, v\theta]$ (see Quaternion::ln()). If the power
	  is a function we have
	  
	  $	\frac{d}{dt}q^{f(t)} = f'(t)q^{f(t)}log(q)	$
	  
	  If $q$ is a function of time and the power is differentiable function of time we have
	  
	  $	\frac{d}{dt}(q(t))^{f(t)} = f'(t)(q(t))^{f(t)}log(q) + f(t)(q(t))^{f(t)-1}q'(t)	$
	  
	  Using these last three equations Squad derivative can be define. Let
	  $U(t)=Slerp(p,q,t)$, $V(t)=Slerp(q,b,t)$, $W(t)=U(t)^{-1}V(t)$. We then
	  have $Squad(p,a,b,q,t)=Slerp(U(t),V(t),2t(1-t))=U(t)W(t)^{2t(1-t)}$
	  
	  $
		Squad'(p,a,b,q,t) = \frac{d}{dt}\Big [ UW^{2t(1-t)}\Big ]
		Squad'(p,a,b,q,t) = U\frac{d}{dt}\Big [ W^{2t(1-t)}\Big ] + U'\Big [W^{2t(1-t)}\Big]
	    Squad'(p,a,b,q,t) = U\Big[(2-4t)W^{2t(1-t)}log(W)+2t(1-t)W^{2t(1-t)-1}W'\Big]
	    + U'\Big[W^{2t(1-t)} \Big]
	  $

	  where $U'=Ulog(p^{-1}q)$, $V'=Vlog(a^{-1},b)$, $W'=U^{-1}V'-U^{-2}U'V$
	  
	  The result is not necessarily a unit quaternion even if all the input quaternions are unit.
	*/
	{
		if ((t < 0) || (t > 1))
			cerr << "Squad_prime(p,a,b,q, t): t < 0 or t > 1. t is set to 0." << endl;

		Quaternion q_squad,
			U = qinterp(p, q, t),
			V = qinterp(a, b, t),
			W = U.qinv()*V,
			U_prime = U*(p.qinv()*q).ln(),
			V_prime = V*(a.qinv()*b).ln(),
			W_prime = U.qinv()*V_prime - U.pow(-2)*U_prime*V;

		q_squad = U*(W.pow(2 * t*(1 - t))*W.ln()*(2 - 4 * t) + W.pow(2 * t*(1 - t) - 1)*W_prime * 2 * t*(1 - t))
			+ U_prime*(W.pow(2 * t*(1 - t)));

		return q_squad;
	}
	//-------------------------------------------------------------------------------------------
	Real Quaternion::dot(const Quaternion& q2) const
	/*!
	  @brief Quaternion dot product. q1.dot(q2) = s1*s2 + dot(v1, v2)
	*/
	{
		return (_s*q2.s()) + (x()*q2.x()) + (y()*q2.y()) + (z()*q2.z());
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::dq_dt(const ColumnVector & w, const short sign)const
	/*!
	  @brief Quaternion time derivative.
	  
	  The quaternion time derivative, quaternion propagation equation, is
	  $
		\dot{s} = - \frac{1}{2}v^Tw_{_0}
		\\
		\dot{v} = \frac{1}{2}E(s,v)w_{_0}
		\\
		E = sI - S(v)
	  $

	  where $w_{_0}$ is the angular velocity vector expressed in the base
	  frame. If the vector is expressed in the object frame, $w_{_b}$, the
	  time derivative becomes

	  $
		\dot{s} = - \frac{1}{2}v^Tw_{_b}
		\\
		\dot{v} = \frac{1}{2}E(s,v)w_{_b}
		\\
		E = sI + S(v)
	  $
	*/
	{
	/*	Quaternion q;
		Real tmp;

		tmp = -0.5* _v*w;
		q.s_ = tmp;
		q.v_ = 0.5*E(sign)*w;
*/
		return Quaternion(-0.5* _v*w, ~(0.5*E(sign)*w));
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::dq_dt(const PositionVector3D & w, const short sign)const{
		return dq_dt(ColumnVector(w(1), w(2), w(3)), sign);
	}
	//-------------------------------------------------------------------------------------------
	Matrix Quaternion::E(const short sign)const
		/*!
		@brief Matrix E.

		See Quaternion::dot for explanation.
		*/
	{
		Matrix E(3, 3), I(3, 3);
		I.eye();

		if (sign == BODY_FRAME)
			E = _s*I + x_prod_matrix(_v);
		else
			E = _s*I - x_prod_matrix(_v);

		return E;
	}
	//-------------------------------------------------------------------------------------------
	ColumnVector Quaternion::Omega(const Quaternion & q, const Quaternion & q_dot)
		/*!
		@brief Return angular velocity from a quaternion and it's time derivative

		See Quaternion::dot for explanation.
		*/
	{
		Matrix A, B(3,1);
		Matrix U;
		ColumnVector w(3), M;
		A = 0.5*q.E(BASE_FRAME);
		B(1, 1) = q_dot.x();
		B(2, 1) = q_dot.y();
		B(3, 1) = q_dot.z();
		if (A.determinant()){
			QRZ(A, U);             //QR decomposition
			QRZ(A, B, M);
			w = U.InverseLU()*M;
		}
		else
			w = 0;

		return w;
	}
	//-------------------------------------------------------------------------------------------
	short Quaternion::Integ_quat(Quaternion & dquat_present, Quaternion & dquat_past,
		Quaternion & quat, const Real dt)
	//! @brief Trapezoidal quaternion integration.
	{
		if (dt < 0)
		{
			cerr << "Integ_Trap(quat1, quat2, dt): dt < 0. dt is set to 0." << endl;
			return -1;
		}

		// Quaternion algebraic constraint
		//  Real Klambda = 0.5*(1 - quat.norm_sqr());

		dquat_present.s(dquat_present.s());//+ Klambda*quat.s());
		dquat_present.v(dquat_present.v()); //+ Klambda*quat.v());

		quat.s(quat.s() + Integ_Trap_quat_s(dquat_present, dquat_past, dt));
		quat.v(quat.v() + Integ_Trap_quat_v(dquat_present, dquat_past, dt));

		dquat_past.s(dquat_present.s());
		dquat_past.v(dquat_present.v());

		quat.Normalize();

		return 0;
	}
	//-------------------------------------------------------------------------------------------
	Real Quaternion::Integ_Trap_quat_s(const Quaternion & present, Quaternion & past,
		const Real dt)
		//! @brief Trapezoidal quaternion scalar part integration.
	{
		Real integ = 0.5*(present.s() + past.s())*dt;
		past.s(present.s());
		return integ;
	}
	//-------------------------------------------------------------------------------------------
	RowVector Quaternion::Integ_Trap_quat_v(const Quaternion & present, Quaternion & past,
		const Real dt)
		//! @brief Trapezoidal quaternion vector part integration.
	{
		RowVector integ = 0.5*(present.v() + past.v())*dt;
		past.v(present.v());
		return integ;
	}
	//-------------------------------------------------------------------------------------------
	math::Matrix Quaternion::SkewSymmetricMatrix()
	/*!
	  @brief Creates skew symmetric matrix from this quaternion.
	*/
	{
		math::Matrix m(3);
		m << 0,		-z(),		 y(),
			z(),	 0,			-x(),
			-y(),	 x(),		 0;
		return m;
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::exp() const
	/*!
	  @brief The exponential of quaternion. 
	  $
		e^\mathbf{q} = e^{\mathbf{q}_s}*\left ( \cos \left \| \mathbf{v} \right \| 
		+ \frac{\mathbf{v}}{\left \| \mathbf{v} \right \|} 
		* \sin \left \| \mathbf{v} \right \| \right )
	  $
	*/
	{
		Real theta = _v.norm();
	
		if ( fabs(sin(theta)) > EPSilon)
			return Quaternion(::exp(_s)*cos(theta), ::exp(_s) * _v * sin(theta)/theta);
	
		return Quaternion(::exp(_s)*cos(theta), ::exp(_s) * _v);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::ln() const
	/*!
	  @brief The logarithm of quaternion

	  $
		\ln \left ( \mathbf{q} \right ) = 
		\ln \left ( \left \| \mathbf{q} \right \| \right )
		+\frac{\mathbf{v}}{\left \| \mathbf{v} \right \|} 
		* \arccos \left ( \frac{\mathbf{q}_s}{\left \| \mathbf{q} \right \|}  \right )
	  $

	  For unit quaternion, the logarithm simplifies to:

	  $
		\ln \left ( \mathbf{q} \right ) = 0+\mathbf{v}*
		\frac{\arccos \left ( \frac{\mathbf{q}_s}{\left \| \mathbf{q} \right \|}  \right )}
		{\sin \left (\arccos \left ( \frac{\mathbf{q}_s}{\left \| 
		\mathbf{q} \right \|}  \right )\right )}
	  $
	*/
	{
		Real n = norm();
		Real phi = acos( _s/n );
		if (n == ONE)
			return Quaternion(ZERO, _v * phi / sin(phi));

		if ( fabs( sin( phi ) ) > EPSilon)
			return Quaternion(::log(n), _v * phi / _v.norm());
	
		return Quaternion(log(n), _v);
	}
	//-------------------------------------------------------------------------------------------
	Quaternion Quaternion::pow(const Real p) const
	/*!
	  @brief Quaternion power.
	  $
		\mathbf{q}^p = e^{\left ( \ln \left ( \mathbf{q} \right ) * p \right )}
	  $
	*/
	{
		return (ln() * p).exp();
	}
	
	//-------------------------------------------------------------------------------------------
	PositionVector3D Quaternion::OMEGA(Quaternion q0, Quaternion q1){
		Quaternion dq = q1*q0.qinv();
		RowVector wi = (RowVector &)(static_cast<Real>(2) * dq.v() / RowVector::norm(dq.v()) * 
								 acos(dq.s()) / static_cast<Real>(0.01));
		return PositionVector3D(wi(1),wi(2),wi(3));
	}
	Real distance(const Quaternion &Q1, const Quaternion &Q2){
		Real theta = Q1.dot(Q2);
		theta = 2 * acos( (theta > 1 ? 1 : (theta < -1 ? -1 : theta)) );
		return abs(theta);
	}
//-------------------------------------------------------------------------------------------
	Quaternion R2Q(const math::RotationMatrix &R)
	/*!
	  @brief Converts a rotation matrix to a Quaternion and saves the Quaternion in this.
	  
	  @param A 3x3 rotation matrix $ \mathbf{R} $
	  
	  $
	   \begin{array}{c}
		 q_x\\ q_y\\ q_z\\ q_w
	   \end{array}
	   =
	   \left[
	     \begin{array}{c}
	       \\
	       \\
	  
	     \end{array}
	   \right]
	  $
	*/
	{
		Quaternion Q;
		Real s,x,y,z;
		Real ss,n,x1,y1,z1;
		bool c;
		s = R(1,1) + R(2,2) + R(3,3) + ONE;
		if (s < ZERO)
			s = ZERO;
		else
			s = static_cast<Real>(0.5) * sqrt( s );
		x = R(3,2) - R(2,3);
		y = R(1,3) - R(3,1);
		z = R(2,1) - R(1,2);
		if ( (R(1,1) >= R(2,2)) & (R(1,1) >= R(3,3)) ){
			x1 = R(1,1) - R(2,2) - R(3,3) + ONE;
			y1 = R(2,1) + R(1,2);
			z1 = R(3,1) + R(1,3);
			c = (x >= 0);
		}else if ( R(2,2) >= R(3,3) ){
			x1 = R(2,1) + R(1,2);
			y1 = R(2,2) - R(1,1) - R(3,3) + ONE;
			z1 = R(3,2) + R(2,3);
			c = (y >= 0);
		}else{
			x1 = R(3,1) + R(1,3);
			y1 = R(3,2) + R(2,3);
			z1 = R(3,3) - R(1,1) - R(2,2) + ONE;
			c = (z >= 0);
		}
		if (c){
			x += x1;
			y += y1;
			z += z1;
		}else{
			x -= x1;
			y -= y1;
			z -= z1;
		}
		n = sqrt( x*x + y*y + z*z );
		if (n == 0){
			math::RowVector v(3);
			Q.setScalar(ONE);
			v(1) = 0;		v(2) = 0;		v(3) = 0;
			Q.setVector(v);
		}else{
			math::RowVector v(3);
			Q.setScalar(s);
			ss = sqrt( 1 - s*s )/n;
			v(1)=ss*x;		v(2)=ss*y;		v(3)=ss*z;
			Q.setVector(v);
		}
		return Q;
	}
	//-------------------------------------------------------------------------------------------
	math::RotationMatrix Q2R(const Quaternion &QQ)
	/*!
	  @brief Calculates the $ 3\times 3 $ Rotation matrix from a quaternion
	  
	  @return A 3x3 rotation matrix $ \mathbf{R} $

	  $
	  \mathbf{R} =
	   \left[
	    \begin{array}{ccc}
	       1-2(q_y^2-q_z^2) & 2(q_x\ q_y+q_z\ q_w)& 2(q_x\ q_z-q_y\ q_w) \\
	       2(q_x\ q_y-q_z\ q_w) & 1-2(q_x^2-q_z^2) & 2(q_y\ q_z+q_x\ q_w)\\
	       2(q_x\ q_z+q_y\ q_w) & 2(q_y\ q_z-q_x\ q_z) & 1-2(q_x^2-q_y^2)
	     \end{array}
	   \right]
	  $
	*/
	{
		math::RotationMatrix R;
		Real xx,xy,xz,yy,yz,zz,xw,yw,zw;
		Quaternion Q;
		Q = QQ.Normalize();
		Real W = Q.s(), X = Q.x(), Y = Q.y(), Z = Q.z();
	
		xx = X * X;
		xy = X * Y;
		xz = X * Z;
	 	xw = X * W;
	
		yy = Y * Y;
		yz = Y * Z;
		yw = Y * W;
	
		zz = Z * Z;
		zw = Z * W;
	
		R(1,1) = 1 - 2 * (yy + zz);
		R(1,2) = 2 * (xy - zw);
		R(1,3) = 2 * (xz + yw);
	
		R(2,1) = 2 * (xy + zw);
		R(2,2) = 1 - 2 * (xx + zz);
		R(2,3) = 2 * (yz - xw);
	
		R(3,1) = 2 * (xz - yw);
		R(3,2) = 2 * (yz + xw);
		R(3,3) = 1 - 2 * (xx + yy);
		return R;
	}
	//-------------------------------------------------------------------------------------------
	math::Transform3D Q2T(const Quaternion &QQ)
	/*!
	  @brief Calculates the $ 4\times 4 $ Homogeneous matrix from a quaternion
	*/
	{
		return math::Transform3D(Q2R(QQ));
	}
	Quaternion T2Q(const math::Transform3D &Ti) {return R2Q(Ti.R());}
};

