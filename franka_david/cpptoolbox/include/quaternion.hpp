#ifndef ___QUATERNION__H
#define ___QUATERNION__H

#include <string.h>
#include <iostream>
#include <iomanip>
#include "Frame3D.h"

#ifdef E
#undef E
#endif

using namespace std;

namespace math{
#define BASE_FRAME 0
#define BODY_FRAME 1

	class Quaternion;
	Real distance(const Quaternion &Q1, const Quaternion &Q2);
	Quaternion operator *(Real d, const Quaternion & q);
	Quaternion R2Q(const math::RotationMatrix &R);
	Quaternion T2Q(const math::Transform3D &Ti);
	RotationMatrix Q2R(const Quaternion &Q);
	Transform3D Q2T(const Quaternion &Q);

	/*!
	 * @brief A Quaternion $ \mathbf{q}\in \mathbf{R}^4 $ a Complex
	 * number used to describe rotations in 3-dimensional space.
	 * $ q_w+{\bf i}\ q_x+ {\bf j} q_y+ {\bf k}\ q_z $
	 *
	 * Quaternions can be added and multiplied in a similar way as usual
	 * algebraic numbers. Though there are differences. Quaternion
	 * multiplication is not commutative which means
	 * ${\bf Q}\cdot {\bf P} \neq {\bf P}\cdot {\bf Q} $
	 */
	
	class Quaternion{
	private:
		Real _s;
		math::RowVector _v;
	public:
		//!++++++++++++++++++++++++++
		//! Constructors
		//!++++++++++++++++++++++++++
		Quaternion();
		Quaternion(const math::ColumnVector & v);
		Quaternion(const math::PositionVector3D & v);
		Quaternion(const math::RowVector & r);
		Quaternion(const Quaternion & q);
		Quaternion(const math::RotationMatrix & R);
		Quaternion(const math::Transform3D & Ti);
		Quaternion(const Real s, const math::RowVector v);
		Quaternion(const math::RowVector axis, const Real angle);
		Quaternion(const math::ColumnVector axis, const Real angle);
		Quaternion(const math::PositionVector3D axis, const Real angle);
		Quaternion(const Real s, const Real x, const Real y, const Real z);
		Quaternion(const Real rzPsi, const Real ryTheta, const Real rxPhi);
		~Quaternion(){}
		//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//! ++++++++++++++++++++++++ Operators overload ++++++++++++++++++++++++
		//! Equal and comparison operators
		Quaternion & operator =(const Quaternion & Q);
		Quaternion & operator =(const math::RotationMatrix & R);
		inline bool operator ==(const Quaternion& Q) const{ return (_s == Q._s) && (_v == Q._v);}
		inline bool operator !=(const Quaternion& Q) const{ return !(*this == Q);}
		//! Arithmetic operators (quaternion and scalar)
		const Quaternion operator  /(Real n)const;
		const Quaternion operator  *(Real n)const	{ return Quaternion(_s*n, _v*n);	}
		const Quaternion&operator /=(Real n)		{ return (*this) = (*this) / n;		}
		const Quaternion&operator *=(Real n)		{ return (*this) = (*this) * n;		}
		friend Quaternion operator * (Real n, const Quaternion & q){	return q * n;	}
		//! Arithmetic operators (quaternion and quaternion)
		const Quaternion operator  +(const Quaternion &Q)const;
		const Quaternion operator  -(const Quaternion &Q)const;
		const Quaternion operator  *(const Quaternion &Q)const;
		const Quaternion operator  /(const Quaternion &r)const;
		const Quaternion&operator +=(const Quaternion &Q)	{ return (*this) = (*this) + Q; }
		const Quaternion&operator -=(const Quaternion &Q)	{ return (*this) = (*this) + Q; }
		const Quaternion&operator *=(const Quaternion &r)	{ return (*this) = (*this) * r; }
		const Quaternion&operator /=(const Quaternion &r)	{ return (*this) = (*this) / r; }
		const Quaternion operator  -()const					{ return *this * -1;			}
		friend ostream & operator << (ostream & o, const Quaternion & Q);
		//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//! ++++++++++++++++++++++++++ Get functions ++++++++++++++++++++++++++
		const math::RowVector & v()	const	 {	return _v;	}   //! Read 3D vector  of the quaternion
		math::RowVector & v()	         	 {	return _v;	}   //! Access 3D vector  of the quaternion
		const math::RowVector & rVec()	const{	return _v;	}   //! Read 3D vector  of the quaternion
		math::RowVector & rVec()	         {	return _v;	}   //! Access 3D vector  of the quaternion
		const math::PositionVector3D cVec()const{	return PositionVector3D(_v(1),_v(2),_v(3));	}   //! Read 3D vector  of the quaternion
		inline const Real & s()const		 {	return _s;	}	//! Read scalar of the quaternion
		inline Real & s()					 {	return _s;	}	//! Access scalar of the quaternion
		inline const Real & x()const		 {	return _v(1);	}	//! Read x coordinate of the quaternion
		inline Real & x()					 {	return _v(1);	}	//! Access x coordinate of the quaternion
		inline const Real & y()const		 {	return _v(2);	}	//! Read y coordinate of the quaternion
		inline Real & y()					 {	return _v(2);	}	//! Access y coordinate of the quaternion
		inline const Real & z()const		 {	return _v(3);	}	//! Read z coordinate of the quaternion
		inline Real & z()					 {	return _v(3);	}	//! Access z coordinate of the quaternion
		//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//! ++++++++++++++++++++++++++ Set functions ++++++++++++++++++++++++++
		void setQuaternion(const Real scalar, const math::RowVector &v);
		void setScalar(const Real scalar)				{	_s = scalar;	}
		void setVector(const math::RowVector & v);
		void v(const math::RowVector & v)	{ setVector(v); }
		void s(const Real s)	{ _s = s; }	//! Set s coordinate of the quaternion
		void x(const Real x)	{	_v(1) = x;	}	//! Set x coordinate of the quaternion
		void y(const Real y)	{	_v(2) = y;	}	//! Set y coordinate of the quaternion
		void z(const Real z)	{	_v(3) = z;	}	//! Set z coordinate of the quaternion
		//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//! +++++++++++++++++++++++++++++ Methods +++++++++++++++++++++++++++++
		bool is_nan()const;
		bool is_neg_inf() const;
		bool is_pos_inf() const;
		bool is_zero() const;
		bool is_inf() const;
		Real normSquared()const;
		Real normSquared(const Quaternion & Q)			{	return Q.normSquared();			}

		Real norm()const;
		Real norm(const Quaternion & Q)					{	return sqrt(Q.normSquared());	}

		Quaternion & Normalize();
		const Quaternion  Normalize()const				{	return (*this)/norm();			}
		static Quaternion Normalize(Quaternion &Q)		{	return Q.Normalize();			}
		static Quaternion Normalize(const Quaternion &Q){	return Q.Normalize();			}

		Quaternion conjugate()const;
		static Quaternion conjugate(const Quaternion &Q){return Q.conjugate();			}

		Quaternion qinv()const;
		static Quaternion qinv(const Quaternion &Q)	{ return Q.qinv();}

		math::RotationMatrix R(){return Q2R(*this);}

		void Angle2Quat(Real rzPsi, Real ryTheta, Real rxPhi, char *type=0);
		static void Quat2Angle(const Quaternion &q, Real &rzPsi, Real &ryTheta, Real &rxPhi);

		Real dot(const Quaternion& q2) const;
		static inline Real dot(const Quaternion& q1, const Quaternion& q2){return q1.dot(q2);}

		//!	qinterp: Quaternion Spherical Linear Interpolation.
		static Quaternion  qinterp(Quaternion Q1, Quaternion Q2, Real r);
		//	lerp: Quaternion Linear Interpolation
		static Quaternion lerp(const Quaternion &q1, const Quaternion &q2, Real r)
			{ return (q1*(1-r) + q2*r).Normalize(); }
		static Quaternion Slerp_prime(const Quaternion & q0, const Quaternion & q1, const Real t);
		static Quaternion Squad(const Quaternion & p, const Quaternion & a, const Quaternion & b,
			const Quaternion & q, const Real t);
		static Quaternion Squad_prime(const Quaternion & p, const Quaternion & a, const Quaternion & b,
			const Quaternion & q, const Real t);

		Quaternion dq_dt(const ColumnVector & w, const short sign)const;
		Quaternion dq_dt(const PositionVector3D & w, const short sign)const;
		Matrix E(const short sign)const;
		static ColumnVector Omega(const Quaternion & q, const Quaternion & q_dot);
		static short Integ_quat(Quaternion & dquat_present, Quaternion & dquat_past, Quaternion & quat, const Real dt);
		static Real Integ_Trap_quat_s(const Quaternion & present, Quaternion & past, const Real dt);
		static RowVector Integ_Trap_quat_v(const Quaternion & present, Quaternion & past, const Real dt);

		void print(char const* ch = 0, int prec = 4, int w = 8) const;

		math::Matrix SkewSymmetricMatrix();

		Quaternion exp() const;
		Quaternion ln() const;
		Quaternion pow(const Real p) const;
		static math::PositionVector3D OMEGA(Quaternion q0, Quaternion q1);
	};
};

#endif
