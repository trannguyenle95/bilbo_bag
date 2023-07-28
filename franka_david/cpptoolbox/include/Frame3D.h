#ifndef ___FRAME_3D_H
#define ___FRAME_3D_H


#include <cmath>
#include "matrix.hpp"
#include "utility_io.h"

namespace math{
	class PositionVector3D;
	class RotationMatrix;
	class Transform3D;
	class Twist;
	class Wrench;
/*********************************************************************************/
/*********************************************************************************/
/******************************  Position Vector3D  ******************************/
/*********************************************************************************/
/*********************************************************************************/
	/*!
	  @brief A 3D vector $ \mathbf{v}\in \mathbf{R}^3 $
	  
	  $ ^{j}\mathbf{v}_i = \left[
		\begin{array}{c}
			v_x \\
			v_y \\
			v_z
		\end{array}
	  \right]
	  $

	  Usage example:
	  const PositionVector3D v1(1.0, 2.0, 3.0);
	  const PositionVector3D v2(6.0, 7.0, 8.0);
	  const PositionVector3D v3 = cross(v1, v2);
	  const double d = dot(v1, v2);
	  const PositionVector3D v4 = v2 - v1;
	*/
	class PositionVector3D{
	private:
		Real _v[3];
	public:
		//!Constructors
		PositionVector3D(Real x=0, Real y=0, Real z=0);
		PositionVector3D(const std::string &str) { set(str);}
		PositionVector3D(const math::ColumnVector &v);
		PositionVector3D(const math::RowVector &v);
		PositionVector3D(const PositionVector3D &p);
		~PositionVector3D(){}

		//!Set functions
		void set(const std::string &str);
		PositionVector3D & set(const PositionVector3D &p)	{*this = PositionVector3D(p); return *this;		}
		PositionVector3D & set(Real x, Real y, Real z){*this = PositionVector3D(x,y,z); return *this;	}

		//!operators
		Real operator () (int i) const;
		Real & operator () (int i);
		Real operator [] (int i) const;
		inline Real & operator [] (int i);
		PositionVector3D & operator =(const PositionVector3D & m);
		PositionVector3D & operator =(const math::ColumnVector & m);
		PositionVector3D & operator =(const math::RowVector & m);
		PositionVector3D & operator =(Real x);
		bool operator ==(const PositionVector3D &m)const;
		bool operator !=(const PositionVector3D &m)const{return !(*this == m);}
		const PositionVector3D operator / ( Real s) const;
		const PositionVector3D operator * ( Real s) const;
		const Real operator * ( const PositionVector3D &s) const;
		const PositionVector3D operator -(const PositionVector3D &m)const;
		const PositionVector3D operator +(const PositionVector3D &m)const;
		PositionVector3D operator *=(Real s)    { return *this = *this * s;	}
		PositionVector3D operator /=(Real s)    { return *this = *this / s;	}
		const PositionVector3D operator +=(const PositionVector3D &p){ return *this = *this + p;	}
		const PositionVector3D operator -=(const PositionVector3D &p){ return *this = *this - p;	}
		const PositionVector3D operator-() const;
		const PositionVector3D   operator / (const PositionVector3D &v)const{return ElementWiseDiv(v);	}
		const PositionVector3D & operator /=(const PositionVector3D &v)		{return *this = *this / v;	}

		friend const PositionVector3D operator * (Real s, const PositionVector3D& v);
		friend const PositionVector3D operator * (const Matrix &s, const PositionVector3D& v);
		friend const Real operator * (const RowVector &s, const PositionVector3D& v);
		friend ostream & operator << (ostream & o, const PositionVector3D & m);
		PositionVector3D & operator >> (const Real *val );

		//!Methods
		void reverseSign();
		void print (char const* ch = 0, int prec=4, int w=8) const;
		double Normalize(double eps=EPSilon);
		Real norm2() const;
		Real norm() const{return norm2();}
		Real norm1() const;
		Real normInf() const;
		bool equal(const PositionVector3D& rot, Real precision=0.000001)const;
		bool equal(const PositionVector3D& rot1, const PositionVector3D& rot2, Real precision=0.000001);
		inline friend void SetToZero(PositionVector3D& v){v = PositionVector3D::zero();}
		inline static PositionVector3D zero(){ return PositionVector3D(0,0,0); }
		static PositionVector3D ElementWiseMul(const PositionVector3D &p1, const PositionVector3D &p2);
		const PositionVector3D ElementWiseDiv(const PositionVector3D &v)const;
		static PositionVector3D ElementWiseDiv(const PositionVector3D &v1, const PositionVector3D &v2){return v1.ElementWiseDiv(v2);}
		Real x() const	{	return _v[0];	}
		Real y() const	{	return _v[1];	}
		Real z() const	{	return _v[2];	}
		Real x(Real x_) {	return _v[0] = x_;	}
		Real y(Real y_) {	return _v[1] = y_;	}
		Real z(Real z_) {	return _v[2] = z_;	}
		inline Real * getData()	{return _v;}
		inline const Real * getData()const	{return _v;}
		void setElements(Real *a, int n)	{_v[0] = a[0];	_v[1] = a[1];	_v[2] = a[2];}
		static Real SumSquare(const PositionVector3D &m)	{ return m.SumSquare();	}
		Real SumSquare()const;
		Matrix skewMatrix()const;
		static Matrix skewMatrix(const PositionVector3D &v)	{return v.skewMatrix();	}
		friend std::istream& operator >> (std::istream& is, PositionVector3D& v){
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
	}
	};
	void random(PositionVector3D &v);
	bool equal(const PositionVector3D& rot1, const PositionVector3D& rot2, Real precision);
	//! Distance Between Two ColumnVectos
	Real distance(const PositionVector3D & v1, const PositionVector3D & v2);
	//! Return the 3D vector cross product
	const PositionVector3D cross(const PositionVector3D& v1, const PositionVector3D& v2);
	void cross(const PositionVector3D& v1, const PositionVector3D& v2, PositionVector3D& v3);
	//! Return the dot product
	Real dot(const PositionVector3D& v1, const PositionVector3D& v2);
	//! Returns the normalized vector
	const PositionVector3D normalize(const PositionVector3D& v);
	//! Calculates the angle from v1 to v2 around the axis defined by v1 x v2 with n determining the sign.
	Real angle(const PositionVector3D& v1, const PositionVector3D& v2, const PositionVector3D& n);
	//! Calculates the angle from v1 to v2 around the axis defined by v1 x v2
	Real angle(const PositionVector3D& v1, const PositionVector3D& v2);
	//! Return the Unit Vector for a given Vector v.
	PositionVector3D UnitVector(const PositionVector3D & v);

/*********************************************************************************/
/*********************************************************************************/
/******************************  Rotation Matrix3D  ******************************/
/*********************************************************************************/
/*********************************************************************************/
	/*!
	  @brief A 3x3 rotation matrix $ \mathbf{R}\in SO(3) $
	  $
	  \mathbf{R}=
	  \left[
		\begin{array}{ccc}
		{}^A\hat{X}_B & {}^A\hat{Y}_B & {}^A\hat{Z}_B
		\end{array}
	  \right]
	  =
	  \left[
		\begin{array}{ccc}
			r_{11} & r_{12} & r_{13} \\
			r_{21} & r_{22} & r_{23} \\
			r_{31} & r_{32} & r_{33}
		\end{array}
	  \right]
	  $
	*/
	
	class RotationMatrix{
	private:
		Real _R[9];
	public:
		//!Default constructor, creates identity rotation matrix
		RotationMatrix();
		RotationMatrix(	Real r11, Real r12, Real r13,
						Real r21, Real r22, Real r23,
						Real r31, Real r32, Real r33);
		RotationMatrix(	const PositionVector3D& i,
						const PositionVector3D& j,
						const PositionVector3D& k);
		RotationMatrix(Real *x);
		RotationMatrix(const RotationMatrix &m);
		~RotationMatrix(){}

		// Operators
		RotationMatrix & operator =(const RotationMatrix & m);
		bool operator ==(const RotationMatrix & m)const;
		bool operator !=(const RotationMatrix & m)const{return !(*this == m);}
		const Real & operator () (int i, int j) const;
		Real & operator () (int i, int j);
		const Real & operator [] (int i) const;
		Real & operator [] (int i);
		const RotationMatrix operator *(const RotationMatrix &m)const;
		const PositionVector3D operator *(const PositionVector3D &v)const;
		const RotationMatrix operator -(const RotationMatrix &m)const;
		const RotationMatrix operator -()const	{ return (*this) * -1;}
		const RotationMatrix operator +(const RotationMatrix &m)const;
		const RotationMatrix operator / ( Real s) const;
		const RotationMatrix operator * ( Real s) const;
		const RotationMatrix operator +=(const RotationMatrix &m);
		const RotationMatrix operator -=(const RotationMatrix &m);
		RotationMatrix operator *=(Real s);
		RotationMatrix operator /=(Real s);
		RotationMatrix operator ~()const;
		friend const RotationMatrix operator * (Real s, const RotationMatrix& v);
		friend const RowVector operator * (const PositionVector3D &v, const RotationMatrix& s);
		friend ostream & operator << (ostream & o, const RotationMatrix & m);

		inline Twist operator * (const Twist& arg) const;
		inline Wrench operator * (const Wrench& arg) const;
		//! Methods
		void print (char const* ch = 0, int prec=4, int w=8) const;
		//! The data array
		inline Real * getData()				{ return _R;	}
		inline const Real * getData()const	{ return _R;	}
		static RotationMatrix identity();
		static RotationMatrix I(){return RotationMatrix::identity();}
		RotationMatrix & normalizeR();
		const PositionVector3D getCol(int col) const;
		const PositionVector3D ii() const{ return PositionVector3D(_R[0],_R[3],_R[6]);}
		const PositionVector3D jj() const{ return PositionVector3D(_R[1],_R[4],_R[7]);}
		const PositionVector3D kk() const{ return PositionVector3D(_R[2],_R[5],_R[8]);}
		bool equal(const RotationMatrix& rot, Real precision=0.000001)const;
		friend bool equal(const RotationMatrix& rot1, const RotationMatrix& rot2, Real precision);
		void zeros ();
		static RotationMatrix skew(const PositionVector3D& v);
		RotationMatrix t()const{return ~(*this);}
		RotationMatrix Transpose()const{return ~(*this);}
		RotationMatrix Inverse()const;
		RotationMatrix i()const{return RotationMatrix::Inverse();}
		PositionVector3D Inverse(const PositionVector3D& v) const;
		inline void SetInverse();
		inline Twist Inverse(const Twist& arg) const;
		inline Wrench Inverse(const Wrench& arg) const;
		//! Rotating methods.
		static PositionVector3D r2rpy(RotationMatrix R);
		PositionVector3D r2rpy(){return r2rpy(*this);}
		static RotationMatrix rpy2r(Real roll, Real pitch, Real yaw);
		void getRPY(Real& roll, Real& pitch, Real& yaw) const;

		void doRotX(Real t);
		void doRotY(Real t);
		void doRotZ(Real t);
		static RotationMatrix roty(Real t);
		static RotationMatrix rotx(Real t);
		static RotationMatrix rotz(Real t);

		static RotationMatrix Rot(const PositionVector3D& rv,Real angle);
		//! Along an arbitrary axes.  rotvec should be normalized.
		static RotationMatrix Rot2(const PositionVector3D& rotvec, Real angle);

		static RotationMatrix eulerzyz(Real ALPHA, Real BETA, Real GAMA);
		void getEulerZYZ(Real& ALPHA, Real& BETA, Real& GAMA) const;
		static RotationMatrix eulerzyx(Real ALPHA, Real BETA, Real GAMA);
		void getEulerZYX(Real& ALPHA, Real& BETA, Real& GAMA) const;

		static ColumnVector irotk(const RotationMatrix & R);
		double GetRotAngle(PositionVector3D& axis, double eps = 0.000001) const;
		PositionVector3D GetRot() const;

		Matrix subMatrix() const;
		//! Access to the underlying unitvectors of the rotation matrix
		inline PositionVector3D UnitX() const		{ return PositionVector3D(_R[0],_R[3],_R[6]);	}
		inline PositionVector3D UnitY() const		{ return PositionVector3D(_R[1],_R[4],_R[7]);	}
		inline PositionVector3D UnitZ() const		{ return PositionVector3D(_R[2],_R[5],_R[8]);	}
		inline void UnitX(const PositionVector3D& X){ _R[0] = X(0); _R[3] = X(1); _R[6] = X(2);		}
		inline void UnitY(const PositionVector3D& X){ _R[1] = X(0); _R[4] = X(1); _R[7] = X(2);		}
		inline void UnitZ(const PositionVector3D& X){ _R[2] = X(0); _R[5] = X(1); _R[8] = X(2);		}

		static RotationMatrix uQuaternion(double x,double y,double z, double w);
		Matrix mat()const;
		friend std::istream& operator >> (std::istream& is, RotationMatrix& r){
		IOTrace("Stream input Rotation (Matrix or EULERZYX, EULERZYZ,RPY, ROT, IDENTITY)");
		char storage[10];
		EatWord(is, "[]", storage, 10);
		if (strlen(storage) == 0) {
			Eat(is, '[');
			for (int i = 1; i<=3; i++) {
				is >> r(i, 1);
				Eat(is, ',');
				is >> r(i, 2);
				Eat(is, ',');
				is >> r(i, 3);
				if (i<3)
					Eat(is, ';');
				else
					EatEnd(is, ']');
			}
			IOTracePop();
			return is;
		}
		PositionVector3D v;
		const double deg2rad = 0.01745329251994329576923690768488;
		const double rad2deg = 57.2957795130823208767981548141052;
		if (strcmp(storage, "EULERZYX") == 0) {
			is >> v;
			v = v*deg2rad;
			r = RotationMatrix::eulerzyx(v(1), v(2), v(3));
			IOTracePop();
			return is;
		}
		if (strcmp(storage, "EULERZYZ") == 0) {
			is >> v;
			v = v*deg2rad;
			r = RotationMatrix::eulerzyz(v(1), v(2), v(3));
			IOTracePop();
			return is;
		}
		if (strcmp(storage, "RPY") == 0) {
			is >> v;
			v = v*deg2rad;
			r = RotationMatrix::rpy2r(v(1), v(2), v(3));
			IOTracePop();
			return is;
		}
		if (strcmp(storage, "ROT") == 0) {
			is >> v;
			double angle;
			Eat(is, '[');
			is >> angle;
			EatEnd(is, ']');
			r = RotationMatrix::Rot(v, angle*deg2rad);
			IOTracePop();
			return is;
		}
		if (strcmp(storage, "IDENTITY") == 0) {
			r = RotationMatrix();
			IOTracePop();
			return is;
		}
		throw Error_Frame_Rotation_Unexpected_id();
		return is;
	}
	};
	void random(RotationMatrix &r);
	bool equal(const RotationMatrix& rot1, const RotationMatrix& rot2, Real precision=0.000001);

/*********************************************************************************/
/*********************************************************************************/
/***************************  Transformation Matrix3D  ***************************/
/*********************************************************************************/
/*********************************************************************************/
/*!
  @brief A 4x4 homogeneous transform matrix @f$ \mathbf{T}\in SE(3) @f$
		 $
		 \mathbf{T} =
		 \left[
			\begin{array}{cc}
			\mathbf{R} & \mathbf{d} \\
			\begin{array}{ccc}0 & 0 & 0\end{array} & 1
			\end{array}
		 \right]
		 $
*/
	class Transform3D{
	private:
		PositionVector3D _P;
		RotationMatrix _R;
    public:
		//++++++++++++++++++++++++++
		// Constructors
		//++++++++++++++++++++++++++
		Transform3D();
		Transform3D(const PositionVector3D& P, const RotationMatrix& R);
		Transform3D(const Real * CartPos);
		Transform3D(const RotationMatrix& R);
		Transform3D(const PositionVector3D& P);
		Transform3D(const Transform3D &m) : _P(m.P()), _R(m.R()) {}
		~Transform3D(){}
		//++++++++++++++++++++++++++
		// Get functions
		//++++++++++++++++++++++++++
		const RotationMatrix & R()const;
		RotationMatrix & R();
		const PositionVector3D & P()const;
		PositionVector3D & P();
		//++++++++++++++++++++++++++
		// Operators overload
		//++++++++++++++++++++++++++
		Real operator () (int i, int j) const;
		Real & operator () (int i, int j);

		Transform3D & operator =(const Transform3D & m);
		Transform3D & operator =(const RotationMatrix & m);
		Transform3D & operator =(const PositionVector3D & m);

		const Transform3D operator *(const Transform3D & tm)const;
		const Transform3D operator -(const Transform3D & tm)const;
		const Transform3D operator +(const Transform3D & tm)const;
		const Transform3D operator -()const{ return (*this) * -1;}

		friend const Transform3D operator * (Real s, const Transform3D& v);

		const Transform3D operator / ( Real s) const;
		const Transform3D operator * ( Real s) const;
		const Transform3D & operator +=(const Transform3D &m);
		const Transform3D & operator -=(const Transform3D &m);

		const Transform3D & operator *=(Real s);
		const Transform3D & operator /=(Real s);

		const PositionVector3D operator*( const PositionVector3D& bP) const;
		friend ostream & operator << (ostream & o, const Transform3D & m);

		inline Wrench operator * (const Wrench& arg) const;
		inline Twist operator * (const Twist& arg) const;
		//++++++++++++++++++++++++++
		// Methods
		//++++++++++++++++++++++++++
		void print (char const* ch = 0, int prec=4, int w=8) const;
		//computes the inverse of t1 and multiplies it with t2. t3 = inv(t1) * t2
		static Transform3D invMult(const Transform3D& t1, const Transform3D& t2);

		Transform3D Inverse()const;
		static Transform3D Inverse(Transform3D Ti){return Ti.Inverse();}
		Transform3D i()const{return Transform3D::Inverse();}
		//! The same as p2=R.Inverse()*p but more efficient.
		inline PositionVector3D Inverse(const PositionVector3D& p) const{return _R.Inverse(p-_P);}
		//! The same as p2=R.Inverse()*p but more efficient.
		inline Wrench Inverse(const Wrench& arg) const;
		//! The same as p2=R.Inverse()*p but more efficient.
		inline Twist  Inverse(const Twist& arg) const;
		//Constructs a homogeneous transform using the original Denavit-Hartenberg notation
		static const Transform3D DH(Real alpha, Real a, Real d, Real theta);
		//Constructs a homogeneous transform using the Craig (modified) Denavit-Hartenberg notation
		static const Transform3D modifiedDH(Real alpha, Real a, Real d, Real theta);
		//Constructs the identity transform
		static Transform3D identity();
		static Transform3D I(){return Transform3D::identity();}
		//Performs an element wise comparison. Two elements are considered equal if the difference
		//are less than precision.
		bool equal(const Transform3D& T, Real precision=epsiLON);
		friend bool equal(const Transform3D& a,const Transform3D& b,double precision);

		static Transform3D MakeT(RotationMatrix &R, PositionVector3D X);
		static Transform3D MakeT(Transform3D &T0, PositionVector3D X);
		// Homogeneous transform for rotation about Y-axis, X-axis, and Z-axis
		static Transform3D roty(Real t);
		static Transform3D rotx(Real t);
		static Transform3D rotz(Real t);

		static Transform3D rotk(Real theta, const PositionVector3D & k);
		static Transform3D rotd(Real theta, const PositionVector3D & k1,const PositionVector3D & k2);
		static Transform3D rpy(const PositionVector3D & a);
		static Transform3D eulzxz(const PositionVector3D & a);
		static ColumnVector irotk(const Transform3D & R);
		static PositionVector3D irpy(const Transform3D & R);
		static PositionVector3D ieulzxz(const Transform3D & R);
		// Euler angle to homogeneous transform
		static Transform3D eul2tr(Real phi , Real theta, Real psi);
		// Roll/pitch/yaw angles to homogeneous transform
		static Transform3D rpy2tr(Real roll, Real pitch, Real yaw);
		// Homogeneous transform to roll/pitch/yaw angles
		static PositionVector3D tr2rpy(Transform3D Ti);
		PositionVector3D tr2rpy(){return tr2rpy(*this);}
		// Homogeneous transform to Euler angles
		static PositionVector3D tr2eul(Transform3D Ti);
		PositionVector3D tr2eul(){return tr2eul(*this);}
		// Normalize a homogeneous transform
		static Transform3D normalizeT(Transform3D Ti);
		Transform3D normalizeT(){return *this = normalizeT((*this));}

		//static Transform3D *  ctraj(const Transform3D & t0, const Transform3D & t1, Real n);
		//static Transform3D trinterp(const Transform3D & T0, const Transform3D & T1, Real r);

		Real * getCartPos()const;
		void getCartPos(Real CartPos[]);
		static Transform3D CV2T(const math::ColumnVector & v);
		static math::ColumnVector T2CV(const Transform3D &Ti);

		//! Get Matrix(r1:r2, c1:c2);
		Matrix subMatrix() const;
		//! Get Matrix as ColumnVector
		//ColumnVector getMatrixAsColumn() const;
		inline void Integrate(const Twist& t_this,double frequency);

		friend std::istream& operator >> (std::istream& is, Transform3D& T){
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
	}
	};
	void random(Transform3D &);
	void jtraj(const math::ColumnVector &q0, const ColumnVector &q1, const ColumnVector &tv,
		Matrix &qt, Matrix &qdt,Matrix &qddt);
	bool equal(const Transform3D& a,const Transform3D& b,double precision=epsiLON);
	void distance(const Transform3D &T0, const Transform3D &T1, Real &dist, Real &angle);
	Transform3D * ctraj(const Transform3D & T0, const Transform3D & T1, Real n);
	Transform3D trinterp(const Transform3D & T0, const Transform3D & T1, Real r);

/*********************************************************************************/
/*********************************************************************************/
/************************************  Twist  ************************************/
/*********************************************************************************/
/*********************************************************************************/
	/**
	 * \brief represents both translational and rotational velocities.
	 *
	 * This class represents a twist.  A twist is the combination of translational
	 * velocity and rotational velocity applied at one point.
	*/
	class Twist {
	private:
		PositionVector3D _vel; //!< The velocity of that point
		PositionVector3D _rot; //!< The rotational velocity of that point.
	public:
		
		//! The default constructor initialises to Zero via the constructor of Vector.
		Twist():_vel(),_rot() {};
		
		Twist(const PositionVector3D& vel,const PositionVector3D& rot):_vel(vel),_rot(rot) {};
		
		inline Twist& operator-=(const Twist& arg);
		inline Twist& operator+=(const Twist& arg);
		//! index-based access to components, first vel(0..2), then rot(3..5)
		double& operator()(int i);
		
		PositionVector3D Vel()const {return _vel;}
		PositionVector3D & Vel() {return _vel;}
		PositionVector3D Rot()const {return _rot;}
		PositionVector3D & Rot() {return _rot;}

		Real Vel(int i)const { return _vel(i); }
		Real & Vel(int i) { return _vel(i); }
		Real Rot(int i)const { return _rot(i); }
		Real & Rot(int i) { return _rot(i); }
		//! index-based access to components, first vel(0..2), then rot(3..5)
		//! For use with a const Twist
		double operator()(int i) const;
		
		inline double operator[] ( int index ) const
		  {
		return this->operator() ( index );
		  }
		
		inline double& operator[] ( int index )
		  {
		return this->operator() ( index );
		  }
		
		//! @return a zero Twist : Twist(Vector::Zero(),Vector::Zero())
		static Twist Zero();

		friend Twist operator*(const Twist& lhs,double rhs);
		friend Twist operator*(double lhs,const Twist& rhs);
		friend Twist operator/(const Twist& lhs,double rhs);
		friend Twist operator+(const Twist& lhs,const Twist& rhs);
		friend Twist operator-(const Twist& lhs,const Twist& rhs);
		friend Twist operator-(const Twist& arg);
		friend double dot(const Twist& lhs,const Wrench& rhs);
		friend double dot(const Wrench& rhs,const Twist& lhs);
		friend void SetToZero(Twist& v);
		/// Spatial cross product for 6d motion vectors, beware all of them have to be expressed in the same reference frame/point
		friend Twist operator*(const Twist& lhs,const Twist& rhs);
		/// Spatial cross product for 6d force vectors, beware all of them have to be expressed in the same reference frame/point
		friend Wrench operator*(const Twist& lhs,const Wrench& rhs);
		
		
		//! Reverses the sign of the twist
		void ReverseSign();
		
		//! Changes the reference point of the twist.
		//! The vector v_base_AB is expressed in the same base as the twist
		//! The vector v_base_AB is a vector from the old point to
		//! the new point.
		//!
		//! Complexity : 6M+6A
		Twist RefPoint(const PositionVector3D& v_base_AB) const;
		
		
		//! do not use operator == because the definition of Equal(.,.) is slightly
		//! different.  It compares whether the 2 arguments are equal in an eps-interval
		friend bool equal(const Twist& a,const Twist& b,double eps);
		
		//! The literal equality operator==(), also identical.
		friend bool operator==(const Twist& a,const Twist& b);
		//! The literal inequality operator!=().
		friend bool operator!=(const Twist& a,const Twist& b);

		void print(char const* ch = 0, int prec = 4, int w = 8) const;

		friend std::istream& operator >> (std::istream& is, Twist& v){
		IOTrace("Stream input Twist");
		Eat(is, '[');
		is >> v(1);
		Eat(is, ',');
		is >> v(2);
		Eat(is, ',');
		is >> v(3);
		Eat(is, ',');
		is >> v(4);
		Eat(is, ',');
		is >> v(5);
		Eat(is, ',');
		is >> v(6);
		EatEnd(is, ']');
		IOTracePop();
		return is;
	}
		friend ostream & operator << (ostream & o, const Twist & m);
	};
	void random(Twist &);
	inline bool equal(const Twist& a,const Twist& b,double eps=epsiLON);

/*********************************************************************************/
/*********************************************************************************/
/************************************  Wrench  ***********************************/
/*********************************************************************************/
/*********************************************************************************/
	/**
	 * \brief represents the combination of a force and a torque.
	 *
	 * This class represents a Wrench.  A Wrench is the force and torque applied at a point
	 */
	class Wrench
	{
	private:
		PositionVector3D _force;       //!< Force that is applied at the origin of the current ref frame
		PositionVector3D _torque;      //!< Torque that is applied at the origin of the current ref frame
	public:
	
		//! Does  initialise force and torque to zero via the underlying constructor of Vector
		Wrench():_force(),_torque() {};
		Wrench(const PositionVector3D& force,const PositionVector3D& torque):_force(force),_torque(torque) {};
	
		// = Operators
		inline Wrench& operator-=(const Wrench& arg);
		inline Wrench& operator+=(const Wrench& arg);
		
		//! index-based access to components, first force(0..2), then torque(3..5)
		double& operator()(int i);
		
		PositionVector3D Force()const {return _force;}
		PositionVector3D & Force() {return _force;}
		PositionVector3D Torque()const {return _torque;}
		PositionVector3D & Torque() {return _torque;}
		//! index-based access to components, first force(0..2), then torque(3..5)
		//! for use with a const Wrench
		double operator()(int i) const;
		
		double operator[] ( int index ) const
		  {
		return this->operator() ( index );
		  }
		
		double& operator[] ( int index )
		  {
		return this->operator() ( index );
		  }
		
		//! Scalar multiplication
		inline friend Wrench operator*(const Wrench& lhs,double rhs);
		//! Scalar multiplication
		inline friend Wrench operator*(double lhs,const Wrench& rhs);
		//! Scalar division
		inline friend Wrench operator/(const Wrench& lhs,double rhs);
		
		inline friend Wrench operator+(const Wrench& lhs,const Wrench& rhs);
		inline friend Wrench operator-(const Wrench& lhs,const Wrench& rhs);
		
		//! An unary - operator
		inline friend Wrench operator-(const Wrench& arg);
		
		//! Sets the Wrench to Zero, to have a uniform function that sets an object or
		//! double to zero.
		inline friend void SetToZero(Wrench& v);
		
		//! @return a zero Wrench
		static inline Wrench Zero();
		
		//! Reverses the sign of the current Wrench
		inline void ReverseSign();
		
		//! Changes the reference point of the wrench.
		//! The vector v_base_AB is expressed in the same base as the twist
		//! The vector v_base_AB is a vector from the old point to
		//! the new point.
		//!
		//! Complexity : 6M+6A
		inline Wrench RefPoint(const PositionVector3D& v_base_AB) const;
		
		
		//! do not use operator == because the definition of Equal(.,.) is slightly
		//! different.  It compares whether the 2 arguments are equal in an eps-interval
		inline friend bool equal(const Wrench& a,const Wrench& b,double eps);
		
		//! The literal equality operator==(), also identical.
		inline friend bool operator==(const Wrench& a,const Wrench& b);
		//! The literal inequality operator!=().
		inline friend bool operator!=(const Wrench& a,const Wrench& b);

		void print(char const* ch = 0, int prec = 4, int w = 8) const;

		friend std::istream& operator >> (std::istream& is, Wrench& v){
		IOTrace("Stream input Wrench");
		Eat(is, '[');
		is >> v(1);
		Eat(is, ',');
		is >> v(2);
		Eat(is, ',');
		is >> v(3);
		Eat(is, ',');
		is >> v(4);
		Eat(is, ',');
		is >> v(5);
		Eat(is, ',');
		is >> v(6);
		EatEnd(is, ']');
		IOTracePop();
		return is;
	}
	};
	void random(Wrench &);
	bool equal(const Wrench& a,const Wrench& b,double eps=epsiLON);

	// Return Twist in ColumnVector
	void Twist2CV(const Twist& t, ColumnVector& e);
	ColumnVector Twist2CV(const Twist& t);
	PositionVector3D diff(const PositionVector3D& a, const PositionVector3D& b, double dt=1);

	PositionVector3D diff(const RotationMatrix& R_a_b1, const RotationMatrix& R_a_b2, double dt=1);

	Twist diff(const Transform3D& F_a_b1, const Transform3D& F_a_b2, double dt = 1);
	PositionVector3D addDelta(const PositionVector3D& p_w_a,const PositionVector3D& p_w_da,double dt=1);

	RotationMatrix addDelta(const RotationMatrix& R_w_a,const PositionVector3D& da_w,double dt=1);

	Transform3D addDelta(const Transform3D& F_w_a,const Twist& da_w,double dt=1);

	Twist addDelta(const Twist& a,const Twist&da,double dt=1);

	Wrench addDelta(const Wrench& a,const Wrench&da,double dt=1);
	RotationMatrix Rot(const PositionVector3D& axis_a_b);
};

#endif
