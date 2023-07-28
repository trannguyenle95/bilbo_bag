//=============================================================================
//==============================================================================

//	file:	RotationMatrix.cpp

//	author:	Fares J. Abu-Dakka
//	date:	May. 2012

//	Description: Implements a 3x3 rotation matrix
//				 It represents the transformation internally in a 3x3 Square Matrix.
//==============================================================================

#include <math.h>
#include "Frame3D.h"

namespace math{


	//-------------------------------------------------------------------------------------------
	RotationMatrix::RotationMatrix()
	/*!
	  @brief Default constructor. 
	  Constructs identity rotation matrix
	  $
	   \mathbf{R} =
	   \left[
		\begin{array}{ccc}
			1 & 0 & 0 \\
			0 & 1 & 0 \\
			0 & 0 & 1
		\end{array}
	   \right]
	  $
	*/
	{
		_R[0] = 1;		_R[1] = 0;		_R[2] = 0;
		_R[3] = 0;		_R[4] = 1;		_R[5] = 0;
		_R[6] = 0;		_R[7] = 0;		_R[8] = 1;
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix::RotationMatrix(	Real r11, Real r12, Real r13,
									Real r21, Real r22, Real r23,
									Real r31, Real r32, Real r33)
	/*!
	  @brief Constructs an initialized 3x3 rotation matrix
	  @param $ r_{11} $
	  @param $ r_{12} $
	  @param $ r_{13} $
	  @param $ r_{21} $
	  @param $ r_{22} $
	  @param $ r_{23} $
	  @param $ r_{31} $
	  @param $ r_{32} $
	  @param $ r_{33} $
	  
	  $
	   \mathbf{R} =
	   \left[
		\begin{array}{ccc}
			r_{11} & r_{12} & r_{13} \\
			r_{21} & r_{22} & r_{23} \\
			r_{31} & r_{32} & r_{33}
		\end{array}
	   \right]
	  $
	*/
	{
		_R[0] = r11;	_R[1] = r12;	_R[2] = r13;
		_R[3] = r21;	_R[4] = r22;	_R[5] = r23;
		_R[6] = r31;	_R[7] = r32;	_R[8] = r33;
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix::RotationMatrix(	const PositionVector3D& i,
									const PositionVector3D& j,
									const PositionVector3D& k)
	/*!
	  @brief Constructs an initialized 3x3 rotation matrix
	  $ ^{a}\mathbf{R}_b =
	  \left[
	   \begin{array}{ccc}
	    ^{a}\mathbf{i}_b & ^{a}\mathbf{j}_b & ^{a}\mathbf{k}_b
	   \end{array}
	  \right]
	  $
	  
	  @param i $ ^{a}\mathbf{i}_b $
	  @param j $ ^{a}\mathbf{j}_b $
	  @param k $ ^{a}\mathbf{k}_b $
	*/
	{
		_R[0] = i(1);	_R[1] = j(1);	_R[2] = k(1);
		_R[3] = i(2);	_R[4] = j(2);	_R[5] = k(2);
		_R[6] = i(3);	_R[7] = j(3);	_R[8] = k(3);
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix::RotationMatrix(Real *x)
	/*!
	  @brief Constructs an initialized 3x3 rotation matrix from a given array
	  @param 9 element array
	*/
	{
		_R[0] = x[0];	_R[1] = x[1];	_R[2] = x[2];
		_R[3] = x[3];	_R[4] = x[4];	_R[5] = x[5];
		_R[6] = x[6];	_R[7] = x[7];	_R[8] = x[8];
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix::RotationMatrix(const RotationMatrix &m)
	/*!
	  @brief Copy constructor.
	*/
	{
		_R[0] = m(1,1);	_R[1] = m(1,2);	_R[2] = m(1,3);
		_R[3] = m(2,1);	_R[4] = m(2,2);	_R[5] = m(2,3);
		_R[6] = m(3,1);	_R[7] = m(3,2);	_R[8] = m(3,3);
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix & RotationMatrix::operator =(const RotationMatrix & m)
	/*!
	  @brief Assignment operator.
	*/
	{
		_R[0] = m(1,1);	_R[1] = m(1,2);	_R[2] = m(1,3);
		_R[3] = m(2,1);	_R[4] = m(2,2);	_R[5] = m(2,3);
		_R[6] = m(3,1);	_R[7] = m(3,2);	_R[8] = m(3,3);
		return (*this);
	}
	//-------------------------------------------------------------------------------------------
	bool RotationMatrix::operator ==(const RotationMatrix & m)const
	/*!
	  @brief Comparison operator.
	  The comparison operator makes a element wise comparison with the precision of T.
	  Returns true only if all elements are equal.
	  @param m Rotation to compare with
	  @return True if equal.
	*/
	{
		for ( int i=0; i<9; i++ )
			if ( _R[i] != m._R[i] )
				return false;
		return true;
	}
	//-------------------------------------------------------------------------------------------
	const Real & RotationMatrix::operator () (int i, int j) const
	/*!
	  @brief Returns const reference to matrix element
	  @param row
	  @param column
	  @return const reference to the element
	*/
	{
		if ( i<1 || i>3 || j<1 || j>3 )
			throw Exception("RotationMatrix::operator () (int i, int j): Invalid Index");
		return _R[(i-1) * 3 + (j-1)];
	}
	//-------------------------------------------------------------------------------------------
	Real & RotationMatrix::operator () (int i, int j)
	/*!
	  @brief Returns reference to matrix element
	  @param row
	  @param column
	  @return reference to the element
	*/
	{
		if ( i<1 || i>3 || j<1 || j>3 )
			throw Exception("RotationMatrix::operator () (int i, int j): Invalid Index");
		return _R[(i-1) * 3 + (j-1)];
	}
	//-------------------------------------------------------------------------------------------
	const Real & RotationMatrix::operator [] (int i) const
	/*!
	  @brief Returns const reference to matrix element
	  @param row
	  @param column
	  @return const reference to the element
	*/
	{
		if ( i<0 || i>8)
			throw Exception("RotationMatrix::operator [] (int i): i must be from 0 to 8");
		return _R[i];
	}
	//-------------------------------------------------------------------------------------------
	Real & RotationMatrix::operator [] (int i)
	/*!
	  @brief Returns reference to matrix element
	  @param row
	  @param column
	  @return reference to the element
	*/
	{
		if ( i<0 || i>8)
			throw Exception("RotationMatrix::operator [] (int i): i must be from 0 to 8");
		return _R[i];
	}
	//-------------------------------------------------------------------------------------------
	const RotationMatrix RotationMatrix::operator *(const RotationMatrix &m)const
	/*!
	  @brief Calculates $ ^{a}\mathbf{R}_c =\, ^{a}\mathbf{R}_b \, ^{b}\mathbf{R}_c $
	  @param $ ^{a}\mathbf{R}_b $
	  @param $ ^{b}\mathbf{R}_c $
	  @return $ ^{a}\mathbf{R}_c $
	*/
	{
		RotationMatrix R;
		Real *pC = R.getData();
		for ( int i = 0; i < 3; i++ )
			for ( int j = 0; j < 3; j++ )
			{
				Real sigma = 0;
				const Real *pA = _R + i * 3;
				const Real *pB = m.getData() + j;

				for ( int k = 0; k < 3; k++ )
				{
					sigma += *pA * *pB;
					pA++; pB += 3;
				}
				*pC++ = sigma;
			}
			return R;
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D RotationMatrix::operator *(const PositionVector3D &v)const
	/*!
	  @brief Calculates $ ^{a}\mathbf{v}_c =\, ^{a}\mathbf{R}_b \, ^{b}\mathbf{v}_c $
	  @param $ ^{a}\mathbf{R}_b $
	  @param $ ^{b}\mathbf{v}_c $
	  @return $ ^{a}\mathbf{v}_c $
	*/
	{
		PositionVector3D pos;
		Real *pC = pos.getData();
		for ( int i = 0; i < 3; i++ )
			for ( int j = 0; j < 1; j++ )
			{
				Real sigma = 0;
				const Real *pA = _R + i * 3;
				const Real *pB = v.getData() + j;

				for ( int k = 0; k < 3; k++ )
				{
					sigma += *pA * *pB;
					pA++; pB += 1;
				}
				*pC++ = sigma;
			}
			return pos;
	}
	//-------------------------------------------------------------------------------------------
	const RotationMatrix RotationMatrix::operator / ( Real s) const
	/*!
	  @brief Scalar division
	*/
	{
		return RotationMatrix(	
			_R[0] / s, _R[1] / s, _R[2] / s,
			_R[3] / s, _R[4] / s, _R[5] / s,
			_R[6] / s, _R[7] / s, _R[8] / s		);
	}
	//-------------------------------------------------------------------------------------------
	const RotationMatrix RotationMatrix::operator * ( Real s) const
	/*!
	  @brief Scalar multiplication
	*/
	{
		return RotationMatrix(	
			_R[0] * s, _R[1] * s, _R[2] * s,
			_R[3] * s, _R[4] * s, _R[5] * s,
			_R[6] * s, _R[7] * s, _R[8] * s		);
	}
	//-------------------------------------------------------------------------------------------
	const RotationMatrix operator * (Real s, const RotationMatrix& v)
	/*!
	  @brief Scalar multiplication
	*/
	{
		return RotationMatrix(	
			v(1,1) * s, v(1,2) * s, v(1,3) * s,
			v(2,1) * s, v(2,2) * s, v(2,3) * s,
			v(3,1) * s, v(3,2) * s, v(3,3) * s	);
	}
	//-------------------------------------------------------------------------------------------
	const RowVector operator * (const PositionVector3D &v, const RotationMatrix& s)
	/*!
	  @brief Calculates $ \mathbf{v}^T =\, \mathbf{v}^T \, \mathbf{R} $
	  @param $ \mathbf{v}^T $
	  @param $ \mathbf{R} $
	  @return $ \mathbf{v}^T $
	*/
	{
		return RowVector(
			s(1,1)*v(1) + s(2,1)*v(2) + s(3,1)*v(3),
			s(1,2)*v(1) + s(2,2)*v(2) + s(3,2)*v(3),
			s(1,3)*v(1) + s(2,3)*v(2) + s(3,3)*v(3)	);
	}
	//-------------------------------------------------------------------------------------------
	const RotationMatrix RotationMatrix::operator -(const RotationMatrix &m)const
	/*!
	  @brief Rotations subtraction
	*/
	{
		RotationMatrix R;
		const Real *pA = (*this).getData(), *pB = m.getData();
		Real *pC = R.getData();
		for ( int i=0; i<9; i++ ){
			*pC = *pA - *pB;
			pA++;	pB++;	pC++;
		}
		return R;
	}
	//-------------------------------------------------------------------------------------------
	const RotationMatrix RotationMatrix::operator +(const RotationMatrix &m)const
	/*!
	  @brief Rotations Addition
	*/
	{
		RotationMatrix R;
		const Real *pA = (*this).getData(), *pB = m.getData();
		Real *pC = R.getData();
		for ( int i=0; i<9; i++ ){
			*pC = *pA + *pB;
			pA++;	pB++;	pC++;
		}
		return R;
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::operator *=(Real s)
	/*!
	  @brief Scalar multiplication
	*/
	{ 
		return (*this) = (*this) * s;	
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::operator /=(Real s)
	/*!
	  @brief Scalar division
	*/
	{ 
		return (*this) = (*this) / s;	
	}
	//-------------------------------------------------------------------------------------------
	const RotationMatrix RotationMatrix::operator +=(const RotationMatrix &m)
	/*!
	  @brief Rotations addition
	*/
	{
		return *this = *this + m;
	}
	//-------------------------------------------------------------------------------------------
	const RotationMatrix RotationMatrix::operator -=(const RotationMatrix &m)
	/*!
	  @brief Rotations subtraction
	*/
	{
		return *this = *this - m;
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::operator ~()const
	/*
	  @brief Calculates the transpose $ ^{b}\mathbf{R}_a =\, 
			 ^{a}\mathbf{R}_b^{T} $ of a rotation matrix

	  @param the rotation matrix $ ^{a}\mathbf{R}_b $

	  @return the matrix transpose $ ^{b}\mathbf{R}_a =\, ^{a}\mathbf{R}_b^{T} $

	  $ ^{b}\mathbf{R}_a =\, ^{a}\mathbf{R}_b^{T} =\, ^{a}\mathbf{R}_b^{-1} $	
	*/
	{
		RotationMatrix R(
			_R[0], _R[3], _R[6],
			_R[1], _R[4], _R[7],
			_R[2], _R[5], _R[8]	);
		R.normalizeR();
		return R;
	}
	//-------------------------------------------------------------------------------------------
	ostream & operator << (ostream & o, const RotationMatrix & m)
	/*
      @brief Writes rotation matrix to stream
      
      @param os output stream to use
      @param r rotation matrix to print
      @return the updated output stream
	*/
	{
		for (int i=1; i<=3; i++){
			for (int j=1; j<=3; j++)
				o << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << m(i,j) << " ";
			o << endl;
		}
		o << endl;
		return o;
	}
	//-------------------------------------------------------------------------------------------
	void RotationMatrix::print (char const* ch, int prec, int w) const
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
			for (int j=1; j<=3; j++)
				cout << std::right << std::fixed<< std::setprecision(4) << std::setw(8) << (*this)(i,j) << " ";
			cout << endl << std::setw(n+2) << " ";
		}
		cout << endl;
	}
	//-------------------------------------------------------------------------------------------
	Twist RotationMatrix::operator * (const Twist& arg) const
	// Transformation of the base to which the twist is expressed.
	// look at Frame*Twist for a transformation that also transforms
	// the velocity reference point.
	// Complexity : 18M+12A
	{
		return Twist((*this)*arg.Vel(),(*this)*arg.Rot());
	}
	//-------------------------------------------------------------------------------------------
	Wrench RotationMatrix::operator * (const Wrench& arg) const
	// Transformation of the base to which the wrench is expressed.
	// look at Frame*Twist for a transformation that also transforms
	// the force reference point.
	{
		return Wrench((*this)*arg.Force(),(*this)*arg.Torque());
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::identity()
	/*!
	  @brief Constructs a 3x3 rotation matrix set to identity

	  @return a 3x3 identity rotation matrix
	  
	  $
	   \mathbf{R} =
	   \left[
	    \begin{array}{ccc}
	      1 & 0 & 0 \\
	      0 & 1 & 0 \\
	      0 & 0 & 1
	    \end{array}
	   \right]
	  $
	*/
	{
		return RotationMatrix(1,0,0,0,1,0,0,0,1);
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix & RotationMatrix::normalizeR()
	/*!
	  @brief Normalizes the rotation matrix to satisfy SO(3).
	  Makes a normalization of the rotation matrix such that the columns
	  are normalized and orthogonal s.t. it belongs to SO(3).
	*/
	{
		Real eps00,eps01,eps02,eps11,eps12,eps22,prod0,prod1,prod2,prod;
		prod0 =  _R[0] * _R[0] + _R[3] * _R[3] + _R[6] * _R[6];
		eps00 = (ONE-prod0)/prod0;

		prod1 =  _R[1] * _R[1] + _R[4] * _R[4] + _R[7] * _R[7];
		eps11 = (ONE-prod1)/prod1;

		prod2 =  _R[2] * _R[2] + _R[5] * _R[5] + _R[8] * _R[8];
		eps22 = (ONE-prod2)/prod2;

		prod  =  _R[0] * _R[1] + _R[3] * _R[4] + _R[6] * _R[7];
		eps01 = -prod/(prod0+prod1);

		prod  =  _R[0] * _R[2] + _R[3] * _R[5] + _R[6] * _R[8];
		eps02 = -prod/(prod0+prod2);

		prod  =  _R[1] * _R[2] + _R[4] * _R[5] + _R[7] * _R[8];
		eps12 = -prod/(prod1+prod2);

		_R[0] += eps00*_R[0] + eps01*_R[1] + eps02*_R[2];
		_R[3] += eps00*_R[3] + eps01*_R[4] + eps02*_R[5];
		_R[6] += eps00*_R[6] + eps01*_R[7] + eps02*_R[8];
		_R[1] += eps01*_R[0] + eps11*_R[1] + eps12*_R[2];
		_R[4] += eps01*_R[3] + eps11*_R[4] + eps12*_R[5];
		_R[7] += eps01*_R[6] + eps11*_R[7] + eps12*_R[8];
		_R[2] += eps02*_R[0] + eps12*_R[1] + eps22*_R[2];
		_R[5] += eps02*_R[3] + eps12*_R[4] + eps22*_R[5];
		_R[8] += eps02*_R[6] + eps12*_R[7] + eps22*_R[8];
		return *this;
	}
	//-------------------------------------------------------------------------------------------
	const PositionVector3D RotationMatrix::getCol(int col) const
	/*!
	  @brief Returns a Column i from this rotation matrix object.
	*/
	{
		if(col > 3 || col < 1)
			throw Exception("error: RotationMatrix::getCol: Invalid Index");
		return PositionVector3D((*this)(1,col),(*this)(2,col),(*this)(3,col));
	}
	//-------------------------------------------------------------------------------------------
	bool RotationMatrix::equal(const RotationMatrix& rot, Real precision)const
	/*!
	  @brief Compares rotations with a given precision
	  @param \b Rot: Rotation to compare with
	  @param \b precision: The precision to use for testing
	  @return True if all elements are less than \b precision apart.
	  
	  Performs an element wise comparison. Two elements are considered equal if the difference
	  are less than \b precision.
	*/
	{
		return (EQUAL(_R[0],rot[0],precision) &&
                EQUAL(_R[1],rot[1],precision) &&
                EQUAL(_R[2],rot[2],precision) &&
                EQUAL(_R[3],rot[3],precision) &&
                EQUAL(_R[4],rot[4],precision) &&
                EQUAL(_R[5],rot[5],precision) &&
                EQUAL(_R[6],rot[6],precision) &&
                EQUAL(_R[7],rot[7],precision) &&
                EQUAL(_R[8],rot[8],precision)    );
	}
	//-------------------------------------------------------------------------------------------
	bool equal(const RotationMatrix& rot1, const RotationMatrix& rot2, Real precision)
	/*!
	  @brief Compares rotations with a given precision
	  @param \b Rot: Rotation to compare with
	  @param \b precision: The precision to use for testing
	  @return True if all elements are less than \b precision apart.
	  
	  Performs an element wise comparison. Two elements are considered equal if the difference
	  are less than \b precision.
	*/
	{
		return rot1.equal(rot2,precision);
	}
	//-------------------------------------------------------------------------------------------
	void RotationMatrix::zeros()
	/*!
	  @brief Sets all elements to zero.
	*/
	{
		// Replaces all values in the Matrix with zero.
		for ( int i=0; i<9; i++ )
			_R[i] = 0;
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::skew(const PositionVector3D& v)
	/*!
	  @brief Creates a skew symmetric matrix from a PositionVector3D. Also
	  known as the cross product matrix of v.
	  
	  @param v vector to create Skew matrix from
	*/
	{
		return RotationMatrix (0, -v(3), v(2), v(3), 0, -v(1), -v(2), v(1), 0);
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::Inverse()const
	/*
	  @brief Calculates the inverse $ ^{b}\mathbf{R}_a =\, ^{a}\mathbf{R}_b^{-1} $ of a rotation matrix

	  @param the rotation matrix $ ^{a}\mathbf{R}_b $

	  @return the matrix inverse $ ^{b}\mathbf{R}_a =\, ^{a}\mathbf{R}_b^{-1} $

	  $ ^{b}\mathbf{R}_a =\, ^{a}\mathbf{R}_b^{-1} =\, ^{a}\mathbf{R}_b^T $	
	*/
	{
		return ~(*this);
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D RotationMatrix::Inverse(const PositionVector3D& v)const
	/*
	  @brief The same as R.Inverse()*v but more efficient.
	*/
	{
		return PositionVector3D(
		_R[0]*v[0] + _R[3]*v[1] + _R[6]*v[2],
		_R[1]*v[0] + _R[4]*v[1] + _R[7]*v[2],
		_R[2]*v[0] + _R[5]*v[1] + _R[8]*v[2] );
	}
	//-------------------------------------------------------------------------------------------
	void RotationMatrix::SetInverse()
	/*
	  @brief Sets the value of *this to its inverse.
	*/
	{
		double tmp;
		tmp = _R[1];	_R[1]=_R[3];	_R[3]=tmp;
		tmp = _R[2];	_R[2]=_R[6];	_R[6]=tmp;
		tmp = _R[5];	_R[5]=_R[7];	_R[7]=tmp;
	}
	//-------------------------------------------------------------------------------------------
	Wrench RotationMatrix::Inverse(const Wrench& arg) const
	/*
	  @brief The same as R.Inverse()*arg but more efficient.
	*/
	{
		return Wrench(Inverse(arg.Force()),Inverse(arg.Torque()));
	}
	//-------------------------------------------------------------------------------------------
	Twist RotationMatrix::Inverse(const Twist& arg) const
	/*
	  @brief The same as R.Inverse()*arg but more efficient.
	*/
	{
		return Twist(Inverse(arg.Vel()),Inverse(arg.Rot()));
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D RotationMatrix::r2rpy(RotationMatrix R)
	/*!
	  @brief Converts a rotation matrix to roll/pitch/yaw angles. RPY = r2rpy(R)
	  @param the rotation matrix $ \mathbf{R} $
	  @Returns $ \mathbf{RPY} $ a vector of roll/pitch/yaw angles corresponding to the rotational
	  		 part of the Rotation Matrix R.
	*/
	{
		PositionVector3D rpy;
		// 	Real eps = ::pow(2.0, -22.0);
		if (::abs(R(1,1)) < EPSilon && ::abs(R(2,1)) < EPSilon){
			// singularity
			rpy(1) = 0;
			rpy(2) = atan2(-R(3,1), R(1,1));
			rpy(3) = atan2(-R(2,3), R(2,2));
		}else{
			rpy(1) = atan2(R(2,1), R(1,1));
			Real sp = sin(rpy(1));
			Real cp = cos(rpy(1));
			rpy(2) = atan2(-R(3,1), cp * R(1,1) + sp * R(2,1));
			rpy(3) = atan2(sp * R(1,3) - cp * R(2,3), cp*R(2,2) - sp*R(1,2));
		}
		return rpy;
	}
	//***************************************************************************************************
	RotationMatrix RotationMatrix::rpy2r(Real roll, Real pitch, Real yaw)
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
	void RotationMatrix::getRPY(Real& roll, Real& pitch, Real& yaw) const
	/*!
	  @brief  Calculates rpy for this RotationMatrix object
      @param  -PI <= roll <= PI
      @param  -PI <= Yaw  <= PI
      @param  -PI/2 <= pitch <= PI/2
	  @Return a homogeneous transformation for the specified roll/pitch/yaw angles.
				These correspond to rotations about the Z, Y, X axes respectively.
	*/
	{
		pitch = atan2(-_R[6], sqrt( SQR(_R[0]) +SQR(_R[3]) )  );
		if ( fabs(pitch) > (M_PI_2 - EPSilon) ) {
			yaw = atan2(	-_R[1], _R[4]);
			roll  = ZERO ;
		} else {
			roll  = atan2(_R[7], _R[8]);
			yaw   = atan2(_R[3], _R[0]);
		}
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::roty(Real t)
	/*!
	  @brief Rotates about Y axis. R = roty(theta)
	  @param angel t
	  @Returns $ \mathbf{R} $ a rotation matrix representing a rotation of theta about the Y axis.
	*/
	{
		Real ct = cos(t), st = sin(t);
		return RotationMatrix(ct,0,st,0,1,0,-st,0,ct);
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::rotx(Real t)
	/*!
	  @brief Rotates about X axis. R = roty(theta)
	  @param angel t
	  @Returns $ \mathbf{R} $ a rotation matrix representing a rotation of theta about the X axis.
	*/
	{
		Real ct = cos(t), st = sin(t);
		return RotationMatrix(1,0,0,0,ct,-st,0,st,ct);
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::rotz(Real t)
	/*!
	  @brief Rotates about Z axis. R = roty(theta)
	  @param angel t
	  @Returns $ \mathbf{R} $ a rotation matrix representing a rotation of theta about the Z axis.
	*/
	{
		Real ct = cos(t) , st = sin(t);
		return RotationMatrix(ct,-st,0,st,ct,0,0,0,1);
	}
	//-------------------------------------------------------------------------------------------
	void RotationMatrix::doRotX(Real t)
	/*!
	  @brief Rotates this object about X
			 *this = *this * rotx(t)
	*/
	{
		Real ct = cos(t), st = sin(t), x1, x2, x3;
		x1  = ct* _R[1] + st* _R[2];
		x2  = ct* _R[4] + st* _R[5];
		x3  = ct* _R[7] + st* _R[8];
		_R[2] = -st* _R[1] + ct* _R[2];
		_R[5] = -st* _R[4] + ct* _R[5];
		_R[8] = -st* _R[7] + ct* _R[8];
		_R[1] = x1;
		_R[4] = x2;
		_R[7] = x3;
	}
	//-------------------------------------------------------------------------------------------
	void RotationMatrix::doRotY(Real t)
	/*!
	  @brief Rotates this object about Y
			 *this = *this * roty(t)
	*/
	{
		Real ct = cos(t), st = sin(t), x1, x2, x3;
		x1  = ct* _R[0] - st* _R[2];
		x2  = ct* _R[3] - st* _R[5];
		x3  = ct* _R[6] - st* _R[8];
		_R[2] = st* _R[0] + ct* _R[2];
		_R[5] = st* _R[3] + ct* _R[5];
		_R[8] = st* _R[6] + ct* _R[8];
		_R[0] = x1;
		_R[3] = x2;
		_R[6] = x3;
	}
	//-------------------------------------------------------------------------------------------
	void RotationMatrix::doRotZ(Real t)
	/*!
	  @brief Rotates this object about Z
			 *this = *this * rotz(t)
	*/
	{
		Real ct = cos(t), st = sin(t), x1, x2, x3;
		x1  = ct* _R[0] + st* _R[1];
		x2  = ct* _R[3] + st* _R[4];
		x3  = ct* _R[6] + st* _R[7];
		_R[1] = -st* _R[0] + ct* _R[1];
		_R[4] = -st* _R[3] + ct* _R[4];
		_R[7] = -st* _R[6] + ct* _R[7];
		_R[0] = x1;
		_R[3] = x2;
		_R[6] = x3;
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::Rot(const PositionVector3D& ra,Real angle){
	/*!
	  @brief Rotates along an arbitrary axes
			 V.(V.tr) + st*[V x] + ct*(I-V.(V.tr))
	  @Return identity rotation matrix in the case that the norm of rv is to small to be used.
	*/
		PositionVector3D rv = ra;
		rv.Normalize();
		return RotationMatrix::Rot2(rv,angle);
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::Rot2(const PositionVector3D& rotvec, Real angle) {
		// rotvec should be normalized !
		// The formula is
		// V.(V.tr) + st*[V x] + ct*(I-V.(V.tr))
		// can be found by multiplying it with an arbitrary vector p
		// and noting that this vector is rotated.
		Real ct = cos(angle);
		Real st = sin(angle);
		Real vt = 1 - ct;
		Real m_vt_0 = vt*rotvec(1);
		Real m_vt_1 = vt*rotvec(2);
		Real m_vt_2 = vt*rotvec(3);
		Real m_st_0 = rotvec(1)*st;
		Real m_st_1 = rotvec(2)*st;
		Real m_st_2 = rotvec(3)*st;
		Real m_vt_0_1 = m_vt_0*rotvec(2);
		Real m_vt_0_2 = m_vt_0*rotvec(3);
		Real m_vt_1_2 = m_vt_1*rotvec(3);
		return RotationMatrix(
			ct + m_vt_0*rotvec(1),
			-m_st_2 + m_vt_0_1,
			m_st_1 + m_vt_0_2,
			m_st_2 + m_vt_0_1,
			ct + m_vt_1*rotvec(2),
			-m_st_0 + m_vt_1_2,
			-m_st_1 + m_vt_0_2,
			m_st_0 + m_vt_1_2,
			ct + m_vt_2*rotvec(3)
			);
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::eulerzyz(Real ALPHA, Real BETA, Real GAMA)
	/*!
	  @brief Converts Euler angels into a rotation matrix
			 by rotating around Z with ALPHA then around Y with BETA then around Z with GAMA
	  @param angels ALPHA, BETA, GAMA
	  @Returns $ \mathbf{R} $ a rotation matrix.
	*/
	{
		Real ca = cos(ALPHA), sa = sin(ALPHA);
		Real cb = cos(BETA ), sb = sin(BETA );
		Real cg = cos(GAMA ), sg = sin(GAMA );
		return RotationMatrix(
			 ca * cb * cg - sa * sg, -ca * cb * sg - sa * cg,  ca * sb,
			 sa * cb * cg + ca * sg, -sa * cb * sg + ca * cg,  sa * sb,
			-sb * cg			   ,  sb * sg				,  cb		);
	}
	//-------------------------------------------------------------------------------------------
	void RotationMatrix::getEulerZYZ(Real& ALPHA, Real& BETA, Real& GAMA)const
	/*!
	  @brief Gives back the EulerZYZ convention description of the rotation matrix
			 rotating around Z with ALPHA then around Y with BETA then around Z with GAMA
	  @param (-PI <  ALPHA  <= PI)
	  @param (0   <= BETA   <= PI),
	  @param (-PI <  GAMA   <= PI)

	  @Returns $ \mathbf{R} $ a rotation matrix.
	*/
	{
		if (fabs(_R[8]) > 1-EPSilon  ) {
			GAMA=ZERO;
			if (_R[8]>0) {
				BETA = ZERO;
				ALPHA= atan2(_R[3],_R[0]);
			} else {
				BETA = PI;
				ALPHA= atan2(-_R[3],-_R[0]);
			}
		} else {
			ALPHA=atan2(_R[5], _R[2]);
			BETA=atan2(sqrt( SQR(_R[6]) +SQR(_R[7]) ),_R[8]);
			GAMA=atan2(_R[7], -_R[6]);
		}
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::eulerzyx(Real ALPHA, Real BETA, Real GAMA)
	/*!
	  @brief Converts Euler angels into a rotation matrix
			 by rotating around Z with ALPHA then around Y with BETA then around X with GAMA
	  @param angels ALPHA, BETA, GAMA
	  @Returns $ \mathbf{R} $ a rotation matrix.
	*/
	{
		return rpy2r(GAMA,BETA,ALPHA);
	}
	//-------------------------------------------------------------------------------------------
	void RotationMatrix::getEulerZYX(Real& ALPHA, Real& BETA, Real& GAMA)const
	/*!
	  @brief Gives back the EulerZYX convention description of the rotation matrix
			 rotating around Z with ALPHA then around Y with BETA then around X with GAMA
	  @param (-PI   <  ALPHA  <= PI  )
	  @param (-PI/2 <= BETA   <= PI/2),
	  @param (-PI   <  GAMA   <= PI  )

	  @Returns $ \mathbf{R} $ a rotation matrix.
	*/
	{
		getRPY(GAMA,BETA,ALPHA);;
	}
	//-------------------------------------------------------------------------------------------
	ColumnVector RotationMatrix::irotk(const RotationMatrix & R)
	/*
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

		return k;
	}
	//-------------------------------------------------------------------------------------------
	Matrix RotationMatrix::subMatrix() const
	/*
	  @brief Returns the rotation matrix as a matrix
	*/
	{
		Matrix mat(3, 3);
		for (int i=1; i<=3; i++)
			for (int j=1; j<=3; j++)
				mat(i,j) = (*this)(i,j);
		return mat;
	}
	//-------------------------------------------------------------------------------------------
	double RotationMatrix::GetRotAngle(PositionVector3D& axis, double eps) const
	/*
	  @brief Returns the rotation angle around the equiv. axis
	  @param axis the rotation axis is returned in this variable
	  @param eps :  in the case of angle == 0 : rot axis is undefined and choosen
			 to be the Z-axis in the case of angle == PI : 2 solutions, positive Z-component
			 of the axis is choosen.
	  @Returns the rotation angle (between [0..PI] )
	*/
	{
		Real ca = (_R[0] + _R[4] + _R[8] - ONE) / TWO;
		Real t = eps*eps / TWO;
		if (ca>1 - t) {
			// undefined choose the Z-axis, and angle 0
			axis = PositionVector3D(ZERO, ZERO, ONE);
			return 0;
		}
		if (ca < -1 + t) {
			// The case of angles consisting of multiples of M_PI:
			// two solutions, choose a positive Z-component of the axis
			Real x = sqrt((_R[0] + ONE) / 2);
			Real y = sqrt((_R[4] + ONE) / 2);
			Real z = sqrt((_R[8] + ONE) / 2);
			if (_R[2] < 0) x = -x;
			if (_R[7] < 0) y = -y;
			if (x*y*_R[1] < 0) x = -x;  // this last line can be necessary when z is 0
			// z always >= 0 
			// if z equal to zero 
			axis = PositionVector3D(x, y, z);
			return PI;
		}
		Real angle;
		Real mod_axis;
		Real axisx, axisy, axisz;
		axisx = _R[7] - _R[5];
		axisy = _R[2] - _R[6];
		axisz = _R[3] - _R[1];
		mod_axis = sqrt(axisx*axisx + axisy*axisy + axisz*axisz);
		axis = PositionVector3D(axisx / mod_axis, axisy / mod_axis,	axisz / mod_axis);
		angle = atan2(mod_axis / 2, ca);
		return angle;
	}
	//-------------------------------------------------------------------------------------------
	PositionVector3D RotationMatrix::GetRot() const
	/*!
	  @brief Returns a vector with the direction of the equiv. axis and its norm is angle
	*/
	{
		PositionVector3D axis;
		double angle;
		angle = RotationMatrix::GetRotAngle(axis);
		return axis * angle;
	}
	//-------------------------------------------------------------------------------------------
	RotationMatrix RotationMatrix::uQuaternion(double x,double y,double z, double w)
    {
        double x2, y2, z2, w2;
        double norm = sqrt(x*x+y*y+z*z+w*w);
        x /= norm;y /= norm;z /= norm;w /= norm;
        x2 = x*x;  y2 = y*y; z2 = z*z;  w2 = w*w;
        return RotationMatrix(w2+x2-y2-z2, 2*x*y-2*w*z, 2*x*z+2*w*y,
                        2*x*y+2*w*z, w2-x2+y2-z2, 2*y*z-2*w*x,
                        2*x*z-2*w*y, 2*y*z+2*w*x, w2-x2-y2+z2);
    }
	Matrix RotationMatrix::mat()const{
		Matrix m(3);
		m(1,1) = _R[0];	m(1,2) = _R[1];	m(1,3) = _R[2];
		m(2,1) = _R[3];	m(2,2) = _R[4];	m(2,3) = _R[5];
		m(3,1) = _R[6];	m(3,2) = _R[7];	m(3,3) = _R[8];
		return m;
	}
	//-------------------------------------------------------------------------------------------
	/*std::istream& operator >> (std::istream& is, RotationMatrix& r)
	{
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
	}*/
	void random(RotationMatrix& r) {
		double alfa, beta, gamma;
		math::random(alfa);
		math::random(beta);
		math::random(gamma);
		r = RotationMatrix::eulerzyx(alfa,beta,gamma);
	}
};
