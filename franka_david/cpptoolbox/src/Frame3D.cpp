/*
 * Frame3D.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: fares
 */

#include <cmath>
#include "Frame3D.h"
namespace math{
	void Twist2CV(const Twist& t, ColumnVector& e) {
		if (e.getnRows() != 6)
			e = ColumnVector(6);
		e(1) = t(1);
		e(2) = t(2);
		e(3) = t(3);
		e(4) = t(4);
		e(5) = t(5);
		e(6) = t(6);
	}
	ColumnVector Twist2CV(const Twist& t) {
		ColumnVector e(6);
		e(1) = t(1);
		e(2) = t(2);
		e(3) = t(3);
		e(4) = t(4);
		e(5) = t(5);
		e(6) = t(6);
		return e;
	}
	RotationMatrix Rot(const PositionVector3D& axis_a_b) {
	    // The formula is
	    // V.(V.tr) + st*[V x] + ct*(I-V.(V.tr))
	    // can be found by multiplying it with an arbitrary vector p
	    // and noting that this vector is rotated.
		PositionVector3D rotvec = axis_a_b;
		double angle = rotvec.Normalize(1E-10);
	    double ct = ::cos(angle);
	    double st = ::sin(angle);
	    double vt = 1-ct;
	    return RotationMatrix(
	        ct            +  vt*rotvec(0)*rotvec(0),
	        -rotvec(2)*st +  vt*rotvec(0)*rotvec(1),
	        rotvec(1)*st  +  vt*rotvec(0)*rotvec(2),
	        rotvec(2)*st  +  vt*rotvec(1)*rotvec(0),
	        ct            +  vt*rotvec(1)*rotvec(1),
	        -rotvec(0)*st +  vt*rotvec(1)*rotvec(2),
	        -rotvec(1)*st +  vt*rotvec(2)*rotvec(0),
	        rotvec(0)*st  +  vt*rotvec(2)*rotvec(1),
	        ct            +  vt*rotvec(2)*rotvec(2)
	        );
	}
	PositionVector3D diff(const PositionVector3D& a, const PositionVector3D& b, double dt)
	/*!
	 * Determines the difference between 2 vectors
	 */
	{
		return (b - a) / dt;
	}

	PositionVector3D diff(const RotationMatrix& r1, const RotationMatrix& r2, double dt)
	/*!
	 * \brief  Determines the scaled rotation axis necessary to rotate from
	 * \warning - The result is not a rotational vector, i.e. it is not a mathematical vector.
	 *          (no communitative addition).
	 *
	 * \warning - When used in the context of numerical differentiation, with the frames b1 and b2 very
	 *           close to each other, the semantics correspond to the twist, scaled by the time.
	 *
	 * \warning - For angles equal to \f$ \pi \f$, The negative of the
	 *          return value is equally valid.
	 */
	{
		RotationMatrix R(r1.Inverse()*r2);
		return r1 * R.GetRot() / dt;
	}

	Twist diff(const Transform3D& F_a_b1, const Transform3D& F_a_b2, double dt)
	/*!
	 * \brief  determines the rotation axis necessary to rotate the frame b1 to the same orientation as frame b2 and the vector
	 * necessary to translate the origin of b1 to the origin of b2, and stores the result in a Twist datastructure.
	 * \param F_a_b1 frame b1 expressed with respect to some frame a.
	 * \param F_a_b2 frame b2 expressed with respect to some frame a.
	 * \warning The result is not a Twist!
	 * see diff() for Rotation and Vector arguments for further detail on the semantics.
	 */
	{
		/*Twist m;
		PositionVector3D p = diff(F_a_b1.P(), F_a_b2.P(), dt);
		PositionVector3D r = diff(F_a_b1.R(), F_a_b2.R(), dt);
		m(1) = p(1);		m(2) = p(2);		m(3) = p(3);
		m(4) = r(1);		m(5) = r(2);		m(6) = r(3);*/
		return Twist(diff(F_a_b1.P(), F_a_b2.P(), dt), diff(F_a_b1.R(), F_a_b2.R(), dt));
	}
	PositionVector3D addDelta(const PositionVector3D& a,const PositionVector3D& da,double dt)
	/**
	 * returns the rotation matrix resulting from the rotation of frame a by the axis and angle
	 * specified with da_w.
	 *
	 * see also the corresponding diff() routine.
	 *
	 * \param R_w_a Rotation matrix of frame a expressed to some frame w.
	 * \param da_w  axis and angle of the rotation expressed to some frame w.
	 * \returns the rotation matrix resulting from the rotation of frame a by the axis and angle
	 *          specified with da.   The resulting rotation matrix is expressed with respect to
	 *          frame w.
	 */
	{
		return a+da*dt;
	}
	RotationMatrix addDelta(const RotationMatrix& a,const PositionVector3D& da,double dt)
	/**
	 * returns the frame resulting from the rotation of frame a by the axis and angle
	 * specified in da_w and the translation of the origin (also specified in da_w).
	 *
	 * see also the corresponding diff() routine.
	 * \param R_w_a Rotation matrix of frame a expressed to some frame w.
	 * \param da_w  axis and angle of the rotation (da_w.rot), together with a displacement vector for the origin (da_w.vel),  expressed to some frame w.
	 * \returns the frame resulting from the rotation of frame a by the axis and angle
	 *          specified with da.rot, and the translation of the origin da_w.vel .  The resulting frame is expressed with respect to frame w.
	 */
	{
		return a*Rot(a.Inverse(da)*dt);
	}
	Transform3D addDelta(const Transform3D& a,const Twist& da,double dt)
	/**
	 * \brief adds the twist da to the twist a.
	 * see also the corresponding diff() routine.
	 * \param a a twist wrt some frame
	 * \param da a twist difference wrt some frame
	 * \returns The twist (a+da) wrt the corresponding frame.
	 */
	{
		return Transform3D(
					addDelta(a.P(),da.Vel(),dt),
					addDelta(a.R(),da.Rot(),dt)
				  );
	}
	Twist addDelta(const Twist& a,const Twist&da,double dt)
	/**
	 * \brief adds the wrench da to the wrench w.
	 * see also the corresponding diff() routine.
	 * see also the corresponding diff() routine.
	 * \param a a wrench wrt some frame
	 * \param da a wrench difference wrt some frame
	 * \returns the wrench (a+da) wrt the corresponding frame.
	 */
	{
		return Twist(addDelta(a.Vel(),da.Vel(),dt),addDelta(a.Rot(),da.Rot(),dt));
	}
	Wrench addDelta(const Wrench& a,const Wrench&da,double dt){
		return Wrench(addDelta(a.Force(),da.Force(),dt),addDelta(a.Torque(),da.Torque(),dt));
	}

};
