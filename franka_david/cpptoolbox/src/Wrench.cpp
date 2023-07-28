/**============================================================================
//==============================================================================

//	file:	Wrench.cpp

//	author:	Fares J. Abu-Dakka
//	date:	2015

//	Description: Adapted from KDL library
//==============================================================================*/

#include <math.h>
#include "Frame3D.h"

namespace math{
	Wrench& Wrench::operator-=(const Wrench& arg){
		_torque -= arg.Torque();
		_force  -= arg.Force();
		return *this;
	}
	Wrench& Wrench::operator+=(const Wrench& arg){
		_torque += arg.Torque();
		_force  += arg.Force();
		return *this;
	}
	double& Wrench::operator()(int i){
		if (i>0 && i<=3) 
			return _force(i);
		else if(i>3 && i<=6)
			return _torque(i-3);
		else
			throw Exception("index i must be from 1 to 6 in Twist");
	}
	double Wrench::operator()(int i)const{
		if (i>0 && i<=3) 
			return _force(i);
		else if(i>3 && i<=6)
			return _torque(i-3);
		else
			throw Exception("index i must be from 1 to 6 in Twist");
	}
	Wrench operator*(const Wrench& lhs,double rhs){
		return Wrench(lhs.Force()*rhs,lhs.Torque()*rhs);
	}
	Wrench operator*(double lhs,const Wrench& rhs){
		return Wrench(lhs*rhs.Force(),lhs*rhs.Torque());
	}
	Wrench operator/(const Wrench& lhs,double rhs){
		return Wrench(lhs.Force()/rhs,lhs.Torque()/rhs);
	}
	Wrench operator+(const Wrench& lhs,const Wrench& rhs){
		return Wrench(lhs.Force()+rhs.Force(),lhs.Torque()+rhs.Torque());
	}
	Wrench operator-(const Wrench& lhs,const Wrench& rhs){
		return Wrench(lhs.Force()-rhs.Force(),lhs.Torque()-rhs.Torque());
	}
	Wrench operator-(const Wrench& arg) {
	return Wrench(-arg.Force(),-arg.Torque());
	}
	void SetToZero(Wrench& v) {
		v.Force()  = PositionVector3D();
		v.Torque() = PositionVector3D();
	}
	Wrench Wrench::Zero(){
		return Wrench(PositionVector3D::zero(),PositionVector3D::zero());
	}
	void Wrench::ReverseSign(){   
		_force.reverseSign();
		_torque.reverseSign();
	}
	Wrench Wrench::RefPoint(const PositionVector3D& v_base_AB) const
	// Changes the reference point of the Wrench.
	// The vector v_base_AB is expressed in the same base as the twist
	// The vector v_base_AB is a vector from the old point to
	// the new point.
	// Complexity : 6M+6A
	{
		return Wrench(this->Force(), this->Torque()+this->Force()*v_base_AB);
	}
	bool Equal(const Wrench& a,const Wrench& b,double eps) {
		return (equal(a.Force(),b.Force(),eps)&&
				equal(a.Torque(),b.Torque(),eps)  );
	}
	bool operator==(const Wrench& a,const Wrench& b ) {
		return (a.Force()==b.Force() &&
				a.Torque()==b.Torque() );
	}
	bool operator!=(const Wrench& a,const Wrench& b) {
		return !operator==(a,b);
	}
	void Wrench::print(char const* ch, int prec, int w)const{
		int n = 0;
		if (ch){
			std::string s(ch);
			std::cout << ch << "= ";
			n = s.size();
		}
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _force(1) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _force(2) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _force(3) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _torque(1) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _torque(2) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _torque(3) << " ";
		cout << endl;
	}

	/*std::istream& operator >> (std::istream& is, Wrench& v)
	{
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
	}*/
	void random(Wrench& a) {
	   	random(a.Torque());
		random(a.Force());
	}
};
