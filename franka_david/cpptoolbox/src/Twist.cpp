/**============================================================================
//==============================================================================

//	file:	Twist.cpp

//	author:	Fares J. Abu-Dakka
//	date:	2015

//	Description: Adapted from KDL library
//==============================================================================*/

#include <cmath>
#include "Frame3D.h"

namespace math{
	Twist& Twist::operator-=(const Twist& arg){
		_vel -= arg.Vel();
		_rot -= arg.Rot();
		return *this;
	}
	Twist& Twist::operator+=(const Twist& arg){
		_vel += arg.Vel();
		_rot += arg.Rot();
		return *this;
	}
	double& Twist::operator()(int i){
		if (i>0 && i<=3) 
			return _vel(i);
		else if(i>3 && i<=6)
			return _rot(i-3);
		else
			throw Exception("index i must be from 1 to 6 in Twist");
	}
	double Twist::operator()(int i)const{
		if (i>0 && i<=3) 
			return _vel(i);
		else if(i>3 && i<=6)
			return _rot(i-3);
		else
			throw Exception("index i must be from 1 to 6 in Twist");
	}
	Twist operator*(const Twist& lhs,double rhs){
		return Twist(lhs.Vel()*rhs,lhs.Rot()*rhs);
	}
	Twist operator*(double lhs,const Twist& rhs){
		return Twist(lhs*rhs.Vel(),lhs*rhs.Rot());
	}
	Twist operator/(const Twist& lhs,double rhs){
		return Twist(lhs.Vel()/rhs,lhs.Rot()/rhs);
	}
	Twist operator+(const Twist& lhs,const Twist& rhs){
		return Twist(lhs.Vel()+rhs.Vel(),lhs.Rot()+rhs.Rot());
	}
	Twist operator-(const Twist& lhs,const Twist& rhs){
		return Twist(lhs.Vel()-rhs.Vel(),lhs.Rot()-rhs.Rot());
	}
	Twist operator-(const Twist& arg){
		return Twist(-arg.Vel(),-arg.Rot());
	}
	double dot(const Twist& lhs,const Wrench& rhs) {
		return dot(lhs.Vel(),rhs.Force())+dot(lhs.Rot(),rhs.Torque());
	}
	double dot(const Wrench& rhs,const Twist& lhs) {
		return dot(lhs.Vel(),rhs.Force())+dot(lhs.Rot(),rhs.Torque());
	}
	void SetToZero(Twist& v) {
		v.Vel() = PositionVector3D();
		v.Rot() = PositionVector3D();
	}
	Twist operator*(const Twist& lhs,const Twist& rhs){
		return Twist(lhs.Rot()*rhs.Vel()+lhs.Vel()*rhs.Rot(),lhs.Rot()*rhs.Rot());
	}
	Wrench operator*(const Twist& lhs,const Wrench& rhs){
		return Wrench(lhs.Rot()*rhs.Force(),lhs.Rot()*rhs.Torque()+lhs.Vel()*rhs.Force());
	}
	Twist Twist::Zero(){
		return Twist(PositionVector3D::zero(),PositionVector3D::zero());
	}
	void Twist::ReverseSign(){
		_vel.reverseSign();
		_rot.reverseSign();
	}
	Twist Twist::RefPoint(const PositionVector3D& v_base_AB) const
	// Changes the reference point of the twist.
	// The vector v_base_AB is expressed in the same base as the twist
	// The vector v_base_AB is a vector from the old point to
	// the new point.
	// Complexity : 6M+6A
	{
		return Twist(this->Vel()+this->Rot()*v_base_AB,this->Rot());
	}
	bool Equal(const Twist& a,const Twist& b,double eps) {
		return (equal(a.Rot(),b.Rot(),eps)&&
				equal(a.Vel(),b.Vel(),eps)  );
	}
	bool operator==(const Twist& a,const Twist& b) {
		return (a.Rot()==b.Rot() &&
				a.Vel()==b.Vel()  );
	}
	bool operator!=(const Twist& a,const Twist& b) {
		return !operator==(a,b);
	}
	void Twist::print(char const* ch, int prec, int w)const{
		int n = 0;
		if (ch){
			std::string s(ch);
			std::cout << ch << "= ";
			n = s.size();
		}
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _vel(1) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _vel(2) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _vel(3) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _rot(1) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _rot(2) << " ";
		cout << std::right << std::fixed << std::setprecision(4) << std::setw(8) << _rot(3) << " ";
		cout << endl;
	}
	ostream & operator << (ostream & o, const Twist & m){
		int n = 0;

		o << std::right << std::fixed << std::setprecision(4) << std::setw(8) << m.Vel(1) << " ";
		o << std::right << std::fixed << std::setprecision(4) << std::setw(8) << m.Vel(2) << " ";
		o << std::right << std::fixed << std::setprecision(4) << std::setw(8) << m.Vel(3) << " ";
		o << std::right << std::fixed << std::setprecision(4) << std::setw(8) << m.Rot(1) << " ";
		o << std::right << std::fixed << std::setprecision(4) << std::setw(8) << m.Rot(2) << " ";
		o << std::right << std::fixed << std::setprecision(4) << std::setw(8) << m.Rot(3) << " ";
		o << endl;
		return o;
	}
	/*std::istream& operator >> (std::istream& is, Twist& v)
	{
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
	}*/
	void random(Twist& a) {
		random(a.Rot());
		random(a.Vel());
	}
};
