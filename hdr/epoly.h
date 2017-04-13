#ifndef EPOLY_H
#define EPOLY_H

#include "poly.h"
#include <iostream>
#include <utility>

using namespace std;

class EPoly : public Poly{
	//gamma * s^beta * exp(alpha * s)
public:
	EPoly();

	EPoly(double b, double a);
	EPoly(const Poly &p, double a);
	EPoly(const EPoly & ep);
	void MultByEPoly(const EPoly &ep2);
	void MultByConstant(double c);
	double Eval(double t);
	EPoly Integral();
	pair<EPoly, double> Integral_0_t();
	friend std::ostream& operator << (std::ostream &o, const EPoly &epoly) {
		
		o << static_cast<const Poly&>(epoly);
		o << "exp component : " << epoly.alpha_ << "\n";;
		return o;
	}
	bool operator == (const EPoly ep2) {
		Poly pp2(static_cast<const Poly>(ep2));
	//	if (this->IsZeroPoly() == (static_cast<const Poly>(ep2)).IsZeroPoly()) return true;
if (this->IsZeroPoly() == pp2.IsZeroPoly()) return true;
Poly cur_poly(static_cast<const Poly>(*this));
return alpha_ == ep2.alpha_ && cur_poly == pp2;
//		return alpha_ == ep2.alpha_ && static_cast<const Poly>(*this) == static_cast<const Poly>(ep2);
	}

	EPoly operator * (const EPoly &p2);
	//bool IsZero();
protected:
	double EIntegral();
	double EDerivitive();
	
	//Poly poly_;
	double alpha_;
	bool IsAlphaZero();
	double sum_eps_;
	EPoly PolyIntegrate();
	EPoly FullIntegrate();
};


#endif
