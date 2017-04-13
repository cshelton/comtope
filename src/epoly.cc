#include "epoly.h"
#include "poly.h"

#include <vector>
#include <iostream>
#include <utility>
#include <cmath>

using namespace std;

// create 0 * exp(s * 0) = 0;
EPoly::EPoly() : Poly() {
	//poly_ = new Poly();
	//Poly::Poly();
	alpha_ = 0;
}

EPoly::EPoly(double b, double a) : Poly() {
	co_[0] = b;
	alpha_ = a;
}

EPoly::EPoly(const Poly &pp, double a) : Poly(pp){
	//EPoly();
	//*poly_ = pp;
	//poly_ = pp;
	alpha_ = a;
	//Poly::operator = (pp);
	//co_ = pp.GetCoeffs();
	//co_.push_back(100);
}

double EPoly::EDerivitive() {
	return alpha_;
}

double EPoly::Eval(double t) {
	return Poly::Eval(t) * std::exp(alpha_ * t);
}

void EPoly::MultByConstant(double c) {
	for (unsigned int i = 0; i < co_.size(); i++)
		co_[i] *= c;
}

double EPoly::EIntegral(){
	//if (abs(alpha_) < 1e-12)
		//return 1e10;
	if (alpha_ == 0)
		return 1e12;
	return 1.0 / alpha_;
}

bool EPoly::IsAlphaZero() {
	return abs(alpha_) < 1e-12;
}

EPoly EPoly::Integral() {
	if (IsZeroPoly())
		return EPoly({ 0 }, 1);
	if (IsAlphaZero())
		return PolyIntegrate();
	else
		return FullIntegrate();
}

pair<EPoly, double> EPoly::Integral_0_t() {
	EPoly integ(this->Integral());
	double eval_part = integ.Eval(0);
	return{ integ, eval_part };
}
EPoly EPoly::FullIntegrate() {
	vector<EPoly> ret;

	Poly cur((co_));// , alpha_);
	Poly cur2;
	int i = 1;
	Poly the_end_poly;// (static_cast<const Poly>(*this));
	double sign = 1.0;
	while (!cur.IsZeroPoly()) {
		EPoly temp((cur), alpha_);
		the_end_poly += cur * (sign * 1.0 /  pow(alpha_, i++));
		sign = -sign;
		cur2 = cur.Derivative();
		cur = cur2;
	}
	EPoly res(the_end_poly, alpha_);
	
	return res;
}


EPoly EPoly::PolyIntegrate() {
	vector<EPoly> ret;
	
	EPoly ep(Poly::Integral(), 0);
	ret.push_back(ep);
	
	return ret[0];
}

EPoly::EPoly(const EPoly & ep) : Poly(static_cast<const Poly&>(ep)) {
	alpha_ = ep.alpha_;
}

EPoly EPoly::operator * (const EPoly &p2) {
	Poly cur(static_cast<const Poly>(*this));
	Poly pp2(static_cast<const Poly &>(p2));
	Poly p_res = cur * pp2;
	//Poly p_res = static_cast<const Poly>(*this) * static_cast<const Poly &>(p2);
	double alpha = this->alpha_ + p2.alpha_;
	EPoly res(p_res, alpha);
	return res;
}
