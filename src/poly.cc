#include "poly.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
//template <class T>

Poly::Poly() {
	co_.push_back(0);
}

Poly::Poly(double b) {
	co_.push_back(b);
}

Poly::Poly(const vector<double> &c) {
	if (c.size() == 0)
		Poly();
	else {
		for (unsigned int i = 0; i < c.size(); i++)
			co_.push_back(c[i]);
	}
}

Poly Poly::Derivative() {
	vector<double> c;
	if (co_.size() == 1) c.push_back(0);
	for (unsigned int i = 1; i < co_.size();  i++) {
		c.push_back(co_[i] * i);
	}
	Poly p(c);
	return p;
}

Poly Poly::Integral() {
	vector<double> ret;
	ret.resize(co_.size() + 1);
	ret[0] = 0;
	for (int i = (int)co_.size() - 1; i >= 0; i--) {
		ret[i + 1] = 1.0 / (i + 1) * co_[i];
	}
	Poly p(ret);
	return p;
}

bool Poly::IsZeroPoly() {
	return (co_.size() == 1 && co_[0] == 0);//abs(co_[0]) < 1e-16);// || IsAlphaZero();
}

Poly Poly::Mult(const Poly &p) {
	vector<double> d;
	d.resize(co_.size() + p.co_.size());
	for (unsigned int i = 0; i < co_.size(); i++) {
		for (unsigned int j = 0; j < p.co_.size(); j++) {
			d[i + j] += co_[i] * p.co_[j];
		}
	}
	Poly p_ret(d);
	return p_ret;
}

bool Poly::IsConstantPoly() {
	return (/*co_[0] == 0 &&*/ co_.size() == 1);
}

Poly Poly::operator=(const Poly &p) {
	if (this == &p) return *this;
	co_ = p.co_;
	return *this;
}

bool Poly::operator== (const Poly &p1) {
	if (p1.co_.size() != co_.size()) return false;
	return IsEqualPoly(p1.co_);
}



bool Poly::IsEqualPoly(const vector<double> &v) {
	for (unsigned int i = 0; i < v.size(); i++) {
		if (!IsEqualConstant(v[i], co_[i])) return false;
	}
	return true;
}

Poly Poly::operator * (double multiplier) {
	vector<double> c =  co_;
	for (unsigned int i = 0; i < c.size(); i++)
		c[i] *= multiplier;
	
	Poly ret(c);
	return ret;
}

bool Poly::IsEqualConstant(const double & d1, const double &d2) {
	return abs(d1 - d2) <= eps_;
}

double Poly::Eval(double t) {
	double ret = co_[0];
	if (co_.size() == 1) return ret;
	double p = t;
	for (unsigned int i = 1; i < co_.size(); i++) {
		ret +=  p * co_[i];
		p *= t;
	}

	return ret;
}

vector<double> Poly::GetCoeffs() {
	return co_;
}

Poly operator +(Poly lhs, const Poly &rhs) {
	lhs += rhs;
	return lhs;
}

void Poly::Negate() {
	for (unsigned int i = 0; i < co_.size(); i++) {
		co_[i] *= -1.0;
	}
}

Poly Poly::Negation() {
	Poly p(*this);
	p.Negate();
	return p;
}

Poly Poly::operator +=(const Poly &rhs) {
	vector<double> new_co, small_co;
	new_co.resize(fmax(co_.size(), rhs.co_.size()) );
	if (co_.size() >= rhs.co_.size()) {
		new_co = co_; small_co = rhs.co_;
	}
	else {
		new_co = rhs.co_; small_co = co_;
	}
	for (int i = 0; i < (int)small_co.size(); i++)
		new_co[i] += small_co[i];
	co_ = new_co;
	return *this;
}

void Poly::NicePrint() {
	if (IsZeroPoly()) cout << "0\n";
	else {
		bool some_print = false;
		for (int i = (int)co_.size() - 1; i >= 0; i--) {
			if (co_[i] != 0) {
				if (co_[i] > 0) cout << " + "; else cout << " - ";
				if (abs(co_[i]) != 1) cout << abs(co_[i]);
			
				if (i != 0) cout << "x^" <<  i;
				//cout << " ";
			}
		}
		cout << "\n";
	}
}

Poly Poly::operator * (const Poly &p) {
	vector<double> co;
	co.resize(co_.size() + p.co_.size() - 1);
	for (int i = 0; i < (int)co_.size(); i++) {
		for (int j = 0; j < (int)p.co_.size(); j++) {
			co[i + j] += co_[i] * p.co_[j];
		}
	}
	Poly ret(co);
	return ret;
}

