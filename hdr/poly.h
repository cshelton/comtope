#ifndef POLY_H
#define POLY_H

#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;


class Poly {
public:
	Poly();
	Poly(double b);
	Poly(const vector<double> &c);

	Poly Mult(const Poly &p);
	void Mult(double c);
	
	Poly Derivative();
	Poly Integral();
	
	bool IsConstantPoly();
	bool IsZeroPoly();

	double Eval(double t);
	
	Poly operator = (const Poly &p);
	bool operator == (const Poly &p1);
	
	Poly operator +=(const Poly &rhs);
	Poly operator * (double multiplier); //Poly& won't work in cout << p * 2 because it creates a reference
	Poly operator * (const Poly &p);
	Poly operator -(){
		return this->Negation();
	}
	//friend Poly operator +( Poly lhs, const Poly &rhs);
	Poly operator +(const Poly &rhs) {
		Poly pp(*this);
		return pp += rhs;
	}
	
	vector<double> GetCoeffs();
	
	void Negate();
	Poly Negation();
	
	
	
	void NicePrint();

protected:
	bool IsEqualPoly(const vector<double> &v);
	bool IsEqualConstant(const double & d1, const double &d2);
	double eps_ = 1e-12;

	vector<double> co_;

public: 
	friend ostream& operator<< (ostream &outp, const Poly &p) {
		for (unsigned int i = 0; i < p.co_.size(); i++) {
			if (i) outp << " , ";
			outp << p.co_[i];// << ", ";
		}
		outp << endl;
		return outp;
	}
};

#endif
