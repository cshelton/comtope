/* Continuous Time Bayesian Network Reasoning and Learning Engine
 * Copyright (C) 2009 The Regents of the University of California
 *
 * see docs/AUTHORS for contributor list
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CTBNRLE_CONTFUNCTION_H_
#define CTBNRLE_CONTFUNCTION_H_
#include <map>
#include <vector>

namespace ctbn {

// This class allows a function of the form scalar -> any numerical type
// to be represented by linear interpolation (scalar defaults to double)

template <class VALUE,class S=double>
class ContFunction {
public:
	typedef std::map<S,VALUE> BaseT;

	ContFunction();
	virtual ~ContFunction();

	// Returns the value at time t
	VALUE GetVal(S t) const;

	// Calls op(f,wt) twice f1*wt1 + f2*wt2 is the value that would
	// be returned by GetVal
	template<class OP>
	void ValOp(S t, OP &op) const;

	// Returns the derivative at time t (not defined at defined t points)
	VALUE GetDeriv(S t);

	// Returns the y-intercept at time t (not defined at defined t points)
	VALUE GetIntercept(S t);

	// Adds a value to the function at time t
	void AddVal(S t, const VALUE& val);

	// Combines the values stored in fun to this function
	void Combine(const ContFunction<VALUE,S> &fun);

	// Erases the values from t0 to t1, inclusively
	void EraseRange(S t0, S t1);

	// Erases the value at time t
	void EraseVal(S t);

	// Replaces the values in the range of fun with the values of fun
	void Replace(const ContFunction<VALUE,S> &fun);

	S GetStart() const {return start;}
	S GetEnd() const {return end;}

	// Returns a vector with all the times that a value is stored
	std::vector<S> GetCutPoints() const;

	// Changes this in place to have cut points at every cut point of
	// this or f.  Then calls op(a,b) on each cut point where a is the
	// value from this and b is the value from f.  op should take the
	// first argument by reference and change it.  Returns reference to *this
	// [Note: this does *NOT* do any sort of analysis to show that the
	//  resulting cut points are suitable for a linearly interpolated fn]
	template<class OpType, class V2>
	ContFunction<VALUE,S> &ApplyOp(const ContFunction<V2,S> &f, OpType op);

	// same, but for a unary op
	template<class OpType>
	ContFunction<VALUE,S> &ApplyOp(OpType op);

	void Print() const;
	void Clear();

private:

	// adds the points from f to this function (without using the values
	// from f)
	template<class V2>
	void Extend(const ContFunction<V2,S> &f);

	// lbound needs to be an iterator that points to the result of
	// calling lower_bound(t) on the map
	VALUE GetVal(S t, const typename BaseT::const_iterator &lbound) const;

	template<class OP>
	void ValOp(S t, OP &op,
			const typename BaseT::const_iterator &lbound) const;

	BaseT values;
	S start, end;
	bool empty;
};


} // end of name ctbn namespace

//#include "contfunction.tcc"
#endif
