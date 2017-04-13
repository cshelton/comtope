/* Continuous Time Bayesian Network Reasoning and Learning Engine
 * Copyright (C) 2009 The Regents of the University of California
 *
 * Authored by Yu Fan, William Lam, Joon Lee, Christian Shelton, and Jing Xu
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
#ifndef CTBNRLE_UNIFORM_FACTORED_INF_H
#define CTBNRLE_UNIFORM_FACTORED_INF_H

#include "markov.h"
#include "ctbndyn.h"
#include "dyncomp.h"
#include "bn.h"
#include "rvcomp.h"
#include "inference.h"
#include "fbinf.h"
#include "factoredvector.h"
#include "factoredmatrix.h"

#include <vector>


namespace ctbn {

// This class performs exact inference by using the conditional distributions
// supplied by the process's (assumed to be a Markov process) dynamics's
// method "Cond."  For compact structures (like a CTBN) this will probably
// require generation of the full "flattened" dynamics which is exponentially
// (in the number of variables) large.  However, for some processes (like a
// fully factored process), it might be possible to return a conditional
// object that can represent this more compactly.
//
// The dynamics's "cond" methods should return RVCondSimple objects that
// are compatible with the RVSimple object's "MakeSimple" method (in that
// the RVCondSimple can accept the output of the MakeSimple method).

class UniformlyFactoredInf : public FBInf {
public:
    UniformlyFactoredInf();
	virtual ~UniformlyFactoredInf() throw();
	virtual UniformlyFactoredInf *Clone() const;
	/*	
	virtual void SetProcess(const Process *p);
	virtual double Smooth(const Instantiation &x, double t,
				bool log=false);
	*/
	void SetL(int l){ lval = l; }
	void SetTheta(double th) { thval = th; }

private:

//  now intervalpropagator
//  RVCondSimple * Cond(double t0, double t1, const Instantiation &x);
//  now pointtransitionpropagator (transition=true) or pointchangeevidencepropagator (transition=false)
//  RVCondSimple * Cond(double t0, const Instantiation &from, const Instantiation &to, bool transition=true);
	RVSimple * Convert(RVCondSimple *&rvcond, RVSimple *&rv);
	void Restrict(FactoredVector * & rv, Instantiation const & x) const;
//  void PrintRVSimple(RVSimple *x);
//  void PrintRVCond(RVCondSimple *x);
	
	
	int lval;
	double thval;
	/*
	double FilterReal(const Instantiation &x, double t,
			bool log=false, bool normalize=true);

	void GetProp();
	*/

public:

	virtual void Restrict(RVSimple *, Instantiation const & variable_assignment) const;
    
	virtual RVSimple * MakeSimple(RV * prior, Instantiation const &, bool normalize = true) const;
    virtual bool IsAlphaElementValid(RVSimple * const alpha_sub_x);

	virtual FactoredMatrix * IntervalPropagator(Dynamics *, double, double, Instantiation const &) const;
	virtual FactoredMatrix * PointTransitionPropagator(Dynamics *, double, Instantiation const &, Instantiation const &) const;
	virtual FactoredMatrix * PointChangeEvidencePropagator(Dynamics *, double, Instantiation const &, Instantiation const &) const;

	void AddExpectedSufficientStatistics(Dynamics const *, unsigned int, SS *, double) const { throw "Not yet implemented."; }
	void AddExpectedTransitionSufficientStatistics(Dynamics const *, unsigned int, SS *, double) const { throw "Not yet implemented."; }
	void AddExpectedInitialSufficientStatistics(Trajectory const *, RVSimple *, const RV *, SS *, double) const { throw "Not yet implemented."; }
};

} // end of ctbn namespace

#endif

