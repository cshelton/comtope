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
#include "fwbase.h"
#include "nullptr03.h"


namespace ctbn {

using namespace std;

void FwBase::SetProcess(const Process *pr) {
	p = pr;
	const CTBNDyn *ctbndyn = 
		dynamic_cast<const CTBNDyn*>(dynamic_cast<const Markov *>(pr)->GetDynamics());
	ctbndyn->GetStructure(str);
	int n = ctbndyn->NumofNodes();
	nodes.clear();
	for (int i=0; i<n; i++)
		nodes.push_back(dynamic_cast<const MarkovDyn*>(ctbndyn->Node(i)));
}

double FwBase::SampleTime(const Instantiation &curri, int id,
		int parindex, int nodeval,
		double t, Random &rand) {
	double ret = 0.0;
	const matrix &Q = (*nodes[id])(parindex)->Intensity();
	//no upcoming evidence
	ret = t + method->SampleTime(Q, id, nodeval, curri, t, evid, -1, -1, rand);
	return ret;
}

double FwBase::SampleTimeWeight(const Instantiation &curri, int nodeid, int transitionid,
		int parindex, int nodeval,
		double currt, double nextt) {
	const matrix &Q = (*nodes[nodeid])(parindex)->Intensity();
	
	return method->TimeWeight(Q, curri, nodeid, transitionid, 
				nodeval, currt, nextt-currt, evid, -1, -1);
	
}

double FwBase::SampleInitial(Instantiation &initval, Random &rand) {
	const BN* bn = dynamic_cast<const BN*>(dynamic_cast<const Markov*>(p)->GetStartDist());
	return bn->ImportanceSample(initval, true,rand);
}

int FwBase::SampleTransition(Instantiation &inst, int id,
		int parindex, int nodeval,
		double t, double &weight, Random &rand) {
	forceflip = 0;
	int retval;
	const matrix &Q = (*nodes[id])(parindex)->Intensity();

		//there is no upcoming evidence
		//matrix Q = (*nodes[id])(parindex)->Intensity();
	retval = method->SampleTransition(Q, id, nodeval, inst, t, evid, -1, -1, weight, rand);
	//update the instantiation
	if (nodeval != retval) // below should probably be left to caller, as std.Node2Var(id)
			// might already have been computed there
		inst.SetVal(str.Node2Var(id), retval);
	return retval;
}

/*
double FwBase::SampleTransitionWeight(const Instantiation &nextinst, int transitionid, int oldval,
		int parindex, int nodeval,
		int evidencetype, Trajectory::Index &nextindex, int e,
		double currt, double nextt) {
	double retweight = 0.0;
	if (evidencetype==2) {//force a transition according to the evidence
		int val = oldval;
		int nextval = nodeval;
		matrix Q = (*nodes[transitionid])(parindex)->Intensity();
		retweight = log(Q[val][nextval]) - log(-(Q[val][val]));
	} else if (evidencetype==1&&e==-1)
		//evidence  change to known value to unknown
		retweight = 0;
	else if (evidencetype==1&&e!=-1) {
		//retweight = method->TransitionWeight();
		if (forceflip==1) {
			int val = oldval;
			int nextval = nodeval;
			matrix Q = (*nodes[transitionid])(parindex)->Intensity();
			retweight = log(Q[oldval][e]) - log(-(Q[val][val]));
			forceflip = 0;
			//cout << "force flip weight: " << retweight << endl;
		} else {
			double te = nextindex.Time();
			matrix Q = (*nodes[transitionid])(parindex)->Intensity();
			retweight = method->TransitionWeight(Q, transitionid, oldval, nextinst, nextinst, currt, nextt-currt, evid, te, e);
		}
	} else {
		retweight = 0;
	}
	return retweight;
}
*/

void FwBase::UndefineTime(int id, vector<int> &undefinelist) {
	undefinelist.clear();
	if (id!=-1) {
		changeindex[id] = 0;
		undefinelist = str.GetChildren(id);
		undefinelist.push_back(id);
		for (unsigned int i=0; i<undefinelist.size(); i++) {
			changeindex[undefinelist[i]] = 0;
		}
	} else {
		//conddomain changes
		Context conddomain = dynamic_cast<const Markov*>(p)->GetDynamics()->CondDomain();
		vector<int> varlist = conddomain.VarList();
		vector<bool> clist(nodes.size(), false);
		for (unsigned int i=0; i<varlist.size(); i++) {
			vector<int> c = str.GetChildren(str.Var2Node(varlist[i]));
			for (unsigned int j=0; j<c.size(); j++) {
				clist[c[j]] = true;
			}
		}
		for (unsigned int i=0; i<clist.size(); i++)
			if (clist[i]) 
				undefinelist.push_back(i);
	}
}



} // end of ctbn namespace
