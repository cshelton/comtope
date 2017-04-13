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
 *\
#include "exactmarkovinf.h"
#include <algorithm>


#include "multisimple.h"

namespace ctbn {

using namespace std;

ExactMarkovInf::ExactMarkovInf() : ForwBackInf() {
	\*
	p = NULL;
	tr = NULL;
	jointP0SS = NULL;
	jointDSS = NULL;
	*\

}

ExactMarkovInf::~ExactMarkovInf() {
	// don't delete anything... we don't own it (?)
	ClearCache();
}


ExactMarkovInf *ExactMarkovInf::Clone() const {
	// a little dangerous, but no more so than having this
	// kind of object in the first place
	ExactMarkovInf* cloneInf = new ExactMarkovInf(*this);
	

	for(unsigned int i=0;i<alpha.size();i++)
		cloneInf->alpha[i] = alpha[i]->Clone();
	for(unsigned int i=0;i<beta.size();i++)
		cloneInf->beta[i] = beta[i]->Clone();
	for(unsigned int i=0;i<prop.size();i++)
		cloneInf->prop[i] = prop[i]->Clone();
	
	if(cloneInf->jointP0SS)
		cloneInf->jointP0SS = jointP0SS->Clone();
	if(cloneInf->jointDSS)
		cloneInf->jointDSS = jointDSS->Clone();
	return cloneInf;
}

RVCondSimple* ExactMarkovInf::Cond(double t0, double t1, const Instantiation &x){
	return p->d->Cond(t0, t1, x);
}

RVCondSimple* ExactMarkovInf::Cond(double t0, const Instantiation &from, const Instantiation &to, bool transition){
	return p->d->Cond(t0, from, to, transition);
}

RVSimple* ExactMarkovInf::MakeSimple(const Instantiation &x, bool normalize){
	return p->p0->MakeSimple(x, normalize);
}
RVSimple* ExactMarkovInf::Convert(RVCondSimple *&rvcond, RVSimple *&rv){
	return rvcond->Convert(rv);
}

void ExactMarkovInf::Restrict(RVSimple *&rv, const Instantiation &x){
	vector<int> ind;
	p->d->Domain().ConsistentIndexes(ind,x);
	rv->Restrict(ind);
}

void ExactMarkovInf::PrintRVCond(RVCondSimple *x){

	if (dynamic_cast<CTBNDyn*>(x) != NULL)
		dynamic_cast<CTBNDyn*>(x)->JointMatrix().niceprint(cout);
		
	else if (dynamic_cast<matrix*>(x) != NULL)
		dynamic_cast<matrix*>(x)->niceprint(cout);
		
}
	
void ExactMarkovInf::PrintRVSimple(RVSimple *x){
//[added to rvsimple.h, rvsimple.cc]
	cout << "0 "; 
	PrintRVMarginal(x);
	//PrintRVJoint(x);
}

void ExactMarkovInf::PrintRVJoint(RVSimple *x){
//[added to rvsimple.h, rvsimple.cc]
	vectr v(p->p0->Domain().Size(), 0.0);
	double logf;
	x->GetDist(v, logf);
	
	
	//Printing scaled to when max = 1
	vectr scaled = v/v.max();
	cout <<"\nJoint Distribution:\t";
	for(int i=0; i<scaled.length(); i++){
		cout <<"["<<i<<"]="<<scaled[i]<<"\t";
	}
	cout << endl;
	//
	
	//Open this up once you're done.
	//cout << "Exact Logf: " << logf << endl;
	//v.niceprint(std::cout);
}

void ExactMarkovInf::PrintRVMarginal(RVSimple *x){
//[added to rvsimple.h, rvsimple.cc]
	vectr v(p->p0->Domain().Size(), 0.0);
	double logf;
	x->GetDist(v, logf);
	//cout << "Exact Logf: " << logf << endl;
	int nn = dynamic_cast<BN*>(p->p0)->NumofNodes();
	int div = 1;
	int ds, ids;
	vector<vectr> mdist;
	mdist.resize(nn);
	for(int i=0; i<nn; i++){
		ds = dynamic_cast<BN*>(p->p0)->NodeByVar(i)->Domain().Size();
		mdist.at(i) = vectr(ds, 0.0);
	}
	for (int ind=0; ind<v.length(); ind++){
		div = 1;
		for (int i=0; i<nn; i++){
			ds = mdist.at(i).length();
			ids = (ind/div)%ds;
			(mdist.at(i))[ids] = (mdist.at(i))[ids] + v[ind];
			div *= ds;
		}
	}
	for(int i=0; i<nn; i++){
		mdist.at(i).normalize();
		(mdist.at(i)).niceprint(std::cout);
	}
}
\*
void ExactMarkovInf::SetProcess(const Process *p) {
	ClearCache();
        this->p=NULL;
	this->p = dynamic_cast<const Markov *>(p);
}

void ExactMarkovInf::SetTrajectory(const Trajectory *tr) {
	ClearCache();
	this->tr = tr;
}

double ExactMarkovInf::FilterReal(const Instantiation &x, double t, bool log,
		bool normalize) {
	if (p==NULL || tr==NULL) return log ? -INFINITY: 0.0;
	MakeCache();

	int i = lower_bound(times.begin(),times.end(),t)-times.begin();
	RVSimple *a;
	if (times[i]==t) a = alpha[i]->Clone();
	else {
		if (i==0) return log ? -INFINITY : 0.0;
		a = alpha[i-1]->Clone();
		RVCondSimple *trans = p->d->Cond(times[i-1],t,
			tr->Values(p->d->Domain()+p->d->CondDomain(),t));
		trans->Mult(a);
		delete trans;
	}

	if (normalize) a->Normalize();
	double ret;
	if (x.NumVars()==0) 
		ret = a->Sum(log);
	else {
		vector<int> ind;
		p->d->Domain().ConsistentIndexes(ind,x);
		a->Restrict(ind);
		ret = a->Sum(log);
	}
	delete a;
	return ret;
}

double ExactMarkovInf::Filter(const Instantiation &x, double t, bool log) {
	return FilterReal(x,t,log);
}

double ExactMarkovInf::Smooth(const Instantiation &x, double t, bool log) {
	if (p==NULL || tr==NULL) return log ? -INFINITY : 0.0;
	MakeCache();
	int i = lower_bound(times.begin(),times.end(),t)-times.begin();
	RVSimple *a;
	if (times[i]==t) {
		a = alpha[i]->Clone();
		a->MultBy(beta[i]);
	} else {
		if (i==0 || i>=times.size()-1) return log ? -INFINITY : 0.0;
		a = alpha[i-1]->Clone();
		Instantiation inst = 
			tr->Values(p->d->Domain()+p->d->CondDomain(),t);

		RVCondSimple *trans = p->d->Cond(times[i-1],t,inst);
		trans->Mult(a);
		delete trans;

		RVSimple *b = beta[i]->Clone();
		trans = p->d->Cond(t,times[i],inst);
		trans->RMult(b);
		delete trans;
		a->MultBy(b);
		delete b;
	}

	a->Normalize();
	double ret;
	if (x.NumVars()==0) 
		ret = a->Sum(log);
	else {
		vector<int> ind;
		p->d->Domain().ConsistentIndexes(ind,x);
		a->Restrict(ind);
		ret = a->Sum(log);
	}
	delete a;
	return ret;
}
double ExactMarkovInf::Prob(double t, bool log) {
	return FilterReal(Instantiation(),t,log,false);
}

void ExactMarkovInf::ClearCache() {
	times.resize(0);
	for(int i=0;i<alpha.size();i++)
		delete alpha[i];
	alpha.resize(0);
	for(int i=0;i<beta.size();i++)
		delete beta[i];
	beta.resize(0);
	for(int i=0;i<prop.size();i++)
		delete prop[i];
	prop.resize(0);
	condi.resize(0);
	transtype.resize(0);
	changevars.resize(0);
	if(jointP0SS) {
		delete jointP0SS;
		jointP0SS = NULL;
	}
	if(jointDSS) {
		delete jointDSS;
		jointDSS = NULL;
	}
}

void ExactMarkovInf::MakeCache() {
	if (times.size()>0) return;
        //cout << "MAKE CACHE" << endl;
	GetProp();
	AlphaPass();
	BetaPass();
}

// This method goes and retrieves all of the sequence of conditional
// distributions (given the evidence in the trajectory) transitioning
// from one event to the next and also "across" events that signify a
// transition.
void ExactMarkovInf::GetProp() {
	times.push_back(tr->TimeBegin());
	Trajectory::Index i = tr->Begin(p->d->Domain()+p->d->CondDomain());
	while(!i.Done()) {
		RVCondSimple *pr = p->d->Cond(i.Time(),
				i.Time()+i.DeltaT(),i.Values());
		prop.push_back(pr);
		transtype.push_back(0);
		changevars.push_back(Context());
		condi.push_back(p->d->CondDomain().Index(i.Values()));


		Instantiation oldv = i.Values();
		int varchange = i.TestInc(p->d->Domain());
		times.push_back(i.Time());

		if (varchange>1 && !i.Done()) {
			pr = p->d->Cond(i.Time(),oldv,i.Values());
			prop.push_back(pr);
			condi.push_back(p->d->CondDomain().Index(i.Values()));
			times.push_back(i.Time());
			transtype.push_back(2);
			changevars.push_back(Context(oldv,i.Values()));
		} else if (varchange==1 && !i.Done()) {
			pr = p->d->Cond(i.Time(),oldv,i.Values(),false);
			prop.push_back(pr);
			condi.push_back(p->d->CondDomain().Index(i.Values()));
			transtype.push_back(1);
			changevars.push_back(Context());
			times.push_back(i.Time());
		} // else a change in conditioning set only...
	}

}

// This method propagates the classic "alpha pass" starting with
// the initial distribution indicated by the Markov process's starting
// distribution and running it forward through the conditional distributions
// set up by GetProp.  They are saved for easy access later.
void ExactMarkovInf::AlphaPass() {
	int n = prop.size();
	alpha.resize(n+1);

	RVSimple *curr = p->p0->MakeSimple(
				tr->Values(
				p->p0->Domain()+p->p0->CondDomain(), 
				tr->TimeBegin()));

	if (n>0) {
		alpha[0] = prop[0]->Convert(curr);
		delete curr;
		curr = alpha[0]->Clone();
	} else {
		alpha[0] = curr->Clone();
	}
	if (!(dynamic_cast<SparseMultiZSimple *>(alpha[0]))->Dist().isvalid())
		cerr << "alpha[0] is invalid" << endl;
	for(int i=0;i<n;i++) {
		prop[i]->Mult(curr);
		alpha[i+1] = curr->Clone();
		if (!(dynamic_cast<SparseMultiZSimple *>(alpha[i+1]))->Dist().isvalid())
			cerr << "alpha[" << i+1 << "] is invalid" << endl;
	}

	delete curr;
}

// This is the same as AlphaPass, but running the distribution (starting
// with a uniform one) backward through the conditional distributions.
void ExactMarkovInf::BetaPass() {
	int i = prop.size();
	beta.resize(i+1);
	RVSimple *curr = alpha[i]->Clone();
	curr->MakeUniform(); // set to all 1s
	beta[i] = curr->Clone();
	for(i--;i>=0;i--) {
		prop[i]->RMult(curr);
		beta[i] = curr->Clone();
	}
	delete curr;
}

void ExactMarkovInf::AddExpSuffStats(SS *ss, double w) {

	MarkovSS *mss = dynamic_cast<MarkovSS *>(ss);
	MakeCache();
	int n = prop.size();
	for(int i=0;i<n;i++) {
		// we could construct the correct full RV and call
		// AddExpSS on p->d, but that's a lot of extra unneeded
		// work given that we know the underlying structure
		if (transtype[i]==0) 
			p->d->AddExpSS(condi[i],alpha[i],beta[i+1],
					times[i],times[i+1]-times[i],
					mss->dss,w);
	   
		else if (transtype[i]==2) 
			p->d->AddExpTransSS(condi[i],alpha[i],beta[i+1],
						changevars[i],times[i],
						mss->dss,w);
	}
	if (alpha.size()>0) {
		// add suffstats for initial distribution
		RVSimple *init = alpha[0]->Clone();
		init->MultBy(beta[0]);
		init->Normalize();
		p->p0->AddExpSS(tr->
				Values(p->p0->CondDomain(),tr->TimeBegin()),
				init,mss->p0ss,w);
		delete init;
	}
}

void ExactMarkovInf::AddExpSuffStats(const Dynamics *dyn, SS *ss, double w) {
	const CTBNDyn* cdyn = dynamic_cast<const CTBNDyn *>(p->d);
	const DynComp* jcdyn = 
			dynamic_cast<const DynComp *>(cdyn->GetJoint());
	if(!jointDSS) {
		DynCompSS* dcss = dynamic_cast<DynCompSS *>(jcdyn->BlankSS());
		MakeCache();
		int n = prop.size();
		for(int i=0;i<n;i++) {
			if (transtype[i]==0) 
				jcdyn->AddExpSS(condi[i],alpha[i],beta[i+1],
						times[i],times[i+1]-times[i],
								dcss,w);
			else if (transtype[i]==2) 
				jcdyn->AddExpTransSS(condi[i],alpha[i],
						beta[i+1],changevars[i],
							times[i],dcss,w);

		}
		jointDSS = dcss;
	}
	dyn->AddSS(jointDSS,jcdyn,ss,w);
}

void ExactMarkovInf::AddExpSuffStats(const RV *p0, SS *ss, double w) {
	const BN* bn = dynamic_cast<const BN *>(p->p0);		
	RVComp *jbn = dynamic_cast<RVComp*>(bn->Joint(
			Instantiation(bn->Domain()+bn->CondDomain())));
	if(!jointP0SS) {
		RVCSCompSS* rvcscss = 
			dynamic_cast<RVCSCompSS *>(jbn->BlankSS());
		MakeCache();
		if (alpha.size()>0) {
			RVSimple *init = alpha[0]->Clone();
			init->MultBy(beta[0]);
			init->Normalize();
			jbn->AddExpSS(tr->Values(jbn->CondDomain(),
						tr->TimeBegin()),
						init,rvcscss,w);
			delete init;
		}
		jointP0SS = rvcscss;
	}
	p0->AddSS(jointP0SS,jbn,ss,w);
}

double ExactMarkovInf::CalcQuery(QueryCalculator &calc) {
	SS *ss = p->BlankSS();
	ClearCache();
	AddExpSuffStats(ss);
	double ret = calc.Calculate(*ss);
	delete ss;
	return ret;
}
*\
} // end of ctbn namespace
*/
