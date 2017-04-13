
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>

#include "bn.h"
#include "em.h"
#include "context.h"
#include "ctbn.h"
#include "ctbndyn.h"
#include "markov.h"
#include "markovdyn.h"
#include "multirv.h"
#include "trajectory.h"

#include "params.h"

#include "comttop.h"

using namespace std;
using namespace ctbn;

double tao = 1;
double beta = 0.5;

void FillIntensityIsing(MarkovDyn &md, Instantiation &i, int pa1, int pa2){
	int pv1, pv2;
	for (int v1=0; v1<2; v1++){
		i.SetVal(pa1, v1);
		for (int v2=0; v2<2; v2++){
			i.SetVal(pa2, v2);
			pv1 = (v1==0)?-1:1;
			pv2 = (v2==0)?-1:1;

			md(i)->Intensity()[1][0] = tao*pow((1+exp(2*beta*(pv1+pv2))), -1.0); 
			md(i)->Intensity()[0][1] = tao*pow((1+exp(-2*beta*(pv1+pv2))), -1.0);	
			md(i)->Intensity()[0][0] = -1*(md(i)->Intensity()[0][1]);
			md(i)->Intensity()[1][1] = -1*(md(i)->Intensity()[1][0]);
		}
	}

}

void FillIntensityIsingOp(MarkovDyn &md, Instantiation &i, int pa1, int pa2){
	int pv1, pv2;
	for (int v1=0; v1<2; v1++){
		i.SetVal(pa1, v1);
		for (int v2=0; v2<2; v2++){
			i.SetVal(pa2, v2);
			pv1 = (v1==0)?-1:1;
			pv2 = (v2==0)?-1:1;

			md(i)->Intensity()[1][0] = tao*pow((1+exp(-2*beta*(pv1+pv2))), -1.0); 
			md(i)->Intensity()[0][1] = tao*pow((1+exp(2*beta*(pv1+pv2))), -1.0);	
			md(i)->Intensity()[0][0] = -1*(md(i)->Intensity()[0][1]);
			md(i)->Intensity()[1][1] = -1*(md(i)->Intensity()[1][0]);
		}
	}
}

void FillIntensityIsing(MarkovDyn &md, Instantiation &i, int pa1){

	int pv1;
	for (int v1=0; v1<2; v1++){
		i.SetVal(pa1, v1);
			pv1 = (v1==0)?-1:1;
			
			md(i)->Intensity()[1][0] = tao*pow((1+exp(2*beta*(pv1))), -1.0); 
			md(i)->Intensity()[0][1] = tao*pow((1+exp(-2*beta*(pv1))), -1.0);	
			md(i)->Intensity()[0][0] = -1*(md(i)->Intensity()[0][1]);
			md(i)->Intensity()[1][1] = -1*(md(i)->Intensity()[1][0]);

	}

}


void FillIntensityIsing(MarkovDyn &md){
	double pv1 = 0;
		//xo and pv = {-1, +1}
		//tao*(1+exp(-2*xo*beta*sum(pv)))^-1
		md(0)->Intensity()[1][0] = tao*pow((1+exp(2*beta*(pv1))), -1.0); 
		md(0)->Intensity()[0][1] = tao*pow((1+exp(-2*beta*(pv1))), -1.0);	
		md(0)->Intensity()[0][0] = -1*(md(0)->Intensity()[0][1]);
		md(0)->Intensity()[1][1] = -1*(md(0)->Intensity()[1][0]);
		
	

}

Markov makesimple(kprod<vectr> &initp) {
	Context x0,x1,x2;
	Context Null;
	x0.AddVar(0,2);
	x1.AddVar(1,2);
	x2.AddVar(2,2);
	MarkovDyn PX0(x0,Context(Null));
	PX0(0)->Intensity()[0][0] = -2;
	PX0(0)->Intensity()[0][1] = 2;
	PX0(0)->Intensity()[1][0] = 1;
	PX0(0)->Intensity()[1][1] = -1;

	MarkovDyn PX1(x1,Context(x0));
	PX1(0)->Intensity()[0][0] = -1;
	PX1(0)->Intensity()[0][1] = 1;
	PX1(0)->Intensity()[1][0] = 3;
	PX1(0)->Intensity()[1][1] = -3;
	PX1(1)->Intensity()[0][0] = -3;
	PX1(1)->Intensity()[0][1] = 3;
	PX1(1)->Intensity()[1][0] = 1;
	PX1(1)->Intensity()[1][1] = -1;

	MarkovDyn PX2(x2,Context(x1));
	PX2(0)->Intensity()[0][0] = -1;
	PX2(0)->Intensity()[0][1] = 1;
	PX2(0)->Intensity()[1][0] = 3;
	PX2(0)->Intensity()[1][1] = -3;
	PX2(1)->Intensity()[0][0] = -3;
	PX2(1)->Intensity()[0][1] = 3;
	PX2(1)->Intensity()[1][0] = 1;
	PX2(1)->Intensity()[1][1] = -1;

	MultiRV P0x0(x0,Null);
	vectr p0(2,0.0);
	p0[0] = 1;  p0[1] = 0;
	P0x0[0].SetDist(p0);
	initp.x.push_back(p0.transpose());

	MultiRV P0x1(x1,Null);
	vectr p1(2,0.0);
	p1[0] = 1;  p1[1] = 0;
	P0x1[0].SetDist(p1);
	initp.x.push_back(p1.transpose());

	MultiRV P0x2(x2,Null);
	vectr p2(2,0.0);
	p2[0] = 1;  p2[1] = 0;
	P0x2[0].SetDist(p2);
	initp.x.push_back(p2.transpose());

	CTBNDyn ctbndyn;
	ctbndyn.AddNode(PX0.Clone());
	ctbndyn.AddNode(PX1.Clone());
	ctbndyn.AddNode(PX2.Clone());

	BN bn;
	bn.AddNode(P0x0.Clone());
	bn.AddNode(P0x1.Clone());
	bn.AddNode(P0x2.Clone());

	Markov ctbn(bn.Clone(),ctbndyn.Clone());
	Context context = ctbndyn.Domain() + ctbndyn.CondDomain();

	return ctbn;
}


Markov maketoroid(int argc, char **argv, kprod<vectr> &initp, kprod<vectr> &endp) {
	if (argc > 2){
		tao = atof(argv[1]);
		beta = atof(argv[2]);
	}

	const int nr = 3;
	const int nc = 3;
	Context Null;
	std::vector<Context> xs;
	for(int i=0;i<nr*nc;i++) {
		Context x;
		x.AddVar(i,2);
		xs.emplace_back(std::move(x));
	}
	CTBNDyn ctbndyn;
	BN bn;
	for(int r=0,k=0;r<nr;r++) {
		for(int c=0;c<nc;c++,k++) {
			int pa1 = (r+nr-1)%nr + c*nr;
			int pa2 = r + ((c+nc-1)%nc)*nr;
			Context pa(xs[pa1],xs[pa2]);
			MarkovDyn X(xs[k],pa);
			Instantiation i(pa);
			FillIntensityIsing(X,i,pa1,pa2);

			MultiRV P0(xs[k],Null);
			vectr p0(2,0.0);
			vectr pT(2,0.0);
			if (r==0) {
				cout << "01 ";
				p0[0] = 1; p0[1] = 0;
				pT[0] = 0; pT[1] = 1;
			} else if (r==nr-1 || c>nc/2) {
				cout << "10 ";
				p0[0] = 0; p0[1] = 1; 
				pT[0] = 1; pT[1] = 0;
			} else {
				cout << "00 ";
				p0[0] = 1; p0[1] = 0; 
				pT[0] = 1; pT[1] = 0;
			}
			/*
			if (k<nr*nc/2) {
				cout << "01 ";
				p0[0] = 1; p0[1] = 0;
				pT[0] = 0; pT[1] = 1;
			} else {
				cout << "10 ";
				p0[0] = 0; p0[1] = 1; 
				pT[0] = 1; pT[1] = 0;
			}
			*/
			P0[0].SetDist(p0);
			initp.x.push_back(p0.transpose());
			endp.x.push_back(pT.transpose());

			ctbndyn.AddNode(X.Clone());
			bn.AddNode(P0.Clone());
		}
		cout << endl;
	}

	Markov ctbn(bn.Clone(),ctbndyn.Clone());
	Context context = ctbndyn.Domain() + ctbndyn.CondDomain();

	return ctbn;

}

template<typename MT>
void dumpdecomp(const ctbncomdecomp<MT> &decomp, double T) {
	Context c;
	for(int i=0;i<decomp.A.x.size();i++)
		c.AddVar(i,decomp.A.x[i].getm());
	cout << "decomp: " << endl;
	cout << "A: " << endl;
	cout << decomp.A << endl;
	cout << "A*" << T << ", flattened: " << endl;
	matrix A = decomp.A.flatten();
	cout << A << endl;
	matrix eAT;
	expmt(eAT,A,1.0);
	cout << "e^{A*" << T << "} = " << eAT << endl;
	cout << "comseqs: " << endl;
	for(int i=0;i<decomp.comseqs.size();i++) {
		cout << "sequence " << i << endl;
		cout << "\tA = " << decomp.comseqs[i].A << endl;
		cout << "\tA, flat = " << decomp.comseqs[i].A.flatten(c) << endl;
		cout << "\tC = " << endl;
		for(int j=0;j<decomp.comseqs[i].C.size();j++) {
			cout << "\t\t" << j << ": " << decomp.comseqs[i].C[j] << endl;
			cout << "\tB, flat = " << decomp.comseqs[i].C[j].flatten(c) << endl;
		}
	}
	/*
	cout << "derivs: " << endl;
	for(int i=0;i<decomp.derivs.size();i++) {
		cout << "deriv #" << i << endl;
		cout << "\tB small = " << decomp.derivs[i].B << endl;
		cout << "\tB flat =  " << decomp.derivs[i].B.flatten(c) << endl;
		cout << "\td small = " << decomp.derivs[i].d << endl;
		cout << "\td flat =  " << decomp.derivs[i].d.flatten(c) << endl;
	}
	*/
}

void dumpexact(const Markov &ctbn, double T, const kprod<vectr> &p0) {
	ofstream exactf("exact.txt");
	const CTBNDyn &dyn = dynamic_cast<const CTBNDyn &>(*ctbn.GetDynamics());
	//const BN &p0 = dynamic_cast<const BN &>(*ctbn.GetStartDist());
	matrix Q = dyn.JointMatrix()*T;
	cout << "Q*" << T << " = " << Q << endl;
	exactf << Q << endl;
	matrix eQT;
	expmt(eQT,Q,1.0);
	cout << "e^{Q*" << T << "} = " << eQT << endl;
	exactf << eQT << endl;
	vectr v0 = flatten(p0);
	cout << "v0 = " << v0 << endl;
	vexpmt(v0,Q,1.0);
	cout << "v0e^{Q*" << T << "} = " << v0 << endl;
	exactf << v0 << endl;
	
}

int main (int argc, char **argv) {
	const double T = 1.0, t=0.5;
	kprod<vectr> p0,pT;
	p0.scalar = 1.0;
	Markov toroid = maketoroid(argc,argv,p0,pT);
	//Markov toroid = makesimple(p0);
	//ctbncomdecomp<matrix> decomp(dynamic_cast<const CTBNDyn&>(*toroid.GetDynamics()));

	//dumpexact(toroid,T,p0);

	/*
	while(1) {
		int i,k;
		cout << "[node #] [commutator power]:";
		cin >> i >> k;
		const partialmatrix<matrix> &m = decomp.comseqs[i][k];
		cout << "full matrix: " << m << endl;
		cout << "as kproduct: " << endl;
		auto pieces = m.tokprod();
		for(int i=0;i<pieces.size();i++)
			cout << "\tpart " << i << ": " << pieces[i] << endl;
	}
	*/

	commexp<matrix,vectr> ce(p0,dynamic_cast<const CTBNDyn&>(*toroid.GetDynamics()),p0,t);
	//smoother<matrix,vectr> ce(p0,pT,dynamic_cast<const CTBNDyn&>(*toroid.GetDynamics()),T,t);

	//dumpdecomp(ce.decomp,T); 

	run(ce,1000000);

	// trajectory
	Trajectory tr;
	tr.SetBeginTime(0.0);
	tr.SetEndTime(1.0);
	tr.AddTransition(0,0.5,1);
	tr.AddTransition(1,0.5,0);

	// query
	//cout << "0 2 0 1.0 0.8983018202" << endl;

	return 0;
}
