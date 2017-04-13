#include <iomanip>
#include <ios>
#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>

#include "bn.h"
#include "em.h"
#include "context.h"
#include "ctbn.h"
#include "ctbndyn.h"
#include "factoredmatrix.h"
#include "factoredvector.h"
#include "markov.h"
#include "markovdyn.h"
#include "multirv.h"
#include "meanfieldinf.h"
#include "exactmarkovfbinf.h"
#include "streamextra.h"
#include "utils.h"
#include "ensurectbn.h"

#include "params.h"
#include "tope.h"
#include "rk.h" //"matexp.h"

using namespace std;
using namespace ctbn;


void fvprint(FactoredVector &fv){
	fv.Print(cout,NULL);
}

void vrprint(vectr &v){
	v.niceprint(cout);
}


int main(int argc, char * argv[]){

    cout.precision(10);
    
    addnan(cout);
    InitParams(argc, argv);
    double evideps = 0.00001;
    
    //CTBN file
	//string inputfile = ParamStr("TOPECTBN", "ctbn/toroid9_t2b1.ctbn");
    string inputfile = ParamStr("TOPECTBN", "ctbn/simple5.ctbn");
    

    //Evidence file
    //string trajfile = ParamStr("TOPETraj", "ctbn/toroid9.evid");
     string trajfile = ParamStr("TOPETraj", "ctbn/simple.evid");

	//TOPE params	
    double stime = 0.0;
    double etime = ParamDouble("TOPEetime",1.0);
    double qtime = ParamDouble("TOPEqtime",0.5);   
    double sseps = ParamDouble("TOPEsseps", 1e-3); //tolerance splitting children
    double maxt = ParamDouble("TOPEmaxt", 1); //max runtime
	
	//FactoredMatrix params
	//0: A is min
    //1: A is min but both A and B sum to zero
    //2: A is the avg according to initial distribution 
    //3: A is zero
    int fmratedist = ParamInt("FMRateDist", 1);
    int fmscale = ParamInt("FMScale", 0);
   	
 	//Load CTBN file
 	ifstream fin(inputfile.c_str());
	Markov m;
	try {
		m.Load(fin);
	} catch (const serial::streamexception &e) {
		cerr << e.what() << endl;
		exit(-1);
	} 	
   	fin.close();

	const CTBNDyn *ctbndyn = dynamic_cast<const CTBNDyn *>(m.GetDynamics());
	const BN *bn = dynamic_cast<const BN *>(m.GetStartDist());
	Context context = ctbndyn->Domain() + ctbndyn->CondDomain();

	//Load evidence file
	Trajectory evid;
	ifstream tfin(trajfile.c_str());
	evid.Load(tfin);
	tfin.close();

    //FactoredVector initialize
	RV *rvbn = (bn->Clone());
	FactoredVector fv(rvbn);
	
	//fv.Print(cout, NULL);
	
	FactoredVector efv(fv);
	
	FactoredVector unifv(fv);
	
	vectr p01(2,0.0);
	p01[0] = 0;  p01[1] = 1;
	
	vectr p10(2,0.0);
	p10[0] = 1;  p10[1] = 0;
	
	vectr p11(2,0.0);
	p11[0] = 0.5;  p11[1] = 0.5;
	
	vectr p001(3, 0.0);
	p001[0] = 0;  p001[1] = 0; p001[2] = 1;
	
	vectr p010(3, 0.0);
	p010[0] = 0;  p010[1] = 1; p010[2] = 0;  
	
	vectr p100(3, 0.0);
	p100[0] = 1;  p100[1] = 0; p100[2] = 0;
	
	vectr p111(3, 1.0);
	p111.normalize();
	
	/*
    for(int i=0; i<ctbndyn->NumofNodes(); i++){
		if(i == 0 || i == 2)
			efv.SetDist(p01, i, 0);
		else if (i==1 || i==3)
			efv.SetDist(p010, i, 0);
		else
		    efv.SetDist(p11, i, 0);
	}
	
	//efv.Print(cout, NULL);
	*/
	for(int i=0; i<ctbndyn->NumofNodes(); i++){
		if(evid.Value(i,evid.TimeBegin()) == 1){
		    fv.SetDist(p01, i, 0);
			/*if(ctbndyn->NodeByVar(i).Domain.Size() == 2)
				fv.SetDist(p01, i, 0);
			else if(ctbndyn->NodeByVar(i).Domain.Size() == 3)
				fv.SetDist(p01, i, 0);
			else
				cerr << "Could not load evidence pattern" <<  endl;
				*/
		}
		else
			fv.SetDist(p10, i, 0);
	}
	//fv.Print(cout,NULL);
	for(int i=0; i<ctbndyn->NumofNodes(); i++){
		if(evid.Value(i,evid.TimeEnd()-evideps) == 1)
			efv.SetDist(p01, i, 0);
		else
			efv.SetDist(p10, i, 0);
	}
	//efv.Print(cout, NULL);
	
	for(int i=0; i<ctbndyn->NumofNodes(); i++){
		unifv.SetDist(p11, i, 0);
	}
	
	//FactoredMatrix initialize
	Dynamics *dyn = (ctbndyn->Clone());
	//vectr avg(2, 0.5);
	FactoredVector avgfv(fv);
	for(int i=0; i<ctbndyn->NumofNodes(); i++){
        //avgfv.SetDist(avg, i, 0);	
        avgfv.SetDist((fv.Dist(i)+efv.Dist(i))/2, i, 0);	
	}

	FactoredMatrix fm(dyn, &avgfv, fmratedist, fmscale);	
	fm.SetL(20);
	fm.SetTheta(1);

    bool transpose = true;
    FactoredMatrix bwfm(dyn, &avgfv, fmratedist, fmscale, transpose);	//avgfv instead of efv
	bwfm.SetL(20);
	bwfm.SetTheta(1);

    // Spectral test
    TOPETree tree2 = TOPETree(fv, fm, qtime);

    //TOPE inference
	clock_t startTime = clock();
	cout << "Tope Result " << "CTBN " << inputfile << " sseps " << sseps << " ";
    //TOPETree tree = TOPETree(etime, qtime, sseps, maxt, fv, efv, fm, bwfm); //nonadaptive version
    TOPETree tree = TOPETree(etime, qtime, maxt, fv, efv, fm, bwfm);
    //(tree.treeFV).Print(cout,NULL);
    cout << endl;
    double tmptime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC ;
    cout << tmptime << endl; // << endl; 
   
    return 0;
}

