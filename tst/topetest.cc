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

    addnan(cout);
    InitParams(argc, argv);
    
    //const char *input_file = "chain.ctbn"; //"chain4.ctbn";// "drug.ctbn"; //"cycle.ctbn";
	//string inputfile = ParamStr("TOPECTBN", "bhpsnew.ctbn");
	//ifstream fin(inputfile.c_str());
	//CTBN file
	//string inputfile = ParamStr("TOPECTBN", "ctbn/toroid9_t2b1.ctbn");
    string inputfile = ParamStr("TOPECTBN", "ctbn/simple5.ctbn");
    

    //Evidence file
     //string trajfile = ParamStr("TOPETraj", "ctbn/toroid9.evid");
     string trajfile = ParamStr("TOPETraj", "ctbn/simple.evid");

	 //string trajfile = ParamStr("TOPETraj", "ctbn/simple.evid");
    ifstream fin(inputfile.c_str());

    double stime = 0.0;
    double etime = ParamDouble("TOPEetime",30.0);
    double eps = ParamDouble("TOPEeps",1e-8);
    int level = ParamInt("TOPElvl",1000); //now used for maxnodes
    int numofintervals = ParamInt("TOPEnuminterval",2); //not used
    
    //0: A is min
    //1: A is min but both A and B sum to zero
    //2: A is the avg according to initial distribution 
    //3: A is zero
    int fmratedist = ParamInt("FMRateDist", 2);
    int fmscale = ParamInt("FMScale", 0);
    
    double sseps = ParamDouble("TOPEsseps", 1e-3); //tolerance splitting children
    
    //stepsize
    double inith = ParamDouble("TOPEinith", 0.1);
    
    //max runtime
    double maxt = ParamDouble("TOPEmaxt", 0.1);
	
	cout.precision(10);
 	
	Markov m;//(fin);
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
		
    	Trajectory evid;
	ifstream tfin(trajfile.c_str());
	evid.Load(tfin);
	tfin.close();


    //evid.SetBeginTime(stime);
	//evid.SetEndTime(etime);
	
	
	ExactMarkovFBInf emfbinf;
	emfbinf.SetProcess(&m);
	emfbinf.SetTrajectory(&evid);
    etime  = 0.5;	
	if (1 == 1 || ParamInt("TOPEdoexact",1)) {
		cout << "Exact Result " << endl;
		clock_t startTimeEx = clock();
		for(int i=0; i<ctbndyn->NumofNodes(); i++){
		    for(int j=0; j<ctbndyn->Node(i)->Domain().Size(); j++){
			   Instantiation x(context, -1);
			   x.SetVal(i, j);
		
			   double filter_result = emfbinf.Filter(x,etime);        
			   cout <<  filter_result << " ";
		    }
		    cout << endl;
		}
		
		double endTimeEx = double( clock() - startTimeEx ) / (double)CLOCKS_PER_SEC ;
		/*
		for(int i=0; i<context.Size(); i++){		       
                	double filter_result = emfbinf.Filter(context.Index(i), etime);        
			//context.Index(i).PrintVal(cout);
			cout << filter_result << " ";
		}
		cout << endl;
		*/
		cout << endl;

		cout << endTimeEx << endl; 

	}
    exit(1);
	RV *rvbn = (bn->Clone());
	FactoredVector fv(rvbn);
/*	
	Dynamics *dyn = (ctbndyn->Clone());
	FactoredMatrix fm(dyn, &fv, fmratedist, fmscale);
	
	fm.SetL(20);
	fm.SetTheta(1);
   
    	
	clock_t startTime = clock();
	cout << "Tope Result " << "CTBN " << inputfile << " ABDist " << fmratedist << " ";
    TOPETree tree = TOPETree(etime, qtime, eps, numofintervals, level, fv, fm, sseps, inith, maxt); 
    (tree.treeFV).Print(cout,NULL);
    //joint
    //tree.trsum->niceprint(cout);
    cout << endl;
    double tmptime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC ;
    cout << tmptime << endl; // << endl; 
  */ 

}

