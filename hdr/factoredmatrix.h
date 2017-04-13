#ifndef CTBNRLE_FACTORED_MATRIX_H
#define CTBNRLE_FACTORED_MATRIX_H

#include "ctbndyn.h"
#include "factoredvector.h"
// cshelton: added
#include "expmcache.h"
//#include "matexp.h"
//#include "rk.h"
#include "eigendecompose.h"
#include <vector>

namespace ctbn {
using namespace std;
#define RKEPS 1e-8

	class FactoredMatrix;

	// private inheritance, b/c this does not really implement everything 
	// from FactoredVector (it just uses one and wants access to dist
	class LagFactoredVector : private FactoredVector {
		friend class FactoredMatrix;
		public:
		LagFactoredVector() : FactoredVector() {
		}

		LagFactoredVector(const FactoredVector &fv, double startt)
			: FactoredVector(fv), realt(dists.size(),startt) {
				t = startt;
			}
		void Assign(const FactoredVector &fv, double startt) {
			FactoredVector::operator=(fv);
			realt.assign(dists.size(),startt);
			t = startt;
		}

		double Diff(const FactoredVector &fv) const {
			return FactoredVector::Diff(fv);
		}
		double Diff2(const FactoredVector &fv) const {
			return FactoredVector::Diff2(fv);
		}


		double AddTo(FactoredVector &fv, double wt, vectr *jointv) const {

			//joint
			//(*jointv) += wt*(fv.Joint(*this));		
			return fv.AddMult(*this,wt);
		}


		double Mag() const {

			double ret = 1.0;
			for(uint i=0;i<dists.size();i++)
				ret *= dists[i].absmax();
			return ret;

			/*
			double ret = dists[0].absmax();//1.0;
			for(uint i=1;i<dists.size();i++){
			if(ret < dists[i].absmax())
			ret = dists[i].absmax();
			}

			return ret;
			*/

		}

		void DotStar(LagFactoredVector &fv, LagFactoredVector &out){
			for(unsigned int i=0; i<dists.size(); i++)
				out.dists[i] = dists[i].dotstar(fv.dists[i]);
			out.logf = logf + fv.logf;
		}

		void Swap(LagFactoredVector &fv) {
			FactoredVector::Swap(fv);
			realt.swap(fv.realt);
			double tt = fv.t;
			fv.t = t;
			t = tt;
		}

		void LFVPrint(std::ostream & out, RV * ignored = nullptr03) const{
			this->Print(out, ignored);
		}

		private:
		std::vector<double> realt;
		double t;
	};

	class FactoredMatrix :public RVCondSimple {

		friend class FactoredVector;
		friend class TOPETree;
		friend class TSPEC;

		public:
		FactoredMatrix();
		FactoredMatrix(Dynamics *dyn, FactoredVector *fv=NULL, int rdist = 0, int scd = 0, 
				bool transpose = false);
		FactoredMatrix(const FactoredMatrix &m);
		~FactoredMatrix();

 

		RVCondSimple *Clone() const;

		double FindMinVal(Dynamics *node) const;
		double FindMinRate() const;
		void Cond(double t0, double t1, const Instantiation &x);
		void Cond(double t0, const Instantiation &from, const Instantiation &to, bool transition=true);
		void CondMultByVec(FactoredVector &v, bool transpose, double &maxval) const;
		void MultByVec(FactoredVector &v, bool transpose, double &maxval) const;
		void TransMultByVec(FactoredVector &v, bool transpose) const;
		void TransMult(RVSimple *&x, bool transpose) const;
		void Mult(RVSimple *&x) const;
		void Mult(RVSimple *&x, bool transpose) const;
		double FindAvgRate(FactoredVector &v, const Dynamics *node, int from, int to,
				const std::vector<int> &cind, bool transpose, double &maxval) const;

		void RMult(RVSimple *&x) const;
		double Unifvexpmt(FactoredVector &v, bool conditioned, bool transition) const;

		//tope
		// cshelton: added
		double Diff(LagFactoredVector &fv1, LagFactoredVector &fv2) const {
			for(uint i=0;i<fv1.dists.size();i++)
				if (fv1.realt[i]!=fv2.realt[i]) {
					Advance(fv1,i); Advance(fv2,i);
				}
			return fv1.Diff(fv2);
		}
		double Diff2(LagFactoredVector &fv1, LagFactoredVector &fv2) const {
			for(uint i=0;i<fv1.dists.size();i++)
				if (fv1.realt[i]!=fv2.realt[i]) {
					Advance(fv1,i); Advance(fv2,i);
				}
			return fv1.Diff2(fv2);
		}
		/*
		   double DiffTrans(LagFactoredVector &fv1, LagFactoredVector &fv2) const {
		   for(uint i=0;i<fv1.dists.size();i++)
		   if (fv1.realt[i]!=fv2.realt[i]) {
		   AdvanceTrans(fv1,i); AdvanceTrans(fv2,i);
		   }
		   return fv1.Diff(fv2);
		   }
		   */
		void initRKcache(double maxtime, double eps=1e-12) const;
		void ForAndBack(FactoredVector &v, double time, int Bid) const;
		void Forward(FactoredVector &v, double time) const;
		void MultByB(FactoredVector &v, int Bid) const;
		int nBs() const;

		typedef std::pair<int,Instantiation> Binfo;
		const Binfo &getBinfo(int Bi) const { return Blist[Bi]; } 
		
		double RKvexpmtA(FactoredVector & v, double time, bool conditioned, bool transpose,
				std::vector<std::vector<double> > *timesteps=NULL, double higheps = 1e-5, double inith = -1.0) const;
		//void BMultByVec(FactoredVector &v, int varid, int bno, bool transpose) const;
		void BMultByVec(FactoredVector &v, int varid,
				const Instantiation &pinst, bool transpose) const;

		double RKvexpmtA(LagFactoredVector & v, double time, bool forceupdate=false) const;
		//double RKvexpmtATrans(LagFactoredVector & v, double time, bool forceupdate=false) const;

		double RKvexpmtA(LagFactoredVector & v, double time, int varid,
				const Instantiation &pinst, const LagFactoredVector &vend) const;
		//double RKvexpmtATrans(LagFactoredVector & v, double time, int varid,
		//	const Instantiation &pinst, const LagFactoredVector &vend) const;

		void BMultByVec(LagFactoredVector &v, int varid,
				const Instantiation &pinst) const;

		//void BMultByVecTrans(LagFactoredVector &v, int varid,
		//	const Instantiation &pinst) const;
		//
		//
		matrix *MakeAB(Dynamics *bnode, FactoredVector *fv, double &asc);
		matrix *MakeABMin(Dynamics *bnode, double &asc);
		matrix *MakeABMinSZ(Dynamics *bnode, double &asc);

		matrix *MakeABInitAvg(Dynamics *bnode, FactoredVector *fv, double &asc);
		matrix *MakeABAvg(Dynamics *bnode, double &asc);
		matrix *MakeAZero(Dynamics *bnode, double &asc);

		void SetL(int terms_in_taylor_expansion);
		void SetTheta(double theta);

		void Load(std::istream &is);
		void Save(std::ostream &os) const;


		// returns the probability of an individual conditional event
		double Prob(int ind, int cind, bool log=false) const;

		// Makes sum to 1 for each conditioning value
		void Normalize();

		// returns an object representing the conditioning
		// on that value
		RVSimple *Condition(int ind) const;

		// Sets to zero all indexes that are not listed in ind and cind
		// That is, for conditions mentioned in cind, the measure is
		// restricted to the indexes in ind.  For conditions not in cind,
		// the measure is made to be zero everywhere.
		void Restrict(const std::vector<int> &ind, 
				const std::vector<int> &cind);


		// Marginalizes, reindexes, or expands depending on the
		// argument
		// n = # new values  nc = # new conditioning values
		// ind[<i,j>] is a mapping that needs to be applied to old 
		//   conditioning index j and then added to i
		// (see Reindex for RVSimple above to see how a mapping is
		//  applied)
		// (Really, this is probably too confusing to be "simple.")
		typedef std::map<std::pair<int,int>, std::vector<std::vector<int> > > RemapT;
		void Reindex(int n, int nc, const RemapT &ind);

		// Multiplies by another measure (point-wise)
		void MultBy(const RVCondSimple *x);

		// return suitable sufficient statistics object
		SS *BlankSS() const;

		// add independent draw of (cind,ind) (conditioning value, value) to ss
		void AddSS(int cind, int ind, SS *ss, double w=1.0) const;
		// add average of independent draws to ss
		void AddExpSS(int cind, SS *ss, double w=1.0) const;
		void AddExpSS(int cind, const RVSimple *rvs,
				SS *ss, double w) const;

		void AddSS(const SS *toadd, const RVCondSimple* rv,
				const std::vector<std::vector<int> > &mapping,
				int mycondi, int rvccondi,
				SS *ss, double w=1.0) const;

		// returns a sample
		int Sample(int cx, Random &rand = randomizer) const;

		// sets parameters to ML estimate
		void Maximize(const SS *ss);

		// Returns a new object that can be multiplied into this
		// object (as per Mult and RMult)
		RVSimple * Convert(const RVSimple *rvs) const;

		// see above class for description
		void Scramble(double alpha=1.0, double degree=1.0, Random &rand=randomizer); //by Yu
		double LLH(const SS *ss) const;
		double GetScore(double numTrans, const SS *ss)const;
		///////

		private:
		void Advance(LagFactoredVector &v) const;
		//void AdvanceTrans(LagFactoredVector &v) const;
		void Advance(LagFactoredVector &v, int id) const;
		//void AdvanceTrans(LagFactoredVector &v, int id) const;



		void Uniformize(Dynamics *node);

		double minrate;
		double t;
		bool trans;

		//Determines how we find A matrices 
		//0: min rates, 1: avg rates
		int ratedist; 

		//Determines if A is scaled or not 
		int scaled; //0: not scaled, 1:scaled
		std::vector<double> ascale;

		int taylor_expansion_term_count;
		double theta_value;

		//holds nodes by varids
		std::vector<Dynamics *> nodes;

		//for tope
		std::vector<Dynamics *> bnodes;
		std::vector<matrix *> amatrices;

		//typedef std::pair<int,Instantiation> Binfo;
		std::vector<Binfo> Blist;
		
		std::vector<vectr > eigenval;
        std::vector<matrix > Us;
        std::vector<matrix > Vs;

		// cshelton: added
		mutable double cachemaxt;
		mutable std::vector<ExpMCache> rkcache;
		//

		//consistent parent indices of the variables that are not in the evidence
		std::map<int, std::vector <int> > condpars;
		std::map<int, std::vector <int> > children_in_evidence;
		///consistent parent indices of the variables that are in the evidence
		std::map<int, std::vector <int> > condpars_evid;
		std::map<int, int> evid;
		std::map<int, std::vector<int> > trans_evid; //for from to conditioning
		std::map<int, int> notrans_evid; //for evidence that went from uninst to inst (when there is no transition)	

	};


}
#endif
