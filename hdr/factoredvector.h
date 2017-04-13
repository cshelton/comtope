/*
 *
 */

#ifndef CTBNRLE_FACTORED_VECTOR_H
#define CTBNRLE_FACTORED_VECTOR_H

//#include "ctbndyn.h"
#include "bn.h"
#include "rvsimple.h"

#include "matrix.h"
#include "ctbn.h"
#include "ctbndyn.h"
#include "context.h"
#include "multirv.h"
#include "notyetimplementederror.h"
#include "nullptr03.h"

#include <string>
#include <sstream>

#include <iostream>
#include <vector>
#include <map>


namespace ctbn {

	class FactoredVector : public RVSimple {

		public:
			FactoredVector();
			//FactoredVector(const FactoredVector & v);
			FactoredVector(RV *rv);
			// constructs ith term in difference btwn v1 and v2
			FactoredVector(const FactoredVector &v1,
				const FactoredVector &v2, int i);
			

			//~FactoredVector();
            friend class TOPETree;
			//copy constructor
			FactoredVector *Clone() const;

			void Swap(FactoredVector &v);

			void makenonneg() {
				for(unsigned int i=0;i<dists.size();i++)
					for(int j=0;j<dists[i].getm();j++)
						if (dists[i][j]<0) dists[i][j] = 0.0;
			}

			// constructs ith term in difference btwn v1 and v2
			void makediff(const FactoredVector &v1,
				const FactoredVector &v2, int i);

			bool isdiff(const FactoredVector &v1, int i) const;

			void Restrict(const std::vector<int>&);
			void RestrictById(int id, int val);
			void RestrictById(int id, int val, double x);
			void Condition(const std::map<int, int> &ev);
			double Normalize(); //(std::vector<double> &normf);
			double absnormalize(vectr &v);
			double maxnormalize(vectr &v);
			double NormalizeNeg();
			int nDist() const { return dists.size(); }
			void GetDist(vectr &v, int varid, double &logfactor);
			// returns ref to dist (without notion of logf)
			vectr &Dist(int varid) { return dists[varid]; }
			const vectr &Dist(int varid) const { return dists[varid]; }
			// cshelton: made first arg const below
			void SetDist(const vectr &v, int varid, double logfactor=0.0);
			void SetLogf(double logfactor);
			double GetLogf() const;


			double GetDistInstNorm(int varid, int index);
			double GetDistInst(int varid, int index);
			void SetDistInst(int varid, int index, double value);
			int GetDistSize(int varid);
			vectr Joint(const FactoredVector &fv);
			int Size();
			int JointSize();
			std::ostream & Print(std::ostream & out, RV * const) const;
			std::ostream & Print(std::ostream & out, const FactoredVector &v2) const;
			std::ostream & PrintSimple(std::ostream & out) const;

			double CompareDist(FactoredVector &fv, int varid) const;

			double GetMargMin();
			double GetJointMin();
			double GetMin(int varid);
			double GetJointMax();
			double GetJointAbsMax();
			double GetMargMax();
			double GetMargAbsMax();
			double GetMax(int varid);
			bool iszero() const;

			void Load(std::istream &is);
			void Save(std::ostream &os) const;
			RVSimple *Condition(int ind) const;
			void GetDist(vectr &d, double &logfactor) const;
			void SetDist(const vectr &d, double logfactor=0.0);

			// returns the probability of an individual event (index)
			double Prob(int ind, bool log=false) const;

			// Returns the absolute sum (or log thereof)
			double Sum(bool log=false) const{
				return ProdofAbsSum(log);
			};

			double ProdofAbsSum(bool log=false) const;

			double MaxofAbsSum(bool log=false) const;

			// Returns the product of the max of absolute sum
			double ProdofMaxAbs(bool log = false) const; 

			double MaxofMaxAbs(bool log = false) const; 

			// Marginalizes, reindexes, or expands depending on the argument
			// ind is a list, for each new index, of the list of old indexes
			// which should be summed
			void Reindex(const std::vector<std::vector<int> > &ind);

			// Multiplies by another measure (point-wise)
			void MultBy(const RVSimple *x);

			void DotStar(const FactoredVector &v, FactoredVector &out) const;

			// returns "amount of change"
			double AddMult(const FactoredVector & v, double x);

			// list AddMult, but for "dot-star"
			double AddMultDotStar(const FactoredVector &v1,
					const FactoredVector &v2, double x);
			// Adds in another measure
			void Add(const RVSimple *x);
			void Add(const RVSimple *x, double w);
            double SumDist(int i) {
                double ret = 0;
                for (int k = 0; k < (int) dists[i].size(); k++)
                    ret += dists[i][k];
                return ret;
            }
            double SumDistOther(int i) {
                double ret = 0;
                for (int k = 0; k < (int) dists.size(); k++) {
                    if (k == i) continue;
                    ret += SumDist(k);
                }
                return ret;
            }

			// Multiplies measure by constant
			void Mult(double x);
			
			void Zero();

			//Scale by a vector
			void Scale(std::vector<double> &sc);

			// return suitable sufficient statistics object
			SS *BlankSS() const;

			// add independent draw of x to ss
			void AddSS(int x, SS *ss, double w=1.0) const;
			// add average of independent draws of x to ss
			void AddExpSS(SS *ss, double w=1.0) const;

			void AddSS(const SS *toadd, const RVSimple* rvs,
					const std::vector<std::vector<int> > &mapping, 
					SS *ss, double w=1.0) const;

			// returns a sample
			int Sample(Random &rand = randomizer) const;

			// sets parameters to ML estimate
			void Maximize(const SS *ss);

			// set measure to be 1 everywhere it has support
			void MakeUniform();

			// Draws parameters from Dirichlet distribution
			// with all alphas the same (equal to alpha)
			// If degree equals 1, all parameters are replaced
			// If degree is between 0 and 1, the old parameters are
			// mixed with the new one (linear combination, degree controls
			// the mixing).  
			void Scramble(double alpha=1.0, double degree=1.0, 
					Random &rand=randomizer); 
			// Returns the log-likelihood of a set of sufficient statistics
			double LLH(const SS *ss) const;

			double GetScore(double numTrans, const SS* ss) const;

			double Diff(const FactoredVector &v) const;
			double Diff2(const FactoredVector &v) const;

			double PairSize(const FactoredVector &v) const;
            void operator += (const FactoredVector &v) {
                for (unsigned int i = 0; i < dists.size(); i++) {
                    for (unsigned p = 0; p < dists[i].size(); p++) {
                        dists[i][p] += v.dists[i][p];
                    }
                }

            }

		protected:
			std::vector<vectr> dists;

			double logf;

	};

	/*
	   inline FactoredVector operator*(const FactoredVector &v, const vector<double> &s) {
	   FactoredVector newv = v;
	   newv *= s;
	   return newv;
	   }

	   inline  FactoredVector operator*(const vector<double> &s, const FactoredVector &v) {
	   FactoredVector newv = v;
	   newv *= s;
	   return newv;
	   }

	   inline FactoredVector operator*(const FactoredVector &v, const double &s) {
	   FactoredVector newv = v;
	   newv *= s;
	   return newv;
	   }

	   inline FactoredVector operator*(const double &s, const FactoredVector &v) {
	   FactoredVector newv = v;
	   newv *= s;
	   return newv;
	   }
	   */
}
#endif
