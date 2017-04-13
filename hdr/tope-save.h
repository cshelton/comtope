#ifndef TOPE_H
#define TOPE_H

#include "factoredvector.h"
#include "factoredmatrix.h"
#include <ctime>
#include <queue>
#include <iostream>
#include <map>

#define TOPEDEBUG 0

namespace ctbn {

	class ppoly {
	public:
		// creates constant fn
		ppoly();
		// copies or creates integral of
		ppoly(const ppoly &p, bool integrate=false);
		// truncates at cutt (cutprev=true => p(t<cutt) = 0)
		ppoly(const ppoly &p, double cutt, bool cutprev=false);
		// truncates between t0 and t1
		ppoly(const ppoly &p, double t0, double t1);
		double eval(double t) const;
		double t0() const {
			if (regions.empty()) return 0.0;
			return regions.begin()->first;
		}
		double t1(double endt) const {
			if (regions.empty()) return endt;
			std::map<double,coeff>::const_reverse_iterator e = regions.rbegin();
			if (e->second.iszero()) return e->first;
			else return endt;
		}
			
		double midpt(double endt) const {
			if (regions.empty()) return endt/2.0;
			double st = regions.begin()->first;
			std::map<double,coeff>::const_reverse_iterator e = regions.rbegin();
			if (e->second.iszero()) return (st+e->first)/2.0;
			else return (st+endt)/2.0;
		}
		double dur(double endt) const {
			if (regions.empty()) return endt;
			double st = regions.begin()->first;
			std::map<double,coeff>::const_reverse_iterator e = regions.rbegin();
			if (e->second.iszero()) return e->first-st;
			else return endt-st;
		}
		void print(std::ostream &os) const;
	private:
		struct coeff {
			std::vector<double> c;

			coeff(double k=1) { c.push_back(k); }
			double eval(double t) const {
				double v = 1.0;
				double ret = 0.0;
				for(uint i=0;i<c.size();i++) {
					ret += v*c[i];
					v *= t;
				}
				return ret;
			}
			bool iszero() const {
				return c.empty() || (c.size()==1 && c[0]==0.0);
			}
			void integrate(double k, double t);
		};
		std::map<double,coeff> regions;
	};

	class TOPENode {
		friend class TOPETree;
		public:
		TOPENode(const ppoly &p0, const FactoredVector &w0, int Bi0, 
				const FactoredVector &u0, int lvl0,
				const FactoredMatrix &fm, double myqtime)
					: p(p0), pbar(p0,true), w(w0), u(u0), newu(w) {
			Bi = Bi0;
			// to be filled in later
			//pri = pri0; //w0.ProdofAbsSum();
			lvl = lvl0;
			fm.ForAndBack(newu,p.midpt(myqtime),Bi);
		}
		private:
		ppoly p,pbar;
		FactoredVector w,u,newu;
		int Bi;
		double pri;
		int lvl;
	};

	class TOPETree {

		friend class TOPENode;
		friend class FactoredMatrix;

		public:

		TOPETree(double time, double qtm, double maxruntime,
				FactoredVector &fv, FactoredVector &efv,
				FactoredMatrix &fm, FactoredMatrix &bwfm) 
			: time(time), qtm(qtm), maxrtime(maxruntime), nexp(0) {
			generate(fv, efv, fm, bwfm);
		}

		void generate(FactoredVector &fv, FactoredVector &efv, FactoredMatrix &fm, FactoredMatrix &bwfm);

		~TOPETree(){
			while(!fqueue.empty()){
				TOPENode *tmp = fqueue.top();
				fqueue.pop();
				delete tmp;
			}
			while(!bqueue.empty()){
				TOPENode *tmp = bqueue.top();
				bqueue.pop();
				delete tmp;
			}
		}

		struct CompareTOPENode : public std::binary_function<TOPENode*, TOPENode*, bool> {
			bool operator()(const TOPENode* lhs, const TOPENode* rhs) const {  		
				return lhs->pri < rhs->pri;
			}  
		}; 

		FactoredVector marginals() const { return ans; }

		private:
		clock_t startrtime;
		FactoredVector ans; // computed marginals
		double time; // total time
		double qtm; //query time 
		double maxrtime; // max running time
		int nexp; // num expanded
		std::map<int,int> lvlexp;
		std::map<int,int> lvlpro;

		static void addmap(std::map<int,int> &m, int i, int j) {
			std::map<int,int>::iterator l = m.lower_bound(i);
			if (l==m.end()) m.insert(std::make_pair(i,j));
			if (l->first==i) l->second += j;
			else m.insert(l,std::make_pair(i,j));
		}

		static void printmap(const std::map<int,int> &m) {
			for(std::map<int,int>::const_iterator i = m.begin();i!=m.end();++i)
				std::cout << i->first << ": " << i->second << std::endl;
		}

		typedef std::priority_queue<TOPENode*, std::vector<TOPENode*>, CompareTOPENode> qtype;

		qtype fqueue,bqueue; // forward and backward queues

		// list of vectors and weights
		typedef std::vector<std::pair<FactoredVector,double> > anstype;

		anstype forlist, backlist;

		double processnode(TOPENode *node, anstype &mylist, anstype &opplist,
				FactoredMatrix &fm, qtype &myqueue, double myqtime,
				const FactoredVector &zero);

		double addnode(TOPENode *node, anstype &mylist, anstype &opplist,
				FactoredMatrix &fm, qtype &myqueue, double myqtime);
	};
}

#endif
