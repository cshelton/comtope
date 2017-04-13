#ifndef TOPE_H
#define TOPE_H

#include "factoredvector.h"
#include "factoredmatrix.h"
#include <ctime>
#include <queue>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/binomial_heap.hpp>
#include <boost/heap/priority_queue.hpp>
#include <boost/unordered_map.hpp>
#include <iostream>
#include <map>

#define TOPEDEBUG 0

#define MERGENODES

namespace ctbn {

struct twotimes {
	twotimes(double tt0,double tt1) : t0(tt0), t1(tt1) {}
    double t0;
    double t1;
};

inline bool operator==(twotimes const& p1, twotimes const& p2)
{
    return p1.t0 == p2.t0 && p1.t1 == p2.t1;
}

struct twotimes_hash 
    : std::unary_function<twotimes, std::size_t>
{
    inline std::size_t operator()(twotimes const& p) const
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.t0);
        boost::hash_combine(seed, p.t1);
        return seed;
    }
};

	class TOPENode {

		friend class TOPETree;

		public:
		TOPENode(const LagFactoredVector &fv, double s1, double s2, const std::vector<double> &coeffs, double calcpri, const LagFactoredVector &endres) 
			: fv(fv), toadd(endres), s1(s1), s2(s2), coeffs(coeffs) {
				//ProdofAbsSum
				//MaxofAbsSum 		
				//ProdofMaxAbs
				//MaxofMaxAbs
				//now done in tope.cc as node is inserted (b/c values are known then)
				pri =  calcpri; // fv.ProdofAbsSum()*pterm;
			};

		/*
		   void setchild(TOPENode* nodeptr){
		   nodeptr->parent = this;
		   children.push_back(nodeptr);
		   };*/

		private:
		LagFactoredVector fv,toadd;
		double s1;
		double s2;
		std::vector<double> coeffs;
		double pri;


	};




	class TOPETree {

		friend class TOPENode;
		friend class FactoredMatrix;

		public:
		clock_t startrtime;
		FactoredVector treeFV;
		vectr *trsum;

		TOPETree(double time, double qtm, double sseps, double maxt,
				FactoredVector &fv, FactoredVector &efv, FactoredMatrix &fm, FactoredMatrix &bwfm) 
			: time(time), qtm(qtm), sseps(sseps), maxrtime(maxt), nexp(1), 
			pq(), epq(), 
			qfvlist(), endfvlist()
		{ trsum = NULL; generate(fv, efv, fm, bwfm); }

		void generate(FactoredVector &fv, FactoredVector &efv, FactoredMatrix &fm, FactoredMatrix &bwfm);

		double MultAddToQ(FactoredMatrix &fm, LagFactoredVector &lfv, double wt, vectr *jointv);
		double MultAddToEnd(FactoredMatrix &bwfm, LagFactoredVector &lfv, double wt, vectr *jointv);
		/*void computevaldf( FactoredMatrix &fm, FactoredVector parentfv,  
		  double evaltime, int nid, int Bi, int level,
		  double s1, double s2, std::vector<double> coeffs);
		  */
		~TOPETree(){
			//clear(root);
			if (trsum) delete trsum;
			while(!pq.empty()){
				TOPENode *tmp = pq.top();
				pq.pop();
				delete tmp;
			}

		}

		struct CompareTOPENode : public std::binary_function<TOPENode*, TOPENode*, bool>                                                                                       
		{  
			bool operator()(const TOPENode* lhs, const TOPENode* rhs) const  
			{  		
				//ProdofAbsSum
				//MaxofAbsSum 		
				//ProdofMaxAbs
				//MaxofMaxAbs
				//return lhs->fv.MaxofAbsSum() < rhs->fv.MaxofAbsSum();  
				//return lhs->fv.ProdofAbsSum()*lhs->pterm < rhs->fv.ProdofAbsSum()*rhs->pterm;  
				return lhs->pri < rhs->pri;
			}  
		}; 

		/*
		   void clear(TOPENode * nodeptr){
		   if(nodeptr == NULL)
		   return;

		   for(uint i=0; i<nodeptr->children.size(); i++){
		   clear(nodeptr->children[i]);
		   }
		   delete nodeptr;
		   }


		   void addnode(TOPENode *parentptr, double stime, double deltat, int nid, int Bind, int level);
		   void generate();
		   void computeval(const vectr &v, const matrix &A, std::vector<matrix> &B, 
		   FactoredVector &fv, FactoredMatrix &fm);

		   std::ostream & print(std::ostream & s);
		   std::ostream & nodeprint(std::ostream & s, TOPENode *nodeptr, int level);
		   */

		private:
		//TOPENode *root;
		double time;
		double qtm; //query  time 
		double sseps;
		double maxrtime;
		int nexp;

		//typedef std::priority_queue<TOPENode*, std::vector<TOPENode*>, CompareTOPENode> qtypebase;
		typedef boost::heap::fibonacci_heap<TOPENode*, boost::heap::compare<CompareTOPENode> > qtypebase;
		//typedef boost::heap::binomial_heap<TOPENode*, boost::heap::compare<CompareTOPENode> > qtypebase;
		//typedef boost::heap::priority_queue<TOPENode*, boost::heap::compare<CompareTOPENode> > qtypebase;
		//
		class qtypeabbr {
		public:
			qtypebase q;

			qtypeabbr(): q() {}

			void push(TOPENode *n, const FactoredMatrix &fm) {
				q.push(n);
			}
			TOPENode *top() const {
				return q.top();
			}

			void pop() {
				q.pop();
			}

			bool empty() const {
				return q.empty();
			}
		};
#ifdef MERGENODES
		class qtypefull {
		public:
			qtypebase q;
			typedef std::pair<TOPENode *,qtypebase::handle_type> ltype;
			//typedef std::multimap<double,
			//	std::pair<TOPENode *,qtypebase::handle_type> > mtype;
			typedef boost::unordered_multimap<twotimes,ltype,twotimes_hash> mtype;
			mtype m;


			qtypefull(): q(), m() {
			}

			void push(TOPENode *n, const FactoredMatrix &fm) {
				const double tres = 0;
				const double vres = 1e-2;
				twotimes trange(n->s1,n->s2);
				std::pair<mtype::iterator,mtype::iterator> r = m.equal_range(trange);
				for(mtype::iterator i = r.first;i!=m.end() && i!=r.second;++i) {
					TOPENode *nmerge = i->second.first;
					double dd = fm.Diff2(nmerge->fv,n->fv);
					if (dd<vres) {
						size_t j;
						for(j=0;j<n->coeffs.size()
									&& j<nmerge->coeffs.size();j++)
							nmerge->coeffs[j] += n->coeffs[j];
						for(;j<n->coeffs.size();j++)
							nmerge->coeffs.push_back(n->coeffs[j]);
						nmerge->pri += n->pri;
						q.increase(i->second.second);
						delete n;
						return;
					}
				}
				qtypebase::handle_type h = q.push(n);
				m.insert(std::make_pair(twotimes(n->s1,n->s2),std::make_pair(n,h)));
			}
			TOPENode *top() const {
				//std::cerr << "top pri = " << q.top()->pri << std::endl;
				return q.top();
			}

			void pop() {
				TOPENode *n = q.top();
				twotimes trange(n->s1,n->s2);
				std::pair<mtype::iterator,mtype::iterator> r = m.equal_range(trange);
				for(mtype::iterator i = r.first;
						i!=m.end() && i!=r.second;++i) 
					if (i->second.first==n) { m.erase(i); break; }
				q.pop();
			}

			bool empty() const {
				return q.empty();
			}
		};
#endif

#ifdef MERGENODES
		typedef qtypefull qtype;
#else
		typedef qtypeabbr qtype;
#endif

		qtype pq,epq;

		std::vector<LagFactoredVector> qfvlist;
		std::vector<LagFactoredVector> endfvlist;

		std::vector<double> qwts;
		std::vector<double> endwts;


		double polyeval(const std::vector<double> &coeffs, double s){
			if(TOPEDEBUG)	std::cout << "Polyeval:\tcoeffs: ";
			double ret = 0;//coeffs[0];
			double tmp = s;
			if(TOPEDEBUG)	std::cout << coeffs[0] << " ";
			for(uint i=1; i<coeffs.size(); i++){     
				if(TOPEDEBUG)	std::cout << coeffs[i] << " ";
				ret += coeffs[i]*tmp;
				tmp *= s;
			}
			if(TOPEDEBUG)	std::cout << "\ts: " << s << "\tret: " << ret << std::endl;

			return ret;
		}

		double polyint(const std::vector<double> &coeffs, double s){
			double tmp = s;
			double ret = 0;
			if(TOPEDEBUG)	std::cout << "Polyint:\tcoeffs: ";
			for(uint i=0; i<coeffs.size(); i++){
				if(TOPEDEBUG)	std::cout << coeffs[i] << " ";
				ret += coeffs[i]/(i+1) * tmp;
				tmp *= s;
			}
			if(TOPEDEBUG)	std::cout << "\ts: " << s << "\tret: " << ret << std::endl;
			return ret;
		}

		void polyup(std::vector<double> &outcoeffs, const std::vector<double> &coeffs, 
				double s1, double s2, double news1, double news2){

			outcoeffs.clear();

			if(TOPEDEBUG){
				std::cout << "Polyup:\ts1: " << s1 << "\ts2: " << s2 << "\tnews1: " 
					<< news1 << "\tnews2: " << news2 << "\tcoeffs: ";

				for(uint i=0; i<coeffs.size(); i++){
					std::cout << coeffs[i] << " ";
				}
				std::cout << std::endl;
				std::cout << "\t";
			}

			if(news2 > s2){       	
				outcoeffs.resize(2);
				outcoeffs[0] = 0;
				outcoeffs[1] = polyeval(coeffs, s2) - polyeval(coeffs, s1);
			} 
			else{
				outcoeffs.resize(coeffs.size()+1);
				outcoeffs[0] = 0;
				for(uint i=0; i<coeffs.size(); i++){
					outcoeffs[i+1] = coeffs[i]/(i+1);
				}    
				outcoeffs[1] -= polyeval(coeffs, s1);
			}

			if(TOPEDEBUG){
				std::cout << "\toutcoeffs: ";
				for(uint i=0; i<outcoeffs.size(); i++){
					std::cout << outcoeffs[i] << " ";
				}
				std::cout << std::endl;
			}
		}

		void addchildren(LagFactoredVector &x1, LagFactoredVector &v1,
				LagFactoredVector &x2, LagFactoredVector &v2,
				double t1,double t2,double tend,const FactoredMatrix &fm,
				int i,const Instantiation &pinst,
				TOPENode *curnode,
				qtype &pq); //,double parstr



	};

}


#endif
