#ifndef TOPE_H
#define TOPE_H

#include "factoredvector.h"
#include "factoredmatrix.h"
#include <ctime>
#include <queue>
#include <iostream>
#include <map>
#include <iomanip>
// change below when converting to C++11
#include <boost/unordered_set.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>
#define TOPEDEBUG 0

#define MIDPT

//#define DEBUGSTATS

namespace ctbn {

    class ppoly {
        public:
            // creates constant fn
            ppoly();
            // copies or creates integral of
            ppoly(const ppoly &p, bool integrate=false, double scale=1.0);
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

#define midptfrac (0.5)

            double midpt(double endt) const {
                if (regions.empty()) return endt*midptfrac;
                double st = regions.begin()->first;
                std::map<double,coeff>::const_reverse_iterator e = regions.rbegin();
                if (e->second.iszero()) return st+(e->first-st)*midptfrac;
                else return st+(endt-st)*midptfrac;
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


    class TOPETree {
        friend class FactoredMatrix;
        public:
        TOPETree(FactoredVector &fv, FactoredMatrix &fm, double qtm);
        TOPETree(double time, double qtm, double maxruntime,
                FactoredVector &fv, FactoredVector &efv,
                FactoredMatrix &fm, FactoredMatrix &bwfm);
        ~TOPETree();
        vector<vectr> CalcMargins( vector<FactoredVector> & fv);
         vectr kronv(const vectr &v1, const vectr &v2);
    vectr kron(const vector<vectr> &v);

        FactoredVector marginals() const { return ans; }
void InitMyBs(const FactoredMatrix &fm, vector<vector<matrix> > *my_Bs);

        private:
        class TOPENode;

        //static void mydel_poly(ppoly *p) { delete p; };
        //static void mydel_fv(FactoredVector *p) { delete p; };
        //static void mydel_node(TOPENode *p) { delete p; };
        static void mydel_poly(ppoly *p) { };
        static void mydel_fv(FactoredVector *p) { };
        static void mydel_node(TOPENode *p) { };
        typedef boost::shared_ptr<ppoly> polyptr;
        typedef boost::shared_ptr<FactoredVector> fvptr;
        typedef boost::shared_ptr<TOPENode> nodeptr;

        struct propinfo {
            propinfo(FactoredMatrix &fm0, double qt,
                    fvptr zero0) :
                fm(fm0), zero(zero0) {
                    qtm = qt;
                }

            FactoredMatrix &fm;
            fvptr zero;
            //nodeptr root;
            double qtm;
#ifdef DEBUGSTATS
            std::map<int,int> lvlexp,lvlgen;
            std::map<std::pair<int,int>,int> lvlexp2;
#endif
        };

        class TOPENode {
            friend class TOPETree;
            public:
            // generate root node
            TOPENode(propinfo &pi0, fvptr fv);

            const std::vector<nodeptr > &getchildren();
            fvptr val() const { return newuprop; }
            fvptr baseval() const { return uprop; }
            double wt(double t) const { return pbar->eval(t); }
            double pri() const { return mypri; }

            bool iszero() const { return w->iszero(); }

            private:
            // generate general node
            TOPENode(TOPENode *par, polyptr p0, fvptr w0, int Bi0, 
                    fvptr u0, fvptr uprop0,
                    int lvl0, int divlvl0, propinfo &pi0)
                : p(p0), pbar(new ppoly(*p0,true),mydel_poly), w(w0),
                u(u0), newu(new FactoredVector(*w),mydel_fv),
                uprop(uprop0), pi(pi0) {
                    parent = par;
                    Bi = Bi0;
                    lvl = lvl0;
                    divlvl = divlvl0;
#ifdef MIDPT
                    pi.fm.ForAndBack(*newu,p->midpt(pi.qtm),Bi);
#else
                    pi.fm.MultByB(*newu,Bi);
#endif
                    newuprop = fvptr(new FactoredVector(*newu),mydel_fv);
                    pi.fm.Forward(*newuprop,pi.qtm);
                    setpri(pi0.qtm);
                }
            //void setpri(double t) { mypri = wt(t)*(w->MaxofMaxAbs()+val()->MaxofMaxAbs());}
            void setpri(double t) { mypri = wt(t)*(1e-4+val()->MaxofMaxAbs());}
            //void setpri(double t) { mypri = wt(t)*val()->ProdofMaxAbs(false);}

            // represents function
            // \int_0^t p(s)[w*e^(As)*B*e^(-As) - u) ds * e^(At)
            //
            // as approximation pbar(t)*(newu - u)*e^(At)
            //   plus children (plus other children representing next
            //                  lvl in the expansion)
            // where pbar(t) = \int_0^t p(s) ds
            //       newu = w*e^(Ar)*b*e^(-Ar) where is the mid-point
            //                                 in the support of p(r)
            //
            // uprop = u*e^(At)
            // newuprop = newu*e^(At)

            polyptr p,pbar;
            fvptr w,u,newu;
            fvptr uprop,newuprop;
            int Bi;
            int lvl,divlvl;
            std::vector<nodeptr > children;
            TOPENode *parent;
            propinfo &pi;
            double mypri;
            void printall(std::ostream &os) const {
                if (parent) parent->printall(os);
                os << lvl << ',' << divlvl << ':' << Bi << ' ';
            }
        };

        // vvvvv Amir vvvvv
        std::vector<std::vector<bool> > B_I;
        std::vector<std::vector<vector<pair<int,int>  > > > Ks;
        void fill(std::vector<int> &rr, int cur , int last, int sz);
        // ^^^^^ Amir ^^^^^

        struct nodepair {
            nodepair(nodeptr f, nodeptr b)
                : fnode(f), bnode(b) {
                }

            nodeptr fnode,bnode;
            double pri;
            bool operator<(const nodepair &p) const {
                return pri < p.pri;
            }
            bool operator==(const nodepair &p) const {
                return fnode==p.fnode && bnode==p.bnode;
            }
        };

        struct nodepairhash {
            inline std::size_t operator()(const nodepair &np) const {
                std::size_t ret = 0;
                boost::hash_combine(ret, np.fnode.get());
                boost::hash_combine(ret, np.bnode.get());
                return ret;
            }
        };

        // open queue type
        typedef std::priority_queue<nodepair> qtype;
        //typedef std::queue<nodepair> qtype;
        // closed list type
        // change to std once moved to C++11
        typedef boost::unordered_set<nodepair,nodepairhash> cltype;
        FactoredVector ans; // computed marginals

        void addqueue(nodeptr fnode, nodeptr bnode,
                qtype &q, cltype &cl, propinfo &fi, propinfo &bi);

        static void addmap(std::map<int,int> &m, int i, int j) {
            std::map<int,int>::iterator l = m.lower_bound(i);
            if (l!=m.end() && l->first==i) l->second += j;
            else m.insert(l,std::make_pair(i,j));
        }
        static void addmap(std::map<std::pair<int,int>,int> &m, int i1, int i2, int j) {
            std::map<std::pair<int,int>,int>::iterator l = m.lower_bound(std::make_pair(i1,i2));
            if (l!=m.end() && l->first.first==i1 && l->first.second==i2)
                l->second += j;
            else m.insert(l,std::make_pair(std::make_pair(i1,i2),j));
        }

        static void printmap(const std::map<int,int> &m) {
            for(std::map<int,int>::const_iterator i = m.begin();i!=m.end();++i)
                std::cout << i->first << ": " << i->second << std::endl;
        }
        static void printmap(const std::map<std::pair<int,int>,int> &m) {
            int mlvl = 0,mdiv=0;
            for(std::map<std::pair<int,int>,int>::const_iterator i = m.begin();i!=m.end();++i) {
                if (i->first.first > mlvl) mlvl = i->first.first;
                if (i->first.second > mdiv) mdiv = i->first.second;
            }
            for(int i=0;i<=mlvl;i++) {
                for(int j=0;j<=mdiv;j++) {
                    std::map<std::pair<int,int>,int>::const_iterator l = m.find(std::make_pair(i,j));
                    std::cout << std::setw(8) << (l==m.end() ? 0 : l->second);
                }
                std::cout << std::endl;
            }
        }
    };
}

#endif
