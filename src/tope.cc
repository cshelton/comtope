#include <queue>
#include "tope.h"
#include "matrix.h"
//#include "matexp.h"
#include "rk.h"
#include "params.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <sstream>
#include "poly.h"
#include "epoly.h"
#include <boost/tuple/tuple.hpp>
//#define NOUSEEIGEN
#define USEEXPMCACHE

using namespace std;

namespace ctbn {

    ppoly::ppoly() {
        regions.insert(std::make_pair(0.0,coeff()));
    }

    void ppoly::coeff::integrate(double k, double t) {
        int i = c.size();
        if (c[i-1] != 0.0) c.push_back(c[i-1]/i);
        while(i>1) {
            i--;
            c[i] = c[i-1]/i;
        }
        c[0] = 0;
        c[0] = k-eval(t);
    }

    ppoly::ppoly(const ppoly &p, bool integrate, double scale) : regions(p.regions) {
        if (integrate) {
            double k = 0.0;
            std::map<double,coeff>::iterator i0=regions.begin(),i1;
            i1 = i0;
            ++i1;
            while(true) {
                i0->second.integrate(k,i0->first);
                if (i1!=regions.end()) {
                    k = i0->second.eval(i1->first);
                    i0 = i1;
                    ++i1;
                } else break;
            }
        }
        if (scale!=1.0) {
            for (std::map<double,coeff>::iterator i=regions.begin();i!=regions.end();i++) 
                for(uint j=0;j<i->second.c.size();j++) i->second.c[j] *= scale;
        }
    }

    void ppoly::print(std::ostream &os) const {
        if (regions.empty()) os << "(empty)" << std::endl;
        else {
            for(std::map<double,coeff>::const_iterator i=regions.begin();
                    i!=regions.end();++i) {
                os << i->first << ": ";
                for(std::vector<double>::const_iterator j=i->second.c.begin();
                        j!=i->second.c.end();++j) {
                    os << *j << ' ';
                }	
                os << std::endl;
            }
        }
    }

    ppoly::ppoly(const ppoly &p, double cutt, bool cutprev) {
        if (cutprev) {
            std::map<double,coeff>::const_iterator i0=p.regions.upper_bound(cutt);
            if (i0==p.regions.end())
                i0 = p.regions.begin();
            else if (i0!=p.regions.begin()) --i0;
            while(i0!=p.regions.end()) {
                if (i0->first<cutt)
                    regions.insert(make_pair(cutt,i0->second));
                else regions.insert(*i0);
                ++i0;
            }
        } else {
            std::map<double,coeff>::const_iterator i0=p.regions.begin();
            while(i0!=p.regions.end() && i0->first<cutt) {
                regions.insert(*i0);
                ++i0;
            }
            regions.insert(make_pair(cutt,coeff(0.0)));
        }
    }

    ppoly::ppoly(const ppoly &p, double t0, double t1) {
        std::map<double,coeff>::const_iterator i0=p.regions.upper_bound(t0);
        if (i0==p.regions.end())
            i0 = p.regions.begin();
        else if (i0!=p.regions.begin()) --i0;
        while(i0!=p.regions.end() && i0->first<t1) {
            if (i0->first<t0)
                regions.insert(make_pair(t0,i0->second));
            else regions.insert(*i0);
            ++i0;
        }
        regions.insert(make_pair(t1,coeff(0.0)));
        /*
           p.print(cout);
           cout << "[" << t0 << ',' << t1 << "] =>" << endl;
           print(cout);
           cout << "=======" << endl;
           */
    }

    double ppoly::eval(double t) const {
        std::map<double,coeff>::const_iterator i = regions.lower_bound(t);
        if (i==regions.begin()) return 0.0;
        if (i==regions.end())
            return regions.rbegin()->second.eval(t);
        --i;
        return i->second.eval(t);
    }

#define Nsplits 10
    //#define NEWNEWLVL

    const std::vector<TOPETree::nodeptr> &TOPETree::TOPENode::getchildren() {

        if (children.empty()) {
#ifdef DEBUGSTATS
            addmap(pi.lvlexp,lvl,1);
            addmap(pi.lvlexp2,lvl,divlvl,1);
#endif
            double tt = p->t0();
            double dt = (p->t1(pi.qtm)-tt)/(Nsplits);
            for(int i=0;i<(Nsplits);i++) { 
                TOPENode *nn = 
                    new TOPENode(this,
                            polyptr(new ppoly(*p,tt,tt+dt),mydel_poly),
                            w,Bi,newu,newuprop,lvl,divlvl+1,pi);
                if (!nn->iszero())
                    children.push_back(nodeptr(nn,mydel_node));
                tt += dt;
            }
#ifdef DEBUGSTATS
            addmap(pi.lvlgen,lvl,(Nsplits));
#endif

            int nBi = pi.fm.nBs();
            int c = 0;
#ifdef NEWNEWLVL
            for(int newBi=0;newBi<nBi;newBi++) {
                children.push_back(nodeptr(new TOPENode(
                                this,pbar,newu,
                                newBi,pi.zero,pi.zero,lvl+1,0,pi),mydel_node));
                c++;
                if (u.get()!=pi.zero.get()) {
                    children.push_back(nodeptr(
                                new TOPENode(this,
                                    polyptr(new ppoly(*pbar,false,-1.0)),u,
                                    newBi,pi.zero,pi.zero,lvl+1,0,pi),mydel_node));
                    c++;
                }
            }
#else
            if (u.get()==pi.zero.get()) {
                if (!newu->iszero())
                    for(int newBi = 0;newBi < nBi; newBi++)
                        children.push_back(nodeptr(
                                    new TOPENode(this,pbar,newu,
                                        newBi,pi.zero,pi.zero,lvl+1,0,pi),
                                    mydel_node));
            } else {
                const FactoredMatrix::Binfo binfo = pi.fm.getBinfo(Bi);
                for(int newBi = 0;newBi < nBi; newBi++) {
                    if (u->isdiff(*newu,binfo.first)) {
                        children.push_back(nodeptr(
                                    new TOPENode(this,pbar,
                                        fvptr(new FactoredVector(
                                                *newu,*u,binfo.first),mydel_fv),
                                        newBi,pi.zero,pi.zero,lvl+1,0,pi),mydel_node));
                        c++;
                    }
                    for(Context::iimap::const_iterator ii
                            =binfo.second.begin();
                            ii!=binfo.second.end();++ii) {
                        if (u->isdiff(*newu,ii->first)) {
                            children.push_back(nodeptr(
                                        new TOPENode(this,pbar,
                                            fvptr(new FactoredVector(
                                                    *newu,*u,ii->first),mydel_fv),
                                            newBi,pi.zero,pi.zero,lvl+1,0,pi),
                                        mydel_node));
                            c++;
                        }
                    }
                }
            }
#endif
#ifdef DEBUGSTATS
            addmap(pi.lvlgen,lvl+1,c);
#endif
        }
        return children;
    }

    TOPETree::TOPENode::TOPENode(propinfo &pi0, fvptr fv) 
        :  pbar(new ppoly(),mydel_poly), w(fv), u(pi0.zero),
        newu(fv), uprop(pi0.zero),
        newuprop(new FactoredVector(*fv),mydel_fv), pi(pi0) {
            lvl = 0;
            divlvl = 0;
            Bi = -1;
            parent = NULL;

            pi0.fm.Forward(*newuprop,pi0.qtm);

            int nBi = pi.fm.nBs();
            for(int newBi = 0;newBi < nBi; newBi++) {
                TOPENode *nn = new TOPENode(this, pbar,fv,newBi,
                        pi.zero,pi.zero,1,0,pi);
                if (!nn->iszero()) children.push_back(nodeptr(nn,mydel_node));
            }
            setpri(pi0.qtm);

#ifdef DEBUGSTATS
            addmap(pi.lvlgen,0,1);
            addmap(pi.lvlgen,1,nBi);
            addmap(pi.lvlexp,0,1);
            addmap(pi.lvlexp2,0,0,1);
#endif
        }

#ifdef DEBUGSTATS
    static int pairc[10][10];
    static int maxpairc = 0;
#endif

    void TOPETree::addqueue(nodeptr fnode, nodeptr bnode,
            qtype &q, cltype &cl, propinfo &fi, propinfo &bi) {
        nodepair np(fnode,bnode);
        std::pair<cltype::iterator, bool> insans = cl.insert(np);
        if (!insans.second) return; // wasn't inserted b/c is already exists
        fvptr f = fnode->val();
        fvptr fbase = fnode->baseval();
        double fwt = fnode->wt(fi.qtm);
        fvptr b = bnode->val();
        fvptr bbase = bnode->baseval();
        double bwt = bnode->wt(bi.qtm);

        double wt = fwt*bwt;

        //FactoredVector oldans(ans);
        /*
           FactoredVector tmp(ans);
           f->DotStar(*b,tmp);          ans.AddMult(tmp,wt);
           f->DotStar(*bbase,tmp);      ans.AddMult(tmp,-wt);
           fbase->DotStar(*b,tmp);      ans.AddMult(tmp,-wt);
           fbase->DotStar(*bbase,tmp);  ans.AddMult(tmp,wt);
           */
        bool fzero = f->iszero(), fbzero = fbase->iszero();
        bool bzero = b->iszero(), bbzero = bbase->iszero();
        np.pri = 0.0;
        if (!fzero && !bzero) np.pri += ans.AddMultDotStar(*f,*b,wt);
        if (!fzero && !bbzero) np.pri += ans.AddMultDotStar(*f,*bbase,-wt);
        if (!fbzero && !bzero) np.pri += ans.AddMultDotStar(*fbase,*b,-wt);
        if (!fbzero && !bbzero) np.pri += ans.AddMultDotStar(*fbase,*bbase,wt);
        //np.pri = ans.Diff2(oldans)+fnode->pri()*bnode->pri();
        np.pri += fnode->pri()*bnode->pri();
        //np.pri = f->MaxofMaxAbs()*fwt * b->MaxofMaxAbs()*bwt;
        //np.pri = fnode->pri()*bnode->pri();
        //np.pri = 1000000-(fnode->lvl+fnode->divlvl+bnode->lvl+bnode->divlvl);
        //if (ans.Diff2(oldans)<1e-10) {
        //	cerr << "got here!" << endl;
        //}

#ifdef DEBUGSTATS
        if (fnode->lvl>maxpairc) maxpairc = fnode->lvl;
        if (bnode->lvl>maxpairc) maxpairc = bnode->lvl;
        if (maxpairc>9) maxpairc = 9;
        else pairc[fnode->lvl][bnode->lvl]++;
#endif

        /*
           if (np.pri > 0.12 && bnode->divlvl>100) {
           cout << "--------" << endl;
           cout << fnode->lvl << '/' << fnode->divlvl << ',' << bnode->lvl << '/' << bnode->divlvl << ": " << np.pri << endl;
           oldans.Print(cout,ans);
           }
        //static int pn = 0.0;
        //np.pri = --pn;
        */
        if (np.pri>0.0) { 
            /*
               cout << "\tadding:\n";
               cout << "\t\t"; fnode->printall(cout); cout << endl;
               cout << "\t\t"; bnode->printall(cout); cout << endl;
               cout << "\t\tpri = " << np.pri << endl;
               */
            q.push(np);
        }
    }

    vectr TOPETree::kronv(const vectr &v1, const vectr &v2) {
        vectr res(v1.size() * v2.size(), 0.0);
        for (int i = 0; i < (int) v1.size(); i++) {
            for (int j = 0; j < (int)v2.size(); j++)
                res[i*v2.size() + j] = v1[i] * v2[j];
        }
        return res;
    }

    vectr TOPETree::kron(const vector<vectr> &v) {
        vectr res1(8, 0.0);
        for (int i = 0; i < 2 ; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    res1[i*4 + j * 2 + k * 1] = v[0][i&1] * v[1][j&1] * v[2][k&1];
                }
            }
        }
        return res1;
        vectr res(pow(v[0].size(), v.size()), 0.0);//v.size() * v[0].size(),0.0);
        //vectr vv = v[0];
        for (int i = 0; i < (int)v[0].size(); i++)
            res[i] = v[0][i];
        for (int i = 1; i < (int) v.size() ; i++) {
            res = kronv(res, v[i]);//res[i*2 + j] += v1[i] * v2[j];
        }
        //res = vv;
        return res;
    }
    vector<vectr> TOPETree::CalcMargins( vector<FactoredVector> & fv) {
        vector<vectr> ans(fv[0].Size() , vectr(2,0.0));
        cout << "vs: \n";
        for (auto v : fv) {
            for (auto vv : v.dists)
                cout << vv << endl;
            cout << endl;
        }
        for (int i = 0; i < min(16,(int)fv.size()); i++){
            //if(i) fm.Forward(spec_ans[i], qtm);

            for (int j = 0; j < (int)fv[0].Size(); j++) {
                //cout << spec_ans[i].dists[j] << endl;
                //cout << spec_ans[i].dists[j] << endl;
                double dd =  fv[i].SumDistOther(j); 
                // cout << "dd: " << dd << endl;
                // cout << fv[i].dists[j] << endl;
                const vectr t = fv[i].dists[j];
                ans[j] += t *  dd; 
                //if (j > 3) exit(1);
            }
            //cout << ans[i] << endl;
        }
        cout << "normalized: \n";
        /*for (auto v : ans) {
          double dd = v[0] + v[1];
          cout << v / dd << endl << endl;
          }*/
        // if (fv.size() > 3)
        //   exit(1);
        return ans;
    }

    /*
    void TOPETree::KronStatementGenerator(const vector<matrix> &m, const string &&mat_name, const FactoredVector &v) {
        for (int i = 0 ; i < (int)v.dists.size(); i++)
            cout << "v" << to_string(i) << " := transpose(matrix([" <<   v.dists[i][0] << " , " << v.dists[i][1] << "])):\n";
        for (int i = 0; i < (int)m.size(); i++) {
            cout << mat_name << to_string(i) << " := matrix([";
            for (int j = 0; j < 2; j++) {
                if(j) cout << ", ";
                cout << "[";
                for (int k = 0; k < 2; k++) {
                    if (k) cout << ", ";
                    cout << m[i].row(j)[k];
                }
                cout << "]";
            }
            cout << "]):\n";
        }
        string base = "linalg::kroneckerProduct(";
        string ret = "";
        int sz = (int)m.size();
        for (int i = 0 ; i < sz; i++) {
            if (i == (int)sz - 2) {
                ret +=  base +"v" + to_string(i) + "* exp(" +  mat_name + to_string(i) + "* t), " + "v" + to_string(i+1) + "* exp(" + mat_name + to_string(i+1)+ "*t)";
                break;
            }
            else {
                ret += base + "v" + to_string(i) + "* exp(" + mat_name + to_string(i) + "*t), ";
            }
        }
        ret += string(")", sz);
        cout << ret << endl;
    }

    void TOPETree::PrintJointInfo( const vector<vector<matrix> > my_Bs, int nB) {
        stringstream ssa, ssb;
        for (int j = 0 ; j < nB; j++) {
            for (int i = 0 ; i < (int) my_Bs[0].size(); i++) {
                ssb << "Bs{" << j+1 << "," << i+1 << "} = [";
                for (int m = 0; m < 2; m++) {
                    for (int n = 0 ; n < 2; n++) {
                        if (n) ssb << ", ";

                        ssb << my_Bs[j][i][m][n];    
                    }
                    ssb << "; ";
                }
                ssb << "];\n";
            }
        }

        cout << ssb.str() << endl;

        cout << my_Bs[3][2] << endl;
    }
    */

    void TOPETree::InitMyBs(const FactoredMatrix &fm, vector<vector<matrix> > *my_Bs) {
        // cout << "sz : " << sz << endl;
        cout.precision(12);
        matrix eye(2,2, 0.0);
        eye[0][0] = 1; eye[1][1] = 1;
        cout << "eye: " << eye << endl;
        matrix d11 (2,2,0.0);	d11[0][0] = 1;
        matrix d22 (2,2,0.0);	d22[1][1] = 1;
        vector<matrix> dd;
        dd.push_back(d11);
        dd.push_back(d22);
        int sz = (int)fm.amatrices.size();
        int n_vars = (int)fm.amatrices.size();

        int nBi = fm.nBs();
        for (int i = 0; i < nBi; i++) {
            //cout << "cc : \n";
            vector<matrix> m;
            for (int j = 0; j < sz; j++) {
                //cout << "dd \n";
                m.push_back(eye);
            }
            //cout << "cur: " << m.size() << endl;
            //cout << "ee\n";
            (*my_Bs)[i]  = (m);
        }
        //  cout << "m: " << my_Bs->size() << endl;
        //  cout << "n : " << my_Bs[0]->size() << endl;
        int j_nb = 0;
        int jj = 0;
        for (int i = 0 ; i < n_vars; i++) {

            MarkovDyn *node = static_cast<MarkovDyn*> (fm.nodes[i]);
            vector<int> pars = node->CondDomain().VarList();
            if(pars.empty()) { j_nb++; continue;}
            Context pcontext = node->CondDomain();
            for (int k = 0; k < pcontext.Size(); k++) {
                int varid = fm.Blist[j_nb].first;
                const Instantiation &pinst =  fm.Blist[j_nb].second;
                MarkovDyn *node_b =(MarkovDyn*) (fm.bnodes[varid]);
                vector<bool> rr(sz,true);;
                rr[i] = false;
                (*my_Bs)[jj][i] =(*node_b)(pinst.Index())->Intensity();
                for (int p = 0; p < (int)pars.size(); p++) {
                    rr[pars[p]] = false;
                    (*my_Bs)[jj][pars[p]] = dd[pcontext.Index(k).Value(pars[p])];
                }
                j_nb++;
                jj++;
                B_I.push_back(rr);

            }
        }

        //PrintJointInfo(*my_Bs, B_I.size());        
        for (int j = 0 ; j < (int)B_I.size(); j++) {
            for (int i = 0 ; i < n_vars; i++) {
                cout << "j: " << j << " , i: " << i << endl;
                cout << (*my_Bs)[j][i] << endl;
            }
        }
        //        exit(1);
    }

    TOPETree::TOPETree(FactoredVector &fv, FactoredMatrix &fm, double qtm) {

        const static double rkcacheeps = ParamDouble("RKcacheeps",1e-12);
        fm.initRKcache(qtm, rkcacheeps);
        //bwfn.initRKcache(time-qtm, rkcacheeps);
        int nBi = fm.nBs();
        cout << "Nbi: " << nBi << endl;
        int n_vars = (int)fm.amatrices.size();
        vector<matrix> xx;
        for (auto x : fm.amatrices)
            xx.push_back(*x);
        //KronStatementGenerator(xx, "A",  fv);
        //exit(1);
        //cout << "d11: " << d11 << endl << " d22: " << d22 << endl;
        for (int i = 0; i < (int)fm.amatrices.size(); i++) {
            //       cout << "i: " << i << " , #p: " << ((MarkovDyn*)fm.bnodes[i])->CondDomain().Size() << endl;
            vector<int> pars = ((MarkovDyn*)fm.bnodes[i])->CondDomain().VarList();
            //cout << "var: " << i << " , #pars: " << pars.size() << endl;
            for (int k = 0; k < (int)pars.size(); k++) {
                //cout << "var: " << i << " , par: " << ((MarkovDyn*)fm.bnodes[i])->CondDomain().VarList()[k] << endl;
            }
            MarkovDyn *node = (MarkovDyn *)(fm.nodes[i]);
            Context pcontext = node->CondDomain();
            cout << "As{" << i+1 << "}= [ ";// << endl << *(fm.amatrices[i]) << endl;
            for (int m = 0 ; m < 2; m++) {
                for (int n = 0 ; n < 2; n++) {
                    if (n) cout << ", ";
                    cout << (*(fm.amatrices[i]))[m][n];
                }
                cout << ";";
            }
            cout << "];\n";

            //cout << "context size: " << pcontext.Size() << " , pars.size(): " << pars.size() << endl;
            for (int k = 0; k < pcontext.Size(); k++) {
                for (unsigned l = 0; l < pars.size(); l++) {
                    //      cout << "k : " << k << " , l : " << l << " --> ";
                    //    cout << pcontext.Index(k).Value(pars[l])<< endl;

                }
            }
        }

        vector<vector<matrix> > my_Bs(nBi, vector<matrix>());
        InitMyBs(fm, &my_Bs);       
        /*
           for (int j = 0; j < (int)my_Bs.size(); j++) {
        //            cout << "j: " << j << endl;
        for (int i = 0 ; i < (int)my_Bs[0].size(); i++) {
        cout << "j,i: " << j <<" , " << i << endl;
        cout << my_Bs[j][i] << endl;
        }
        }
        cout << "orig b:\n";
        for(int j = 0;j < nBi; j++) {
        int varid = fm.Blist[j].first;
        const Instantiation &pinst = fm.Blist[j].second;
        //cout << "par ins: " << pinst.second << endl;
        //for (auto v: pinst)
        //{
        MarkovDyn *node = (MarkovDyn *)(fm.bnodes[varid]);
        cout << "i: " << varid << endl << (*node)(pinst.Index())->Intensity() << endl;
        //}
        }

        cout << "B_I\n";
        for (int i = 0 ; i < (int)B_I.size(); i++) {
        for (int j = 0; j < (int)B_I[0].size(); j++) {
        cout << B_I[i][j] << " " ;
        }
        cout << endl;
        }
        */
        pair<int, int>  a1(make_pair(0,0));
        pair<int, int>  a2(make_pair(0,1));        
        pair<int, int>  a3(make_pair(1,0));      
        pair<int, int>  a4(make_pair(1,1));
        //        vector<vector<pair<int,int> > > aa = {a1, a2 ,a3 , a4};
        vector<pair<int, int> > aa;
        Ks.resize(5, vector<vector<pair<int, int> > >()); 
        //        for (int i =0 ; i < 4; i++)
        aa.push_back(a1);            aa.push_back(a2);            aa.push_back(a3);            aa.push_back(a4);

        //        Ks[0] = (aa);
        for (int i = 0; i < 4; i++)
            Ks[0].push_back({aa[i]});
        for (int l = 1; l < 4; l++) {
            for (int cnt = 0; cnt < (int)Ks[l-1].size() ;cnt++ ) {
                for (int i = 0; i < 4; i++) {
                    vector<pair<int, int> >  cur =  Ks[l-1][cnt];
                    cur.push_back(aa[i]);
                    //                    cur.push_back(aa.at(i));
                    Ks[l].push_back(cur);
                    //              cout << "l : " << l << " , sz: " << Ks[l].size() << endl;
                }
            }
        }


        // cout << "filling: \n";
        /*
           for (vector<vector<pair<int,int> > > r : Ks) {
           for (auto c : r) {
           for (auto e : c )
           cout << "<" <<e << "> , ";
           cout << endl;
           }
           cout << endl;
           }
           */
        //        exit(1);

        vector<FactoredVector> spec_ans;
        FactoredVector my_v(fv);
        fm.Forward(my_v, qtm);
        for (int i = 0; i < my_v.Size(); i++) {
            cout << my_v.dists[i] << endl;
        }
        spec_ans.push_back(my_v);
        cout << "eigen values:\n";
        for (auto v : fm.eigenval) {
            cout << v << endl;
        }
        typedef tuple<FactoredVector*, vector<EPoly*>, int, int> tp; //int j, int level
        typedef tuple<int, int, int, int> right_; // i, j, k1i, k2i
        map<right_, matrix> mp;
        for (int j = 0; j < (int)B_I.size(); j++) {
            for (int k1 = 0; k1 < 2; k1++) {
                for (int k2 = 0; k2 < 2; k2++) {
                    for (int i = 0 ; i < (int)B_I[0].size(); i++) {
                        matrix mat = fm.Us[i].col(k1) * fm.Vs[i].row(k1) * my_Bs[j][i] * fm.Us[i].col(k2) * fm.Vs[i].row(k2);
                        //printf("j:%d, i:%d, k1:%d, k2: %d\n", j,i,k1,k2);
                        //
                        //cout << mat << endl;
                        mp.emplace(right_(i,j,k1,k2), mat);
                    }
                }
            }
        }
        queue<tp> q;
        EPoly *e =  new EPoly({1}, 0);
        vector<EPoly*> ve = {e};
        for (unsigned int i = 0; i < B_I.size(); i++)
            q.push(tp(new FactoredVector(fv), {e}, i, 1));
        //q.push(tp(fv, {e}, 1));
        int iter = 0;
        int pushed = 0, cnt_prune = 0;
        vector<vector<FactoredVector> > v_level(7, vector<FactoredVector>());
        v_level[0].push_back(my_v);
        cout << "while\n";
        int prev_level = 0;
        vector<vectr> kron_res(5, vectr((int)pow(2, n_vars), 0.0));
        vectr kron_res2((int)pow(2, n_vars), 0.0);
        FactoredVector margs(my_v), lvl1m;
        margs.Print(cout, NULL); //`exit(1);
        //lvl1m.Zero();
        //margs.Zero();
        fvptr zero(new FactoredVector(fv),mydel_fv);
        zero->Zero();
        lvl1m = *zero;
        //margs = *zero;
        //for (int i = 0; i < 2; i++)
        //  for (int j = 0; j < 2; j++)
        //    margs.dists[i][j] = 0.0;
        margs.Print(cout, NULL); //exit(1);
        cout << " lvl1 margs: \n"; lvl1m.Print(cout, NULL);
        int tot_nodes = 0;
        auto start_time = chrono::system_clock::now();
        int alloted = 87500;
        int print_time = 0;
        while (!q.empty() && ++iter < 50 /*(int)B_I.size() + 1*/) {
            //cout << "front while: q.size() : " << q.size() << endl;
            auto dur = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start_time);
            //cout << "time so far: " << dur.count() << endl;
            if (dur.count() >= alloted) break;

            tp node = q.front();

            q.pop();
            int node_level = get<3>(node);
            if (get<3>(node) == prev_level + 1) {
                auto ss = chrono::system_clock::now();
                cout << "iter: " << iter << endl;
                //cout << "size: " << v_level[prev_level].size() << endl;
                //cout << /*"rr.size: " << rr.size() <<*/ " , ans for level: " << prev_level << endl;
                //auto rr = CalcMargins(v_level[prev_level]);
                //cout << "res: " << kron_res[prev_level] << endl;
                cout << "margs: \n"  ; margs.Print(cout,NULL); cout << endl;
                cout << "tot nodes: " << tot_nodes << endl;
                cout << "\n-------------------\n";
                ++prev_level;
                auto ee = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - ss);
                alloted += ee.count();
            }
            vector<int> nIs;
            FactoredVector *mbase = new FactoredVector(*(get<0>(node)));
            int j = get<2>(node);

            //mbase.logf = 0;
            for (int i = 0; i < n_vars; i++) {
                if (B_I[j][i] == 1) {
                    mbase->dists[i] =get<0>(node)->dists[i];
                }
                else{
                    nIs.push_back(i);
                    //cout << "Not I: " << i << endl;
                }
            }
            //tp next_node;
            int cur_lvl = get<3>(node);
            //FactoredVector *cur_node_res = new FactoredVector(*mbase);
            int pushed_nodes = 0;
            for (vector<pair<int, int> > k: Ks[nIs.size() - 1]) {

                // cout << "inside k \n";
                FactoredVector *cur_v = new FactoredVector(*mbase);
                double new_lambda = 0;

                //bool prune = false;
                bool prune = false;
                //for (auto ttt : k)
                //  cout << "< " << ttt << " > , ";
                //cout << endl;
                int idx = 0;
                for (auto i : nIs) {
                    //cout << "to mult: \n";
                    //printf("j:%d, i:%d, k1:%d, k2:%d\n",j,i,k[idx].first,k[idx].second);
                    //cout << "ki: " << k[idx] << endl;
                    //cout << "mp is :\n" << mp[right_(i, j, k[idx].first, k[idx].second )];
                    //cout << "dist " << i << " \n";
                    //cout << cur_v.dists[i] << endl;
                    //cout << cur_v.dists[i] << endl;
                    //printf("C%d := linalg::col(U%d,%d) * linalg::row(V%d,%d) * B%d%d * linalg::col(U%d,%d) * linalg::row(V%d,%d):\n",idx, i+1, k[idx].first+1, i+1, k[idx].first+1, j+1, i+1, i+1, k[idx].second+1,i+1 ,k[idx].second+1);
                    matrix tmp = mp[right_(i, j, k[idx].first, k[idx].second )]  ;//my_v.dists[i];
                    //cout << "1\n";
                    //cout << "cached part\n" << tmp << endl;
                    cur_v->dists[i] = get<0>(node)->dists[i] * mp[right_(i, j, k[idx].first, k[idx].second )]  ;//my_v.dists[i];

                    //cout << "2\n";
                    if (cur_v->dists[i][0] == 0 && cur_v->dists[i][1] == 0) {
                        prune = true;
                        cnt_prune++;
                        break;
                    }
                    // cout << "i: " << i << ", j: " << j << " , k[" << i << "]: " << k[i].first << " , " << k[i].second << endl;
                    //cout << "after mult:\n " << cur_v.dists[i] << endl;
                    //cout << "3\n";

                    new_lambda += ((fm.eigenval[i]))[k[idx].first] - ((fm.eigenval[i]))[k[idx].second];//(*(fm.eigenval[i]))[k1[0]] - (*(fm.eigenval[i]))[k1[1]] ;//get<1>(node).GetExponent();
                    idx++;
                    //cout << "4\n";
                }
                if (prune) {
                    delete cur_v;
                    continue;
                }
                //printf("CC := linalg::kroneckerProduct(C0 * numeric::expMatrix((1/2) * A1), C1 * numeric::expMatrix((1/2) * A2 ))\n");
                //printf("integ_res := int(exp(%f*s), s = 0 .. t)\n", new_lambda);
                //printf("val1 := float(subs(integ_res, t = qtm))\n");
                // printf("val1 * CC\n");
                //cout << "\nlamb: " << new_lambda << endl;
                //exit(1);

                // EPoly *right = new EPoly({ 1 }, new_lambda);
                vector<EPoly*> v_ep;
                double cur_eval = 0;
                for (int integrand_idx = 0; integrand_idx < (int)get<1>(node).size(); integrand_idx++) {
                    if (iter >= 3) {
                        int pk = 0;
                    }
                    EPoly ep(*(get<1>(node))[integrand_idx]);
                    EPoly l_r({1} , new_lambda);
                    ep = ep * l_r;
                    pair<EPoly, double> integ_res = ep.Integral_0_t();
                    EPoly *ep1 = new EPoly(integ_res.first );//* (*right));
                    EPoly *ep2  = new EPoly({-integ_res.second}, 0);
                    cur_eval += ep1->Eval(qtm) + ep2->Eval(qtm);
                    //ep2.MultByConstant(-integ_res.second);
                        if (!ep1->IsZeroPoly())
                        v_ep.push_back(ep1);
                        else delete ep1;
                    if (!ep2->IsZeroPoly())
                        v_ep.push_back(ep2);
                        else delete ep2;

                }
                //cout << "1- " << cur_v.dists[0] << endl;
                //next_node<0> = cur_v;
                for (int new_j = 0; new_j < (int)B_I.size(); new_j++) {
                    tp next_node(cur_v, v_ep, new_j, cur_lvl + 1);
                    q.push(next_node);
                    //cout << "new_j: " << new_j << " , j: " << j << " , pushed: " << ++pushed << " , q.size: " << q.size() <<endl;
                    ++pushed_nodes;

                }
                //cur_eval = 1;
                //                for (int  i = 0 ;  i < (int)cur_v.dists.size(); i++){
                //cout << "eval: " << cur_eval << " , before eval: " << cur_v.dists[i]  << endl; 
                //cout << "before forward: \n";
                //cout << cur_v.dists[0] << endl;
                //cout << cur_v.dists[1] << endl;

                /*for (int p = 0; p < (int)cur_v.dists[0].size(); p++){
                  cur_v.dists[0][p] *= cur_eval;

                  }*/
                //cout << "after eval : " << cur_v.dists[i]  << endl;

                //              }
                //printf("j:%d,  k1:%d, k2:%d,n",j,i,k[idx].first,k[idx].second);
                fm.Forward(*cur_v, qtm);
                //kron_res[cur_lvl] += kron(cur_v.dists) * cur_eval;
                //kron_res2 += kron(cur_v.dists) * cur_eval;

                //cout << "1\n";
                margs.AddMult(*cur_v, cur_eval);
                //printf("val_res{%d ,%d, %d, %d, %d, %d, %d} = %.10f;\n",j+1, k[0].first+1, k[0].second+1, k[1].first+1, k[1].second+1, k[2].first+1, k[2].second+1,cur_eval);

                //if (node_level == 1)
                //    lvl1m.AddMult(cur_v, cur_eval);
                cur_eval  =1 ;
                //printf("c_res{%d ,%d, %d, %d, %d, %d, %d} = [%.8f %.8f; %.8f %.8f; %.8f %.8f];\n", j+1, k[0].first+1, k[0].second+1, k[1].first+1, k[1].second+1, k[2].first+1, k[2].second+1, cur_v.dists[0][0] * cur_eval, cur_v.dists[0][1] * cur_eval, cur_v.dists[1][0] * cur_eval, cur_v.dists[1][1] * cur_eval , cur_v.dists[2][0] * cur_eval , cur_v.dists[2][1] * cur_eval);
                //cout << "2\n";
                //cout << "eval: " << cur_eval << " , after e^at:\n";
                //cout << cur_v.dists[0] << endl;
                //cout << cur_v.dists[1] << endl;
                //cout << "after multy by cur_eval:\n" << kronv(cur_v.dists[0], cur_v.dists[1]) * cur_eval << endl << endl;

                //for (int p = 0; p < (int)cur_v.dists[0].size(); p++){
                //  cur_v.dists[0][p] *= cur_eval;

                //}
                //                cout << cur_v.dists[0] << endl;
                //              cout << cur_v.dists[1] << endl << endl;

                //cur_node_res += cur_v;
                //spec_ans.push_back(cur_v);
                //v_level[cur_lvl].push_back(cur_v);

            } // end for (k ..)
            //cout << "pushed nodes: "  << pushed_nodes << endl;
            tot_nodes += pushed_nodes;
            //spec_ans.push_back(cur_node_res)
            //;
            delete mbase;


        } // end while
        cout << "lvl1 margs: \n"; lvl1m.Print(cout, NULL); cout <<"\n";
        for (auto tt : kron_res)
            cout << tt << endl;
        cout << "iter: " << iter << endl << "qtime: " << qtm << endl;
          cout << "#pruned: " << cnt_prune << endl;
        //    for (auto v : aans) {
        //      double n = v[0] + v[1];
        //    cout << v / n << endl << "-------------------\n";
        // }
        /*    cout << "rr\n";
              for (auto v : rr)
              cout << v << " , " ;
              cout << endl;
              ofstream outp("res.txt");
              outp << "y" <<  " = [ ";
              outp.precision(11);

              for (int i = 0; i < 4; i++) {
              outp << "[ ";
              for (int j = 0; j < (int)rr2[0].size(); j++){
              if (j) outp << " , ";
              outp << rr2[i][j];
              }
              outp << "];\n";
              outp.flush();
              }
              outp << "];\n";
              outp.close();
              */  cout << "margs: \n";// << margs << endl;
        cout << setprecision(12);
        margs.Print(cout, NULL);
        cout << "time: " << chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start_time).count() << endl;
        /*
           double dd = 0, tot = 0;;
           for (int i = 0 ; i < (int)kron_res2.size(); i++) {
           if (i < kron_res2.size() / 2)    dd += kron_res2[i];
           tot += kron_res2[i];
           }
           cout << dd / tot << endl;
           */
        while(!q.empty()) {
            tp tt = q.front();
            q.pop();
            //FactoredVector *p = get<0>(tt);
            //delete p;
            //for(int i = 0 ; i < (int) get<1>(tt).size(); i++)
            //delete get<1>(tt)[i] ;
        }
    }

    void TOPETree::fill(vector<int> &rr, int cur, int last, int sz) {
        /*       if (cur > last) {
                 vector<int> t;
                 for (auto v : rr)
                 t.push_back(v);
        //for (auto v : rr)
        //  t.push_back(v);
        Ks[cur].push_back(t);
        return;
        }
        for (int i = 0; i < 2; i++) {
        rr.push_back( i);
        fill(rr, cur + 1, last, sz);
        rr.pop_back();
        }*/
    }
    TOPETree::TOPETree(double time, double qtm, double maxruntime,
            FactoredVector &fv, FactoredVector &efv,
            FactoredMatrix &fm, FactoredMatrix &bwfn) {
        clock_t startrtime = clock();

        const static double rkcacheeps = ParamDouble("RKcacheeps",1e-12);
        fm.initRKcache(qtm, rkcacheeps);
        bwfn.initRKcache(time-qtm, rkcacheeps);

        fvptr zero(new FactoredVector(fv),mydel_fv);
        zero->Zero();
        ans = *zero;
        propinfo forinfo(fm,qtm,zero), backinfo(bwfn,time-qtm,zero);

        fvptr ffvptr(new FactoredVector(fv),mydel_fv);
        fvptr bfvptr(new FactoredVector(efv),mydel_fv);
        //forinfo.root = nodeptr(new TOPENode(forinfo,ffvptr),mydel_node);
        //backinfo.root = nodeptr(new TOPENode(backinfo,bfvptr),mydel_node);
        nodeptr forroot = nodeptr(new TOPENode(forinfo,ffvptr),mydel_node);
        nodeptr backroot = nodeptr(new TOPENode(backinfo,bfvptr),mydel_node);

        qtype q;
        cltype cl;

#ifdef DEBUGSTATS
        for(int i=0;i<10;i++) for(int j=0;j<10;j++) pairc[i][j] = 0;
#endif

        //addqueue(forinfo.root,backinfo.root,q,cl,forinfo,backinfo);
        addqueue(forroot,backroot,q,cl,forinfo,backinfo);
        int nprocessed = 0;
        while(!q.empty()
                && (double)(clock()-startrtime)<maxruntime*CLOCKS_PER_SEC) {
            nprocessed++;
            nodepair np = q.top(); //q.front();
            /*
               cout << "priority: " << np.pri << endl;
               cout << np.fnode->lvl << '/' << np.fnode->divlvl << ',' << np.bnode->lvl << '/' << np.bnode->divlvl << ": " << np.pri << endl;
               cout << "forward: ";
               np.fnode->printall(cout);
               cout << endl;
               cout << "backward: ";
               np.bnode->printall(cout);
               cout << endl;
               */
            q.pop();
            const vector<nodeptr > fch = np.fnode->getchildren();
            for(uint i=0;i<fch.size();i++)
                addqueue(fch[i],np.bnode,q,cl,forinfo,backinfo);
            const vector<nodeptr > bch = np.bnode->getchildren();
            for(uint i=0;i<bch.size();i++)
                addqueue(np.fnode,bch[i],q,cl,forinfo,backinfo);
            //if (nprocessed == 500) break;
        }

        cout << " processed nodes: " << nprocessed << endl;
#ifdef DEBUGSTATS
        cout << " forward exp by lvl:" << endl;
        printmap(forinfo.lvlexp);
        cout << " forward gen by lvl:" << endl;
        printmap(forinfo.lvlgen);
        cout << " forward exp by lvl-div: " << endl;
        printmap(forinfo.lvlexp2);
        cout << " backward exp by lvl:" << endl;
        printmap(backinfo.lvlexp);
        cout << " backward gen by lvl:" << endl;
        printmap(backinfo.lvlgen);
        cout << " backward exp by lvl-div: " << endl;
        printmap(backinfo.lvlexp2);
        cout << " processed pairs by lvl: " << endl;
        int ctot = 0;
        for(int i=0;i<maxpairc;i++) {
            for(int j=0;j<maxpairc;j++) {
                cout << setw(8) << pairc[i][j];
                ctot += pairc[i][j];
            }
            cout << endl;
        }
        cout << " total pairs = " << ctot << endl;

        if(fm.scaled == 1) ans.Scale(fm.ascale);
        cout << "------" << endl;

        ans.Print(cout,(RV *)(10));
        cout << "------" << endl;
#endif
        ans.Print(cout,NULL);
        // should clean up tree nodes!
    }

    TOPETree::~TOPETree() {
    }
}
