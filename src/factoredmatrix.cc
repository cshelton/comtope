#include "factoredmatrix.h"
#include "factoredvector.h"
#include "bn.h"
#include "ctbn.h"
#include "context.h"
#include "dynamics.h"
#include "matrix.h"
#include "notyetimplementederror.h"
#include "nullptr03.h"
#include "params.h"
#include "rk.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace ctbn {

    using namespace std;


    static const string the_class_name("FactoredMatrix");

    FactoredMatrix::FactoredMatrix() :
        RVCondSimple(),
        minrate(0.0),
        t(-1.0),
        trans(false),
        ratedist(0),
        scaled(0),
        taylor_expansion_term_count(10),
        theta_value(10),
        cachemaxt(-1.0) {}

    FactoredMatrix::FactoredMatrix(const FactoredMatrix & copy_from) :
        RVCondSimple(),
        minrate(copy_from.minrate),
        t(copy_from.t),
        trans(copy_from.trans),
        ratedist(copy_from.ratedist),
        scaled(copy_from.scaled),
        taylor_expansion_term_count(copy_from.taylor_expansion_term_count),
        theta_value(copy_from.theta_value),
        nodes(),
        Blist(copy_from.Blist),
        cachemaxt(copy_from.cachemaxt),
        rkcache(copy_from.rkcache) {
            unsigned int const nodes_size(copy_from.nodes.size());
            nodes.resize(nodes_size);
            ascale.resize(nodes_size);
            for (unsigned int i(0); i != nodes_size; ++i) {
                // TODO: Ask Busra: should this actually be cloning the elements?
                nodes[i] = copy_from.nodes[i];
                ascale[i] = copy_from.ascale[i];
            }


            //copy maps as well 
            if (copy_from.condpars.size() > 0)
                condpars = copy_from.condpars;
            if (copy_from.children_in_evidence.size() > 0)
                children_in_evidence = copy_from.children_in_evidence;
            if (copy_from.condpars_evid.size() > 0)
                condpars_evid = copy_from.condpars_evid;
            if (copy_from.evid.size() > 0)
                evid = copy_from.evid;
            if (copy_from.trans_evid.size()>0)
                trans_evid = copy_from.trans_evid;
            if (copy_from.notrans_evid.size()>0)
                notrans_evid = copy_from.notrans_evid;

        }


    FactoredMatrix::FactoredMatrix(Dynamics * dyn, FactoredVector *fv, int rdist, int scd, bool transpose) :
        RVCondSimple(),
        minrate(0.0),
        t(-1.0),
        trans(false),
        ratedist(rdist),
        scaled(scd),
        taylor_expansion_term_count(10),
        theta_value(10),
        cachemaxt(-1.0) {
            CTBNDyn *ctbndyn = dynamic_cast<CTBNDyn*>(dyn);
            int n = ctbndyn->NumofNodes(); 
            nodes.clear();
            nodes.resize(n);

            bnodes.clear();
            bnodes.resize(n);

            amatrices.clear();
            amatrices.resize(n);

            ascale.clear();
            ascale.resize(n);

            eigenval.clear();
            eigenval.resize(n);

            Us.clear();
            Us.resize(n);

            Vs.clear();
            Vs.resize(n);
            for (int i=0; i<n; i++){


                nodes[i] = ctbndyn->NodeByVar(i)->Clone();
                minrate += FindMinVal(nodes[i]);		


                bnodes[i] = ctbndyn->NodeByVar(i)->Clone();


                if(transpose){
                    int cds = ((MarkovDyn*)bnodes[i])->CondDomain().Size();
                    for(int j=0; j<cds; j++){
                        ((MarkovDyn*)bnodes[i])->operator()(((MarkovDyn*)bnodes[i])-> \
                            CondDomain().Index(j))->Intensity() =
                            ((MarkovDyn*)nodes[i])->operator()(((MarkovDyn*)nodes[i])-> \
                                CondDomain().Index(j))->Intensity().transpose();
                    }
                }

            }
            for (int i=0; i<n; i++){
                ascale[i] = 0.0;
                //cout << "Calling MakeAB for node: " << i << endl;
                amatrices[i] = MakeAB(bnodes[i], fv, ascale[i]);	 
                //cout << "ascale " << ascale[i] << endl;
                EigenDecompose ev(*amatrices[i]);
                eigenval[i] =   (ev.getEigenValues().real());
                Us[i] = (ev.getRightEigenVectors().real());
                Vs[i] = (ev.getLeftEigenVectors().real().transpose());
//		Vs[i] = new matrix(Us[i]->transpose());
		cout << "right: \n" << Us[i] << endl;
		cout << "left: \n" << Vs[i] << endl;
		cout << "row : \n" << Vs[i].row(1) << endl;

            }
            /*
               for (int i=0; i<n; i++){
               if(transpose){
               int cds = ((MarkovDyn*)bnodes[i])->CondDomain().Size();
               for(int j=0; j<cds; j++){
               cout << "Node " << i << " component " << j << endl;
               ((MarkovDyn*)bnodes[i])->operator()(((MarkovDyn*)bnodes[i])-> \
               CondDomain().Index(j))->Intensity().niceprint(cout);

               matrix m =   ((MarkovDyn*)bnodes[i])->operator()(((MarkovDyn*)bnodes[i])-> \
               CondDomain().Index(j))->Intensity().transpose();

               ((MarkovDyn*)bnodes[i])->operator()(((MarkovDyn*)bnodes[i])-> \
               CondDomain().Index(j))->Intensity() = m;


               cout << "AFter " << endl;
               ((MarkovDyn*)bnodes[i])->operator()(((MarkovDyn*)bnodes[i])-> \
               CondDomain().Index(j))->Intensity().niceprint(cout);
               }


               matrix m = *amatrices[i];
             *amatrices[i] = m.transpose();
             }

             }*/
            for(uint i=0;i<nodes.size();i++)
                forallassign(Instantiation &pinst, nodes[i]->CondDomain())
                    Blist.push_back(std::make_pair(i,pinst));
        }


    FactoredMatrix::~FactoredMatrix(){
        unsigned int const nodes_size(nodes.size());
        for (unsigned int i=0; i != nodes_size; ++i){
            if (nodes[i] != nullptr03) {
                delete nodes[i];
                delete bnodes[i];
                delete amatrices[i];
            }
        }
    }

    int FactoredMatrix::nBs() const { 
        return Blist.size();
    }

    RVCondSimple* FactoredMatrix::Clone() const{
        return new FactoredMatrix(*this);
    }

    void FactoredMatrix::Uniformize(Dynamics *node){
        int domainsize = ((MarkovDyn*)node)->Domain().Size();
        //	matrix eye(domainsize, domainsize, vectr(domainsize, 1.0));
        matrix eye(domainsize, domainsize, 0.0);
        for (int i(0); i != domainsize; ++i) {
            eye[i][i] = 1.0;
        };

        for (int i=0; i<((MarkovDyn*)node)->CondDomain().Size(); i++){
            matrix &m = ((MarkovDyn*)node)->operator()(((MarkovDyn*)node)->CondDomain().Index(i))->Intensity();
            m += eye;
            m = m/minrate;
        }
    }


    double FactoredMatrix::FindMinVal(Dynamics *node) const{
        double min = 1.0;
        double tmp = 1.0;
        for (int i=0; i< ((MarkovDyn*)node)->CondDomain().Size(); i++){
            matrix &m =  ((MarkovDyn*)node)->operator()( ((MarkovDyn*)node)->CondDomain().Index(i))->Intensity();
            tmp = m.min();
            if ( min >= tmp)		
                min = tmp;
        }
        return -min;
    }

    double FactoredMatrix::FindMinRate() const{
        double ret=0.0;
        for (unsigned int i=0; i<nodes.size(); i++){
            ret += FindMinVal(nodes[i]);
        }
        return ret;
    }



    void FactoredMatrix::Cond(double t0, double t1, const Instantiation &x){

        t = t1 - t0;	
        trans = false;

        evid.clear();
        condpars_evid.clear();
        condpars.clear();
        children_in_evidence.clear();		

        vector<int> cind;
        std::map<int, int>::iterator it;
        std::map<int, vector<int> >::iterator itv;	
        for (unsigned int i=0; i<nodes.size(); i++){
            cind.clear();		
            (nodes[i]->Domain()).ConsistentIndexes(cind, x);


            if (cind.size()==1){
                it = evid.begin();
                evid.insert(it, pair<int,int>(i, cind[0]));
                (nodes[i]->CondDomain()).ConsistentIndexes(condpars_evid[i],x);

            }
            //don't include observed parents
            else {
                (nodes[i]->CondDomain()).ConsistentIndexes( condpars[i], x);						
            }

        }
        int id;
        for (it=evid.begin(); it!=evid.end(); it++){
            id = (*it).first;
            MarkovDyn *evnode = dynamic_cast<MarkovDyn*>(nodes[id]);
            vector<int> parlist = evnode->CondDomain().VarList();
            for (unsigned int k=0; k<parlist.size(); k++){
                if (condpars.find(parlist.at(k)) != condpars.end()){
                    children_in_evidence[parlist.at(k)].push_back(id);
                }
            }
        }

    }

    void FactoredMatrix::Cond(double t0, const Instantiation &from, const Instantiation &to, bool transition){
        t = t0;
        trans = true;
        if(transition){
            for (unsigned int i=0; i<nodes.size(); i++){
                vector<int> fromind;
                nodes.at(i)->Domain().ConsistentIndexes(fromind,from);
                vector<int> toind;
                nodes.at(i)->Domain().ConsistentIndexes(toind,to);

                int fromsize = fromind.size();
                int tosize = toind.size();

                if (fromsize == 1 && tosize == 1){
                    if (fromind.at(0) != toind.at(0)){
                        trans_evid[i].push_back(fromind.at(0));
                        trans_evid[i].push_back(toind.at(0));
                        nodes.at(i)->CondDomain().ConsistentIndexes(condpars_evid[i], from);
                    }
                }
            }
        }
        else{
            //No transition
            for (unsigned int i=0; i<nodes.size(); i++){

                vector<int> fromind;
                nodes.at(i)->Domain().ConsistentIndexes(fromind,from);
                vector<int> toind;
                nodes.at(i)->Domain().ConsistentIndexes(toind,to);

                int fromsize = fromind.size();
                int tosize = toind.size();
                if (fromsize > 1 && tosize == 1){
                    notrans_evid[i] = toind.at(0);
                }

            }	
        }
    }

    double FactoredMatrix::FindAvgRate(FactoredVector &v, const Dynamics *node, int from, int to,
            const vector<int> &cind, bool transpose, double &maxval) const{
        double par_sum=0.0;
        double tmp_mult;
        MarkovDyn const *nnode = dynamic_cast<const MarkovDyn*>(node);
        vector<int> parlist = nnode->CondDomain().VarList();

        maxval = 0.0;

        for (unsigned int j=0; j<cind.size(); j++){ 
            //for each consistent parent instantiation
            tmp_mult = 1.0;			
            //mult all v(u) according to inst
            for (unsigned int k=0; k<parlist.size(); k++){
                tmp_mult *= v.GetDistInst(parlist.at(k), nnode->CondDomain().Index(cind.at(j)).Value(parlist.at(k)));
            }
            if (transpose)
                tmp_mult *= nnode->operator()(nnode->CondDomain().Index(cind.at(j)))->Intensity()[to][from];
            else
                tmp_mult *= nnode->operator()(nnode->CondDomain().Index(cind.at(j)))->Intensity()[from][to];

            if (maxval < fabs(tmp_mult))  maxval = tmp_mult;

            par_sum += tmp_mult;
        }
        return par_sum;
    }

    //rewrites on v
    void FactoredMatrix::MultByVec(FactoredVector & v, bool transpose, double &maxval) const{

        FactoredVector newv(v);

        maxval = 0.0;

        double tmpsum = 0.0;
        int numofnodes = nodes.size();
        for (int id=0; id<numofnodes; id++){
            MarkovDyn *node = (MarkovDyn *)(nodes[id]);
            Context pcontext = node->CondDomain();
            vector<int> parlist = pcontext.VarList();
            double node_sum = 0.0;
            double tmp_mult = 1.0;
            double par_sum = 0.0;
            for (int xout=0; xout<node->Domain().Size(); xout++){
                node_sum = 0.0;
                for (int xin=0; xin<node->Domain().Size(); xin++){
                    par_sum = 0;
                    for (int k=0; k<pcontext.Size(); k++){
                        tmp_mult = 1.0;
                        for (unsigned int l=0; l<parlist.size(); l++) {
                            tmp_mult *= v.GetDistInst( parlist.at(l), 
                                    pcontext.Index(k).Value(parlist.at(l)) ); 
                        }
                        if (transpose)
                            tmpsum = tmp_mult*(dynamic_cast<MarkovDyn *>(node)->operator()(pcontext.Index(k))->Intensity()[xout][xin]);	
                        else					
                            tmpsum = tmp_mult*(dynamic_cast<MarkovDyn *>(node)->operator()(pcontext.Index(k))->Intensity()[xin][xout]);

                        if (maxval < fabs(tmpsum) ) maxval = tmpsum;		
                        par_sum += tmpsum;

                    }
                    par_sum = par_sum / minrate;
                    if (xin == xout){ 
                        par_sum += 1.0;
                    }		

                    node_sum += v.GetDistInst(id, xin)*par_sum;	

                }

                newv.SetDistInst(id, xout, node_sum);
            }
        }

        v = newv;

    }



    void FactoredMatrix::TransMultByVec(FactoredVector & v, bool transpose) const{
        FactoredVector newv(v);
        map<int, vector<int> >::const_iterator it;	
        int id, from, to;
        double par_sum, node_sum;
        MarkovDyn const *node;
        if (trans_evid.size() > 0){
            for(it=trans_evid.begin(); it!=trans_evid.end(); it++){
                id = (*it).first;
                from = ((*it).second).at(0);
                to = ((*it).second).at(1);
                node = dynamic_cast<MarkovDyn*>(nodes[id]);
                const vector<int> cind = (*condpars_evid.find(id)).second;
                //par_sum = FindAvgRate(v, node, from, to, transpose);
                double tmp_mult = 1.0;
                double tmp_pardist;
                par_sum = 0.0;
                vector<int> parlist = node->CondDomain().VarList();
                int parid, parval;

                //Initialize parent distributions to zero
                for (unsigned int k=0; k<parlist.size(); k++){
                    parid = parlist.at(k);
                    vectr tmpv(v.GetDistSize(parid), 0.0);
                    newv.SetDist(tmpv, parid, 0.0);
                }

                for(unsigned int i=0; i<cind.size(); i++){
                    tmp_mult = 1.0;	
                    for (unsigned int k=0; k<parlist.size(); k++){
                        parid = parlist.at(k);
                        parval = node->CondDomain().Index(cind.at(i)).Value(parid);
                        tmp_mult *= v.GetDistInst(parid, parval);
                    }
                    if (transpose)
                        tmp_mult *= node->operator()(node->CondDomain().Index(cind.at(i)))->Intensity()[to][from];
                    else
                        tmp_mult *= node->operator()(node->CondDomain().Index(cind.at(i)))->Intensity()[from][to];

                    par_sum += tmp_mult;

                    for (unsigned int k=0; k < parlist.size(); k++){
                        parid = parlist.at(k);
                        parval = node->CondDomain().Index(cind.at(i)).Value(parid);
                        tmp_pardist = newv.GetDistInst(parid, parval);	
                        newv.SetDistInst(parid, parval, (tmp_pardist+tmp_mult));				
                    }			
                }	
                if (transpose){
                    int tmp = from;
                    from = to;
                    to = tmp;
                }
                node_sum = par_sum*v.GetDistInst(id, from);
                newv.RestrictById(id, to, node_sum);
            }
        }
        v = newv;
    }


    void FactoredMatrix::CondMultByVec(FactoredVector & v, bool transpose, double &maxval) const{
        FactoredVector newv(v);
        vector <double> leakp;
        leakp.resize(nodes.size());
        vector<int> cind;
        vector <int> cpind;
        map<int,int>::const_iterator it;
        int id, val;
        double tmp_mult;
        double leaksum = 0.0;
        std::vector<int> parinst;

        maxval = 0.0;

        for (it=evid.begin(); it!=evid.end(); it++){
            id = (*it).first;
            val = (*it).second;

            MarkovDyn *evnode = dynamic_cast<MarkovDyn*>(nodes[id]);
            vector<int> parlist = evnode->CondDomain().VarList();
            leakp[id] = 0.0;
            tmp_mult=1.0;
            parinst.clear();			
            parinst = condpars_evid.find(id)->second;

            for (unsigned int j=0; j<parinst.size(); j++){ //for each consistent parent instantiation
                tmp_mult = evnode->operator()(evnode->CondDomain().Index(parinst.at(j)))->Intensity()[val][val];

                for (unsigned int k=0; k<parlist.size(); k++){
                    if (condpars.find(parlist.at(k)) != condpars.end()){
                        tmp_mult *= v.GetDistInst(parlist.at(k), 
                                evnode->CondDomain().Index(parinst.at(j)).Value(parlist.at(k)));																									
                    }
                }

                leakp[id] += tmp_mult;
            }
            leakp[id] = leakp[id]/minrate;	
            leaksum += leakp[id];				
        }

        map<int, vector<int> >::const_iterator itv;

        for (itv=condpars.begin(); itv!=condpars.end(); itv++){
            id = (*itv).first;
            std::vector<int> condparlist;
            std::vector<int> chie;
            vector<int> parlist;
            vector<int>parinst;
            condparlist = (*itv).second;
            MarkovDyn *nnode = dynamic_cast<MarkovDyn*>(nodes[id]);
            double par_sum, node_sum;

            for (int xout=0; xout<nnode->Domain().Size(); xout++){
                node_sum=0.0;

                for (int xin=0; xin<nnode->Domain().Size(); xin++){
                    if (v.GetDistInst(id, xin) > 0.0){
                        par_sum = FindAvgRate(v, nnode, xin, xout, condparlist, transpose, maxval);
                        par_sum /= minrate;			
                        if (xin == xout){

                            double chsum = 0.0;
                            double depsum = 0.0;
                            double depmult = 1.0;

                            par_sum = par_sum + leaksum + 1.0;
                            chie.clear();

                            map<int, vector<int> >::const_iterator tmpit;
                            tmpit = children_in_evidence.find(id); 

                            if (tmpit != children_in_evidence.end())
                                chie = tmpit->second;

                            if (chie.size() >0 ){
                                chsum = 0.0;
                                depsum = 0.0;
                                for (unsigned int cc=0; cc<chie.size(); cc++){
                                    int ch_id = chie.at(cc);

                                    chsum += leakp[ch_id];
                                    MarkovDyn *evnode = dynamic_cast<MarkovDyn*>(nodes[ch_id]);
                                    parlist.clear();							
                                    parlist = evnode->CondDomain().VarList();
                                    parinst.clear();							
                                    parinst = condpars_evid.find(ch_id)->second;

                                    depmult=1.0;
                                    for (unsigned int ch=0; ch<parinst.size(); ch++){ //for each consistent parent instantiation
                                        if (evnode->CondDomain().Index(parinst.at(ch)).Value(id) == xin){
                                            depmult = evnode->operator()(evnode->CondDomain().Index(parinst.at(ch)))->Intensity()
                                                [evid.find(ch_id)->second][evid.find(ch_id)->second];

                                            for (unsigned int k=0; k<parlist.size(); k++){
                                                if (parlist.at(k) != id ){
                                                    depmult *= v.GetDistInst(parlist.at(k), 
                                                            evnode->CondDomain().Index(parinst.at(ch)).Value(parlist.at(k)) );
                                                }
                                            }
                                            depsum += depmult;
                                        }		
                                    }
                                }
                                depsum = depsum /minrate;
                                par_sum = par_sum - chsum + depsum;
                                if (maxval < (leaksum - chsum + depsum))  maxval = leaksum - chsum + depsum;
                            }			
                        } 

                        node_sum += par_sum*v.GetDistInst(id, xin);
                    }
                }
                newv.SetDistInst(id, xout, node_sum);	
            }	
        }
        v = newv;
    }

    void FactoredMatrix::TransMult(RVSimple *&x, bool transpose) const {
        FactoredVector *v = dynamic_cast<FactoredVector *>(x);

        if (trans_evid.size()>0){ 
            double res = log(v->Normalize());
            TransMultByVec(*v, transpose);
            res += log(v->Normalize());
            v->SetLogf(res);
        }
        else if (notrans_evid.size() >0){
            v->Condition(notrans_evid);
        }
    }


    void FactoredMatrix::Mult(RVSimple *&x, bool transpose) const{

        FactoredVector *v = dynamic_cast<FactoredVector *>(x);
        assert (nullptr03 != v);

        bool conditioned;
        if (evid.size()>0){
            v->Condition(evid);
            conditioned = true;
            // If some of the variables are not observed
            if (condpars.size()>0){

                double res = Unifvexpmt(*v, conditioned, transpose);
                v->SetLogf(res);
            }

        }
        else{

            //No evidence, calling normal multiplication
            conditioned = false;
            double res =  Unifvexpmt(*v, conditioned, transpose);
            v->SetLogf(res);
        }
    }




    void FactoredMatrix::Mult(RVSimple *&x)const{
        bool transpose = false;
        if (t != -1.0){	
            if (trans){ TransMult(x, transpose); }
            else{ Mult(x, transpose); }
        }
        else { 
            cout << "Cannot multiply: Interval t not set" << endl;
        }
    }

    // this is the reverse 
    void FactoredMatrix::RMult(RVSimple *&x) const{
        bool transpose = true;
        if (t != -1.0){	
            if (trans){ TransMult(x, transpose); }
            else{ Mult(x, transpose); }
        }
        else { 
            cout << "Cannot multiply: Interval t not set" << endl;
        }
    }

    double FactoredMatrix::Unifvexpmt(FactoredVector & v, bool conditioned, bool transpose) const {

        double a, theta, s, r, gamma, m, beta;
        double wnormf;
        double tnew;
        double tmpmaxval;

        FactoredVector f, w;

        //Compute norm of A
        a = minrate; 

        //Loop parameters
        theta = theta_value;
        int l = taylor_expansion_term_count; //MAX_L;
        m = ceil(a*theta/l); //ceil(theta/l); //1,10,100  // ceil(a*t/theta);
        tnew = t/m;
        s = a*tnew;
        r = exp(-s);

        w = v;
        wnormf = 0.0;
        wnormf += log(w.Normalize());
        for (int i=1; i<=m; i++){
            f = w;
            gamma = 1.0;
            beta = 1.0;
            for (int k=1; k<=l; k++){
                gamma *= (s/k);
                if (conditioned)
                    CondMultByVec(f, transpose, tmpmaxval);
                else
                    MultByVec(f, transpose, tmpmaxval);

                beta *= f.Normalize();			
                RVSimple *tmp = f.Clone();
                tmp->Mult(exp(log(gamma)+log(beta)));

                w.Add(tmp); 

            }

            w.Mult(r); 		

            //normalize w
            wnormf += log(w.Normalize());

        }

        v = w;	
        return wnormf; 
    }

    // cshelton: added
    void FactoredMatrix::initRKcache(double maxtime, double eps) const {
        rkcache.clear();
        for(uint id=0;id<nodes.size();id++)
            rkcache.push_back(ExpMCache(*amatrices[id],maxtime,eps,true));
        cachemaxt = maxtime;
    }

    void FactoredMatrix::ForAndBack(FactoredVector &v,
            double time, int Bid) const {
        int varid = Blist[Bid].first;
        const Instantiation &pinst = Blist[Bid].second;
        for(Context::iimap::const_iterator i = pinst.begin();
                i!=pinst.end();++i) {
            rkcache[i->first].vexpmt(v.Dist(i->first),time);
            v.RestrictById(i->first,i->second);
            rkcache[i->first].vexpmt(v.Dist(i->first),-time);
        }
        MarkovDyn *node = (MarkovDyn *)(bnodes[varid]);
        rkcache[varid].vexpmt(v.Dist(varid),time);
        v.Dist(varid) = v.Dist(varid)*((*node)(pinst.Index())->Intensity());
        rkcache[varid].vexpmt(v.Dist(varid),-time);
    }

    void FactoredMatrix::MultByB(FactoredVector &v, int Bid) const {
        int varid = Blist[Bid].first;
        const Instantiation &pinst = Blist[Bid].second;
        for(Context::iimap::const_iterator i = pinst.begin();
                i!=pinst.end();++i)
            v.RestrictById(i->first,i->second);
        MarkovDyn *node = (MarkovDyn *)(bnodes[varid]);
        v.Dist(varid) = v.Dist(varid)*((*node)(pinst.Index())->Intensity());
    }

    void FactoredMatrix::Forward(FactoredVector &v, double t) const {
        for(uint id=0;id<nodes.size();id++) 
            rkcache[id].vexpmt(v.Dist(id),t);
    }

    // cshelton:  not sure why this returns a value -- not sure what to return
    double FactoredMatrix::RKvexpmtA(FactoredVector & v, double time, bool conditioned, bool transpose, 
            vector< vector<double> > *timesteps, double higheps, double inith) const{

        if (cachemaxt<0.0 || time>cachemaxt) {
            for(uint id=0; id<nodes.size(); id++){
                if (timesteps != NULL ){
                    timesteps->push_back(vector<double>());
                    vector<double> *ts = &(timesteps->back());
                    vectr tmp; double logfactor;
                    v.GetDist(tmp, id, logfactor);
                    vexpmt(tmp, *amatrices[id], time, higheps, inith, ts);
                }
                double ret = vexpmt(v.Dist(id), *amatrices[id], time);
                v.Dist(id) *= exp(ret);
            } 
        } else { 
            for(uint id=0;id<nodes.size();id++) {
                if (timesteps != NULL){
                    timesteps->push_back(vector<double>());
                    vector<double> *ts = &(timesteps->back());

                    double currtime = 0.0;
                    double deltat = inith;
                    while(currtime<time){
                        ts->push_back(currtime);
                        currtime += deltat;
                    }
                    ts->push_back(currtime);
                }
                rkcache[id].vexpmt(v.Dist(id),time);
            }
        }
        return 0.0;
    }

    double FactoredMatrix::RKvexpmtA(LagFactoredVector & v, double time,
            bool forceupdate) const{
        v.t += time;
        if (forceupdate)
            for(uint id=0;id<nodes.size();id++)
                Advance(v,id);
        return 0.0;
    }



    double FactoredMatrix::RKvexpmtA(LagFactoredVector &v, double time,
            int varid, const Instantiation &pinst,
            const LagFactoredVector &vend) const {
        v.t += time;
        for(int id=0;id<(int)nodes.size();id++)
            if (id==varid || pinst.HasId(id)) Advance(v,id);
            else {
                v.Dist(id) = vend.Dist(id);
                v.realt[id] = v.t;
            }
        return 0.0;
    }



    void FactoredMatrix::Advance(LagFactoredVector &v) const {
        for(int id=0;id<(int)nodes.size();id++){
            Advance(v, id);
        }
    }




    void FactoredMatrix::Advance(LagFactoredVector &v, int id) const {
        if (v.realt[id]>=v.t) return;
        double tm = v.t-v.realt[id];
        if (cachemaxt<0.0 || tm>cachemaxt) {
            double ret = vexpmt(v.Dist(id), *amatrices[id], tm);
            v.Dist(id) *= ret;
        } else rkcache[id].vexpmt(v.Dist(id),tm);
        v.realt[id] = v.t;
    }



    void FactoredMatrix::BMultByVec(LagFactoredVector &v, int varid,
            const Instantiation &pinst) const{

        double logfactor = 0;
        MarkovDyn *node = (MarkovDyn *)(bnodes[varid]);
        const Context &pcontext = node->CondDomain();
        vector<int> parlist = pcontext.VarList();

        for (uint l=0; l<parlist.size(); l++){
            Advance(v,parlist[l]);
            v.RestrictById( parlist[l],pinst.Value(parlist[l]));
        }
        Advance(v,varid);
        v.Dist(varid) = v.Dist(varid)*((*node)(pinst.Index())->Intensity());
        if (ratedist==0) v.makenonneg();
    }



    matrix* FactoredMatrix::MakeAB(Dynamics *bnode, FactoredVector *fv, double &asc){

        if(ratedist == 0) //min
            return MakeABMin(bnode, asc);
        else if(ratedist == 1) //min sum zero
            return MakeABMinSZ(bnode, asc);
        else if(ratedist == 2) //avg according to init dist
            return MakeABInitAvg(bnode, fv, asc);
        else if(ratedist == 3) // A zero
            return MakeAZero(bnode, asc);
        else{
            cerr << "Rate dist. is not implemented in FactoredMatrix" << endl;
            return NULL;
        }

    }

    matrix* FactoredMatrix::MakeABMin(Dynamics *bnode, double &asc){

        int ds = ((MarkovDyn*)bnode)->Domain().Size();
        int cds = ((MarkovDyn*)bnode)->CondDomain().Size();

        //vector<int> parlist = ((MarkovDyn*)bnode)->CondDomain().VarList();
        //vector<double> avgdiags(ds, 0.0);

        matrix *amatrix = NULL;

        //calculate the matrix of minimal values among all rate matrices
        for (int i=0; i<cds; i++){
            if(i==0){
                amatrix = new matrix( ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                            CondDomain().Index(0))->Intensity() );
            }
            else{
                matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                        CondDomain().Index(i))->Intensity();

                for(int j=0; j<ds; j++){ 			
                    for(int k=0; k<ds; k++){
                        if( (*amatrix)[j][k] > m[j][k] )
                            (*amatrix)[j][k] = m[j][k];					 					
                    }			
                }
            }
        }		

        //calculate the remaining values in each rate matrix
        for(int i=0; i<cds; i++){
            matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(i))->Intensity();

            /*
               cout << "Q" << i <<endl;
               m.niceprint(cout);
               cout << "A" << i <<endl;
               (*amatrix).niceprint(cout);
               */
            m -= (*amatrix);
            //m.unfuzz(1e-15);
            /*		
                    cout << "B" << i <<endl;
                    m.niceprint(cout);
                    */
        }

        if(scaled == 1){
            double mindiag = 0.0;
            for(int j=0; j<ds; j++){
                if( mindiag > (*amatrix)[j][j] )
                    mindiag = (*amatrix)[j][j];
            }
            asc = abs(mindiag);
            for(int j=0; j<ds; j++){
                (*amatrix)[j][j] += asc;

            }
            /*cout << "scaled A" <<endl;
              (*amatrix).niceprint(cout);
              */
        }

        return amatrix;

    }

    matrix* FactoredMatrix::MakeABMinSZ(Dynamics *bnode, double &asc){

        int ds = ((MarkovDyn*)bnode)->Domain().Size();
        int cds = ((MarkovDyn*)bnode)->CondDomain().Size();

        //vector<int> parlist = ((MarkovDyn*)bnode)->CondDomain().VarList();
        //vector<double> avgdiags(ds, 0.0);

        matrix *amatrix = NULL;

        //calculate the matrix of minimal values among all rate matrices
        for (int i=0; i<cds; i++){
            if(i==0){
                amatrix = new matrix( ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                            CondDomain().Index(0))->Intensity() );
            }
            else{
                matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                        CondDomain().Index(i))->Intensity();

                for(int j=0; j<ds; j++){ 			
                    for(int k=0; k<ds; k++){
                        if( (*amatrix)[j][k] > m[j][k] )
                            (*amatrix)[j][k] = m[j][k];					 					
                    }			
                }
            }
        }		

        double rsum = 0.0;
        for(int j=0; j<ds; j++){ 			
            rsum = 0.0;
            for(int k=0; k<ds; k++){
                if(j!=k)
                    rsum += (*amatrix)[j][k];					 					
            }			
            (*amatrix)[j][j] = -1*rsum;
        }

        //calculate the remaining values in each rate matrix
        for(int i=0; i<cds; i++){
            matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(i))->Intensity();

            /*		
                    cout << "Q" << i << endl;
                    m.niceprint(cout);
                    cout << "A" << i <<endl;
                    (*amatrix).niceprint(cout);
                    */		
            m -= (*amatrix);
            //m.unfuzz(1e-15);
            /*			
                        cout << "B" << i <<endl;
                        m.niceprint(cout);
                        */
        }

        if(scaled == 1){
            double mindiag = 0.0;
            for(int j=0; j<ds; j++){
                if( mindiag > (*amatrix)[j][j] )
                    mindiag = (*amatrix)[j][j];
            }
            asc = abs(mindiag);
            for(int j=0; j<ds; j++){
                (*amatrix)[j][j] += asc;

            }
            /*cout << "scaled A" <<endl;
              (*amatrix).niceprint(cout);
              */
        }

        return amatrix;
    }




    matrix* FactoredMatrix::MakeABInitAvg(Dynamics *bnode, FactoredVector *fv, double &asc){

        if(fv == NULL){
            cerr << "FV is not defined for MakeABInitAvg in FactoredMatrix"<<endl;
            exit(-1);
        }
        int ds = ((MarkovDyn*)bnode)->Domain().Size();
        int cds = ((MarkovDyn*)bnode)->CondDomain().Size();

        vector<int> parlist = ((MarkovDyn*)bnode)->CondDomain().VarList();
        //vector<double> avgdiags(ds, 0.0);

        matrix *amatrix = NULL;
        /* = new matrix( ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
           CondDomain().Index(0))->Intensity() );
           */
        double parprod;
        int parid, parval;
        //calculate the matrix of avg values among all rate matrices
        for (int i=0; i<cds; i++){
            //cout << "CondDomain: " << i << endl;
            matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(i))->Intensity();

            parprod = 1.0;
            for(uint p=0; p<parlist.size(); p++){
                parid = parlist.at(p);
                parval = (((MarkovDyn*)bnode)-> \
                        CondDomain().Index(i)).Value(parid);

                parprod *= fv->GetDistInstNorm(parid, parval);
                //cout << "parid:" << parid << " parval:" << parval << " parprod:" << parprod << endl;

            }

            if(i == 0) {
                amatrix = new matrix( ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                            CondDomain().Index(i))->Intensity() );
                (*amatrix) *= parprod;
            }
            else{
                (*amatrix) += m*parprod;
            }
            /*for(int j=0; j<ds; j++){
            //double rsum = 0.0;
            for(int k=0; k<ds; k++){
            if(j!=k){
            if(i == 1)
            (*amatrix)[j][k] =  (*amatrix)[j][k]/cds;

            (*amatrix)[j][k] += m[j][k]/cds;		 
            //if(j!=k)	rsum += (*amatrix)[j][k];
            }
            else{
            if( (*amatrix)[j][k] > m[j][k] )
            (*amatrix)[j][k] = m[j][k];
            }
            }
            //(*amatrix)[j][j] = -1*rsum;

            }*/
        }		

        //calculate the remaining values in each rate matrix
        for(int i=0; i<cds; i++){
            matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(i))->Intensity();

            /*cout << "Q" << i <<endl;
              m.niceprint(cout);
              cout << "A" <<endl;
              (*amatrix).niceprint(cout);
              */
            m -= (*amatrix);
            //m.unfuzz(1e-15);

            /*
               cout << "B" << i <<endl;
               m.niceprint(cout);
               */
        }

        if(scaled == 1){
            double mindiag = 0.0;
            for(int j=0; j<ds; j++){
                if( mindiag > (*amatrix)[j][j] )
                    mindiag = (*amatrix)[j][j];
            }
            asc = abs(mindiag);
            for(int j=0; j<ds; j++){
                (*amatrix)[j][j] += asc;

            }
            //cout << "scaled A" <<endl;
            //(*amatrix).niceprint(cout);

        }

        return amatrix;

    }

    matrix* FactoredMatrix::MakeABAvg(Dynamics *bnode, double &asc){

        int ds = ((MarkovDyn*)bnode)->Domain().Size();
        int cds = ((MarkovDyn*)bnode)->CondDomain().Size();

        //vector<int> parlist = ((MarkovDyn*)bnode)->CondDomain().VarList();
        //vector<double> avgdiags(ds, 0.0);

        matrix *amatrix = new matrix( ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(0))->Intensity() );

        //calculate the matrix of avg values among all rate matrices
        for (int i=1; i<cds; i++){
            matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(i))->Intensity();

            for(int j=0; j<ds; j++){
                //double rsum = 0.0;
                for(int k=0; k<ds; k++){
                    if(j!=k){
                        if(i == 1)
                            (*amatrix)[j][k] =  (*amatrix)[j][k]/cds;

                        (*amatrix)[j][k] += m[j][k]/cds;		 
                        //if(j!=k)	rsum += (*amatrix)[j][k];
                    }
                    else{
                        if( (*amatrix)[j][k] > m[j][k] )
                            (*amatrix)[j][k] = m[j][k];
                    }
                }
                //(*amatrix)[j][j] = -1*rsum;

            }
        }		

        //calculate the remaining values in each rate matrix
        for(int i=0; i<cds; i++){
            matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(i))->Intensity();

            /*cout << "Q" << i <<endl;
              m.niceprint(cout);
              cout << "A" <<endl;
              (*amatrix).niceprint(cout);
              */
            m -= (*amatrix);
            //m.unfuzz(1e-15);

            /*
               cout << "B" << i <<endl;
               m.niceprint(cout);
               */
        }

        return amatrix;

    }

    matrix* FactoredMatrix::MakeAZero(Dynamics *bnode, double &asc){

        int ds = ((MarkovDyn*)bnode)->Domain().Size();
        int cds = ((MarkovDyn*)bnode)->CondDomain().Size();

        matrix *amatrix = NULL;

        amatrix = new matrix( ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(0))->Intensity() );
        for(int j=0; j<ds; j++){ 			
            for(int k=0; k<ds; k++){
                (*amatrix)[j][k] = 0.0;					 					
            }			
        }


        //calculate the remaining values in each rate matrix
        for(int i=0; i<cds; i++){
            matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> \
                    CondDomain().Index(i))->Intensity();
            m -= (*amatrix);
        }


        return amatrix;

    }

    void FactoredMatrix::SetL(int terms_in_taylor_expansion) {
        taylor_expansion_term_count = terms_in_taylor_expansion;
    }


    void FactoredMatrix::SetTheta(double theta) {
        theta_value = theta;
    }


    void FactoredMatrix::Load(std::istream &is) {
        throw not_yet_implemented_error(the_class_name, "Load");
    }

    void FactoredMatrix::Save(std::ostream &os) const {
        throw not_yet_implemented_error(the_class_name, "Save");
    }

    double FactoredMatrix::Prob(int ind, int cind, bool log) const {
        throw not_yet_implemented_error(the_class_name, "Prob");
    }

    void FactoredMatrix::Normalize() {
        throw not_yet_implemented_error(the_class_name, "Normalize");
    }

    RVSimple * FactoredMatrix::Condition(int ind) const {
        throw not_yet_implemented_error(the_class_name, "Condition");
    }

    void FactoredMatrix::Restrict(const std::vector<int> &ind, const std::vector<int> &cind) {
        throw not_yet_implemented_error(the_class_name, "Restrict");
    }


    typedef std::map<std::pair<int,int>, std::vector<std::vector<int> > > RemapT;
    void FactoredMatrix::Reindex(int n, int nc, const RemapT &ind) {
        throw not_yet_implemented_error(the_class_name, "Reindex");
    }

    // Multiplies by another measure (point-wise)
    void FactoredMatrix::MultBy(const RVCondSimple *x) {
        throw not_yet_implemented_error(the_class_name, "MultBy");
    }

    // return suitable sufficient statistics object
    SS* FactoredMatrix::BlankSS() const {
        throw not_yet_implemented_error(the_class_name, "BlankSS");
    }

    // add independent draw of (cind,ind) (conditioning value, value) to ss
    void FactoredMatrix::AddSS(int cind, int ind, SS *ss, double w) const{
        throw not_yet_implemented_error(the_class_name, "AddSS");
    }


    void FactoredMatrix::AddSS(const SS *toadd, const RVCondSimple* rv,
            const std::vector<std::vector<int> > &mapping,
            int mycondi, int rvccondi,
            SS *ss, double w) const {
        throw not_yet_implemented_error(the_class_name, "AddSS");
    }


    // add average of independent draws to ss
    void FactoredMatrix::AddExpSS(int cind, SS *ss, double w) const {
        throw not_yet_implemented_error(the_class_name, "AddExpSS");
    }


    void FactoredMatrix::AddExpSS(int cind, const RVSimple *rvs, SS *ss, double w) const {
        throw not_yet_implemented_error(the_class_name, "AddExpSS");
    }


    // returns a sample
    int FactoredMatrix::Sample(int cx, Random &rand ) const {
        throw not_yet_implemented_error(the_class_name, "Sample");
    }


    // sets parameters to ML estimate
    void FactoredMatrix::Maximize(const SS *ss) {
        throw not_yet_implemented_error(the_class_name, "Maximize");
    }


    // Returns a new object that can be multiplied into this
    // object (as per Mult and RMult)
    RVSimple * FactoredMatrix::Convert(const RVSimple *rvs) const {
        //  FactoredVector * v =  new FactoredVector(rvs);
        return rvs->Clone();
    }


    // see above class for description
    void FactoredMatrix::Scramble(double alpha, double degree, Random &rand) {
        throw not_yet_implemented_error(the_class_name, "Scramble");
    }


    double FactoredMatrix::LLH(const SS *ss) const{
        throw not_yet_implemented_error(the_class_name, "LLH");
    }


    double FactoredMatrix::GetScore(double numTrans, const SS *ss) const {
        throw not_yet_implemented_error(the_class_name, "GetScore");
    }
    ////////


} // end of ctbn namespace
