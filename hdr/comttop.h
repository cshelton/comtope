#ifndef COMTTOP_H
#define COMTTOP_H

#include "kron.h"
#include "ctbndyn.h"
#include <algorithm>
#include <map>
#include <iostream> // to be removed
#include <queue>

namespace ctbn {

// commutatorseq[i] = {B,A^i} = [{B,A^(i-1)}, A]/i! with {B,A^0} = B
template<typename T>
struct commutatorseq {
     T A;
     std::vector<T> C;

	commutatorseq(T initA, T initB) : A(std::move(initA)) {
		C.emplace_back(std::move(initB));
	}

     const T &operator[](std::size_t i) {
          for(std::size_t j=C.size();j<=i;j++)
               C.emplace_back((C[j-1]*A - A*C[j-1])/j);
          return C[i];
     }
};

template<typename T>
struct deriv {
	T d;
	T B;

	deriv(T initA, T initB) : d(initB), B(initB) {
		commutatorseq<T> comseq(std::move(initA),std::move(initB));
		int k;
		for(k=1;k<100;k++) {
			const auto &M = comseq[k];
			d += M;
			if (M.M.absmax()<1e-6) break;
		}
		std::cout << "deriv maxk = " << k << std::endl;
	}
};

template<typename MT>
struct ctbncomdecomp {
	ksum<MT> A;
	std::vector<commutatorseq<partialmatrix<MT>>> comseqs;
	//std::vector<deriv<partialmatrix<MT>>> derivs;

	ctbncomdecomp(const CTBNDyn &dyn, double scale=1.0, bool transpose=false) {
		build(dyn, [dyn](int i) {
				return getA(
						*(dynamic_cast<const MarkovDyn *>(dyn.Node(i))),
						static_cast<const double &(*)(const double &,const double &)>(std::max),
						static_cast<const double &(*)(const double &,const double &)>(std::min)
					);
			},
			scale, transpose);
	}

	ctbncomdecomp(const CTBNDyn &dyn, const kprod<vectr> &p0, double scale=1.0, bool transpose=false) {
		build(dyn, [dyn,p0](int i) {
				return getA(
						*(dynamic_cast<const MarkovDyn *>(dyn.Node(i))),
						p0.x[i]
					);
			},
			scale, transpose);
	}

	template<typename F>
	void build(const CTBNDyn &dyn, F f, double scale=1.0, bool transpose=false) {
		std::map<int,int> var2node; // rebuilding from dyn... seems silly
		int nn = dyn.NumofNodes();
		for(int i=0;i<nn;i++) {
			// assuming 1 var/node (otherwise, this is doesn't work)
			var2node[dyn.Node(i)->Domain().VarList()[0]] = i;
			if (!transpose) A.x.emplace_back(f(i));
			else A.x.emplace_back(f(i).transpose());
		}
		for(int i=0;i<nn;i++) {
			partialmatrix<MT> Ai((dyn.Node(i)->Domain()+dyn.Node(i)->CondDomain()).VarList(),
					A,i,var2node);
			partialmatrix<MT> Bi(dynamic_cast<const MarkovDyn &>(*dyn.Node(i)),var2node,A.x[i]);
			Bi.M = Bi.M.transpose();
			comseqs.emplace_back(Ai*scale,Bi*scale);
			//derivs.emplace_back(Ai*scale,Bi*scale);
		}
		for(int i=0;i<nn;i++)
			A.x[i] *= scale;
	}

	template<typename DiagOp, typename OffDiagOp>
	static MT getA(const MarkovDyn &node, DiagOp dop, OffDiagOp odop) {
		int ds = node.Domain().Size();
		int cds = node.CondDomain().Size();
		MT ret = node.operator()(0)->Intensity();
		for(int i=1;i<cds;i++) {
			const MT &nm= node.operator()(i)->Intensity();
			for(int j=0;j<ret.getm();j++)
				for (int k=0;k<ret.getn();k++)
					ret[j][k] = (j==k) ? dop(ret[j][k],nm[j][k])
							: odop(ret[j][k],nm[j][k]);
		}
		return ret;
	}

	static MT getA(const MarkovDyn &node, const vectr &p0) {
		int ds = node.Domain().Size();
		int cds = node.CondDomain().Size();
		MT ret = node.operator()(0)->Intensity()*p0[0];
		for(int i=1;i<cds;i++) {
			const MT &nm= node.operator()(i)->Intensity();
			ret = ret + nm*p0[i];
		}
		return ret;
	}

};

vectr flattenmarg(const kprod<vectr> &X) {
	Context c;
	for(int i=0;i<X.x.size();i++)
		c.AddVar(i,X.x[i].getm());
	int retn = c.Size();
	vectr ret(retn,1.0);
	double scale = 1.0, s1=0.0;
	for(int v=0;v<X.x.size();v++) {
		double s = X.x[v].sum();
		scale *= s;
		s1 += s;
	}
	forallassign(Instantiation &i,c) {
		int vi = 0;
		for(int v = 0; v < X.x.size();v++) {
			ret[i.Index()] *= X.x[vi][i.Value(v)];
			vi++;
		}
	}
	ret *= X.scalar/scale*(s1/X.x.size());
	return ret;
}


template<typename MT, typename VT>
struct commexp {
	struct node {
		double pri;
		double alpha;
		int d,i,k,l;
		kprod<MT> M;
		// for debug
		std::vector<int> ds;
		std::vector<int> is;
		std::vector<int> ks;

		node(double p, double a, int dd, int ii, int kk, int ll,kprod<MT> MM,
				const node &par)
			: pri(p),
			//: pri(a/(dd+kk+1)),
			alpha(a), d(dd), i(ii), k(kk), l(ll), M(std::move(MM)),
	    		ds(par.ds), is(par.is), ks(par.ks) {
			ds.push_back(d);
			is.push_back(i);
			ks.push_back(k);

		}

		node(double p, double a, int dd, int ii, int kk, int ll,kprod<MT> MM)
			: pri(p),
			//: pri(a/(dd+kk+1)),
			alpha(a), d(dd), i(ii), k(kk), l(ll), M(std::move(MM)) { }

		bool operator<(const node &n) const { return pri<n.pri; }
		bool operator<=(const node &n) const { return pri<=n.pri; }
		bool operator>(const node &n) const { return pri>n.pri; }
		bool operator>=(const node &n) const { return pri>=n.pri; }
		bool operator==(const node &n) const { return pri==n.pri; }
		bool operator!=(const node &n) const { return pri!=n.pri; }
	};

	ctbncomdecomp<MT> decomp;
	kprod<VT> veAt;
	std::priority_queue<node> q;
	Context full; // for debug

	template<typename... Ds>
	commexp(const kprod<VT> &v0, const CTBNDyn &dyn, Ds &&...ds)//double delt, bool transpose=false)
				: decomp(dyn,std::forward<Ds>(ds)...),
				veAt(vexpmt(v0,decomp.A,1.0)),
				full(dyn.Domain()) { 
		//currpri = 0;
		kprod<MT> eye;
		for(auto &Ai : decomp.A.x) {
			MT eyei(Ai.getm(),Ai.getn(),0.0);
			for(int i=0;i<eyei.getm() && i<eyei.getn();i++)
				eyei[i][i] = 1.0;
			eye.x.emplace_back(std::move(eyei));
		}
		double factor = 1.0;
		/*
		for(auto &s : decomp.comseqs) {
			auto &B = s.C[0].M;
			double mind = 0.0;
			for(int i=0;i<B.getm();i++)
				if (mind>B[i][i]) mind=B[i][i];
			for(int i=0;i<B.getm();i++)
				B[i][i] -= mind;
			factor *= std::exp(mind);
		}
		*/
		q.emplace(0.0,factor,0,-1,0,0,eye);
	}
	int currpri = 0.0;

	template<typename... Ts>
	void addtoq(double pri, Ts &&...ts) {
		q.emplace(pri,std::forward<Ts>(ts)...);
		//q.emplace(currpri,std::forward<Ts>(ts)...);
		currpri--;
	}

	std::vector<kprod<VT>> list;

	const kprod<VT> &operator[](int i) {
		while(list.size()<=i) {
			auto n = getmore();
			for(auto &x : n)
				list.emplace_back(std::move(x));
		}
		return list[i];
	}

	int curri = 0;

	const kprod<VT> &next() {
		return (*this)[curri++];
	}

	std::vector<kprod<VT>> getmore() {
		std::vector<kprod<VT>> ret;
		if (q.empty()) return ret;
		node n = std::move(q.top()); // no whammies between here and next line!
		q.pop();
		if (n.i<0) {
			double pri = veAt.absmax();
			for(int i=0;i<decomp.A.x.size();i++)
				addtoq(pri,1.0,0,i,0,1,n.M,n);
			ret.emplace_back(veAt*n.alpha);
		} else {
			double scalar = n.alpha/(n.d+n.k+1);
			//std::cout << "d = " << n.ds << " & i = " << n.is << " & k = " << n.ks << std::endl;
			//std::cout << "n.pri = " << n.pri << std::endl;

			auto &C = decomp.comseqs[n.i][n.k];
			//auto &C = decomp.derivs[n.i].d;
			//std::cout << "M = " << flatten(n.M) << std::endl;
			//std::cout << "M = " << C.M << std::endl;
			//std::cout << "veAt = " << flatten(veAt) << std::endl;
			//std::cout << "C = " << C.flatten(full) << std::endl;
			//std::cout << "veAt*C = " << flatten(veAt)*C.flatten(full) << std::endl;
			std::vector<sparsekprod<MT>> Ms = C.tokprod(); //decomp.comseqs[n.i].C[n.k].tokprod();
			double maxpri = 0.0;
			for(auto &Mi : Ms) {
				//std::cout << "Mi = " << Mi.flatten(full) << std::endl;
				kprod<MT> newM = Mi*n.M;
				kprod<VT> r(veAt);
				r.postmultby(newM,scalar);
				ret.emplace_back(r);
				double pri = newM.absmax()*scalar;
				if (pri==0.0) continue;
				//maxpri = std::max(maxpri,pri);
				maxpri = maxpri+pri;
				//std::cout << "pri=" << pri << std::endl;
				//std::cout << "for " << flatten(r) << std::endl;
				for(int i=0;i<decomp.A.x.size();i++)
					addtoq(pri,n.alpha/(n.d+n.k+1),n.d+n.k+1,i,0,n.l+1,newM,n);
			}

			maxpri = C.M.absmax()*scalar; 
			if (maxpri>0)
				addtoq(maxpri,n.alpha,n.d,n.i,n.k+1,n.l,std::move(n.M),n);
		}
		return ret;
	}

};

template<typename VT>
void addtomarg(kprod<VT> &x, const kprod<VT> &y) {
	double s = y.scalar;
	int zindex = -1;
	std::vector<double> ss(y.x.size());
	for(int i=0;i<y.x.size();i++) {
		double t = y.x[i].sum();
		ss[i] = t;
		if (t==0) {
			if (zindex>=0) return;
			zindex=i;
		} else s *= t;
	}
	if (zindex>=0)
		x.x[zindex] += y.x[zindex]*s;
	else
		for(int i=0;i<x.x.size();i++)
			x.x[i] += y.x[i]*s/ss[i];
}

template<typename T>
auto run(T &lister, int maxitt) {
	auto toadd = lister.next(); //getmore()[0];
	decltype(toadd) ret;
	for(int i=0;i<toadd.x.size();i++)
		ret.x.emplace_back(toadd.x[i].getm(),0.0);
	addtomarg(ret,toadd);
	for(int i=1;i<maxitt;i++) {
		if (i%1000==0)
			std::cout << i << "th iteration: " << std::endl;

		/*
		   auto res = getmore();
		   kprod<VT> toadd;
		   for(int i=0;i<ret.x.size();i++)
		   toadd.x.emplace_back(ret.x[i].getm(),0.0);
		   for(auto &v : res) addtomarg(toadd,v);
		   std::cout << "added: " << std::endl;
		   std::cout << toadd << std::endl;
		   std::cout << flattenmarg(toadd) << std::endl;
		   if (!res.empty()) {
		   vectr toadd2 = flatten(res[0]);
		   for(int i=1;i<res.size();i++)
		   toadd2 += flatten(res[i]);
		   std::cout << toadd2 << std::endl;
		   }
		   for(auto &v : res) addtomarg(ret,v);
		   */
		addtomarg(ret,lister.next());

		if (i%1000==0) {
			std::cout << "running sum: " << std::endl;
			std::cout << ret << std::endl;
		}
		//std::cout << flattenmarg(ret) << std::endl;
	}
	return ret;
}

template<typename MT, typename VT>
struct smoother {

	commexp<MT,VT> alpha,beta;

	smoother(const kprod<VT> &v0, const kprod<VT> &vT, const CTBNDyn &dyn, double T, double t)
		//: alpha(v0,dyn,t), beta(vT,dyn,T-t,true) {
		: alpha(v0,dyn,v0,t), beta(vT,dyn,vT,T-t,true) {
		ai=0;
		bi=0;
	}

	auto operator[](int i) {
		int d = floor((sqrt(8.0*i+1)-1)/2);
		int r = i-(d*(d+1))/2;
		return get(d-r,r);
	}

	int ai=0,bi=0;
	auto next() {
		auto ret = get(ai,bi);
		ai--; bi++;
		if (ai<0) {
			ai = bi;
			bi = 0;
		}
		return ret;
	}

	auto get(int a,int b) {
		return alpha[a].hadamard(beta[b]);
	}
};

}

#endif
