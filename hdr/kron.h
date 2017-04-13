#ifndef KRON_H
#define KRON_H

#include <vector>
#include <utility>
#include "context.h"
#include "rk.h"
#include "markovdyn.h"
#include <map>

namespace ctbn {
template<typename T>
struct kprod {
	std::vector<T> x;
	double scalar=1.0;

	template<typename S>
	kprod<decltype(std::declval<T>()*std::declval<S>())>
	operator*(const kprod<S> &y) const {
		kprod<decltype(std::declval<T>()*std::declval<S>())> ret;
		// assert(x.size()==y.x.size());
		for(int i=0;i<x.size();i++)
			ret.x.emplace_back(x[i]*y.x[i]);
		ret.scalar = scalar*y.scalar;
		return ret;
	}

	template<typename S>
	void postmultby(const kprod<S> &y, double f) {
		for(int i=0;i<x.size();i++)
			x[i] = x[i]*y.x[i];
		scalar *= f;
	}

	/*
	template<typename S>
	kprod<decltype(std::declval<T>()*std::declval<S>())>
	operator*(const S &y) const {
		kprod<decltype(std::declval<T>()*std::declval<S>())> ret;
		for(int i=0;i<x.size();i++)
			ret.x.emplace_back(x[i]*y);
		return ret;
	}
	*/
	kprod<T> operator*(double y) const {
		kprod<T> ret(*this);
		ret.scalar *= y;
		return ret;
	}

	double absmax() const {
		double ret = x[0].absmax();
		for(int i=1;i<x.size();i++)
			ret *= x[i].absmax();
		return ret*scalar;
	}

	kprod<T> hadamard(const kprod<T> &p) const {
		kprod<T> ret;
		ret.scalar = scalar*p.scalar;
		for(int i=0;i<x.size();i++)
			ret.x.emplace_back(x[i].dotstar(p.x[i]));
		return ret;
	}

};

matrix flatten(const kprod<matrix> &X) { //, const Context &c) {
	Context c1,c2;
	for(int i=0;i<X.x.size();i++) {
		c1.AddVar(i,X.x[i].getm());
		c2.AddVar(i,X.x[i].getn());
	}
	matrix ret(c1.Size(),c2.Size(),1.0);
	forallassign(Instantiation &i,c1)
		forallassign(Instantiation &j,c2) {
			int vi = 0;
			for(int v = 0; v < X.x.size();v++) {
				ret[i.Index()][j.Index()] *= X.x[vi][i.Value(v)][j.Value(v)];
				vi++;
			}
		}
	ret *= X.scalar;
	return ret;
}

vectr flatten(const kprod<vectr> &X) { //, const Context &c) {
	Context c;
	for(int i=0;i<X.x.size();i++)
		c.AddVar(i,X.x[i].getm());
	int retn = c.Size();
	vectr ret(retn,1.0);
	forallassign(Instantiation &i,c) {
		int vi = 0;
		for(int v = 0; v < X.x.size();v++) {
			ret[i.Index()] *= X.x[vi][i.Value(v)];
			vi++;
		}
	}
	ret *= X.scalar;
	return ret;
}

/*
template<typename S, typename T>
kprod<decltype(std::declval<S>()*std::declval<T>())>
operator*(const S &x, const kprod<T> &y) {
	kprod<decltype(std::declval<S>()*std::declval<T>())> ret(y);
	ret.scalar *=x;
	return ret;
}
*/

template<typename T>
kprod<T> operator*(double &x, const kprod<T> &y) {
	kprod<T> ret(y);
	ret.scalar *=x;
	return ret;
}

template<typename T>
struct ksum {
	std::vector<T> x;

	template<typename S>
	ksum<decltype(std::declval<T>()*std::declval<S>())>
	operator+(const ksum<S> &y) const {
		ksum<decltype(std::declval<T>()*std::declval<S>())> ret;
		// assert(x.size()==y.x.size());
		for(int i=0;i<x.size();i++)
			ret.x.emplace_back(x[i]+y.x[i]);
		return ret;
	}

	T flatten() const {
		Context c;
		for(int i=0;i<x.size();i++)
			c.AddVar(i,x[i].getm());
		T ret(c.Size(),c.Size(),0.0);
		forallassign(Instantiation &i,c) {
			foralldiffone_diff(Instantiation &j, int dvar, i)
				ret[i.Index()][j.Index()] += x[dvar][i.Value(dvar)][j.Value(dvar)];
			for(int dvar : c.Vars())
				ret[i.Index()][i.Index()] += x[dvar][i.Value(dvar)][i.Value(dvar)];
		}
		return ret;
	}

	T flatten(const Context &c) const {
		T ret(c.Size(),c.Size(),0.0);
		forallassign(Instantiation &i,c) {
			foralldiffone_diff(Instantiation &j, int dvar, i)
				ret[i.Index()][j.Index()] += x[dvar][i.Value(dvar)][j.Value(dvar)];
			for(int dvar : c.Vars())
				ret[i.Index()][i.Index()] += x[dvar][i.Value(dvar)][i.Value(dvar)];
		}
		return ret;
	}
};



template<typename V>
kprod<V> vexpmt(const kprod<V> &v, const ksum<matrix> &A, double t,
		double eps=-1.0, double h=-1.0) {
	kprod<V> ret;
	for(std::size_t i=0;i<v.x.size();i++) {
		V a = v.x[i];
		vexpmt(a,A.x[i],t,eps,h);
		ret.x.emplace_back(std::move(a));
	}
	ret.scalar = v.scalar;
	return ret;
}

template<typename V>
kprod<V> expmtv(const kprod<V> &v, const ksum<matrix> &A, double t,
		double eps=-1.0, double h=-1.0) {
	kprod<V> ret;
	for(std::size_t i=0;i<v.x.size();i++) {
		V b = v.x[i];
		expmtv(b,A.x[i],t,eps,h);
		ret.x.emplace_back(std::move(b));
	}
	ret.scalar = v.scalar;
	return ret;
}

template<typename T>
struct sparsekprod {
	std::vector<T> x;
	std::vector<std::size_t> i;
	double scalar=1.0;

	sparsekprod<T> operator*(const sparsekprod<T> &y) const {
		sparsekprod<T> ret;
		for(int j=0,k=0;j<i.size() && k<y.i.size();) {
			if (j<i.size() && k<y.i.size() && i[j]==y.i[k]) {
				ret.x.emplace_back(x[j]*y.x[k]);
				ret.i.emplace_back(i[j]);
				j++; k++;
			} else if (k>=y.i.size() || i[j]<y.i[k]) {
				ret.x.emplace_back(x[j]);
				ret.i.emplace_back(i[j]);
				j++;
			} else {
				ret.x.emplace_back(y.x[k]);
				ret.i.emplace_back(y.i[k]);
				k++;
			}
		}
		ret.scalar = scalar*y.scalar;
		return ret;
	}

	template<typename S>
	kprod<S> operator*(const kprod<S> &y) const {
		kprod<S> ret;
		for(int j=0,k=0;j<i.size() || k<y.x.size();k++) {
			if (j<i.size() && i[j]==k) {
				ret.x.emplace_back(x[j]*y.x[k]);
				j++;
			} else ret.x.emplace_back(y.x[k]);
		}
		ret.scalar = scalar*y.scalar;
		return ret;
	}

	// assumes square matrix (at least currently)
	T flatten(const Context &full) const {
		T ret(full.Size(),full.Size(),0.0);
		Context vars;
		for(int l=0;l<i.size();l++)
			vars.AddVar(i[l],x[l].getm());
		Context notvars(full,vars,Context::DIFFERENCE);
		forallassign(Instantiation &j, vars)
			forallassign(Instantiation &k, vars)
				forallassign_except(Instantiation &jj, full, j) {
					Instantiation kk(jj);
					for(auto v : vars.Vars())
						kk.SetVal(v,k.Value(v));
					double val = 1.0;
					for(int l=0;l<i.size();l++)
						val *= x[l][j.Value(l)][k.Value(l)];
					ret[jj.Index()][kk.Index()] = val;
				}
		return ret;
	}

};

template<typename T, typename S>
kprod<S> operator*(const kprod<S> &x, const sparsekprod<T> &y) {
	kprod<S> ret;
	for(int j=0,k=0;j<y.i.size() || k<x.x.size();k++) {
		if (j<y.i.size() && y.i[j]==k) {
			ret.x.emplace_back(x.x[k]*y.x[j]);
			j++;
		} else ret.x.emplace_back(x.x[k]);
	}
	ret.scalar = x.scalar*y.scalar;
	return ret;
}

// A matrix over a subset of context (in vars)
// with an implied Kronecker product with I for the other variables
template<typename MT>
struct partialmatrix {
	Context vars;
	int primary;
	MT M;

	partialmatrix(const partialmatrix &pm, bool copycache) : vars(pm.vars), primary(pm.primary), M(pm.M) {
		if (copycache) sumkprod = pm.sumkprod;
	}

	partialmatrix(const Context &v, const ksum<MT> &A, int prime) : vars(v), primary(prime), M(v.Size(),v.Size(),0.0) {
		
		forallassign(Instantiation &i,vars)
			foralldiffone_incself_diff(Instantiation &j, int dvar, i)
				M[i.Index()][j.Index()] += A.x[dvar][i.Value(dvar)][j.Value(dvar)];
	}

	template<typename I>
	static Context makecontext(const std::vector<I> &v, const ksum<MT> &A, const std::map<int,int> &varmap) {
		Context ret;
		for(auto &i : v) ret.AddVar(varmap.at(i),A.x[varmap.at(i)].getm());
		return ret;
	}

	template<typename I>
	static Context makecontext(const std::vector<I> &v, const ksum<MT> &A) {
		Context ret;
		for(auto &i : v) ret.AddVar(i,A.x[i].getm());
		return ret;
	}
		
	template<typename I>
	partialmatrix(const std::vector<I> &v, const ksum<MT> &A,
			int prime, const std::map<int,int> &varmap)
			: vars(makecontext(v,A,varmap)), primary(varmap.at(prime)),
			M(A.flatten(vars)) {
	}

	template<typename I>
	partialmatrix(const std::vector<I> &v, const ksum<MT> &A, int prime)
			: vars(makecontext(v,A)), primary(prime),
			M(vars.Size(),vars.Size(),0.0) {
		forallassign(Instantiation &i,vars)
			foralldiffone_incself_diff(Instantiation &j, int dvar, i)
				M[i.Index()][j.Index()] += A.x[dvar][i.Value(dvar)][j.Value(dvar)];
	}

	partialmatrix(const MarkovDyn &node, const std::map<int,int> &varmap, const MT &A) {
		int primaryorig = *(node.Domain().Vars().begin());
		primary = varmap.at(primaryorig);
		vars.AddVar(primary,node.Domain().Cardinality(primaryorig));
		Context notprime;
		for(auto &vorig : node.CondDomain().Vars()) {
			int v = varmap.at(vorig);
			vars.AddVar(v,node.CondDomain().Cardinality(vorig));
			notprime.AddVar(v,node.CondDomain().Cardinality(vorig));
		}

		M = MT(vars.Size(),vars.Size(),0.0);
		Instantiation myi(notprime);
		forallassign(Instantiation &i, node.CondDomain()) {
			for(auto vpair : i) // it would be great to only set those that changed
				myi.SetVal(varmap.at(vpair.first),vpair.second);
			const MT &Qcond = node.operator()(i)->Intensity();
			forallassign_except(Instantiation &ii, vars, myi)
				forallassign_except(Instantiation &jj, vars, myi)
					M[ii.Index()][jj.Index()]
						= Qcond[ii.Value(primary)][jj.Value(primary)] - A[ii.Value(primary)][jj.Value(primary)];
		}
	}

	partialmatrix(const MarkovDyn &node, const MT &A) {
		primary = *(node.Domain().Vars().begin());
		Context notprime = node.CondDomain();
		vars = node.CondDomain()+notprime;

		M = MT(vars.Size(),vars.Size(),0.0);
		forallassign(Instantiation &i, notprime) {
			const MT &Qcond = node.operator()(i)->Intensity();
			forallassign_except(Instantiation &ii, vars, i)
				forallassign_except(Instantiation &jj, vars, i)
					M[ii.Index()][jj.Index()]
						= Qcond[ii.Value(primary)][jj.Value(primary)] - A[ii.Value(primary)][jj.Value(primary)];
		}
	}

	partialmatrix operator*(const partialmatrix &m) const {
		partialmatrix ret(*this,false);
		ret.sumkprod.clear();
		ret.M *= m.M;
		return ret;
	}
	partialmatrix operator*(double f) const {
		partialmatrix ret(*this,false);
		ret.sumkprod.clear();
		ret.M *= f;
		return ret;
	}
	partialmatrix operator/(double f) const {
		partialmatrix ret(*this,false);
		ret.sumkprod.clear();
		ret.M /= f;
		return ret;
	}
	partialmatrix operator+(const partialmatrix &m) const {
		partialmatrix ret(*this,false);
		ret.sumkprod.clear();
		ret.M += m.M;
		return ret;
	}
	partialmatrix &operator+=(const partialmatrix &m) {
		M += m.M;
		sumkprod.clear();
		return *this;
	}
	partialmatrix operator-(const partialmatrix &m) const {
		partialmatrix ret(*this,false);
		ret.sumkprod.clear();
		ret.M -= m.M;
		return ret;
	}

	mutable std::vector<sparsekprod<MT>> sumkprod;

	const std::vector<sparsekprod<MT>> &tokprod() const {
		if (sumkprod.empty()) {
		int np = vars.Cardinality(primary);
		Context primec;
		primec.AddVar(primary,np);
		Context notprime(vars,primec,Context::DIFFERENCE);
		forallassign(Instantiation &i, notprime)
			forallassign(Instantiation &j, notprime) {
				MT m(np,np,0.0);
				bool nonzero=false;
				forallassign_except(Instantiation &ii, vars, i)
					forallassign_except(Instantiation &jj, vars, j) {
						double v = M[ii.Index()][jj.Index()];
						if (v!=0.0) {
							nonzero = true;
							m[ii.Value(primary)][jj.Value(primary)] = v;
						}
					}
				if (nonzero) {
					sumkprod.emplace_back();
					sparsekprod<MT> &km = sumkprod.back();
					for(auto &iv : vars.Vars()) {
						km.i.emplace_back(iv);
						if (iv==primary) km.x.emplace_back(m);
						else {
							km.x.emplace_back(notprime.Cardinality(iv),notprime.Cardinality(iv),0.0);
							km.x.back()[i.Value(iv)][j.Value(iv)] = 1.0;
						}
					}
				}
			}
		}
		return sumkprod;
	}

	MT flatten(const Context &full) const {
		MT ret(full.Size(),full.Size(),0.0);
		Context notvars(full,vars,Context::DIFFERENCE);
		forallassign(Instantiation &i, vars)
			forallassign(Instantiation &j, vars)
				forallassign_except(Instantiation &ii, full, i) {
					Instantiation jj(ii);
					for(auto v : vars.Vars())
						jj.SetVal(v,j.Value(v));
					ret[ii.Index()][jj.Index()] = M[i.Index()][j.Index()];
				}
		return ret;
	}

};

template<typename T>
std::ostream &operator<<(std::ostream &os, const kprod<T> &kp) {
	bool first=true;
	for(auto &m : kp.x) {
		if (!first) os << "(*)";
		os << m;
		first=false;
	}
	os << " (scalar=" << kp.scalar << ")";
	return os;
}
template<typename T>
std::ostream &operator<<(std::ostream &os, const ksum<T> &kp) {
	bool first=true;
	for(auto &m : kp.x) {
		if (!first) os << "(+)";
		os << m;
		first=false;
	}
	return os;
}
template<typename T>
std::ostream &operator<<(std::ostream &os, const sparsekprod<T> &kp) {
	for(int i=0;i<kp.x.size();i++) {
		if (!i) os << "(*)";
		os << '[' << kp.i[i] << ']' << ' ' << kp.x[i];
	}
	os << " (scalar=" << kp.scalar << ")";
	return os;
}

template<typename MT>
std::ostream &operator<<(std::ostream &os, const partialmatrix<MT> &pm) {
	bool first=true;
	for(int v : pm.vars.Vars()) {
		if (!first) os << ',';
		if (v==pm.primary) os << '{' << v << '}';
		else os << v;
		first = false;
	}
	os << ' ';
	os << pm.M;
	return os;
}

}

#endif
