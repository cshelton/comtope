#include "factoredvector.h"



namespace ctbn {

using namespace std;


static const string the_class_name("FactoredMatrix");


FactoredVector::FactoredVector() : RVSimple()
{
	logf = 0.0;
}

/*
FactoredVector::FactoredVector(const FactoredVector &v) : RVSimple(), dists(v.dists), logf(v.logf) {
	//unsigned int n = v.dists.size();
	//dists.resize(n);
	//for (unsigned int i = 0; i != n; ++i) {
//		dists[i] = v.dists[i];
//	}
//	logf = v.logf;
}
*/

FactoredVector::FactoredVector(RV *rv) : RVSimple() {
	BN *bn = dynamic_cast<BN*>(rv); 
	int n = bn->NumofNodes();
	dists.resize(n);
	int domainsize; 
	logf = 0.0;
	vector<double> sums(n,0.0);
	for (int i=0; i<n; i++){
		MultiRV *newrv = (MultiRV*) bn->NodeByVar(i);
		domainsize = newrv->Domain().Size();
		vectr d(domainsize, -1.0);
		double logfactor;
		newrv->operator[](0).GetDist(d, logfactor);
		dists[i] = d*exp(logfactor);
		sums[i] = dists[i].sum();
		//logf += logfactor;

	}
	
	for(unsigned int i=0;i<dists.size();i++) {
	    
	    
          double s = 1.0;
          for(unsigned int j=0;j<dists.size();j++)
               if (j!=i) s *= sums[j];
          dists[i] = dists[i]*s;
          
     }
    
	logf = 0.0;
}

FactoredVector::FactoredVector(const FactoredVector &v1,
		const FactoredVector &v2, int i) : dists(v1.dists) {
	uint ii = (uint)i;
	if (v1.logf==v2.logf) {
		logf = v1.logf;
		for(uint j=0;j<ii;j++)
			dists[j] = v2.dists[j];
		dists[ii] = v1.dists[ii] - v2.dists[ii];
/*		for(uint j=ii+1;j<dists.size();j++)
			dists[j] = v1.dists[j]; */
	} else if (v1.logf>v2.logf) {
		logf = v1.logf;
		double m = exp(v2.logf-v1.logf);
		for(uint j=0;j<ii;j++)
			dists[j] = v2.dists[j];
		dists[ii] = v1.dists[ii] - v2.dists[ii]*m;
		for(uint j=ii+1;j<dists.size();j++)
			dists[j] *= m;
	} else {
		logf = v2.logf;
		double m = exp(v1.logf-v2.logf);
		for(uint j=0;j<ii;j++)
			dists[j] = v2.dists[j]*m;
		dists[ii] = v1.dists[ii]*m - v2.dists[ii];
/*		for(uint j=ii+1;j<dists.size();j++)
			dists[j] = v1.dists[j]; */
	}
	
	
}

bool FactoredVector::isdiff(const FactoredVector &v1, int i) const {
	return dists[i]!=v1.dists[i];
}

void FactoredVector::makediff(const FactoredVector &v1,
		const FactoredVector &v2, int i) {
	uint ii = (uint)i;
	if (v1.logf==v2.logf) {
		logf = v1.logf;
		for(uint j=0;j<ii;j++)
			dists[j] = v2.dists[j];
		dists[ii] = v1.dists[ii] - v2.dists[ii];
		for(uint j=ii+1;j<dists.size();j++)
			dists[j] = v1.dists[j]; 
	} else if (v1.logf>v2.logf) {
		logf = v1.logf;
		double m = exp(v2.logf-v1.logf);
		for(uint j=0;j<ii;j++)
			dists[j] = v2.dists[j];
		dists[ii] = v1.dists[ii] - v2.dists[ii]*m;
		for(uint j=ii+1;j<dists.size();j++)
			dists[j] *= m;
	} else {
		logf = v2.logf;
		double m = exp(v1.logf-v2.logf);
		for(uint j=0;j<ii;j++)
			dists[j] = v2.dists[j]*m;
		dists[ii] = v1.dists[ii]*m - v2.dists[ii];
		for(uint j=ii+1;j<dists.size();j++)
			dists[j] = v1.dists[j];
	}
}

/*
FactoredVector::~FactoredVector() {

}
*/

FactoredVector * FactoredVector::Clone() const {

	return new FactoredVector(*this);
}

void FactoredVector::Swap(FactoredVector &v) {

	dists.swap(v.dists);
}

void FactoredVector::Condition(const map<int, int> &ev) {
	map<int, int>::const_iterator it;
	int id, val;
	double tmp;
	for (it = (ev).begin(); it != ev.end(); it++){
		id = (*it).first;
		val = (*it).second;
		RestrictById(id, val);
	}
	
}

double FactoredVector::Diff2(const FactoredVector &v) const {
	double ret = 0.0;
	for(unsigned int i=0;i<dists.size();i++) {
		//double s = 0.0;
		//for(int j=0;j<dists[i].getm();j++)
		//	s += abs(dists[i][j]);
		for(int j=0;j<dists[i].getm();j++) {
			double val = abs((dists[i][j]-v.dists[i][j]));///s;
			if (val>ret) ret = val;
			//ret += val;
		}
	}
	return ret;
	
/*
	double ret = 0.0;
	for(unsigned int i=0;i<dists.size();i++) {
		static const double eps = 1e-10,eps2 = 1e-10;
		double d1=0.0,d2=0.0;
		const vectr &v1 = dists[i];
		const vectr &v2 = v.dists[i];
		int m = v1.getm();
		for(int j=0;j<m;j++)
			d1 += (v1[j]<eps) ? eps2 : v1[j];
		for(int j=0;j<m;j++)
			d2 += (v2[j]<eps) ? eps2 : v2[j];
		double r = d1/d2;
		double s = 0.0;
		for(int j=0;j<m;j++) {
			if (v1[j]<eps) {
				if (v2[j]>=eps)
					s += v2[j]*::logf(v2[j]/eps2);
				else
					s += v2[j]*::logf(eps2/eps2);
			} else {
				if (v2[j]>=eps)
					s += v2[j]*::logf(v2[j]*r/v1[j]);
				else s += eps2*::logf(eps2*r/v1[j]);
			}
		}
		ret += s/d2;
	}
	ret /= dists.size();
	return ret;
	*/
}

double FactoredVector::PairSize(const FactoredVector &v) const {
	double ret = 0.0;
	double m1 = logf!=0.0 ? exp(logf) : 1.0;
	double m2 = v.logf!=0.0 ? exp(v.logf) : 1.0;
     vector<double> s1(dists.size()),s2(dists.size());
     for(unsigned int i=0;i<v.dists.size();i++) {
          s1[i] = dists[i].sum();
          s2[i] = v.dists[i].sum();
	}
	for(unsigned int i=0;i<dists.size();i++) {
		double p1=m1,p2=m2;
          for(unsigned int j=0;j<dists.size();j++)
               if (j!=i) {
				p1 *= s1[j];
				p2 *= s2[j];
			}
		for(int j=0;j<dists[i].getm();j++) {
			double val = abs(dists[i][j]*p1 - v.dists[i][j]*p2);
			if (ret<val) ret = val;
		}
	}
	return ret;
}

double FactoredVector::Diff(const FactoredVector &v) const {
	
	double ret = 0.0;
	double m1 = logf!=0.0 ? exp(logf) : 1.0;
	double m2 = v.logf!=0.0 ? exp(v.logf) : 1.0;
     vector<double> s1(dists.size()),s2(dists.size());
     for(unsigned int i=0;i<v.dists.size();i++) {
          s1[i] = dists[i].absmax();
          s2[i] = v.dists[i].absmax();
	}
	for(unsigned int i=0;i<dists.size();i++) {
		double p1=m1,p2=m2;
          for(unsigned int j=0;j<dists.size();j++)
               if (j!=i) {
				p1 *= s1[j];
				p2 *= s2[j];
			}
		for(int j=0;j<dists[i].getm();j++) {
			double val = abs(dists[i][j]*p1 - v.dists[i][j]*p2);
			if (ret<val) ret = val;
		}
	}
	return ret;

}

void FactoredVector::Restrict(const vector<int> &) {
	cout << "FactoredVector::Restrict not implemented" << endl;
}

void FactoredVector::RestrictById(int id, int val) {
	double tmp = dists[id][val];
	dists[id] = 0.0;
	dists[id][val] = tmp;
}

void FactoredVector::RestrictById(int id, int val, double x) {
	dists[id] = 0.0;
	dists[id][val] = x;
}

void FactoredVector::MultBy(const RVSimple *x) {
	const FactoredVector *v = dynamic_cast<const FactoredVector *>(x);
	for(unsigned int i=0; i<dists.size(); i++) {
		dists[i].multby(v->dists[i]);
	}
}

void FactoredVector::Zero() {
	for(uint i=0;i<dists.size();i++)
		dists[i] = 0.0;
	logf = 0.0;
}

void FactoredVector::DotStar(const FactoredVector &v, FactoredVector &out) const{
	for(unsigned int i=0; i<dists.size(); i++)
		out.dists[i] = dists[i].dotstar(v.dists[i]);
	out.logf = logf+v.logf;
}

void FactoredVector::MakeUniform() {
	for (unsigned int i=0; i<dists.size(); i++){
		dists[i] = 1.0;
	}
	logf = 0.0;
	
}

double FactoredVector::Normalize() {

	double mofsum, tmp, diff;
	double eps = 0.00001;

	tmp = dists[0].normalize();
	mofsum = tmp;
	for (unsigned int i=1; i<dists.size(); i++){
		tmp = dists[i].normalize();

		diff = tmp - mofsum;
		if (!((diff < eps)  && (diff > - eps))){
			mofsum *= tmp;

		}
	
	}
	double ret= logf!=0.0 ? exp(logf)*mofsum: mofsum;
	logf = 0.0;
	return ret;
	
	/*
	double tmp = 1.0;	
	for (unsigned int i=0; i<dists.size(); i++){
		tmp *= dists[i].normalize();	
	}
	
	double ret = exp(logf)*tmp;
	logf = 0.0;
	return ret;
	*/
}

double FactoredVector::absnormalize(vectr &v) {
	double tmpsum = 0.0;
	for (int i=0; i<v.length(); i++){
		tmpsum += fabs(v[i]);
	}
	for (int i=0; i<v.length(); i++){
		v[i] = fabs(v[i]) / tmpsum;
	}
	return tmpsum;
}

double FactoredVector::maxnormalize(vectr &v) {
	double max = v[0];
	double min = v[0];
	for (int i=0; i<v.length(); i++){
		if (max < v[i])
			max = v[i];
		if (min > v[i])
			min = v[i];
	}
	for (int i=0; i<v.length(); i++){
		v[i] = (v[i] - min)/ (max-min);
	}
	return max-min;
}

double FactoredVector::NormalizeNeg(){
	double mofsum, diff;
//  double tmp;
	double eps = 0.00001;
	mofsum = 1.0;

	for (unsigned int i=0; i<dists.size(); i++) {
		mofsum *= dists[i].sum();
		for (int j=0; j<(dists[i]).length(); j++) {
			if ((dists[i])[j] < 0.0){

				dists[i][j] = 0.0;
			}
		}
		
//		tmp = dists[i].normalize();	
	}			


	double ret= logf!=0.0 ? exp(logf)*mofsum: mofsum;

	logf = 0.0;
	return ret;
	
}


double FactoredVector::GetLogf() const {
	return logf;
}


void FactoredVector::SetLogf(double logfactor) {
	logf = logfactor;
}
	

double FactoredVector::ProdofAbsSum(bool islog) const {
	double nc = 0.0;
	double tmps;
	//cout << "logf: " << logf << endl;
	for (unsigned int i=0; i<dists.size(); i++) {
		//nc += log(dists[i].norm2()); l2 square
		tmps = 0.0;
	    for(int j=0; j<dists[i].getm(); j++){
	        tmps += ( abs((dists[i])[j]) );
            //cout <<"d:"<< (dists[i])[j] <<" absd:"<< abs((dists[i])[j]) << " log:"<<log( abs((dists[i])[j]) ) << " nc:"<<nc<<endl;
	    }
	    nc += log(tmps);
		//cout << "node:" << i << " " << abs(dists[i].sum()) << endl;
	}
	//cout << "nc " << nc <<endl;
	if (islog) {
	    return (logf+nc);
	} else {
		return (exp(logf+nc));
	}
	
}

double FactoredVector::MaxofAbsSum(bool islog) const {
	double nc = 0.0;
	double tmps;
	for (unsigned int i=0; i<dists.size(); i++) {
		tmps = 0.0;
	    for(int j=0; j<dists[i].getm(); j++){
	        tmps += ( abs((dists[i])[j]) );
            //cout <<"d:"<< (dists[i])[j] <<" absd:"<< abs((dists[i])[j]) << " log:"<<log( abs((dists[i])[j]) ) << " tmps:"<<tmps<<endl;
	    }
	    if(nc < tmps){
	    	nc = tmps;
	    }
		//cout << "node:" << i << " " << nc << endl;
	}
	if (islog) {
	    return (logf+log(nc));
	} else {
		return (exp(logf)*nc);
	}
	
}

double FactoredVector::ProdofMaxAbs(bool islog) const {
	double nc = islog ? 0.0 : 1.0;
	double tmps;
	for (unsigned int i=0; i<dists.size(); i++) {
		tmps = 0.0;
	    for(int j=0; j<dists[i].getm(); j++){
	    	if (tmps < abs((dists[i])[j]) )
	    		tmps = abs((dists[i])[j]);
            	//cout <<"d:"<< (dists[i])[j] <<" absd:"<< abs((dists[i])[j]) << " log:"<<log( abs((dists[i])[j]) ) << " tmps:"<<tmps<<endl;
	    }
		if (islog) nc += log(tmps);
		else nc *= tmps;
		//cout << "node:" << i << " " << nc << endl;
	}
	if (islog) {
	    return (logf+nc);
	} else {
		return logf==0.0 ? nc : exp(logf)*nc;
		//return (exp(logf+nc));
	}
	
}

double FactoredVector::MaxofMaxAbs(bool islog) const {
	double nc = 0.0;
	double tmps;
	for (unsigned int i=0; i<dists.size(); i++) {
		tmps = 0.0;
	    for(int j=0; j<dists[i].getm(); j++){
	    	if (tmps < abs((dists[i])[j]) )
	    		tmps = abs((dists[i])[j]);
            	//cout <<"d:"<< (dists[i])[j] <<" absd:"<< abs((dists[i])[j]) << " log:"<<log( abs((dists[i])[j]) ) << " tmps:"<<tmps<<endl;
	    }
	    if(nc < tmps){
	    	nc = tmps;
	    }
		//cout << "node:" << i << " " << nc << endl;
	}
	if (islog) {
	    return (logf+log(nc));
	} else {
		return logf==0.0 ? nc : (exp(logf)*nc);
	}
	
}

/*
double FactoredVector::L2Sum(bool islog) const {
	double nc = 0.0;
	for (unsigned int i=0; i<dists.size(); i++) {
		nc += log(dists[i].sum());
	}
	if (islog) {
	    return (logf+nc);
	} else {
		return (exp(logf+nc));
	}
	
}
*/

double FactoredVector::AddMultDotStar(const FactoredVector &v1,
			const FactoredVector &v2, double x) {

	if (logf!=0.0) {
          double m = exp(logf);
          for(unsigned int i=0;i<dists.size();i++)
               dists[i] *= m;
		logf = 0.0;
     }
     double m = v1.logf!=0.0 ? exp(v1.logf)*x : x;
	if (v2.logf!=0.0) m *= exp(v2.logf);
     vector<double> sums(dists.size(),0.0);
	double prod = 1.0;
	uint skip = 0;
	double ret = 0.0;
     for(unsigned int i=0;i<dists.size();i++) {
		const vectr &vv1 = v1.dists[i];
		const vectr &vv2 = v2.dists[i];
		for(int j=0;j<vv1.size();j++)
			sums[i] += vv1[j]*vv2[j];
		if (sums[i]==0.0) {
			if (skip>0) return 0.0;
			skip = i+1;
		} else prod *= sums[i];
	}
	if (skip>0) {
		const vectr &vv1 = v1.dists[skip-1];
		const vectr &vv2 = v2.dists[skip-1];
		vectr &d = dists[skip-1];
		for(int j=0;j<vv1.size();j++) {
			double v = vv1[j]*vv2[j]*m*prod;
			if (ret<abs(v)) ret = abs(v);
			d[j] += v;
		}
	} else {
		for(unsigned int i=0;i<dists.size();i++) {
			// hopefully accurate enough
			double sprod = m*prod/sums[i];
				
			const vectr &vv1 = v1.dists[i];
			const vectr &vv2 = v2.dists[i];
			vectr &d = dists[i];
			for(int j=0;j<vv1.size();j++) {
				double v = vv1[j]*vv2[j]*sprod;
				if (ret<abs(v)) ret = abs(v);
				d[j] += vv1[j]*vv2[j]*sprod;
			}
		}
	}
	return ret;
}
	
//Assumes both vectors have the same length
double FactoredVector::AddMult(const FactoredVector & v, double x) {

	double ret;
	//double ret2=0.0;
    //cout << "--0\n";
	if (logf!=0.0) {
      //      cout << "-- 0-1\n";
          double m = exp(logf);
          for(unsigned int i=0;i<dists.size();i++)
               dists[i] *= m;
		logf = 0.0;
     }
     double m = v.logf!=0.0 ? exp(v.logf)*x : x;
     ret = m;
     vector<double> sums(dists.size());
	double prod = 1.0;
	uint skip = 0;
    //cout << "--1\n";
     for(unsigned int i=0;i<v.dists.size();i++) {
      //   cout << ".. 1\n";
		double isum = v.dists[i].sum();
        //cout << ".. 2\n";
        //cout << "i: " << i << " , sumsize: " << sums.size() << endl;
		sums[i] = isum;
        //cout << ".. 3\n";
		if (isum==0.0) {
			if (skip>0) return 0.0;
			skip = i+1;
		} else prod *= isum;
	}
    //cout << "-- 2\n";
	if (skip>0)
		dists[skip-1] += v.dists[skip-1]*m*prod;
	else {
		for(unsigned int i=0;i<dists.size();i++) {
/*
			double sprod = m;
			for(unsigned int j=0;j<dists.size();j++)
				if (j!=i) sprod *= sums[j];
*/
			// hopefully accurate enough
			double sprod = m*prod/sums[i];
				
			dists[i] += v.dists[i]*sprod;
			//cout << "sprod" << sprod << endl;
			//cout<<(v.dists[i])<<endl;
			
			double delta = v.dists[i].absmax();
			//double delta2 = delta*sprod;
			//if (delta2>ret2) ret2 = delta2;
			ret *= delta;
				
		}
	}
	//cout << ret << endl;
	return ret;
	/*
	const FactoredVector * a = dynamic_cast<const FactoredVector *>(&v);
	double difflogf, maxlogf;
	double alogf = a->logf + log(x);
	if(logf < alogf){
		for (unsigned int i=0; i<dists.size(); i++) {
		dists[i] = (this->dists[i])/expf + a->dists[i];
		dists[i] = dists[i]/(1+1/expf);
		}
	}
	else{
		for (unsigned int i=0; i<dists.size(); i++) {
		dists[i] = (minfv->dists[i])/expf + maxfv->dists[i];
		dists[i] = dists[i]/(1+1/expf);
		}
	}
	*/
	/*const FactoredVector * a = dynamic_cast<const FactoredVector *>(&v);
	
	const FactoredVector *minfv, *maxfv; 
	double difflogf, maxlogf;
	double alogf = a->logf + log(x);
	if(logf < alogf){
		minfv = this;
		maxfv = a;
		difflogf = alogf - logf;
		maxlogf = alogf;
	}
	else{
		minfv = a;
		maxfv = this;
		difflogf = logf - alogf;
		maxlogf = logf;
	}

	double expf = exp(difflogf);
	for (unsigned int i=0; i<dists.size(); i++) {
		dists[i] = (minfv->dists[i])/expf + maxfv->dists[i];
		dists[i] = dists[i]/(1+1/expf);
	}
	logf = maxlogf + log(1+1/expf);
	*/
}
	
void FactoredVector::Add(const RVSimple *rv) {
	const FactoredVector * v = dynamic_cast<const FactoredVector *>(rv);
	if (logf!=0.0) {
          double m = exp(logf);
          for(unsigned int i=0;i<dists.size();i++)
               dists[i] *= m;
     }
     double m = v->logf!=0.0 ? exp(v->logf) : 1.0;
     vector<double> sums(dists.size());
     for(unsigned int i=0;i<v->dists.size();i++)
          sums[i] = v->dists[i].sum();
     for(unsigned int i=0;i<dists.size();i++) {
          double sprod = m;
          for(unsigned int j=0;j<dists.size();j++)
               if (j!=i) sprod *= sums[j];
          dists[i] += v->dists[i]*sprod;
     }
     logf = 0.0;
	/*const FactoredVector * a = dynamic_cast<const FactoredVector *>(x);
	
	const FactoredVector *minfv, *maxfv; 
	double difflogf, maxlogf;
	if(logf < a->logf){
		minfv = this;
		maxfv = a;
		difflogf = a->logf - logf;
		maxlogf = a->logf;
		
	}
	else{
		minfv = a;
		maxfv = this;
		difflogf = logf - a->logf;
		maxlogf = logf;
	}

	double expf = exp(difflogf);
	for (unsigned int i=0; i<dists.size(); i++) {
		dists[i] = (minfv->dists[i])/expf + maxfv->dists[i];
		dists[i] = dists[i]/(1+1/expf);
	}
	logf = maxlogf + log(1+1/expf);
	
	*/
	
}

void FactoredVector::Mult(double x) {

   	if(x > 0){
	    logf += log(x);
	    /*
	   	for (unsigned int i=0; i<dists.size(); i++) {
		    dists[i] *= x;
	    }
	    */
	}
	else if( x == 0){
		logf = 0;
	 	for (unsigned int i=0; i<dists.size(); i++) {
		    dists[i] *= 0;
	    }
	}	
	else {
		cout << "Multiplying with < 0 !!!" << endl;
		cout << __LINE__ << " logf " << logf << endl;
		/*
	    for (unsigned int i=0; i<dists.size(); i++) {
		    dists[i] *= x;
	    }
	    */
	}

}

void FactoredVector::Scale(vector<double> &sc){
    
    double scale = 1.0;
    for (unsigned int i=0; i<sc.size(); i++) {
	    scale += sc[i];
	    cout << "sc " << i << ": " << sc[i] << endl;
	}
	cout << "scale: "<< scale<<endl;
    for (unsigned int i=0; i<dists.size(); i++) {
	    cout << "Before ";
	    dists[i].niceprint(cout);
	    dists[i] /= exp(scale);
	    cout << "\nAfter ";
	    dists[i].niceprint(cout);
	    
	}
}

double FactoredVector::CompareDist(FactoredVector &fv, int varid) const {
    unsigned varid_as_uint(varid);
	if (/*varid_as_uint >= 0 &&*/ varid_as_uint < dists.size()) {
		vectr tmp = fv.dists[varid_as_uint] - dists[varid_as_uint];
		return tmp.norm();
	}
	else{
		cerr << "Index out of bounds at FactoredVector" << endl;
		return -1;
	}
	
}

void FactoredVector::GetDist(vectr &v, int varid, double &logfactor) {
    unsigned varid_as_uint(varid);
	if (/*varid_as_uint >= 0 &&*/ varid_as_uint < dists.size()) {
		v = dists[varid_as_uint];
		logfactor = logf;
	}
}

	
void FactoredVector::SetDist(const vectr &v, int varid, double logfactor) {
    unsigned varid_as_uint(varid);
    if ( varid_as_uint < dists.size()){
		dists[varid_as_uint] = v;	
		logf = logfactor;
	}
//	cout << __LINE__ << " logf " << logf << endl;
    /*
	v.niceprint(cout);
	double vsum = v.sum();
	bool allzero = false;
	if (vsum == 0){
		allzero = true;
		cout << "sum zero" << endl;
		for(int j=0; j<v.getm(); j++){
			if(v[j] > 0)	allzero = false;
		}
	}
		
		
	if ( varid_as_uint < dists.size()) {//varid_as_uint >= 0 &&
		if (allzero){	
			dists[varid_as_uint] = v+1e-10;	
			vsum = 2*(1e-10);
		}
		else
			dists[varid_as_uint] = v;
			
		logf = logfactor;
		
		if(vsum > 0){
			for(uint i=0; i<dists.size(); i++){
				if(i != varid_as_uint)	dists[i] = dists[i]*vsum; 
			}
		}
	}
	*/
		
}

double FactoredVector::GetDistInstNorm(int varid, int index) {
    double dsum = dists[varid].sum();
    if(dsum > 0)
        return (dists[varid])[index]/dsum;
    else
	    return (dists[varid])[index];
}

double FactoredVector::GetDistInst(int varid, int index) {
	return (dists[varid])[index];
}

void FactoredVector::SetDistInst(int varid, int index, double value) {

	(dists[varid])[index] = value;
}

double FactoredVector::GetMargMin() {
	double minval = GetMin(0);
	double tmp = 0.0;
	for(unsigned int i=1; i<dists.size(); i++) {
		tmp = GetMin(i);
		if (minval > tmp) minval = tmp;
	}
	return minval;
} 

double FactoredVector::GetJointMin() {
	double minval = 1.0;
	for (unsigned int i=0; i<dists.size(); i++) {
		minval *= GetMin(i);
	}
	return minval;
}

double FactoredVector::GetMin(int varid) {
	double minval = fabs((dists[varid])[0]);
	for (int i=0; i<dists[varid].length(); i++) {
		if (minval < fabs((dists[varid])[i]))
			minval = (dists[varid])[i];
	}
	return minval;
}

double FactoredVector::GetMargAbsMax() {
	double maxval = dists[0].absmax();
	double tmp = 0.0;
	for(unsigned int i=1; i<dists.size(); i++) {
		tmp = dists[i].absmax();
		if (maxval < tmp) maxval = tmp;
	}
	return maxval;
} 

double FactoredVector::GetMargMax() {
	double maxval = GetMax(0);
	double tmp = 0.0;
	for(unsigned int i=1; i<dists.size(); i++) {
		tmp = GetMax(i);
		if (maxval < tmp) maxval = tmp;
	}
	return maxval;
} 

bool FactoredVector::iszero() const {
	for(unsigned int i=1; i<dists.size(); i++) {
		bool r = true;
		for(int j=0;j<dists[i].length();j++)
			if (abs(dists[i][j])>1e-16) { r = false; break; }
		if (r) return true;
	}
	return false;
}

double FactoredVector::GetJointAbsMax() {
	double maxval = 1.0;
	for(unsigned int i=0; i<dists.size(); i++)
		maxval *= dists[i].absmax();
	return maxval;
}

double FactoredVector::GetJointMax() {
	double maxval = 1.0;
	for(unsigned int i=0; i<dists.size(); i++) {
		maxval *= GetMax(i);
	}
	return maxval;
}

double FactoredVector::GetMax(int varid) {
	double maxval = fabs((dists[varid])[0]);
	for (int i=0; i<dists[varid].length(); i++) {
		if (maxval > fabs((dists[varid])[i]))
			maxval = (dists[varid])[i];
	}
	return maxval;
}

int FactoredVector::GetDistSize(int varid) {
	return (dists[varid]).length();
}

vectr FactoredVector::Joint(const FactoredVector &fv){
	vectr ret(fv.dists[0]);
	for(uint i=1; i<fv.dists.size(); i++){
		ret = outer(ret,fv.dists[i]);
	}
	return ret;
}

int FactoredVector::Size() { return dists.size(); }

int FactoredVector::JointSize() { 
   int jsize = 1;
   for(uint i=0; i<dists.size(); i++){
       jsize *= dists[i].getm();     
   }
   return jsize;
}

ostream & FactoredVector::Print(ostream & out, const FactoredVector &v2) const {
	for (unsigned int i=0; i<dists.size(); i++){
		stringstream ss;
		vectr tmp = dists[i]*exp(logf);
		/*tmp.niceprint(ss);
		ss.seekp(-1,ios_base::end);
		for(int j=ss.str().length();j<35;j++)
			ss << ' ';
		*/
		tmp.niceprint(ss);
		string sss = ss.str();
		sss.erase(sss.size()-1,1);
		out << sss;
		for(uint j=0;j<35-sss.size();j++) out << ' ';
		tmp = v2.dists[i]*exp(v2.logf);
		tmp.niceprint(out);

		////tmp.normalize();
		//tmp.niceprint(out);	
	}
	//out << endl;
	return out;
}

ostream & FactoredVector::Print(ostream & out, RV * ignored = nullptr03) const {
	for (unsigned int i=0; i<dists.size(); i++){
		stringstream ss;
		vectr tmp = dists[i]*exp(logf);
		/*tmp.niceprint(ss);
		ss.seekp(-1,ios_base::end);
		for(int j=ss.str().length();j<35;j++)
			ss << ' ';
		*/
		if (ignored==nullptr03) tmp.normalize();
		tmp.niceprint(ss);
		string sss = ss.str();
		sss.erase(sss.size()-1,1);
		out << sss << endl; //" (" << dists[i].sum() << ")" << endl;
		

		////tmp.normalize();
		//tmp.niceprint(out);	
	}
	//out << endl;
	return out;
}


ostream & FactoredVector::PrintSimple(ostream & out) const {
    double logf = this->GetLogf();
    out << "FactoredVector logf: " << logf << endl;
    this->Print(out);
    return out;
}


//From RVSimple
void FactoredVector::GetDist(vectr &d, double &logfactor) const {
    throw not_yet_implemented_error(the_class_name, "GetDist");
}


void FactoredVector::SetDist(const vectr &d, double logfactor) {
    throw not_yet_implemented_error(the_class_name, "SetDist");
}


double FactoredVector::Prob(int ind, bool log) const{
    throw not_yet_implemented_error(the_class_name, "Prob");
}


void FactoredVector::Load(istream &is) {
    throw not_yet_implemented_error(the_class_name, "Load");
}

void FactoredVector::Save(ostream &os) const {
    throw not_yet_implemented_error(the_class_name, "Save");
}

void FactoredVector::Reindex(const vector<vector<int> > &ind) {
    throw not_yet_implemented_error(the_class_name, "Reindex");
}

void FactoredVector::Add(const RVSimple *x, double w) {
    throw not_yet_implemented_error(the_class_name, "Add");
}

SS *FactoredVector::BlankSS() const {
    throw not_yet_implemented_error(the_class_name, "BlankSS");
}

void FactoredVector::AddSS(int x, SS *ss, double w) const {
    throw not_yet_implemented_error(the_class_name, "AddSS");
}

void FactoredVector::AddExpSS(SS *ss, double w) const {
    throw not_yet_implemented_error(the_class_name, "AddExpSS");
}

void FactoredVector::AddSS(const SS *toadd, const RVSimple* rvs,
			const vector<vector<int> > &mapping, 
			SS *ss, double w) const {
    throw not_yet_implemented_error(the_class_name, "AddSS");
}

int FactoredVector::Sample(Random &rand) const {
    throw not_yet_implemented_error(the_class_name, "Sample");
}

void FactoredVector::Maximize(const SS *ss) {
    throw not_yet_implemented_error(the_class_name, "Maximize");
}

void FactoredVector::Scramble(double alpha, double degree, Random &rand) {
    throw not_yet_implemented_error(the_class_name, "Scramble");
} 

double FactoredVector::LLH(const SS *ss) const {
    throw not_yet_implemented_error(the_class_name, "LLH");
}

double FactoredVector::GetScore(double numTrans, const SS* ss) const {
    throw not_yet_implemented_error(the_class_name, "GetScore");
}

}

