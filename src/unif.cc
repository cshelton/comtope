/* Continuous Time Bayesian Network Reasoning and Learning Engine
 * Copyright (C) 2009 The Regents of the University of California
 *
 * Authored by Yu Fan, William Lam, Joon Lee, Christian Shelton, and Jing Xu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "unif.h"
#include "matrix.h"
#include <map>
#include <cmath>
#include <vector>

namespace ctbn{

using namespace std;
double iuvexpmt(vectr &a, vectr &ind, const matrix &Q, double t, int l, int limit){

	bool print = false;
	double alpha, theta, s, r, gamma, m;
	//int  l= 10; //20, 40.. //MAX_L;
	double wnorm = a.normalize();
	vectr f, w, indf, indw;

	//Compute norm of Q
	alpha = fabs(Q.min()); 
	cout << "ALPHA" << alpha << endl;
	//int n = Q.getm();
	//matrix I(n,n,vectr(n, 1.0));
	//matrix P = (Q / alpha) + I;

	//Loop parameters
	theta = 1;
	l = 5;
	//if (alpha*t <= theta)
	//	theta = 1;

	m = ceil(alpha*t*theta/l);
	//bu sekilde orjinalinde ama biz degistirdik: m = ceil(alpha*t/theta);
	//originally produced in this way, but we ???: m = ceil (alpha * t / theta);
	t = t/m;
	s = alpha*t;
	r = exp(-s);

	//Q.niceprint(cout);
	//a.niceprint(cout);
	//ind.niceprint(cout);
	if (a.length() > 1){
	
		//sort ind into ascending order and change a accordingly - not necessary since ind is already sorted?
		/*
		int maxind = (int)ind.max();
		vectr full(maxind+1, 0.0);
		for(int i=0; i<ind.length(); i++){
			full[(int)ind[i]] = a[i];
			if (ind[i] > maxind) 
				maxind = (int)ind[i];
		}
		int j=0;
		for(int i=0; i<maxind+1; i++){
			if (full[i] > 0){
				
				ind[j] = i;
				a[j] = full[i];
				j++;
			}
		}
		*/
/* requires matrix-noeigen
		if (print){
			cout <<"\nInitial distribution"<<endl;
			(a/a.max()).niceprint(cout);
		}
*/
		w = a;
		indw = ind;
//		double ret=0;
		wnorm += log(w.normalize());
		for (int i=1; i<=m; i++){
	
			if (print){
				cout <<"Uniformization: m= "<<i<<endl;
			}
			f = w;
			indf = indw;
			gamma = 1;		
			for (int k=1; k<=l; k++){
				if (print)	cout <<"\tUniformization: l= "<<k<<endl;
				
				gamma = gamma * s/k;

/* requires matrix-noeigen
				if (print){	
					cout <<"\n\tF before multiplication"<<endl;
					(f/f.max()).niceprint(cout);
				}
*/
//				ret = imult(f, Q, limit);
				
				if (print){
/* requires matrix-noeigen
					cout <<"\n\tF after multiplication"<<endl;
					(f/f.max()).niceprint(cout);
					cout <<"\n\tW before addition"<<endl;
					(w/w.max()).niceprint(cout);
*/
					cout << "GAMMA " << gamma << endl;
				}
				addinplace(w, vectr(f*gamma));
				
/* requires matrix-noeigen
				if (print){
					cout <<"\n\tW after addition"<<endl;
					(w/w.max()).niceprint(cout);
				}
*/
			}
				
			w *= r; //= r*w;
			//normalize w
			wnorm += log(w.normalize());

/* requires matrix-noeigen
			if (print){
				cout <<"\nW after the loop and normalization"<<endl;
				(w/w.max()).niceprint(cout);
			}
*/
		}
		a = w;
		ind = indw;
	}

	return wnorm;
}

double iuexpmtv(vectr &a, vectr &ind, const matrix &Q, double t, int l, int limit){
	const matrix Qtrans(Q, true);
	return iuvexpmt(a, ind, Qtrans, t, l, limit);
}

double imult(vectr &a, const matrix &Q, int limit){

	bool print = false;
	
	//double eps = 0.0000001;
	int n = Q.getm();
	double alpha = fabs(Q.min());
	
	vectr tmpa(n, 0.0);
	double tmp = 0.0;
	double rate = 0.0;
	int k=0;
	
	double eps = 0;
	int ln= a.length();
	if (limit < ln) {	
		vector<double> sorted;
		sorted.resize(ln);
		for (int i=0; i<ln; i++){
			sorted[i] = a[i];
		}
		sort(sorted.begin(), sorted.end());
		eps = sorted[ln-limit];
		//eger birden fazla ayni degerli eleman varsa olmuyor - ama zaten double oldugundan herhalde pek olmaz
		// If you have more than one element if it is not the same value - but many probably do not already have double
	}

	for (int j=0; j<n; j++){
		tmp = 0.0;
		for (int i=0; i<a.length(); i++){
			if ( (a[i] > 0) && (a[i] >= eps) && (fabs(Q[i][j]) > 0.0) ){
				if (j==i){
					tmp += a[i]*(Q[i][j]/alpha + 1.0);
					
					if(print){	
						cout << "HERE Multiplyying the diagonal " << i << "->" << i << " : " << a[i] << "*Q(" <<  Q[i][j] << "/" << alpha << "+ 1)"<< endl; 
						cout << "= " << a[i]*(Q[i][j]/alpha + 1.0) << endl;
					}
				}
				else
					tmp += a[i]*(Q[i][j]/alpha);
			}
		
		}
		if (tmp > 0.0){
			tmpa[j] = tmp;
		}
	}
	
	a = tmpa;
    // scaling factor was unchanged
    return 1.0;
}

//su an calismiyor, ayrica eps'yi hesaplamayi eklemek lazim limitten.
// Currently not working, the limit also necessary to add the EPS calculation.
double imult(vectr &a, vectr &ind, const matrix &Q, int limit){
	//double eps = 0.0000001;
	int n = Q.getm();
	double alpha = fabs(Q.min());

	double eps = 0;
	vector<double> tmpa;
	vector<int> tmpind;
	double tmp = 0.0;
	int k=0;
	for (int j=0; j<n; j++){
		tmp = 0.0;
		for (int i=0; i<ind.length(); i++){
			if ( a[i] > eps && fabs(Q[(int)ind[i]][j]) > 0.0){
				if (j== ind[i])
					tmp += a[i]*(Q[(int)ind[i]][j]/alpha + 1.0);
				else
					tmp += a[i]*(Q[(int)ind[i]][j]/alpha);
			}
		
		}
		if (tmp > 0.0){
			tmpa.push_back(tmp);
			tmpind.push_back(j);
			k++;
		}
	}
	
	if (k > 0){
		a.resize(k);
		ind.resize(k);
		for(int i=0; i<k; i++){
			a[i] = tmpa[i];
			ind[i] = tmpind[i];
		}
		return 0;
	}
	else{
		return -1;
	}
}


void addinplace(vectr &d,const vectr &s){
	bool print = false;
	
	int dn = d.length();
	int sn = s.length();

	if (print) {
		cout << "==adding w: " << endl;
			d.niceprint(cout);
		cout << "==adding f*gamma: " << endl;
			s.niceprint(cout);
	}
	
	if (dn == sn) {
		for(int i=0; i<dn; i++){
			d[i] += s[i];
		}
	}
	else { 
		cout << "Cannot add in place" << endl;
		d.niceprint(cout);
		s.niceprint(cout);
	}
	
}
void addinplace(vectr &d, vectr &indd, const vectr &s, const vectr &inds){
	int dn = d.length();
	int sn = s.length();

	vectr tmpd(dn+sn, 0.0);
	vectr tmpind(dn+sn, 0.0);
	int k =0;
	int i=0, j=0; 
	while (i<dn && j<sn){
		if ((int)indd[i] == (int)inds[j]) {
			tmpind[k] = (int)indd[i];
			tmpd[k] = d[i]+s[i];
			 i++;
			 j++;
			 k++;
		}
		else if ((int)indd[i] < (int)inds[j]){
			tmpind[k] = (int)indd[i];
			tmpd[k] = d[i];
			 i++;
			 k++;		
		}
		else{
			tmpind[k] = (int)inds[i];
			tmpd[k] = s[i];
			 j++;
			 k++;	
		}
	}

	tmpind.resize(k);
	tmpd.resize(k);
	d = tmpd;
	indd = tmpind;
	
}
}
