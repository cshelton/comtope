//============================================================================
// Author      : Sepideh Azarnoosh
// Description : Calulating EigenValues and EigenVectors of a matrix!
//============================================================================

#include <iostream>
#include "Eigen/Eigenvalues"
#include "eigendecompose.h"

namespace ctbn {

using namespace std;
using namespace Eigen;

    EigenDecompose::EigenDecompose(MatrixXd M2):M(M2){}

VectorXcd EigenDecompose::getEigenValues() const{
	EigenSolver<MatrixXd> es(M, false);
	return es.eigenvalues();
}
    
MatrixXcd EigenDecompose::getRightEigenVectors() const{
        EigenSolver<MatrixXd> es2(M);
        return es2.eigenvectors();
    }
    
MatrixXcd EigenDecompose::getLeftEigenVectors() const{
        EigenSolver<MatrixXd> es3(M);
        EigenSolver<MatrixXd> est(M.transpose());
        MatrixXcd left = est.eigenvectors();
        MatrixXcd temp = left.transpose() * es3.eigenvectors();
        VectorXcd temp2 = temp.col(0);
        complex<double> mycomp;
        mycomp = 1.0;
        //cout << "mycomp" << mycomp << endl;
        mycomp = mycomp/temp2[0];
        MatrixXcd result = left * mycomp;
        return result;
    }

}

