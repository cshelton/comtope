//============================================================================
// Author      : Sepideh Azarnoosh
// Description : Calulating EigenValues and EigenVectors of a matrix!
//============================================================================
#include <iostream>
#include "Eigen/Eigenvalues"

namespace ctbn {

using namespace std;
using namespace Eigen;
    
class EigenDecompose {
public:
    EigenDecompose(MatrixXd M2);
    VectorXcd getEigenValues() const;
    MatrixXcd getRightEigenVectors() const;
    MatrixXcd getLeftEigenVectors() const;
    
    MatrixXd MyTranspose(MatrixXd M2){return M2.transpose();}

private:
    MatrixXd M;
    
};
}
