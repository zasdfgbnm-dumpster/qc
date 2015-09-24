#ifndef __HPP_QUANTUM_MAT__
#define __HPP_QUANTUM_MAT__

/* This header file defines matrix and its operations
 * 
 * The implementation now is just an adapter of the
 * "Eigen" library. The reason for not directly using
 * "Eigen" is to make sure that these codes are matrix
 * library independent, so that it's easy to switch to
 * other library if we want to.
 */

#include <Eigen/Eigen>
namespace quantum_chem{

/* class matrix:
 * the class for general matrix
 */
class matrix {
protected:
	Eigen::MatrixXd mat;
public:
	matrix(int rows,int cols):mat(rows,cols){}
	virtual int n_rows(){ return mat.rows(); }
	virtual int n_columns(){ return mat.cols(); }
	virtual double &operator()(int i,int j){ return mat(i,j); }
	virtual matrix left_columns(int n){ return mat.leftCols(n); }
	virtual matrix conjugate_transpose(){ return mat.adjoint(); }
	virtual double norm(){ return mat.norm(); }
};

/* class hermitian_matrix:
 * the class for Hermitian matrix
 * Storage may be optimize considering that the matrix
 * is Hermitian.
 */
class hermitian_matrix : public matrix {
public:
	hermitian_matrix(int n):matrix(n,n){}
};

/* eigen problem solver for Hermitian matrix
 * Hv=uSv
 */
void generalized_hermitian_eigen_solver(const hermitian_matrix &H, const hermitian_matrix &S, 
   matrix &vs, array<double> &eigenvalues) {
	Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> es(H,S);
	vs = es.eigenvectors()
	for(int i=0;i<H.n_cols();i++)
		eigenvalues[i] = es.eigenvalues()[i];
	return;
}

/* identical matrix I=diag(1,1,1,...)
 */
hermitian_matrix idmat(int size){
	return Eigen::MatrixXd::Identity(size);
}

}

#endif
