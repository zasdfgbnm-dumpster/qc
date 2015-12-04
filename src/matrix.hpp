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
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <iostream> //delete


namespace quantum_chem{

/* class matrix:
 * the class for general matrix
 */
class matrix {
//protected:
public:
	matrix(Eigen::MatrixXd mat):mat(mat){}
public:
	Eigen::MatrixXd mat;
	matrix(int rows,int cols):mat(rows,cols){}
	matrix(){}
	virtual int n_rows() const { return mat.rows(); }
	virtual int n_columns() const { return mat.cols(); }
	virtual double &operator()(int i,int j){ return mat(i,j); }
	virtual double operator()(int i,int j) const { return mat(i,j); }
	virtual matrix left_columns(int n) const { return matrix(mat.leftCols(n)); }
	virtual matrix conjugate_transpose() const { return matrix(mat.adjoint()); }
	virtual double norm() const { return (mat.lpNorm<Eigen::Infinity>()); }
	matrix operator-(matrix rhs) const {
		matrix ret = *this;
		ret.mat -= rhs.mat;
		return ret;
	}
	matrix operator+(matrix rhs) const {
		matrix ret = *this;
		ret.mat += rhs.mat;
		return ret;
	}
	matrix operator*(matrix rhs) const {
		matrix ret = *this;
		ret.mat *= rhs.mat;
		return ret;
	}
};

/* class hermitian_matrix:
 * the class for Hermitian matrix
 * Storage may be optimize considering that the matrix
 * is Hermitian.
 */
class hermitian_matrix : public matrix {
protected:
	hermitian_matrix(Eigen::MatrixXd mat):matrix(mat){}
public:
	hermitian_matrix(){}
	hermitian_matrix(int n):matrix(n,n){}
	static void lowdin_diagonalization(const hermitian_matrix &H, const hermitian_matrix &S, matrix &vs, std::vector<double> &eigenvalues);
	static void eigen_solver(const hermitian_matrix &A, matrix &vs, std::vector<double> &eigenvalues);
	static hermitian_matrix idmat(int size);
	std::tuple<hermitian_matrix,hermitian_matrix> pmsqrt() const;
	hermitian_matrix operator=(matrix m){
		m.mat += m.conjugate_transpose().mat;
		m.mat /= 2;
		mat = m.mat;
		return *this;
	}
};

/* eigen problem solver for Hermitian matrix
 * Hv=uSv
 */
void hermitian_matrix::lowdin_diagonalization(const hermitian_matrix &H, const hermitian_matrix &S, matrix &vs, std::vector<double> &eigenvalues) {
	// Lowdin diagonalization
	eigenvalues.resize(H.n_columns());

	std::tuple<hermitian_matrix,hermitian_matrix> S_sqrt = S.pmsqrt();
	hermitian_matrix &sqrt = std::get<0>(S_sqrt);
	hermitian_matrix &nsqrt = std::get<1>(S_sqrt);
	hermitian_matrix H2(S.n_columns());
	H2.mat = nsqrt.mat*H.mat*nsqrt.mat;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H2.mat);
	vs = nsqrt*matrix(es.eigenvectors());
	for(int i=0;i<H.n_columns();i++)
		eigenvalues[i] = es.eigenvalues()[i];
}

void eigen_solver(const hermitian_matrix &A, matrix &vs, std::vector<double> &eigenvalues){
	eigenvalues.resize(A.n_columns());
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A.mat);
	vs = matrix(es.eigenvectors());
	for(int i=0;i<A.n_columns();i++)
		eigenvalues[i] = es.eigenvalues()[i];
}

void eigen_solver(const matrix &A, matrix &vs, std::vector<std::complex<double>> &eigenvalues){
	eigenvalues.resize(A.n_columns());
	Eigen::EigenSolver<Eigen::MatrixXd> es(A.mat);
	vs = matrix(A.n_columns(),A.n_columns());
	for(int i=0;i<A.n_columns();i++)
		for(int j=0;j<A.n_columns();j++)
			vs(i,j) = std::real(es.eigenvectors()(i,j));
	for(int i=0;i<A.n_columns();i++)
		eigenvalues[i] = es.eigenvalues()[i];
}

/* identical matrix I=diag(1,1,1,...)
 */
hermitian_matrix hermitian_matrix::idmat(int size){
	return hermitian_matrix(Eigen::MatrixXd::Identity(size,size));
}

/* S^(1/2) and S^(-1/2) */
std::tuple<hermitian_matrix,hermitian_matrix> hermitian_matrix::pmsqrt() const{
	hermitian_matrix sqrt(this->n_rows());
	sqrt.mat = this->mat.sqrt();
	hermitian_matrix nsqrt(this->n_rows());
	nsqrt.mat = sqrt.mat.inverse();
	return std::tuple<hermitian_matrix,hermitian_matrix>(sqrt,nsqrt);
}

// scalar times matrix
hermitian_matrix operator*(double l,const hermitian_matrix &r){
	hermitian_matrix ret = r;
	ret.mat = l*ret.mat;
	return ret;
}

}

#endif
