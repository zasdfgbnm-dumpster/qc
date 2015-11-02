#ifndef __HPP_QUANTUM_RHF__
#define __HPP_QUANTUM_RHF__

#include <string>
#include "matrix.hpp"
#include "hartree_fock.hpp"
#include "utils.hpp"

#ifdef DEBUG
#include <iostream>
#endif

namespace quantum_chem {

/* rhf:
 * the routine that calculates Restricted Hartree Fock
 */
hf_output rhf(int n_electrons, hermitian_matrix h, hermitian_matrix overlap, coulomb_exchange ce, double tol=1e-12){

	int n_bases = h.n_rows();
	int dime = n_electrons/2;

	// check data consistency
	std::string errinfo;
	errinfo = "bad matrix size in parameters of rhf()";
	qc_assert(n_bases==overlap.n_rows(),errinfo);
	qc_assert(n_bases==ce.n_bases(),errinfo);
	errinfo = "Restricted HF only apply to doubly occupied orbitals";
	qc_assert(n_electrons%2==0,errinfo);
	errinfo = "Number of base set must be larger than number of space orbitals";
	qc_assert(dime<=n_bases,errinfo);

	// SCF loop
	hermitian_matrix P = hermitian_matrix::idmat(n_bases);
	hermitian_matrix Fuv(n_bases);
	matrix C(n_bases,n_bases);
	std::vector<double> eigenvalues(n_bases);
	bool converged = false;
	double E = 0;
	#ifdef DEBUG
	int count=0;
	#endif
	do{
		#ifdef DEBUG
		count++;
		std::cout << "Loop number: " << count << std::endl;
		matrix F_old = Fuv;
		matrix C_old = C;
		#endif
		matrix P_old = P;
		// build Fock matrix
		Fuv = h;
		for(int u=0;u<n_bases;u++)
			for(int v=0;v<n_bases;v++)
				for(int a=0;a<n_bases;a++)
					for(int b=0;b<n_bases;b++)
						Fuv(u,v) += P(a,b)*(2*ce(u,b,v,a)-ce(u,b,a,v));
		// calculate new C, P and e
		// FC=SCe, C = [c1,c2,...], e = diag(e1,e2,...)
		// Fc1 = e1*Sc1, Fc2 = e2*Sc2, ......
		hermitian_matrix::generalized_eigen_solver(Fuv,overlap,C,eigenvalues);
		matrix Cocc = C.left_columns(dime);
		P = Cocc*Cocc.conjugate_transpose();

		// calculate energy
		E = 0;
		for(int i=0;i<dime;i++)
			E += eigenvalues[i];
		for(int i=0;i<n_bases;i++)
			for(int j=0;j<n_bases;j++)
				E += P(j,i)*h(i,j);
		// test convergence
		matrix diffP = P-P_old;
		double errorP = diffP.norm();
		converged = ( errorP < tol );
		#ifdef DEBUG
		matrix diffC = C-C_old;
		matrix diffF = Fuv-F_old;
		double errorC = diffC.norm();
		double errorF = diffF.norm(); 
		std::cout << "Energy: " << E << std::endl;
		std::cout << "Error in C: " << errorC << std::endl;
		std::cout << "Error in P: " << errorP << std::endl;
		std::cout << "Error in F: " << errorF << std::endl;
		std::cout << "Eigenvalues:";
		for(int i=0;i<n_bases;i++)
			std::cout << eigenvalues[i] << ' ';
		std::cout << std::endl;
		for(int i=0;i<n_bases;i++)
			for(int j=0;j<n_bases;j++){
				if(abs(diffC(i,j))<tol)
					diffC(i,j) = 0;
				if(abs(diffP(i,j))<tol)
					diffP(i,j) = 0;
				if(abs(diffF(i,j))<tol)
					diffF(i,j) = 0;
			}
		std::cout << "C: " << std::endl << C.mat << std::endl;
		std::cout << "delta C: " << std::endl << diffC.mat << std::endl;
		std::cout << "P: " << std::endl << P.mat << std::endl;
		std::cout << "delta P: " << std::endl << diffP.mat << std::endl;
		std::cout << "F: " << std::endl << Fuv.mat << std::endl;
		std::cout << "delta F: " << std::endl << diffF.mat << std::endl;
		std::cout << "-----------------" << std::endl;
		#endif
	}while(!converged);
	#ifdef DEBUG
	std::cout << "Converged." << std::endl;
	std::cout << "Energy: " << E << std::endl;
	std::cout << "Eigenvalues:" << std::endl;
	for(int i=0;i<n_bases;i++)
		std::cout << eigenvalues[i] << ' ';
	std::cout << std::endl;
	std::cout << "==========================================" << std::endl;
	#endif
	// return
	return make_tuple(E,eigenvalues,C);
}

}

#endif
