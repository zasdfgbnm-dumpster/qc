#ifndef __HPP_QUANTUM_RHF__
#define __HPP_QUANTUM_RHF__

#include <string>
#include "matrix.hpp"
#include "hartree_fock.hpp"
#include "utils.hpp"

namespace quantum_chem {

/* rhf:
 * the routine that calculates Restricted Hartree Fock
 */
hf_output rhf(int n_electrons, hermitian_matrix h, hermitian_matrix overlap, coulomb_exchange ce, double tol){

	int n_bases = h.n_rows();
	int dime = n_electrons/2;

	// check data consistency
	std::string errinfo;
	errinfo = "bad matrix size in parameters of rhf()";
	assert(n_bases==overlap.n_rows(),errinfo);
	assert(n_bases==ce.n_bases(),errinfo);
	errinfo = "Restricted HF only apply to doubly occupied orbitals";
	assert(n_electrons%2==0,errinfo);
	errinfo = "Number of base set must be larger than number of space orbitals";
	assert(dime<=n_bases,errinfo);

	// SCF loop
	hermitian_matrix P = idmat(dime);
	matrix C(n_bases);
	array<double> eigenvalues(n_bases);
	bool converged = false;
	do{
		// build Fock matrix
		hermitian_matrix Fuv = h;
		for(int u=0;u<n_bases;u++)
			for(int v=0;v<n_bases;v++)
				for(int a=0;a<dime;a++)
					for(int b=0;b<dime;b++)
						Fuv(u,v) += P(a,b)*(2*ce(u,v,a,b)-ce(u,b,a,v));
		// calculate new C, P and e
		// FC=SCe, C = [c1,c2,...], e = diag(e1,e2,...)
		// Fc1 = e1*Sc1, Fc2 = e2*Sc2, ......
		matrix C_old = C;
		generalized_hermitian_eigen_solver(F,overlap,C,eigenvalues);
		matrix Cocc = C.left_columns(dime);
		P = Cocc.conjugate_transpose()*Cocc;

		// test convergence
		converged = ( (C-C_old).norm() < tol );
	}while(!converged);
	return make_tuple(eigenvalues,C);
}

#endif
