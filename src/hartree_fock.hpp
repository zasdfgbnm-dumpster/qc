#ifndef __HPP_QUANTUM_HF__
#define __HPP_QUANTUM_HF__

#include <string>
#include "matrix.hpp"
#include "common.hpp"

#ifdef DEBUG
#include <iostream>
#endif

namespace quantum_chem {

/* rhf:
 * The routine that do Hartree Fock calculation, that is, to calculate:
 *  1) Hartree Fock eigenvalues, 
 *  2) Hartree Fock molecular orbitals,
 *  3) 0th, 1th order energy, Hartree Fock energy
 */
void hartree_fock(calculation_data &data, const double tol=1e-12){
	const int n_bases = data.n_baseset;
	const int dime = data.n_paired/2;
	const hermitian_matrix &h = data.ao_h;
	const hermitian_matrix &overlap = data.ao_overlap;
	const coulomb_exchange &ce = data.ao_2eint;
	matrix &C = data.ao2mo;
	std::vector<double> &eigenvalues = data.eigenvalues;

	// SCF loop
	hermitian_matrix P = hermitian_matrix::idmat(n_bases);
	hermitian_matrix Fuv(n_bases);
	C = matrix(n_bases,n_bases);
	eigenvalues.resize(n_bases);
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
		matrix P_old = P;
		#endif

		matrix C_old = data.ao2mo;
		// build Fock matrix
		Fuv = h;
		for(int u=0;u<n_bases;u++)
			for(int v=0;v<n_bases;v++)
				for(int a=0;a<n_bases;a++)
					for(int b=0;b<n_bases;b++)
						Fuv(u,v) += P(b,a)*(2*ce(u,a,v,b)-ce(u,a,b,v));
		// calculate new C, P and e
		// FC=SCe, C = [c1,c2,...], e = diag(e1,e2,...)
		// Fc1 = e1*Sc1, Fc2 = e2*Sc2, ......
		hermitian_matrix::generalized_eigen_solver(Fuv,overlap,C,eigenvalues);
		matrix Cocc = C.left_columns(dime);
		P = Cocc*(Cocc.conjugate_transpose());

		// test convergence
		matrix diffC = C-C_old;
		double errorC = diffC.norm();
		converged = ( errorC < tol );


		#ifdef DEBUG
		// calculate energy
		E = data.nuclear_repulsion;
		for(int i=0;i<dime;i++)
			E += eigenvalues[i];
		for(int i=0;i<n_bases;i++)
			for(int j=0;j<n_bases;j++)
				E += P(j,i)*h(i,j);

		matrix diffP = P-P_old;
		matrix diffF = Fuv-F_old;
		double errorP = diffP.norm();
		double errorF = diffF.norm(); 
		std::cout << "Energy: " << E << std::endl;
		std::cout << "Error in C: " << errorC << std::endl;
		std::cout << "Error in P: " << errorP << std::endl;
		std::cout << "Error in F: " << errorF << std::endl;
		std::cout << "Eigenvalues:";
		for(int i=0;i<n_bases;i++)
			std::cout << eigenvalues[i] << ' ';
		std::cout << std::endl;
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

	// calculate energy
	data.E0 = 0;
	for(int i=0;i<dime;i++)
		data.E0 += 2*data.eigenvalues[i];
	E = data.nuclear_repulsion + data.E0/2;
	data.E0 += data.nuclear_repulsion;
	for(int i=0;i<n_bases;i++)
		for(int j=0;j<n_bases;j++)
			E += P(j,i)*h(i,j);
	for(int i=0;i<dime;i++)
	data.E1 = E - data.E0;
	return;
}

}

#endif
