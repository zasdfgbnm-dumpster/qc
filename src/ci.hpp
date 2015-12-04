#ifndef __HPP_QUANTUM_CI__
#define __HPP_QUANTUM_CI__

#include "matrix.hpp"
#include "common.hpp"
#include <functional>
#include <cmath>

namespace quantum_chem {
	
std::chrono::duration<double> cis(calculation_data &data){
	std::chrono::duration<double> elapsed_seconds(0);
	auto start = std::chrono::steady_clock::now();
	
	// build A matrix
	int sz = (data.n_paired/2)*(data.n_baseset-data.n_paired/2);
	hermitian_matrix epsilon(sz);
	hermitian_matrix J(sz);
	hermitian_matrix K(sz);
	std::function<int(int,int)> idx = [&](int i,int a){
		return (data.n_baseset-data.n_paired/2)*i + (a-data.n_paired/2);
	};
	// build epsilon, J and K
	for(int i=0;i<data.n_paired/2;i++){
		for(int a=data.n_paired/2;a<data.n_baseset;a++){
			int ia = idx(i,a);
			for(int j=0;j<data.n_paired/2;j++){
				for(int b=data.n_paired/2;b<data.n_baseset;b++){
					int jb = idx(j,b);
					epsilon(ia,jb) =  ( i==j&&a==b ? (data.eigenvalues[a]-data.eigenvalues[i]) : 0 );
					J(ia,jb) = data.mo_2eint(a,j,i,b);
					K(ia,jb) = data.mo_2eint(a,j,b,i);
				}
			}
		}
	}
	// build A_singlet and A_triplet
	data.A_singlet = epsilon + 2*J - K;
	data.A_triplet = epsilon - K;
	
	// solve for eigenvalues
	eigen_solver(data.A_singlet, data.ci_singlet_wavefunc, data.ci_singlet_excited);
	eigen_solver(data.A_triplet, data.ci_triplet_wavefunc, data.ci_triplet_excited);
	
	auto end = std::chrono::steady_clock::now();
	elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
	return elapsed_seconds;
}

}

#endif