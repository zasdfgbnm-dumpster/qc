#ifndef __HPP_QUANTUM_TDHF__
#define __HPP_QUANTUM_TDHF__

#include "matrix.hpp"
#include "common.hpp"
#include <functional>
#include <cmath>
#include <iostream> // delete

namespace quantum_chem {
	
std::chrono::duration<double> tdhf(calculation_data &data){
	std::chrono::duration<double> elapsed_seconds(0);
	auto start = std::chrono::steady_clock::now();
	
	// build B matrix
	int sz = (data.n_paired/2)*(data.n_baseset-data.n_paired/2);
	hermitian_matrix J(sz);
	hermitian_matrix K(sz);
	std::function<int(int,int)> idx = [&](int i,int a){
		return (data.n_baseset-data.n_paired/2)*i + (a-data.n_paired/2);
	};
	// build J and K
	for(int i=0;i<data.n_paired/2;i++){
		for(int a=data.n_paired/2;a<data.n_baseset;a++){
			int ia = idx(i,a);
			for(int j=0;j<data.n_paired/2;j++){
				for(int b=data.n_paired/2;b<data.n_baseset;b++){
					int jb = idx(j,b);
					J(ia,jb) = data.mo_2eint(a,b,i,j);
					K(ia,jb) = data.mo_2eint(a,b,j,i);
				}
			}
		}
	}
	// build B_singlet
	data.B_singlet = 2*J - K;
	
	// build [[A B][-B -A]]
	matrix ABBA(sz*2,sz*2);
	ABBA.mat <<  data.A_singlet.mat ,  data.B_singlet.mat,
	            -data.B_singlet.mat , -data.A_singlet.mat;
	
	// solve for eigenvalues
	std::vector<std::complex<double>> tdhf_singlet_excited_complex(2*sz);
	data.tdhf_singlet_excited.resize(2*sz);
	eigen_solver(ABBA, data.tdhf_singlet_wavefunc, tdhf_singlet_excited_complex);
	transform(tdhf_singlet_excited_complex.begin(),tdhf_singlet_excited_complex.end(),
		  data.tdhf_singlet_excited.begin(),[](std::complex<double> c){ return std::real(c); });
					      
	
	auto end = std::chrono::steady_clock::now();
	elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
	return elapsed_seconds;
}

}

#endif