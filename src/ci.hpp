#ifndef __HPP_QUANTUM_CI__
#define __HPP_QUANTUM_CI__

#include "matrix.hpp"
#include "common.hpp"
#include <functional>
#include <cmath>
#include <iostream> // delete

namespace quantum_chem {
	
std::chrono::duration<double> cis(calculation_data &data){
	std::chrono::duration<double> elapsed_seconds(0);
	auto start = std::chrono::steady_clock::now();
	
	// build A matrix
	int sz = (data.n_paired/2)*(data.n_baseset-data.n_paired/2);
	data.A = hermitian_matrix(sz);
	std::function<int(int,int)> idx = [&](int i,int a){
		return (data.n_baseset-data.n_paired/2)*i + (a-data.n_paired/2);
	};
	for(int i=0;i<data.n_paired/2;i++){
		for(int a=data.n_paired/2;a<data.n_baseset;a++){
			int ia = idx(i,a);
			for(int j=0;j<data.n_paired/2;j++){
				for(int b=data.n_paired/2;b<data.n_baseset;b++){
					int jb = idx(j,b);
					double term1 = (i==j&&a==b?(data.eigenvalues[a]-data.eigenvalues[i]):0);
					double term2 = 2*data.mo_2eint(a,i,j,b)-data.mo_2eint(a,b,j,i);
					data.A(ia,jb) = term1 + term2;
				}
			}
		}
	}
	
	// solve for eigenvalues
	eigen_solver(data.A, data.ci_wavefunc, data.ci_excited);
	
	auto end = std::chrono::steady_clock::now();
	elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
	return elapsed_seconds;
}
	
}

#endif