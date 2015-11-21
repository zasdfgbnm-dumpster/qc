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
	int sz = data.n_paired*(data.n_baseset*2-data.n_paired);
	data.A = hermitian_matrix(sz);
	std::function<int(int,int)> idx = [&](int i,int a){
		return (data.n_baseset*2-data.n_paired)*i + (a-data.n_paired);
	};
	for(int i=0;i<data.n_paired;i++){
		for(int a=data.n_paired;a<data.n_baseset*2;a++){
			int ia = idx(i,a);
			for(int j=0;j<data.n_paired;j++){
				for(int b=data.n_paired;b<data.n_baseset*2;b++){
					int jb = idx(j,b);
					double term1 =  ( i==j&&a==b ? (data.eigenvalues[a/2]-data.eigenvalues[i/2]) : 0 );
					double term2 =  ( a%2==i%2 && b%2==j%2 ? data.mo_2eint(a/2,j/2,i/2,b/2)      : 0 );
					double term3 = -( a%2==b%2 && i%2==j%2 ? data.mo_2eint(a/2,j/2,b/2,i/2)      : 0 );
					data.A(ia,jb) = term1 + term2 + term3;
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