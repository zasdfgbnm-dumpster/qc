#ifndef __HPP_QUANTUM_MBPT__
#define __HPP_QUANTUM_MBPT__

#include <vector>
#include "common.hpp"

namespace quantum_chem {

std::chrono::duration<double> mbpt(calculation_data &data){
	std::chrono::duration<double> elapsed_seconds(0);
	auto start = std::chrono::steady_clock::now();

	data.mbpt.resize(1);
	std::vector<double> &e = data.eigenvalues;
	
	int n = data.n_paired/2;
	int m = data.n_baseset;
	double E = 0;

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			for(int a=n;a<m;a++){
				for(int b=n;b<m;b++){
					double abij = data.mo_2eint(a,b,i,j);
					double ijba = data.mo_2eint(i,j,b,a);
					double eijab = e[i] + e[j] - e[a] - e[b];
					E += (2*abij*abij - abij*ijba)/eijab;
				}
			}
		}
	}

	data.mbpt[0] = E;

	auto end = std::chrono::steady_clock::now();
	elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
	return elapsed_seconds;
}

}


#endif
