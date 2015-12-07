#ifndef __HPP_QUANTUM_CPHF__
#define __HPP_QUANTUM_CPHF__

#include "matrix.hpp"
#include "common.hpp"
#include <functional>
#include <cmath>

namespace quantum_chem {
	
std::chrono::duration<double> cphf(calculation_data &data){
	std::chrono::duration<double> elapsed_seconds(0);
	auto start = std::chrono::steady_clock::now();
	
	// build Theta matrix
	
	
	auto end = std::chrono::steady_clock::now();
	elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
	return elapsed_seconds;
}

}

#endif