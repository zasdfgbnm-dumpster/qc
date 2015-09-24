#ifndef __HPP_QUANTUM_HF__
#define __HPP_QUANTUM_HF__

#include <array>
#include <tuple>
#include "matrix.hpp"

namespace quantum_chem {

class coulomb_exchange {
	const int n;
	array<double> itgl;
public:
	coulomb_exchange(int n):
		n(n),itgl(array<double>( n*(n+1)*(n*n+n+2)/8 ){}
	int n_bases(){ return n; }
	double operator()(int i1,int j1, int i2, int j2){
		return 0;
	}
};

using hf_output = tuple<array<double>,matrix>;

}

#endif
