#ifndef __HPP_QUANTUM_HF__
#define __HPP_QUANTUM_HF__

#include <memory>
#include "matrix.hpp"

namespace quantum_chem {

class coulomb_exchange {
	/* 1212 notation */
	int n;
	std::vector<double> itgl;
public:
	coulomb_exchange(int n):n(n),itgl(n*n*n*n){}
	coulomb_exchange(){}
	int n_bases(){ return n; }
	void resize(int new_n) { n=new_n; itgl.resize(n*n*n*n); }
	double &operator()(int i1,int j1, int i2, int j2){
		return itgl[i1*n*n*n+j1*n*n+i2*n+j2];
	}
	double operator()(int i1,int j1, int i2, int j2) const {
		return itgl[i1*n*n*n+j1*n*n+i2*n+j2];
	}
};

using hf_output = std::tuple<double,std::vector<double>,matrix>;


}

#endif
