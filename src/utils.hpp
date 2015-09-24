#ifndef __HPP_QUANTUM_UTILS__
#define __HPP_QUANTUM_UTILS__

#include <string>
#include <iostream>

namespace quantum_chem {

void assert(bool supposed_to_be_true, std::string errinfo){
	if(supposed_to_be_true)
		return;
	std::cerr << errinfo << std::endl;
	throw errinfo;
}

}

#endif
