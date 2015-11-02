//#define DEBUG
#include "aces.hpp"
#include "rhf.hpp"
#include <iostream>
#include <cstdlib>

using namespace std;
using namespace quantum_chem;

int main(int argc, char *argv[]){
	try {
		if(argc!=2){
			cout << "Usage: rhf <filename>" << endl;
			exit(1);
		}

		string filename = argv[1];
		int n_alpha,n_beta;
		hermitian_matrix h,overlap;
		coulomb_exchange ce;
		
		#ifdef DEBUG
		cout << "==========================================" << endl;
		cout << "Reading configuration file" << endl;
		#endif
		read_configuration(filename,n_alpha,n_beta,h,overlap,ce);
		#ifdef DEBUG
		cout << "Reading configuration file done" << endl;
		cout << "==========================================" << endl;
		cout << "Starting SCF loop" << endl;
		#endif

		hf_output output = rhf(n_alpha+n_beta,h,overlap,ce);
	
		cout << get<0>(output) << endl;
		return 0;
	}catch(string s){
		cerr << s << endl;
		return 1;
	}
}
