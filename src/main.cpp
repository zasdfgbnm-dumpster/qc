#define DEBUG
#include "aces.hpp"
#include "hartree_fock.hpp"
#include "common.hpp"
#include <iostream>
#include <cstdlib>

using namespace std;
using namespace quantum_chem;

int main(int argc, char *argv[]){
	string filename = "samples/H2O.txt";
	calculation_data data;

	// read from file
	try { read_configuration(filename,data); }
	catch(string s){ cerr << s << endl; return 1; }
		
	data.nuclear_repulsion = 7.3256933549;
	hartree_fock(data);

	return 0;
}
