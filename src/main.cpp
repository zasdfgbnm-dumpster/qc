#define DEBUG
#include "common.hpp"
#include "aces.hpp"
#include "hartree_fock.hpp"
#include "itgl_transform.hpp"
#include "mbpt.hpp"
#include <iostream>
#include <cstdlib>

using namespace std;
using namespace quantum_chem;

int main(int argc, char *argv[]){
	calculation_data data;

	// read from file
	string filename = "samples/H2O.txt";
	data.nuclear_repulsion = 7.3256933549;
	try { read_configuration(filename,data); }
	catch(string s){ cerr << s << endl; return 1; }
		
	// Hartree Fock
	auto t = hartree_fock(data);
	cout << "---------- Hartree Fock Calculation ----------" << endl;
	cout << "Hartree Fock energy: " << data.Ehf() << endl;
	cout << "Time usage: " << t.count() << " seconds" << endl;

	// integral transform
	t = itgl_transform(data);
	cout << "---------- Integrals Transformation ----------" << endl;
	cout << "Time usage: " << t.count() << " seconds" << endl;

	// MBPT
	t = mbpt(data);
	cout << "------- Many Body Perturbation Theory --------" << endl;
	cout << "MBPT2 energy: " << data.Ehf()+data.mbpt[0] << endl;
	cout << "Time usage: " << t.count() << " seconds" << endl;

	return 0;
}
