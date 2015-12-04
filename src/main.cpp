//#define DEBUG
#include "common.hpp"
#include "aces.hpp"
#include "hartree_fock.hpp"
#include "itgl_transform.hpp"
#include "mbpt.hpp"
#include "ci.hpp"
#include "tdhf.hpp"
#include <iostream>
#include <cstdlib>
#include <algorithm>

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
	cout << "Time usage: " << t.count() << " seconds" << endl << endl;

	// integral transform
	t = itgl_transform(data);
	cout << "---------- Integrals Transformation ----------" << endl;
	cout << "Time usage: " << t.count() << " seconds" << endl << endl;

	// MBPT
	t = mbpt(data);
	cout << "------- Many Body Perturbation Theory --------" << endl;
	cout << "MBPT2 energy: " << data.Ehf()+data.mbpt[0] << endl;
	cout << "Time usage: " << t.count() << " seconds" << endl << endl;
	
	// CIS
	t = cis(data);
	cout << "----------------- CI Single ------------------" << endl;
	cout << "CIS singlet excited states: " << endl;
	for(double i:data.ci_singlet_excited)
		cout << i << " ";
	cout << endl;
	cout << "CIS triplet excited states: " << endl;
	for(double i:data.ci_triplet_excited)
		cout << i << " ";
	cout << endl;
	cout << "Time usage: " << t.count() << " seconds" << endl << endl;

	// CIS
	t = tdhf(data);
	cout << "------------------- TDHF ---------------------" << endl;
	cout << "TDHF singlet excited states: " << endl;
	auto sorted = data.tdhf_singlet_excited;
	sort(sorted.begin(),sorted.end());
	for_each(sorted.begin()+sorted.size()/2,sorted.end(),[](double i){ cout << i << " "; });
	cout << endl;
	cout << "Time usage: " << t.count() << " seconds" << endl << endl;
	
	return 0;
}
