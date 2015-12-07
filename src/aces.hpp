#ifndef __HPP_QUANTUM_ACES__
#define __HPP_QUANTUM_ACES__

#include <string>
#include <fstream>
#include "matrix.hpp"
#include "hartree_fock.hpp"
#include "utils.hpp"

#ifdef DEBUG
#include <iostream>
#endif

namespace quantum_chem {

void read_configuration(std::string filename, calculation_data &data){
	hermitian_matrix &h = data.ao_h;
	hermitian_matrix &overlap = data.ao_overlap;
	dbl_e_itgls &ce = data.ao_2eint;
	std::ifstream in(filename);
	// number of base functions
	int &n_bases = data.n_baseset;
	int n_alpha,n_beta;
	in >> n_bases;
	// number of electrons
	in >> n_alpha;
	in >> n_beta;
	data.n_paired = n_alpha+n_beta;
	#ifdef DEBUG
	std::cout << "n_bases: " << n_bases << std::endl;
	std::cout << "n_alpha: " << n_alpha << std::endl;
	std::cout << "n_beta: " << n_beta << std::endl;
	#endif

	// check next line
	std::string line,errinfo;
	std::getline(in, line);
	std::getline(in, line);
	errinfo = "bad format of input file";
	qc_assert(line == std::string(" one electron h"),errinfo);
	// read one h
	#ifdef DEBUG
	std::cout << "Start reading integrals h..." << std::endl;
	#endif
	h = hermitian_matrix(n_bases);
	for(int i=0;i<n_bases;i++)
		for(int j=0;j<n_bases;j++)
			in >> h(i,j);

	// check next line
	std::getline(in, line);
	std::getline(in, line);
	errinfo = "bad format of input file";
	qc_assert(line == std::string(" overlap integrals"),errinfo);
	// read ovrlp
	#ifdef DEBUG
	std::cout << "Start reading overlap integrals..." << std::endl;
	#endif
	overlap = hermitian_matrix(n_bases);
	for(int i=0;i<n_bases;i++)
		for(int j=0;j<n_bases;j++)
			in >> overlap(i,j);

	// check next line
	std::getline(in, line);
	std::getline(in, line);
	errinfo = "bad format of input file";
	qc_assert(line == std::string(" 2-electron integrals in 1122 notation"),errinfo);
	// read 2e int
	#ifdef DEBUG
	std::cout << "Start reading 2 electron integrals..." << std::endl;
	#endif
	ce = dbl_e_itgls(n_bases);
	for(int i=0;i<n_bases;i++)
		for(int j=0;j<n_bases;j++)
			for(int k=0;k<n_bases;k++)
				for(int l=0;l<n_bases;l++)
					in >> ce(i,k,j,l);

	// check next line
	std::getline(in, line);
	std::getline(in, line);
	errinfo = "bad format of input file";
	qc_assert(line == std::string(" dipole x"),errinfo);
	// read ao_dipole_x
	#ifdef DEBUG
	std::cout << "Start reading dipole x integrals..." << std::endl;
	#endif
	data.ao_dipole_x = hermitian_matrix(n_bases);
	for(int i=0;i<n_bases;i++)
		for(int j=0;j<n_bases;j++)
			in >> data.ao_dipole_x(i,j);
	
	// check next line
	std::getline(in, line);
	std::getline(in, line);
	errinfo = "bad format of input file";
	qc_assert(line == std::string(" dipole y"),errinfo);
	// read ao_dipole_y
	#ifdef DEBUG
	std::cout << "Start reading dipole y integrals..." << std::endl;
	#endif
	data.ao_dipole_y = hermitian_matrix(n_bases);
	for(int i=0;i<n_bases;i++)
		for(int j=0;j<n_bases;j++)
			in >> data.ao_dipole_y(i,j);
		
	// check next line
	std::getline(in, line);
	std::getline(in, line);
	errinfo = "bad format of input file";
	qc_assert(line == std::string(" dipole z"),errinfo);
	// read ao_dipole_z
	#ifdef DEBUG
	std::cout << "Start reading dipole z integrals..." << std::endl;
	#endif
	data.ao_dipole_z = hermitian_matrix(n_bases);
	for(int i=0;i<n_bases;i++)
		for(int j=0;j<n_bases;j++)
			in >> data.ao_dipole_z(i,j);
		
}

}

#endif
