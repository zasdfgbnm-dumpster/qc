#ifndef __HPP_QUANTUM_COMMON__
#define __HPP_QUANTUM_COMMON__

#include <memory>
#include <vector>
#include "matrix.hpp"

namespace quantum_chem {

class dbl_e_itgls {
	/* 2 electron integrals <ij|kl> in 1212 notation
	 * no symmetry implemented yet*/
	int n;
	std::vector<double> itgl;
	int idx(int i1,int j1, int i2, int j2) const {
		if(i1>=n||j1>=n||i2>=n||j2>=n)
			throw "out of range in i1, i2, j1 or j2 in class dbl_e_itgls";
		int idx = i1*n*n*n+j1*n*n+i2*n+j2;
		if(idx>=itgl.size())
			throw "class dbl_e_itgls not correctly initialized"; 
		return idx;
	}
public:
	dbl_e_itgls(int n):n(n),itgl(n*n*n*n){}
	dbl_e_itgls(){}
	int n_bases() const { return n; }
	void resize(int new_n) { n=new_n; itgl.resize(n*n*n*n); }
	double &operator()(int i1,int j1, int i2, int j2){
		return itgl[idx(i1,j1,i2,j2)];
	}
	double operator()(int i1,int j1, int i2, int j2) const {
		return itgl[idx(i1,j1,i2,j2)];
	}
};

class calculation_data {
public:
	int n_baseset;                             //number of atomic base functions
	int n_paired;                              //number of paired electrons
	int n_alpha_electrons;                     //number of unpaired alpha electrons: not implemented yet
	int n_beta_electrons;                      //number of unpaired beta electrons: not implemented yet

	hermitian_matrix ao_h;                     //one electron integrals h in atomic orbitals
	hermitian_matrix ao_overlap;               //overlap matrix of atomic orbitals
	hermitian_matrix ao_dipole_x;              //one electron dipole moment at direction x in atomic orbitals
	hermitian_matrix ao_dipole_y;              //one electron dipole moment at direction y in atomic orbitals
	hermitian_matrix ao_dipole_z;              //one electron dipole moment at direction z in atomic orbitals
	dbl_e_itgls ao_2eint;                      //2 electron integrals <ij|kl> in atomic orbitals
	matrix ao2mo;                              //transmition matrix between atomic orbitals and molecular orbitals
	hermitian_matrix mo_h;                     //one electron integrals h in molecular orbitals
	hermitian_matrix mo_dipole_x;              //one electric dipole moment at direction x in molecular orbitals
	hermitian_matrix mo_dipole_y;              //one electric dipole moment at direction y in molecular orbitals
	hermitian_matrix mo_dipole_z;              //one electric dipole moment at direction z in molecular orbitals
	dbl_e_itgls mo_2eint;                      //2 electron integrals <ij|kl> in molecular orbitals
	
	double nuclear_repulsion;                  //nuclear repulsion
	std::vector<double> eigenvalues;           //Hartree Fock eigenvalues
	double E0;                                 //0th order energy
	double E1;                                 //1st order energy
	double Ehf() {                             //return Hartree Fock energy
		 return E0 + E1;
	}
	std::vector<double> mbpt;                  //MBPT energies, starting from 2nd order
	
	hermitian_matrix A_singlet;                //CIS A matrix for singlet spin states
	hermitian_matrix A_triplet;                //CIS A matrix for triplet spin states
	matrix ci_singlet_wavefunc;                //CIS coefficient for spacial orbitals for singlet spin wavefunction
	matrix ci_triplet_wavefunc;                //CIS coefficient for spacial orbitals for triplet spin wavefunction
	std::vector<double> ci_singlet_excited;    //energy of CIS singlet excited states
	std::vector<double> ci_triplet_excited;    //energy of CIS triplet excited states
	
	hermitian_matrix B_singlet;                //TDHF B matrix for singlet spin states
	std::vector<double> tdhf_singlet_excited;  //energy of TDHF singlet excited states
	matrix tdhf_singlet_wavefunc;              //TDHF coefficient for spacial orbitals for singlet spin wavefunction
};


}

#endif
