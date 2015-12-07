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
	
	int sz = (data.n_paired/2)*(data.n_baseset-data.n_paired/2);
	std::function<int(int,int)> idx = [&](int i,int a){
		return (data.n_baseset-data.n_paired/2)*i + (a-data.n_paired/2);
	};
	
	// build A + B matrix for singlet
	hermitian_matrix AB = data.A_singlet + data.B_singlet;
	hermitian_matrix inv_AB = AB.inverse();
	
	
	// build theta matrix
	matrix theta_x(sz,1);
	matrix theta_y(sz,1);
	matrix theta_z(sz,1);
	for(int i=0;i<data.n_paired/2;i++){
		for(int a=data.n_paired/2;a<data.n_baseset;a++){
			int ia = idx(i,a);
			theta_x(ia,0) = data.mo_dipole_x(i,a);
			theta_y(ia,0) = data.mo_dipole_y(i,a);
			theta_z(ia,0) = data.mo_dipole_z(i,a);
		}
	}
	
	// solve (A+B)C=Theta
	matrix Cx = inv_AB * theta_x;
	matrix Cy = inv_AB * theta_y;
	matrix Cz = inv_AB * theta_z;
	
	// calculate the polarizability
	matrix *Cs[3] = { &Cx,&Cy,&Cz };
	matrix *thetas[3] = { &theta_x,&theta_y,&theta_z };
	for(int xyz1=0;xyz1<3;xyz1++)
		for(int xyz2=0;xyz2<3;xyz2++)
			data.polarizability[xyz1][xyz2] = 0;
	for(int i=0;i<data.n_paired/2;i++){
		for(int a=data.n_paired/2;a<data.n_baseset;a++){
			int ia = idx(i,a);
			for(int xyz1=0;xyz1<3;xyz1++)
				for(int xyz2=0;xyz2<3;xyz2++)
					data.polarizability[xyz1][xyz2] += 4 * (*thetas[xyz1])(ia,0) * (*Cs[xyz2])(ia,0);
		}
	}
	
	auto end = std::chrono::steady_clock::now();
	elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
	return elapsed_seconds;
}

}

#endif