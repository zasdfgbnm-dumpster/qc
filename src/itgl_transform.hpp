#ifndef __HPP_QUANTUM_ITGLTRANS__
#define __HPP_QUANTUM_ITGLTRANS__

#include "matrix.hpp"
#include "common.hpp"

namespace quantum_chem {

std::chrono::duration<double> itgl_transform(calculation_data &data) {
	std::chrono::duration<double> elapsed_seconds(0);
	auto start = std::chrono::steady_clock::now();

	int n = data.n_baseset;
	hermitian_matrix &mm = data.mo_h;
	const hermitian_matrix &aa = data.ao_h;
	dbl_e_itgls &mmmm = data.mo_2eint;
	const dbl_e_itgls &aaaa = data.ao_2eint;
	matrix &ao2mo = data.ao2mo;

	mm = hermitian_matrix(n);
	mmmm = dbl_e_itgls(n);

	// transform h
	hermitian_matrix am(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			am(i,j) = 0;
			for(int k=0;k<n;k++)
				am(i,j) += ao2mo(k,j)*aa(i,k);
		}
	}
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			mm(i,j) = 0;
			for(int k=0;k<n;k++)
				mm(i,j) += ao2mo(k,i)*am(k,j);
			mm(j,i) = mm(i,j);
		}
	}
	
	// transform <ij|kl>
	dbl_e_itgls aaam(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				for(int l=0;l<n;l++){
					aaam(i,j,k,l) = 0;
					for(int p=0;p<n;p++)
						aaam(i,j,k,l) += ao2mo(p,l)*aaaa(i,j,k,p);
				}
			}
		}
	}
	dbl_e_itgls aamm(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				for(int l=0;l<n;l++){
					aamm(i,j,k,l) = 0;
					for(int p=0;p<n;p++)
						aamm(i,j,k,l) += ao2mo(p,k)*aaam(i,j,p,l);
				}
			}
		}
	}
	dbl_e_itgls ammm(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				for(int l=0;l<n;l++){
					ammm(i,j,k,l) = 0;
					for(int p=0;p<n;p++)
						ammm(i,j,k,l) += ao2mo(p,j)*aamm(i,p,k,l);
				}
			}
		}
	}
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				for(int l=0;l<n;l++){
					mmmm(i,j,k,l) = 0;
					for(int p=0;p<n;p++)
						mmmm(i,j,k,l) += ao2mo(p,i)*ammm(p,j,k,l);
				}
			}
		}
	}
	auto end = std::chrono::steady_clock::now();
	elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
	return elapsed_seconds;
}

}

#endif
