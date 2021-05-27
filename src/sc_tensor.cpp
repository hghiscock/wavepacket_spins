#include <math.h>
#include <complex>
#include "param.h"
#include "procedures.h"
#include "hamiltonians.h"
#include "sc_tensor.h"
#include <iostream>

using namespace std;
using namespace N;
using namespace N1;
using namespace N2;

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

spin_corr::spin_corr(param_dat d, hamiltonian h) {

	int i, j, k;

	complex<double> htmp[h.m*h.m*d.ntheta] = {};
	complex<double> ptmp[h.m*h.m*d.ntheta] = {};

	p = new complex<double>[h.m*h.m*d.ntheta]();
	q = new complex<double>[h.m*h.m*d.ntheta]();
	rtmp = new complex<double>[9*d.ntheta]();
	rab = new double[9*d.ntheta]();

	xj = complex<double>(0.0, 1.0);
	m = h.m;
	mprime = m / 2;

	for (i=0; i<d.ntheta; i++) {
		for (j=0; j<m; j++) {
			for (k=0; k<m; k++) {
				htmp[j*m+k] = h.hze[i*m*m+j*m+k] + h.h_td0[j*m+k];
			}
			q[i*m*m+j*m+j] = complex<double>(1.0, 0.0);
		}
		matmul(htmp, q, ptmp, m, 0, i*m*m, 0);
		for (j=0; j<m; j++) {
			for (k=0; k<m; k++) {
				p[i*m*m+j*m+k] = -1.0*xj*ptmp[j*m+k];
			}
		}
	}

};

//---------------------------------------------------------------------------//

void spin_corr::evolve(param_dat d, hamiltonian h){

	int i, j, k;
	int ii, itmp, jtmp;

	complex<double> htmp[m*m] = {};
	complex<double> wtmp[m*m] = {};
	complex<double> ptmp[m*m] = {};
	complex<double> w[m*m] = {};

	for (i=0; i<d.ntheta; i++) {
		for (j=0; j<m; j++) {
			for (k=0; k<m; k++) {
				itmp = i*m*m+j*m+k;
				q[itmp] = q[itmp] + d.a[0]*p[itmp];
			}
		}
	}

	for (ii=1; ii<4; ii++) {
		h.time_step(d);
		for (i=0; i<d.ntheta; i++) {
			for (j=0; j<m; j++) {
				for (k=0; k<m; k++) {
					jtmp = j*m+k;
					itmp = jtmp + i*m*m;
					htmp[jtmp] = h.hze[itmp] + h.h_td[jtmp];
				}
			}
			matmul(htmp, htmp, wtmp, m);
			add_mat(wtmp, h.dhdt, w, xj, m);
			matmul(w, q, ptmp, m, 0, i*m*m, 0);
			for (j=0; j<m; j++) {
				for (k=0; k<m; k++) {
					jtmp = j*m+k;
					itmp = jtmp + i*m*m;
					p[itmp] = p[itmp] - d.b[ii]*ptmp[jtmp];
					q[itmp] = q[itmp] + d.a[ii]*p[itmp];
				}
			}
		}
	}

};

//---------------------------------------------------------------------------//

void spin_corr::evolve_t_ind(param_dat d, hamiltonian h){

	int i, j, k;
	int ii, itmp, jtmp;

	complex<double> ptmp[m*m] = {};
	complex<double> w[m*m] = {};

	for (i=0; i<d.ntheta; i++) {
		for (j=0; j<m; j++) {
			for (k=0; k<m; k++) {
				itmp = i*m*m+j*m+k;
				q[itmp] = q[itmp] + d.a[0]*p[itmp];
			}
		}
	}

	for (ii=1; ii<4; ii++) {
		h.time_step(d);
		for (i=0; i<d.ntheta; i++) {
			matmul(h.hze, h.hze, w, m, i*m*m, i*m*m, 0); 
			matmul(w, q, ptmp, m, 0, i*m*m, 0);
			for (j=0; j<m; j++) {
				for (k=0; k<m; k++) {
					jtmp = j*m+k;
					itmp = jtmp + i*m*m;
					p[itmp] = p[itmp] - d.b[ii]*ptmp[jtmp];
					q[itmp] = q[itmp] + d.a[ii]*p[itmp];
				}
			}
		}
	}

};

//---------------------------------------------------------------------------//

void spin_corr::calc_rab(param_dat d) {

	int i, j, k;
	int jprime, kprime;
	int itmp, itmp1;

	complex<double> rab_tmp[9*d.ntheta] = {};

	for (j=0; j<mprime; j++) {
		jprime = j + mprime;
		for (k=0; k<mprime; k++) {
			kprime = k + mprime;
			for (i=0; i<d.ntheta; i++) {
				itmp1 = i*9;
				itmp = i*m*m;

				rab_tmp[itmp1] = rab_tmp[itmp1] +
					conj(q[itmp+j*m+kprime])*q[itmp+jprime*m+k];
				rab_tmp[itmp1] = rab_tmp[itmp1] +
					conj(q[itmp+j*m+k])*q[itmp+jprime*m+kprime];

				rab_tmp[itmp1+4] = rab_tmp[itmp1+4] -
					conj(q[itmp+jprime*m+k])*q[itmp+j*m+kprime];
				rab_tmp[itmp1+4] = rab_tmp[itmp1+4] +
					conj(q[itmp+jprime*m+kprime])*q[itmp+j*m+k];

				rab_tmp[itmp1+8] = rab_tmp[itmp1+8] +
					conj(q[itmp+j*m+k])*q[itmp+j*m+k];
				rab_tmp[itmp1+8] = rab_tmp[itmp1+8] +
					conj(q[itmp+jprime*m+kprime])*q[itmp+jprime*m+kprime];
				rab_tmp[itmp1+8] = rab_tmp[itmp1+8] -
					conj(q[itmp+jprime*m+k])*q[itmp+jprime*m+k];
				rab_tmp[itmp1+8] = rab_tmp[itmp1+8] -
					conj(q[itmp+j*m+kprime])*q[itmp+j*m+kprime];

				rab_tmp[itmp1+2] = rab_tmp[itmp1+2] +
					conj(q[itmp+j*m+kprime])*q[itmp+j*m+k];
				rab_tmp[itmp1+2] = rab_tmp[itmp1+2] -
					conj(q[itmp+jprime*m+kprime])*q[itmp+jprime*m+k];

				rab_tmp[itmp1+6] = rab_tmp[itmp1+6] +
					conj(q[itmp+j*m+k])*q[itmp+jprime*m+k];
				rab_tmp[itmp1+6] = rab_tmp[itmp1+6] -
					conj(q[itmp+j*m+kprime])*q[itmp+jprime*m+kprime];

			}
		}
	}

	for (i=0; i<d.ntheta*9; i++) {
		rtmp[i] = rab_tmp[i];
	}

};

//---------------------------------------------------------------------------//

void spin_corr::ex_rab(param_dat d) {

	int i, itmp;

	for (i=0; i<d.ntheta; i++) {
		itmp = i*9;
		rab[itmp] = rtmp[itmp].real() / double(m);
		rab[itmp+4] = rtmp[itmp+4].real() / double(m);
		rab[itmp+8] = rtmp[itmp+8].real() * 0.5 / double(m);
		rab[itmp+1] = rtmp[itmp].imag() / double(m);
		rab[itmp+3] = rtmp[itmp+4].imag() / double(m);
		rab[itmp+2] = rtmp[itmp+2].real() / double(m);
		rab[itmp+6] = rtmp[itmp+6].real() / double(m);
		rab[itmp+5] = rtmp[itmp+2].imag() / double(m);
		rab[itmp+7] = rtmp[itmp+6].imag() / double(m);
	}

};

//---------------------------------------------------------------------------//

void spin_corr::calc_sc(param_dat d, hamiltonian h) {

	if (d.t_ind == 1) {
		evolve(d, h);
	}
	else {
		evolve_t_ind(d, h);
	}

	calc_rab(d);
	ex_rab(d);

};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

