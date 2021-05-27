#include <math.h>
#include <complex>
#include "param.h"
#include "procedures.h"
#include "hamiltonians.h"
#include <iostream>

using namespace std;
using namespace N;
using namespace N1;

double const pi = 3.14159265359;

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

hamiltonian::hamiltonian(radical_dat rad) {

	double rttwo = 1.0 / sqrt(2.0);

	s2 = new complex<double>[16];	
	s3 = new complex<double>[36];	

	m = rad.m;
	nspins = rad.nspins;
	inow = 0;

	s2[0] = complex<double>(0.0, 0.0);
	s2[1] = complex<double>(0.5, 0.0);
	s2[2] = complex<double>(0.5, 0.0);
	s2[3] = complex<double>(0.0, 0.0);

	s2[4] = complex<double>(0.0, 0.0);
	s2[5] = complex<double>(0.0, -0.5);
	s2[6] = complex<double>(0.0, 0.5);
	s2[7] = complex<double>(0.0, 0.0);

	s2[8] = complex<double>(0.5, 0.0);
	s2[9] = complex<double>(0.0, 0.0);
	s2[10] = complex<double>(0.0, 0.0);
	s2[11] = complex<double>(-0.5, 0.0);

	s2[12] = complex<double>(1.0, 0.0);
	s2[13] = complex<double>(0.0, 0.0);
	s2[14] = complex<double>(0.0, 0.0);
	s2[15] = complex<double>(1.0, 0.0);

	s3[0] = complex<double>(0.0, 0.0);
	s3[1] = complex<double>(rttwo, 0.0);
	s3[2] = complex<double>(0.0, 0.0);
	s3[3] = complex<double>(rttwo, 0.0);
	s3[4] = complex<double>(0.0, 0.0);
	s3[5] = complex<double>(rttwo, 0.0);
	s3[6] = complex<double>(0.0, 0.0);
	s3[7] = complex<double>(rttwo, 0.0);
	s3[8] = complex<double>(0.0, 0.0);

	s3[9] = complex<double>(0.0, 0.0);
	s3[10] = complex<double> (0.0, -rttwo);
	s3[11] = complex<double> (0.0, 0.0);
	s3[12] = complex<double> (0.0, rttwo);
	s3[13] = complex<double> (0.0, 0.0);
	s3[14] = complex<double> (0.0, -rttwo);
	s3[15] = complex<double> (0.0, 0.0);
	s3[16] = complex<double> (0.0, rttwo);
	s3[17] = complex<double> (0.0, 0.0);

	s3[18] = complex<double>(1.0, 0.0);
	s3[19] = complex<double>(0.0, 0.0);
	s3[20] = complex<double>(0.0, 0.0);
	s3[21] = complex<double>(0.0, 0.0);
	s3[22] = complex<double>(0.0, 0.0);
	s3[23] = complex<double>(0.0, 0.0);
	s3[24] = complex<double>(0.0, 0.0);
	s3[25] = complex<double>(0.0, 0.0);
	s3[26] = complex<double>(-1.0, 0.0);

	s3[27] = complex<double>(1.0, 0.0);
	s3[28] = complex<double>(0.0, 0.0);
	s3[29] = complex<double>(0.0, 0.0);
	s3[30] = complex<double>(0.0, 0.0);
	s3[31] = complex<double>(1.0, 0.0);
	s3[32] = complex<double>(0.0, 0.0);
	s3[33] = complex<double>(0.0, 0.0);
	s3[34] = complex<double>(0.0, 0.0);
	s3[35] = complex<double>(1.0, 0.0);

};

//---------------------------------------------------------------------------//

void hamiltonian::build(param_dat d, radical_dat rad) {

	double a = 1.76e8;
	int i, j, k;
	int dim_tmp;

	complex<double> hxtmp[m*m] = {};
	complex<double> hytmp[m*m] = {};
	complex<double> hztmp[m*m] = {};
	complex<double> hxtmp_hf[m*m] = {};
	complex<double> hytmp_hf[m*m] = {};
	complex<double> hztmp_hf[m*m] = {};
	complex<double> hhftmp[m*m] = {};
	complex<double> hhftmp1[m*m] = {};
	complex<double> hhftmp2[m*m] = {};
	complex<double> hhftmp3[m*m] = {};
	complex<double> hhftmp4[m*m] = {};

	hx = new complex<double>[m*m]();	
	hy = new complex<double>[m*m]();	
	hz = new complex<double>[m*m]();	
	h_hf = new complex<double>[m*m]();	
	hze = new complex<double>[m*m*d.ntheta]();

	dim_tmp = 2;

	cpy_mat(hx, s2, 2, 0);
	cpy_mat(hy, s2, 2, 4);
	cpy_mat(hz, s2, 2, 8);

	for (i=0; i<nspins; i++) {
		cpy_mat(hxtmp, hx, dim_tmp);
		cpy_mat(hytmp, hy, dim_tmp);
		cpy_mat(hztmp, hz, dim_tmp);
		if (rad.m_array[i] == 2) {
			kron(hxtmp, s2, hx, dim_tmp, 2, 12);
			kron(hytmp, s2, hy, dim_tmp, 2, 12);
			kron(hztmp, s2, hz, dim_tmp, 2, 12);
			dim_tmp = dim_tmp * 2;
		}
		else if (rad.m_array[i] == 3) {
			kron(hxtmp, s3, hx, dim_tmp, 3, 27);
			kron(hytmp, s3, hy, dim_tmp, 3, 27);
			kron(hztmp, s3, hz, dim_tmp, 3, 27);
			dim_tmp = dim_tmp * 3;
		}
	}

	for (i=0; i<nspins; i++) {
		dim_tmp = 2;
		cpy_mat(hxtmp_hf, s2, 2, 0);
		cpy_mat(hytmp_hf, s2, 2, 4);
		cpy_mat(hztmp_hf, s2, 2, 8);
		for (j=0; j<i; j++) {
			cpy_mat(hxtmp, hxtmp_hf, m);
			cpy_mat(hytmp, hytmp_hf, m);
			cpy_mat(hztmp, hztmp_hf, m);
			if (rad.m_array[j] == 2) {
				kron(hxtmp, s2, hxtmp_hf, dim_tmp, 2, 12);
				kron(hytmp, s2, hytmp_hf, dim_tmp, 2, 12);
				kron(hztmp, s2, hztmp_hf, dim_tmp, 2, 12);
				dim_tmp = dim_tmp*2;
			}
			else if (rad.m_array[j] == 3) {
				kron(hxtmp, s3, hxtmp_hf, dim_tmp, 3, 27);
				kron(hytmp, s3, hytmp_hf, dim_tmp, 3, 27);
				kron(hztmp, s3, hztmp_hf, dim_tmp, 3, 27);
				dim_tmp = dim_tmp*3;
			}
		}
		if (rad.m_array[i] == 2) {
			kron(hxtmp_hf, s2, hhftmp, dim_tmp, 2, 0);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i], dim_tmp*2);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*2);
			kron(hxtmp_hf, s2, hhftmp, dim_tmp, 2, 4);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+1], dim_tmp*2);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*2);
			kron(hxtmp_hf, s2, hhftmp, dim_tmp, 2, 8);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+2], dim_tmp*2);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*2);

			kron(hytmp_hf, s2, hhftmp, dim_tmp, 2, 0);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+3], dim_tmp*2);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*2);
			kron(hytmp_hf, s2, hhftmp, dim_tmp, 2, 4);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+4], dim_tmp*2);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*2);
			kron(hytmp_hf, s2, hhftmp, dim_tmp, 2, 8);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+5], dim_tmp*2);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*2);

			kron(hztmp_hf, s2, hhftmp, dim_tmp, 2, 0);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+6], dim_tmp*2);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*2);
			kron(hztmp_hf, s2, hhftmp, dim_tmp, 2, 4);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+7], dim_tmp*2);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*2);
			kron(hztmp_hf, s2, hhftmp, dim_tmp, 2, 8);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+8], dim_tmp*2);

			dim_tmp = dim_tmp * 2;
		}
		else if (rad.m_array[i] == 3) {
			kron(hxtmp_hf, s3, hhftmp, dim_tmp, 3, 0);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i], dim_tmp*3);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*3);
			kron(hxtmp_hf, s3, hhftmp, dim_tmp, 3, 9);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+1], dim_tmp*3);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*3);
			kron(hxtmp_hf, s3, hhftmp, dim_tmp, 3, 18);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+2], dim_tmp*3);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*3);

			kron(hytmp_hf, s3, hhftmp, dim_tmp, 3, 0);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+3], dim_tmp*3);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*3);
			kron(hytmp_hf, s3, hhftmp, dim_tmp, 3, 9);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+4], dim_tmp*3);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*3);
			kron(hytmp_hf, s3, hhftmp, dim_tmp, 3, 18);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+5], dim_tmp*3);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*3);

			kron(hztmp_hf, s3, hhftmp, dim_tmp, 3, 0);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+6], dim_tmp*3);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*3);
			kron(hztmp_hf, s3, hhftmp, dim_tmp, 3, 9);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+7], dim_tmp*3);
			cpy_mat(hhftmp1, hhftmp2, dim_tmp*3);
			kron(hztmp_hf, s3, hhftmp, dim_tmp, 3, 18);
			add_mat(hhftmp1, hhftmp, hhftmp2, a*rad.a_tensor[9*i+8], dim_tmp*3);

			dim_tmp = dim_tmp * 3;
		}
		for (j=i+1; j<nspins; j++) {
			if (rad.m_array[j] == 2) {
				cpy_mat(hhftmp, hhftmp2, dim_tmp);
				kron(hhftmp, s2, hhftmp2, dim_tmp, 2, 12);
				dim_tmp = dim_tmp*2;
			}
			if (rad.m_array[j] == 3) {
				cpy_mat(hhftmp, hhftmp2, dim_tmp);
				kron(hhftmp, s2, hhftmp2, dim_tmp, 3, 27);
				dim_tmp = dim_tmp*3;
			}
		}

		add_mat(hhftmp2, hhftmp3, hhftmp4, 1.0, m);
		cpy_mat(hhftmp3, hhftmp4, m);
	}

	cpy_mat(h_hf, hhftmp4, m);

	for (i=0; i<d.ntheta; i++) {
		for (j=0; j<m; j++) {
			for (k=0; k<m; k++) {
				hze[i*m*m+j*m+k] = d.omega_0*(hx[j*m+k]*sin(d.theta[i]) 
						+ hz[j*m+k]*cos(d.theta[i])) + h_hf[j*m+k];
			}
		}
	}

};

//---------------------------------------------------------------------------//

void hamiltonian::noise(param_dat d) {

	int i, j, k;

	h_td = new complex<double>[m*m]();	
	dhdt = new complex<double>[m*m]();	
	h_td0 = new complex<double>[m*m]();	
	hfourier = new complex<double>[m*m*d.nfourier]();	

	for (i=0; i<d.nfourier; i++) {
		for (j=0; j<m; j++) {
			for (k=0; k<m; k++) {
				hfourier[i*m*m+j*m+k] = d.coeff_rf[i]*(sin(d.theta_rf[i])*
					(hx[j*m+k]*cos(d.phi_rf[i]) + hy[j*m+k]*sin(d.phi_rf[i]))
					+ hz[j*m+k]*cos(d.theta_rf[i]));
				h_td0[j*m+k] = h_td0[j*m+k] + hfourier[i*m*m+j*m+k]*sin(d.phase_rf[i]);
			}
		}
	}

};


//---------------------------------------------------------------------------//

void hamiltonian::time_step(param_dat d) {

	int i, j, k;
	double arg;

	complex<double> hi[m*m] = {};
	complex<double> dhdti[m*m] = {};

	for (i=0; i<d.nfourier; i++) {
		arg = d.w_rf[i]*d.t_array[inow] + d.phase_rf[i];
		for (j=0; j<m; j++) {
			for(k=0; k<m; k++) {
				hi[j*m+k] = hi[j*m+k] + sin(arg)*hfourier[i*m*m+j*m+k];
				dhdti[j*m+k] = dhdti[j*m+k] + d.w_rf[i]*cos(arg)*hfourier[i*m*m+j*m+k];
			}
		}
	}

	inow = inow + 1;
	cpy_mat(h_td, hi, m);
	cpy_mat(dhdt, dhdti, m);

};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
