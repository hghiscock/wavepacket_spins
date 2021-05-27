#pragma once
#include <string>

using namespace std;
namespace N
{
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

class param_dat {

    public:
	int ncores;
	double a[4], b[4];
	double omega_0, gamma_e;

	double w_min, w_max;
	int nfourier;
	double *w_rf;
	double *phase_rf;
	double *coeff_rf;
	double *theta_rf;
	double *phi_rf;

	int t_ind;
	bool iden;
	bool rms;

	int nt;
	double time, tau, k;
	double *t_array;

	double theta_min, theta_max;
	double phi_min, phi_max;
	double d_theta, d_phi;
	int ntheta, nphi;
	double phi;
	double *theta;

	string rada, radb;
	int na, nb;

	string outputfile, logfile;
	int start, finish;

	void read_in();
	void initialise();

};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

class radical_dat {

    public:
 	static int nrads;

	int nspins, m;
	int *m_array;
	double *a_tensor;
	double alpha;

	radical_dat() {nrads++; };
	void read_in(param_dat);
};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
}
