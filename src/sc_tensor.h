#pragma once

using namespace std;
using namespace N;
using namespace N1;
namespace N2
{
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

class spin_corr {

    public:
	int m, mprime;
	int inow;
	complex<double> xj;

	double *rab;
	complex<double> *p;
	complex<double> *q;
	complex<double> *rtmp;

	spin_corr(param_dat, hamiltonian);
	void evolve(param_dat, hamiltonian);
	void evolve_t_ind(param_dat, hamiltonian);
	void calc_rab(param_dat);
	void ex_rab(param_dat);
	void calc_sc(param_dat, hamiltonian);
};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
}
