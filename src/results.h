#pragma once
#include "param.h"
#include "procedures.h"
#include "hamiltonians.h"
#include "sc_tensor.h"

using namespace std;
using namespace N;
using namespace N1;
using namespace N2;
namespace N3
{
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

class calc_dat {

    public:
	int jnow;

	double *ps;
	double *phis;
	double *tarray;

	calc_dat(param_dat);
	void integrate_sy(param_dat);
	void calc_ps(param_dat, spin_corr, spin_corr);
	void write_results(param_dat);

};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
}
