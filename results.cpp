#include <math.h>
#include <complex>
#include "param.h"
#include "procedures.h"
#include "hamiltonians.h"
#include "sc_tensor.h"
#include "results.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace N;
using namespace N1;
using namespace N2;
using namespace N3;

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

calc_dat::calc_dat(param_dat d) {

	int i;
	double third;

	phis = new double[d.ntheta]();
	ps = new double[d.ntheta]();
	tarray = new double[d.nt]();

	tarray[0] = 0.0;
	for (i=1; i<d.nt; i++) {
		tarray[i] = tarray[i-1] + d.tau;
	}

	third = 1.0/3.0;
	tarray[0] = third*exp(-1.0*d.k*tarray[0])*d.tau*d.k;
	tarray[d.nt-1] = third*exp(-1.0*d.k*tarray[d.nt-1])*d.tau*d.k;
	for (i=1; i<d.nt-1; i++) {
		tarray[i] = (double(i%2) + 1.0)*2.0*third*
			exp(-1.0*d.k*tarray[i])*d.tau*d.k;
	}

	jnow = 0;
}


//---------------------------------------------------------------------------//

void calc_dat::integrate_sy(param_dat d) {

	int i;

	for (i=0; i<d.ntheta; i++) {
		phis[i] = phis[i] + ps[i]*tarray[jnow];
	}

	jnow = jnow + 1;

}

//---------------------------------------------------------------------------//

void calc_dat::calc_ps(param_dat d, spin_corr sca, spin_corr scb) {

	int i, j;
	int itmp, jtmp;

	for (i=0; i<d.ntheta; i++) {
		ps[i] = 0.25; 
		itmp = i*9;
		for (j=0; j<9; j++) {
			jtmp = itmp+j;
			ps[i] = ps[i] + sca.rab[jtmp]*scb.rab[jtmp];
		}
	}

};

//---------------------------------------------------------------------------//

void calc_dat:: write_results(param_dat d) {

	int i;

	ofstream myfile ("phis.dat");
	if (myfile.is_open()) {

		for (i=0; i<d.ntheta; i++) {
			myfile << phis[i] << endl;
		}
		myfile.close();
	}
	else cout << "Unable to open file";

};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

