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

int main(){

	int i;

	param_dat p1;
	p1.read_in();
	p1.initialise();

	radical_dat r1;
	r1.read_in(p1);

	radical_dat r2;
	r2.read_in(p1);

	hamiltonian h1(r1);
	hamiltonian h2(r2);

	h1.build(p1, r1);
	h2.build(p1, r2);

	h1.noise(p1);
	h2.noise(p1);

	spin_corr sc1(p1, h1);
	spin_corr sc2(p1, h2);

	calc_dat c(p1);

	for (i=0; i<p1.nt; i++) {

		sc1.calc_sc(p1, h1);
		sc2.calc_sc(p1, h2);

		c.calc_ps(p1, sc1, sc2);
		c.integrate_sy(p1);

	}

	cout << c.phis[0] << "\t" << c.phis[1] << endl;

	c.write_results(p1);

	return 0;
};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

