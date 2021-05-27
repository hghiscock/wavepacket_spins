#pragma once

using namespace std;
using namespace N;
namespace N1
{
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

class hamiltonian {

	complex<double> *s2;
	complex<double> *s3;	

    public:
	int m, nspins, inow;
	complex<double> *hx;	
	complex<double> *hy;	
	complex<double> *hz;	
	complex<double> *h_hf;	
	complex<double> *hze;	
	complex<double> *h_td;	
	complex<double> *h_td0;	
	complex<double> *dhdt;
	complex<double> *hfourier;

	hamiltonian(radical_dat);
	void build(param_dat, radical_dat);
	void noise(param_dat);
	void time_step(param_dat);

};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
}
