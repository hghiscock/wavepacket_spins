#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include "param.h"

using namespace std;
using namespace N;

double const pi = 3.14159265359;

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

void param_dat::read_in(){

	string line;
	double b0, rtthree, tmin, tmax, pmin, pmax;

	ifstream indat("input/input.dat");

	getline(indat, line);
	getline(indat, line);
	b0 = stod(line);
	getline(indat, line);
	gamma_e = stod(line);
	getline(indat, line);
	k = stod(line);

	getline(indat, line);
	getline(indat, line);
	tmin = stod(line);
	getline(indat, line);
	tmax = stod(line);
	getline(indat, line);
	ntheta = stoi(line);

	getline(indat, line);
	getline(indat, line);
	time = stod(line);
	getline(indat, line);
	nt = stoi(line);

	getline(indat, line);
	getline(indat, line);
	rada = line;
	getline(indat, line);
	na = stoi(line);

	getline(indat, line);
	getline(indat, line);
	radb = line;
	getline(indat, line);
	nb = stoi(line);

	getline(indat, line);
	getline(indat, line);
	outputfile = line;

	getline(indat, line);
	getline(indat, line);
	ncores = stoi(line);

	indat.close();

	theta_min = tmin * pi;
	theta_max = tmax * pi;
	d_theta = (theta_max - theta_min) / double(ntheta);

	omega_0 = b0 * gamma_e;
	tau = time / double(nt);
	rtthree = 1.0 / sqrt(3.0);

	a[0] = 0.5*(1.0-rtthree)*tau;
	a[1] = rtthree*tau;
	a[2] = -rtthree*tau;
	a[3] = 0.5*(1.0+rtthree)*tau;

	b[0] = 0.0;
	b[1] = 0.5*(0.5+rtthree)*tau;
	b[2] = 0.5*tau;
	b[3] = 0.5*(0.5-rtthree)*tau;

};

//---------------------------------------------------------------------------//

void param_dat::initialise() {

	int i, j;
	double dw, t_tmp;

	string line;
	ifstream indat("input/noise.dat");

	getline(indat, line);
	getline(indat, line);
	nfourier = stoi(line);

	if (nfourier == 0){
		t_ind = 1;
	}
	else {
		t_ind = 0;
	}

	w_rf = new double[nfourier];
	phase_rf = new double[nfourier];
	coeff_rf = new double[nfourier];
	theta_rf = new double[nfourier];
	phi_rf = new double[nfourier];

	getline(indat, line);

	for (i=0; i<nfourier; i++) {
		getline(indat,line);
		w_rf[i] = stod(line);
		getline(indat,line);
		phase_rf[i] = stod(line) * pi;
		getline(indat,line);
		coeff_rf[i] = stod(line) * gamma_e * 1.0e-6;
		phi_rf[i] = 0.0;
		theta_rf[i] = 0.0;
	}

	indat.close();

	t_array = new double[3*nt];

	t_tmp = 0.0;
	for (i=0; i<nt; i++) {
		t_tmp = t_tmp + a[0];
		for (j=1; j<4; j++) {
			t_array[3*i-1+j] = t_tmp;
			t_tmp = t_tmp + a[j];
		}
	}

	theta = new double[ntheta];

	for (i=0; i<ntheta; i++) {
		theta[i] = theta_min + double(i)*d_theta;
	}

	phi = 0.0 * pi;
};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

int radical_dat::nrads=0;

//---------------------------------------------------------------------------//

void radical_dat::read_in (param_dat ptmp){

	int i, j, k, mtmp;
	string line;

	mtmp = 2;
	if (nrads == 1) {

		ifstream indat("input/" + ptmp.rada);
		nspins = ptmp.na;
		m_array = new int[nspins];
		for (i=0; i<nspins; i++) {
			getline(indat, line);
			m_array[i] = stoi(line);
			mtmp = mtmp * m_array[i];
		}
		m = mtmp;

		a_tensor = new double[(nspins+1)*9];
		for (i=0; i<nspins; i++) {
			getline(indat, line);
			for (j=0; j<9; j++) {
				getline(indat, line);
				a_tensor[9*i+j] = stod(line);
			}
		}
		indat.close();
	}
	else {
		ifstream indat("input/" + ptmp.radb);
		nspins = ptmp.nb;
		m_array = new int[nspins];
		for (i=0; i<nspins; i++) {
			getline(indat, line);
			m_array[i] = stoi(line);
			mtmp = mtmp * m_array[i];
		}
		m = mtmp;

		a_tensor = new double[(nspins+1)*9];
		for (i=0; i<nspins; i++) {
			getline(indat, line);
			for (j=0; j<9; j++) {
				getline(indat, line);
				a_tensor[9*i+j] = stod(line);
			}
		}
		indat.close();
	}
};

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
