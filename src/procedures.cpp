#include <math.h>
#include <complex>
#include "procedures.h"
#include <iostream>

using namespace std;

//---------------------------------------------------------------------------//

void kron(complex<double>* mat1, complex<double>* mat2, complex<double>* mat3, int dim1, int dim2, int offset) {

	int i, j, k, l;

	for (i=0; i<dim1; i++) {
		for (j=0; j<dim2; j++) {
			for (k=0; k<dim1; k++) {
				for (l=0; l<dim2; l++) {
					mat3[i*dim1*dim2*dim2 + j*dim1*dim2 + k*dim2 + l] = mat1[i*dim1+k] * mat2[j*dim2+l+offset];
				}
			}
		}
	}	

}

//---------------------------------------------------------------------------//

void cpy_mat(complex<double>* mat1, complex<double>* mat2, int dim, int offset1, int offset2) {

	int i, itmp1, itmp2;

	for (i=0; i < dim*dim; i++) {
		itmp1 = i + offset1;
		itmp2 = i + offset2;
		mat1[itmp1] = mat2[itmp2];
	}
}

//---------------------------------------------------------------------------//

void cpy_mat(complex<double>* mat1, complex<double>* mat2, int dim, int offset) {

	int i, itmp;

	for (i=0; i < dim*dim; i++) {
		itmp= i + offset;
		mat1[i] = mat2[itmp];
	}
}

//---------------------------------------------------------------------------//

void cpy_mat(complex<double>* mat1, complex<double>* mat2, int dim) {

	int i;

	for (i=0; i < dim*dim; i++) {
		mat1[i] = mat2[i];
	}
}

//---------------------------------------------------------------------------//

void add_mat(complex<double>* mat1, complex<double>* mat2, complex<double>* mat3, double factor, int dim) {

	int i;

	for (i=0; i < dim*dim; i++) {
		mat3[i] = mat1[i] + (mat2[i] * factor);
	}
}

//---------------------------------------------------------------------------//

void add_mat(complex<double>* mat1, complex<double>* mat2, complex<double>* mat3, complex<double> factor, int dim) {

	int i;

	for (i=0; i < dim*dim; i++) {
		mat3[i] = mat1[i] + (mat2[i] * factor);
	}
}

//---------------------------------------------------------------------------//

void print_mat(double* mat, int dim, int offset) {

	int i, j;

	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			cout << mat[i*dim + j + offset] << '\t';
		}
		cout << endl;
	}

}

//---------------------------------------------------------------------------//

void print_mat(complex<double>* mat, int dim, int offset) {

	int i, j;

	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			cout << mat[i*dim + j + offset] << '\t';
		}
		cout << endl;
	}

}

//---------------------------------------------------------------------------//

void matmul(complex<double>* mat1, complex<double>* mat2, complex<double>* mat3, int dim, int offset1, int offset2, int offset3) {

	int i, j, k;
	complex<double> tmp;

	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			tmp = complex<double>(0.0, 0.0);
			for (k=0; k<dim; k++) {
				tmp = tmp + mat1[offset1+i*dim+k]*mat2[offset2+j+k*dim];
			}
			mat3[i*dim+j] = tmp;
		}
	}


}

//---------------------------------------------------------------------------//

void matmul(complex<double>* mat1, complex<double>* mat2, complex<double>* mat3, int dim) {

	int i, j, k;
	complex<double> tmp;

	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			tmp = complex<double>(0.0, 0.0);
			for (k=0; k<dim; k++) {
				tmp = tmp + mat1[i*dim+k]*mat2[j+k*dim];
			}
			mat3[i*dim+j] = tmp;
		}
	}


}

//---------------------------------------------------------------------------//

template <class T, class U>
void print_mat_test(T mat, U dim) { 

	int i, j;

	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			cout << mat[i*dim + j] << '\t';
		}
		cout << endl;
	}

}

//---------------------------------------------------------------------------//

int foo(){
	return 1;
}
