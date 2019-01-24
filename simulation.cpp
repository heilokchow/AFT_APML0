#define _CRT_SECURE_NO_DEPRECATE
#define TRUE_PARAMETER 15

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include "cdlasso.h"
#include "function.h"

using namespace std;

int n = 200;
int p = 1000;
int maxit = 50;
int rep = 1;
int s0 = 1;
double lambda = 0.01;
char np = 'Y';

int main(int argc, char *argv[])
{
	double *beta = new double[p];
	ofstream myfile;
	ofstream lasso;

//    ORIGINAL TESTING
//	for (int i = 0; i < TRUE_PARAMETER; i++) {
//		beta[i] = pow(-1.0, i) * 2 * exp(-i / 15.0);
//	}

//	for (int i = 15; i < p; i++) {
//		beta[i] = 0;
//	}

//    ROC CURVE SETTING
	for (int i = 0; i < TRUE_PARAMETER; i++) {
		beta[i] = pow(-1.0, i)* 2 * std::exp(-i / 15.0);
	}
	for (int i = TRUE_PARAMETER; i < 100; i++) {
        beta[i] = 0;
	}
	for (int i = 100; i < p; i++) {
		beta[i] = 0;
	}

	if (np == 'Y' or np == 'y') {

        ALPath pre_lasso(p);
        ALPath pre_apml0(p);

        for (int z = 0; z < rep; z++) {
            XY_old test(n, p, s0 + z, beta);
            XY_new test1(test);
            cv_path(test1, 5, pre_lasso, pre_apml0);
            test.delete_old();
            test1.delete_new();
        }

        pre_lasso.ALprint("_lasso");
        pre_apml0.ALprint("_apml0");

	delete[] beta;
	}
	else {
		for (int z = 0; z < rep; z++) {
			myfile.open("simulation_APML0.txt",ios_base::app);
			lasso.open("simulation_LASSO.txt",ios_base::app);
			XY_old test(n, p, z + s0, beta);
			XY_new test1(test);
			test1.cross_validation(myfile, lasso, maxit);
			test1.delete_new();
			test.delete_old();
			myfile.close();
			lasso.close();
			cout << "-----------" << z << endl;
		}

	myfile.open("simulation_APML0.txt",ios_base::app);
	lasso.open("simulation_LASSO.txt",ios_base::app);
	time_t result = std::time(nullptr);
	myfile << asctime(localtime(&result)) << "\n";
	lasso << asctime(localtime(&result)) << "\n";
	myfile.close();
	lasso.close();
	}

	return 0;
}
