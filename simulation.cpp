#define _CRT_SECURE_NO_DEPRECATE
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
int p = 10000;
int maxit = 50;
int rep = 1;
int s0 = 1;
double lambda = 0.01;
char np = 'Y';

int main(int argc, char *argv[])
{
//	int n, p, maxit, rep, s0;
//	double lambda;
//	char np;

//	n = atoi(argv[1]);
//	p = atoi(argv[2]);
//	maxit = atoi(argv[3]);
//	rep = atoi(argv[4]);
//	s0 = atoi(argv[5]);
//	lambda = atof(argv[6]);
//	np = argv[7][0];

//	cout << "Please provide the sample size:\n";
//	cin >> n;
//	cout << "Please provide the number of parameters (>15):\n";
//	cin >> p;
//	cout << "Single Run (Y) or Simulation (N)?\n";
//	cin >> np;
//	cout << "Please Provide Lambda\n";
//	cin >> lambda;
//	cout << "Please Provide number of Replications\n";
//	cin >> rep;
//	cout << "Please Provide Max number of paramters selected (<1000)\n";
//	cin >> maxit;
//	cout << "Please Provide Initial Seed\n";
//	cin >> s0;

	double *beta = new double[p];
	double t1, t2, t3;
	ofstream myfile;
	ofstream lasso;
	ofstream time;

	time.open("time.txt");

	for (int i = 0; i < 15; i++) {
		beta[i] = pow(-1.0, i) * 2 * exp(-i / 15.0);
	}

	for (int i = 15; i < p; i++) {
		beta[i] = 0;
	}

	if (np == 'Y' or np == 'y') {
		myfile.open("single_test.txt");

		t1 = clock();
		XY_old test(n, p, s0, beta);
		XY_new test1(test);

		t2 = clock();
//		double *beta1 = cdLasso(test1.x, test1.y, (test1.n1 + 1), p, lambda*(test1.n1 + 1));
		test.print_old();

		cv_path(test1, 3);

		t3 = clock();
//		double s = test1.LG_all(beta1, lambda*(test1.n1 + 1));
//		myfile << s << "\n";
//		for (int i = 0; i < p; i++) {
//			myfile << beta1[i] << "\n";
//		}
		myfile.close();

		test.delete_old();
		test1.delete_new();
		delete[] beta;
//		delete[] beta1;

		time << (t2 - t1) / (double)CLOCKS_PER_SEC << endl;
		time << (t3 - t2) / (double)CLOCKS_PER_SEC << endl;

		time.close();
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
