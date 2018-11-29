#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include "cdlasso.h"
#include "function.h"

using namespace std;

int main() {
	int n, p, maxit;
	char np;

	cout << "Use default settings (n = 200, p = 10000)? (Y/N)\n";
	cin >> np;
	if (np == 'N' or np == 'n') {
		cout << "Please provide the sample size:\n";
		cin >> n;
		cout << "Please provide the number of parameters (>15):\n";
		cin >> p;
	}
	else {
		n = 200;
		p = 10000;
	}

	double *beta = new double[p];
	double t1, t2, t3;
	ofstream myfile;
	ofstream lasso;
	//ofstream check;
	ofstream time;

	time.open("time.txt");

	for (int i = 0; i < 15; i++) {
		beta[i] = pow(-1.0, i) * 2 * exp(-i / 15.0);
	}

	for (int i = 15; i < p; i++) {
		beta[i] = 0;
	}

	cout << "Single Run (Y) or Simulation (N)?\n";
	cin >> np;
	if (np == 'Y' or np == 'y') {
		myfile.open("single_test.txt");
		double lambda;
		cout << "Please Provide Lambda\n";
		cin >> lambda;

		t1 = clock();
		XY_old test(n, p, 1, beta);
		XY_new test1(test);

		t2 = clock();
		double *beta1 = cdLasso(test1.x, test1.y, (test1.n1 + 1), p, lambda*(test1.n1 + 1));

		//check.open("check.txt");
		//for (int i = 0; i < (test1.n1 + 1); i++) {
		//	check << test1.y[i] << ",";
		//	for (int j = 0; j < p; j++) {
		//		check << test1.x[i][j] << ",";
		//	}
		//	check << '\n';
		//}
		//check.close();

		t3 = clock();
		double s = test1.LG_all(beta1, lambda*(test1.n1 + 1));
		myfile << s << "\n";
		for (int i = 0; i < p; i++) {
			myfile << beta1[i] << "\n";
		}
		myfile.close();

		test.delete_old();
		test1.delete_new();
		delete[] beta;
		delete[] beta1;

		time << (t2 - t1) / (double)CLOCKS_PER_SEC << endl;
		time << (t3 - t2) / (double)CLOCKS_PER_SEC << endl;

		time.close();
	}
	else {
		cout << "Please Provide number of Replications\n";
		int rep;
		cin >> rep;
		cout << "Please Provide Max number of paramters selected (<1000)\n";
		cin >> maxit;
		cout << "Please Provide Initial Seed\n";
		int s0;
		cin >> s0;

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