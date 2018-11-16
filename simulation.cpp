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

int main()
{
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
		double *beta1 = cdLasso(test1.x, test1.y, (test1.n1 + 1), p, lambda*pow(test1.n, 2));

		t3 = clock();
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
		double lambda[20];
		double *LG = new double[maxit];
		double *stage2_beta = new double[p];
		int fold = floor(n/10);
		for (int i = 0; i < 20; i++) {
			lambda[i] = 0.01 * (i + 1);
		}

		for (int z = 0; z < rep; z++) {
			myfile.open("simulation_APML0.txt",ios_base::app);
			XY_old test(n, p, z + s0, beta);
			XY_new test1(test);
			double *best_beta;
			double min_LG = 1e8, min_lambda = 1;
			int min_k = maxit;
			
			for(int m = 0; m < 20; m++)
			{
				for (int i = 0; i < maxit; i++) {
					LG[i] = 0;
				}

				for (int i = 0; i < 10; i++) {
					int low = fold * i;
					int high = fold * (i + 1) - 1;
					if (i == 9) {
						high = n - 1;
					}

					int d = high - low + 1;
					int *target = new int[d];
					for (int j = 0; j < d; j++) {
						target[j] = low + j;
					}
				
					test1.seperate(target, d);
					double *rep_beta = cdLasso(test1.x_train, test1.y_train, (test1.n1_train + 1), p, lambda[m]*pow(test1.n_train, 2));
					int *rank = beta_rank(rep_beta, p);
					for (int j = 0; j < maxit; j++) {
						for (int j_ = 0; j_ < p; j_++) {
							if (rank[j_] <= j) {
								stage2_beta[j_] = rep_beta[j_];
							}
							else {
								stage2_beta[j_] = 0;
							}
						}
						LG[j] += test1.LG(stage2_beta);
					}
					delete[] rank;
					delete[] rep_beta;
					delete[] target;
					test1.delete_new_cv();
				}

				for (int j = 0; j < maxit; j++) {
					if (min_LG > LG[j]) {
						min_LG = LG[j];
						min_lambda = lambda[m];
						min_k = j + 1;
					}
				}
				cout << m << "\n";
			}

			myfile << min_LG << "," << min_lambda << "," << min_k << ",";
			best_beta = cdLasso(test1.x, test1.y, (test1.n1 + 1), p, min_lambda*pow(test1.n, 2));
			int *rank = beta_rank(best_beta, p);
			for (int j = 0; j < p; j++) {
				if (rank[j] < min_k) {
					stage2_beta[j] = best_beta[j];
					myfile << stage2_beta[j] << ",";
				}
				else {
					stage2_beta[j] = 0;
					myfile << stage2_beta[j] << ",";
				}
			}
			myfile << "\n";
			delete[] best_beta;
			delete[] rank;
			test1.delete_new();
			test.delete_old();
			myfile.close();
			cout << "-----------" << z << endl;
		}
		delete[] LG;
		delete[] stage2_beta;
	myfile.open("simulation_APML0.txt",ios_base::app);
	time_t result = std::time(nullptr);
	myfile << asctime(localtime(&result)) << "\n";
	myfile.close();
	}

	return 0;
}