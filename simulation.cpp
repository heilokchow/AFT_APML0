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
	int n, p;
	char np;
	double lambda;
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

	//double **A = new double *[n];
	//double *B = new double[n];
	double *beta = new double[p];
	double t1, t2, t3;
	ofstream myfile;
	ofstream time;

	time.open("time.txt");
	myfile.open("test0.txt");

	t1 = clock();

	for (int i = 0; i < 15; i++)
	{
		beta[i] = pow(-1.0, i) * 2 * exp(-i / 15.0);
	}

	for (int i = 15; i < p; i++)
	{
		beta[i] = 0;
	}

	XY_old test(n, p, 1, beta);
	XY_new test1(test);
	
	//seed_seq seed{ 1 };
	//mt19937 e(seed);
	//normal_distribution<> normal_dist(0, 1);

	//for (int i = 0; i < n; i++)
	//{
	//	A[i] = new double[p];
	//}

	//for (int i = 0; i < n; i++)
	//{
	//	B[i] = normal_dist(e);
	//	for (int j = 0; j < p; j++)
	//	{
	//		A[i][j] = normal_dist(e);
	//		if (beta[j] != 0)
	//		{
	//			B[i] += A[i][j] * beta[j];
	//		}
	//	}
	//}

	t2 = clock();
	double *beta1 = cdLasso(test1.x, test1.y, n, p, 0.1);

	t3 = clock();

	for (int i = 0; i < p; i++) {
		myfile << beta1[i] << "\n";
	}

	myfile.close();
	//for (int i = 0; i < n; i++)
	//{
	//	delete A[i];
	//}

	//delete[] A;
	//delete[] B;
	test.delete_old();
	test1.delete_new();
	delete[] beta;
	delete[] beta1;
	//delete[] beta1;

	time << (t2 - t1) / (double)CLOCKS_PER_SEC << endl;
	time << (t3 - t2) / (double)CLOCKS_PER_SEC << endl;

	time.close();

	return 0;
}