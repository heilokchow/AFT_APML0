#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <math.h>

#define ITERATOR 10000

using namespace std;

void swap_(double &x, double &y) {
	double temp = x;
	x = y;
	y = temp;
}

int partition(double *x, double *z, int low, int high)
{
	double pivot = z[high];
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++)
	{
		if (z[j] <= pivot)
		{
			i++;
			swap_(z[i], z[j]);
			swap_(x[i], x[j]);
		}
	}
	swap_(z[i + 1], z[high]);
	swap_(x[i + 1], x[high]);

	return (i + 1);
}

void quickSort(double *x, double *z, int low, int high)
{
	if (low < high)
	{
		int pi = partition(x, z, low, high);

		quickSort(x, z, low, pi - 1);
		quickSort(x, z, pi + 1, high);
	}
}

double key_sort(double *x, double *z, int n) {

	double *xx = new double[n + 1];
	double *zz = new double[n + 1];

	for (int i = 0; i < (n + 1); i++)
	{
		xx[i] = x[i];
		zz[i] = z[i];
	}

	quickSort(xx, zz, 0, n);

	double s = 0, s1 = 0;

	for (int i = 0; i < (n + 1); i++)
	{
		s += xx[i];
	}

	int i = 0;
	while (s1 < s / 2)
	{
		s1 += xx[i];
		i++;
	}
	i--;

	double output = zz[i];
	delete[] xx;
	delete[] zz;
	return output;
}

double *cdLasso(double **A, double *B, int n, int p, double lambda) {

	//ofstream myfile;
	//myfile.open("test_LG.txt");

	double *r = new double[n];
//	double **x = new double *[p];
    double *x = new double[n + 1];
	double s = 0, s1 = 0;
	double *z = new double[n + 1];
	double *beta = new double[p];
	double beta_ = 0;
	double sum_ = 0;
	double *fdd = new double[p];
	double *zdd = new double[p];
	double c, c1;
	double *r_ = new double[n];
	int flag1 = 0;

    x[n] = lambda;
	for (int i = 0; i < p; i++)
	{
//      CACHE OTIMIZTION WITH MEMORY TRADE OFF
//		x[i] = new double[n + 1];
		for (int j = 0; j < n; j++)
		{
//			x[i][j] = abs(A[j][i]);
			r[j] = B[j];
		}
//		x[i][n] = lambda;
		beta[i] = 0;
	}

	z[n] = 0;
	int k = 0;
	int flag = 0;

	for (int i = 0; i < p; i++)
	{
		s1 += lambda * abs(beta[i]);
	}

	for (int i = 0; i < p; i++)
	{
		if (beta[i] > 1e-8)
		{
			fdd[i] = lambda;
			zdd[i] = 0;
		}
		else if (beta[i] < -1e-8)
		{
			fdd[i] = -lambda;
			zdd[i] = 0;
		}
		else
		{
			fdd[i] = 0;
			zdd[i] = lambda;
		}
	}

	for (int i = 0; i < n; i++)
	{
		if (r[i] > 1e-8)
		{
			for (int j = 0; j < p; j++)
			{
				fdd[j] += -A[i][j];
			}
		}
		else if (r[i] < -1e-8)
		{
			for (int j = 0; j < p; j++)
			{
				fdd[j] += A[i][j];
			}
		}
		else
		{
			for (int j = 0; j < p; j++)
			{
				zdd[j] += abs(A[i][j]);
			}
		}
	}

	while (flag < ITERATOR)
	{
		c = 10;
		k = 0;
		for (int i = 0; i < p; i++)
		{
			c1 = min(fdd[i] + zdd[i], -fdd[i] + zdd[i]);
			if (c1 < c)
			{
				c = c1;
				k = i;
			}
		}

		//cout << c << endl;

		if (c >= 0) {
			break;
		}

		if (abs(beta[k]) > 1e-8) {
			for (int i = 0; i < n; i++)
			{
				if (A[i][k] != 0)
				{
					z[i] = (r[i] + beta[k] * A[i][k]) / A[i][k];
				}
			}
		}
		else {
			for (int i = 0; i < n; i++)
			{
				if (A[i][k] != 0)
				{
					z[i] = r[i] / A[i][k];
				}
			}
		}

        for (int j = 0; j < n; j++)
		{
			x[j] = abs(A[j][k]);
		}
		beta_ = key_sort(x, z, n);

		s = 0;
		for (int i = 0; i < n; i++)
		{
			r_[i] = r[i] + A[i][k] * (beta[k] - beta_);

			if (i < (n - 1) && r_[i] < 0)
				s += -2 * r_[i];

			flag1 = 0;
			flag1 += (r[i] < -1e-8 && r_[i] < -1e-8) ? 1 : 0;
			flag1 += (r[i] > 1e-8 && r_[i] > 1e-8) ? 1 : 0;
			flag1 += (abs(r[i]) <= 1e-8 && abs(r_[i]) <= 1e-8) ? 1 : 0;

			if (flag1 == 0)
			{
				if (r[i] > 1e-8)
				{
					if (r_[i] < -1e-8)
					{
						for (int j = 0; j < p; j++)
						{
							fdd[j] += 2 * A[i][j];
						}
					}
					else
					{
						for (int j = 0; j < p; j++)
						{
							fdd[j] += A[i][j];
							zdd[j] += abs(A[i][j]);
						}
					}
				}
				else if (r[i] < -1e-8)
				{
					if (r_[i] > 1e-8)
					{
						for (int j = 0; j < p; j++)
						{
							fdd[j] += -2 * A[i][j];
						}
					}
					else
					{
						for (int j = 0; j < p; j++)
						{
							fdd[j] += -A[i][j];
							zdd[j] += abs(A[i][j]);
						}
					}
				}
				else
				{
					if (r_[i] > 1e-8)
					{
						for (int j = 0; j < p; j++)
						{
							fdd[j] += -A[i][j];
							zdd[j] += -abs(A[i][j]);
						}
					}
					else
					{
						for (int j = 0; j < p; j++)
						{
							fdd[j] += A[i][j];
							zdd[j] += -abs(A[i][j]);
						}
					}
				}
			}
			r[i] = r_[i];
		}

		if (beta[k] > 1e-8) {
			fdd[k] += -lambda;
		}
		else if (beta[k] < -1e-8) {
			fdd[k] += lambda;
		}
		else {
			zdd[k] += -lambda;
		}

		if (beta_ > 1e-8) {
			fdd[k] += lambda;
		}
		else if (beta_ < -1e-8) {
			fdd[k] += -lambda;
		}
		else {
			zdd[k] += lambda;
		}

		s1 += lambda * (abs(beta_) - abs(beta[k]));
		sum_ = s + s1;

		beta[k] = beta_;

		//myfile << sum_ << "\n";
		flag++;
	}

//	for (int i = 0; i < p; i++)
//	{
//		delete[] x[i];
//	}

	delete[] x;
	delete[] r;
	delete[] r_;
	delete[] z;
	delete[] fdd;
	delete[] zdd;
	//myfile.close();

	return beta;
}
