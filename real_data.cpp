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
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

using namespace std;

int main()
{
    ifstream y_input("y_input_1.txt");
    ifstream status_input("status_input_1.txt");
    ifstream x_input("x_input_1.txt");
    ofstream out;
	ofstream myfile;
	ofstream lasso;

	out.open("real_out.txt", ios_base::app);
	myfile.open("real_random_cv.txt", ios_base::app);
	lasso.open("real_random_cv_LASSO.txt", ios_base::app);

    int n, p;

    cout << "Please Provide the Sample Size:\n";
    cin >> n;
    cout << "Please Provide number of Parameters\n";
    cin >> p;

    double *y = new double[n];
    double *status = new double[n];
    double **x = new double*[n]; 
    for(int i = 0; i < n; i++)
    {
        x[i] = new double[p];
    }
    
    string y1, x1, status1;
    
    for(int i = 0; i < n; i++)
    {
        getline(y_input, y1, '\n');
        getline(status_input, status1, '\n');
        y[i] = atof(y1.c_str());
        status[i] = atof(status1.c_str());

        for(int j = 0; j < (p - 1); j++)
        {
            getline(x_input, x1, ',');
            x[i][j] = atof(x1.c_str());
        }
        getline(x_input, x1, '\n');
        x[i][p - 1] = atof(x1.c_str());
    }

	char sr;
	cout << "Single Run?(Y/N)\n";
	cin >> sr;

	if (sr == 'Y' || sr == 'y') {
		double *y_pro_lasso = new double[n];
		double *y_pro_apml0 = new double[n];
		ofstream prognostic_y;
		prognostic_y.open("pronostic_y.txt", ios_base::app);

		XY_old test0;
		test0.x = x;
		test0.y = y;
		test0.status = status;
		test0.n = n;
		test0.p = p;
		
		XY_new test00(test0);
		test00.cross_validation(myfile, lasso, 100);
		y_pro_lasso = test00.prognostic_index(test00.lasso_lambda, 0);
		y_pro_apml0 = test00.prognostic_index(test00.apml0_lambda, test00.apml0_k);

		for (int i = 0; i < n; i++) {
			prognostic_y << y_pro_lasso[i] << ',' << y_pro_apml0[i] << '\n';
		}

		//double *beta = cdLasso(test00.x, test00.y, test00.n1 + 1, p, lambda*pow(n,2.0));
		//double *beta = cdLasso(x, y, n, p, lambda);
		//for (int i = 0; i < p; i++)
		//{
		//	myfile << beta[i] << ",";
		//}
		//myfile << "\n";

		//delete[] beta;

		test00.delete_new();
		test00.delete_new_beta();
		test0.delete_old();
		delete[] x;
		delete[] y;
		delete[] status;
		delete[] y_pro_lasso;
		delete[] y_pro_apml0;
		time_t result = std::time(nullptr);
		myfile << asctime(localtime(&result)) << "\n";
		prognostic_y.close();
		myfile.close();
		lasso.close();
		out.close();
	}
	else {
		std::random_device rd;
		std::mt19937 g(rd());
		int *key = new int[n];
		for (int i = 0; i < n; i++)
			key[i] = i;

		int n0 = int(floor(n / 2.0) + 0.5);
		int n1 = n - n0;

		double c = 0, c_LASSO = 0;
		XY_old test0, test1;

		for (int z = 0; z < 5; z++) {
			std::shuffle(&key[0], &key[n - 1], g);
			double *y0 = new double[n0];
			double **x0 = new double*[n0];
			double *status0 = new double[n0];
			double *y1 = new double[n1];
			double **x1 = new double*[n1];
			double *status1 = new double[n1];

			for (int i = 0; i < n0; i++) {
				y0[i] = y[key[i]];
				x0[i] = x[key[i]];
				status0[i] = status[key[i]];
			}

			for (int i = n0; i < n; i++) {
				y1[i - n0] = y[key[i]];
				x1[i - n0] = x[key[i]];
				status1[i - n0] = status[key[i]];
			}

			test0.y = y0;
			test0.x = x0;
			test0.status = status0;
			test0.n = n0;
			test0.p = p;
			test1.y = y1;
			test1.x = x1;
			test1.status = status1;
			test1.n = n1;
			test1.p = p;

			XY_new test10(test0), test11(test1);
			test10.cross_validation(myfile, lasso, 100);
			double *e_beta = test10.best_beta;
			double *e_beta_LASSO = test10.best_beta_LASSO;
			c = test11.c_index(e_beta);
			c_LASSO = test11.c_index(e_beta_LASSO);
			out << c << "," << c_LASSO << '\n';

			test10.delete_new();
			test10.delete_new_beta();
			test11.delete_new();
			delete[] y0;
			delete[] x0;
			delete[] status0;
			delete[] y1;
			delete[] x1;
			delete[] status1;
		}

		time_t result = std::time(nullptr);
		myfile << asctime(localtime(&result)) << "\n";
		lasso << asctime(localtime(&result)) << "\n";
		out << asctime(localtime(&result)) << "\n";
		myfile.close();
		lasso.close();
		out.close();

		for (int i = 0; i < n; i++) {
			delete[] x[i];
		}
		delete[] y;
		delete[] status;
	}
    return 0;
}
