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

using namespace std;

int main()
{
    ifstream y_input("y_input_1.txt");
    ifstream status_input("status_input_1.txt");
    ifstream x_input("x_input_1.txt");
    ofstream out;
    out.open("real_out.txt");

    int n, p;

    cout << "Please Provide the Sample Size:\n";
    cin >> n;
    cout << "Please Provide number of Parameters\n";
    cin >> p;
    double lambda;
	cout << "Please Provide Lambda\n";
	cin >> lambda;

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


    XY_old test;
    test.y = y;
    test.x = x;
    test.status = status;
    test.n = n;
    test.p = p;
    XY_new test1(test);

    double *beta1 = cdLasso(test1.x, test1.y, (test1.n1 + 1), p, lambda*pow(test1.n, 2));
    for (int i = 0; i < p; i++) {
		out << beta1[i] << "\n";
	}
    out.close();
    test1.delete_new();
    for(int i = 0; i < n; i++)
    {
        delete[] x[i];
    }
    delete[] x;
    delete[] y;
    delete[] status;
    delete[] beta1;

    return 0;
}
