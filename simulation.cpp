#define _CRT_SECURE_NO_DEPRECATE
#define TRUE_PARAMETER 15
// #define TEST_MODE

#include <iostream>
#include <fstream>
#include <string>
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

int n = N_SIM;
int p = P_SIM;
int maxit = 50;
int rep = N_REP;
int s0 = 1;
double lambda = 0.01;
char np = 'N';

int main(int argc, char *argv[])
{
#ifndef TEST_MODE
    double *beta = new double[p];
    ofstream myfile;
    ofstream lasso;

//    ORIGINAL TESTING
//    for (int i = 0; i < TRUE_PARAMETER; i++) {
//        beta[i] = pow(-1.0, i) * 2 * exp(-i / 15.0);
//    }

//    for (int i = 15; i < p; i++) {
//        beta[i] = 0;
//    }

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
            cv_path(test1, 50, pre_lasso, pre_apml0);
            test.delete_old();
            test1.delete_new();
        }

        pre_lasso.ALprint("_lasso");
        pre_apml0.ALprint("_apml0");

    delete[] beta;
    }
    else {
        string s1 = "simulation_APML0";
        string s2 = "simulation_LASSO";
        s1 = s1 + "_" + to_string(N_SIM) + "_" + to_string(P_SIM) + ".txt";
        s2 = s2 + "_" + to_string(N_SIM) + "_" + to_string(P_SIM) + ".txt";
        for (int z = 0; z < rep; z++) {
            myfile.open(s1, ios_base::app);
            lasso.open(s2, ios_base::app);
            XY_old test(n, p, z + s0, beta);
            XY_new test1(test);
            test1.cross_validation(myfile, lasso, maxit);
            test1.delete_new();
            test.delete_old();
            myfile.close();
            lasso.close();
            cout << "-----------" << z << endl;
        }

        myfile.open(s1, ios_base::app);
        lasso.open(s2, ios_base::app);
        time_t result = std::time(nullptr);
        myfile << asctime(localtime(&result)) << "\n";
        lasso << asctime(localtime(&result)) << "\n";
        myfile.close();
        lasso.close();
    }
#endif

#ifdef TEST_MODE

    double* beta = new double[p];
    for (int i = 0; i < 15; i++) {
        beta[i] = pow(-1.0, i) * 2 * std::exp(-i / 15.0);
    }
    for (int i = 15; i < p; i++) {
        beta[i] = 0;
    }

    XY_old test(n, p, 0, beta);
    XY_new test1(test);

    cdLasso(test1.x, test1.y, (test1.n1 + 1), test1.p, 0.02 * (test1.n + 1));

    test.delete_old();
    test1.delete_new();

    delete[] beta;

#endif // TEST_MODE
    return 0;
}
