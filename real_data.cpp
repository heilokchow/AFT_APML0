#define _CRT_SECURE_NO_DEPRECATE
#define TRUE_PARAMETER 15
#define REPLI 1

//#define CINDEX
#define FOLDS 4
#define TRAIN_TEST 1
//#define MSE
//#define MULTI_SIMULATION
//#define TEST_MODE
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

using namespace std;

#ifndef ILLUSTRATION
int n = 235;
int p = 7399;
char sr = 'N';
char c_sim = 'N';

#else
int n = 100;
int p = 200;
char sr = 'N';
char c_sim = 'Y';

#ifndef CINDEX
#define CINDEX
#endif // CINDEX

#ifdef FOLDS
#undef FOLDS
#endif // FOLDS
#define FOLDS 4

#ifdef TRAIN_TEST
#undef TRAIN_TEST
#endif // TRAIN_TEST
#define TRAIN_TEST 0.667

#ifndef MSE
#define MSE
#endif // MSE

#ifdef MULTI_SIMULATION
#undef MULTI_SIMULATION
#endif // MULTI_SIMULATION

#undef REPLI
#define REPLI 1

#ifdef TEST_MODE
#undef TEST_MODE
#endif // TEST_MODE
#endif // ILLUSTRATION

#include "cdlasso.h"
#include "function.h"

int main()
{
    ofstream out;
    ofstream myfile;
    ofstream lasso;

    out.open("real_out.txt", ios_base::app);
    myfile.open("real_random_cv.txt", ios_base::app);
    lasso.open("real_random_cv_LASSO.txt", ios_base::app);

    double *y = new double[n];
    double *status = new double[n];
    double **x = new double*[n];

    for(int i = 0; i < n; i++) {
        x[i] = new double[p];
    }

    if (c_sim != 'Y') {
        ifstream y_input("y_input_1.txt");
        ifstream status_input("status_input_1.txt");
        ifstream x_input("x_input_1.txt");
        string y1, x1, status1;

        for(int i = 0; i < n; i++) {
            getline(y_input, y1, '\n');
            getline(status_input, status1, '\n');
            y[i] = atof(y1.c_str());
            status[i] = atof(status1.c_str());

            for(int j = 0; j < (p - 1); j++) {
                getline(x_input, x1, ',');
                x[i][j] = atof(x1.c_str());
            }
            getline(x_input, x1, '\n');
            x[i][p - 1] = atof(x1.c_str());
        }
    }
    else {
        double *beta = new double[p];
        for (int i = 0; i < 15; i++) {
            beta[i] = pow(-1.0, i)* 2 * std::exp(-i / 15.0);
        }
        for (int i = 15; i < p; i++) {
            beta[i] = 0;
        }
        XY_old test(n, p, 0, beta);
        for (int i = 0; i < n; i++) {
            y[i] = test.y[i];
            status[i] = test.status[i];
            for (int j = 0; j < p; j++) {
                x[i][j] = test.x[i][j];
            }
        }
        delete[] beta;
        test.delete_old();
    }

#ifndef TEST_MODE
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
        y_pro_lasso = test00.prognostic_index(test00.get_lasso_lambda(), 0);
        y_pro_apml0 = test00.prognostic_index(test00.get_apml0_lambda(), test00.get_apml0_k());

        for (int i = 0; i < n; i++) {
            prognostic_y << y_pro_lasso[i] << ',' << y_pro_apml0[i] << '\n';
        }

//        double *beta = cdLasso(test00.x, test00.y, test00.n1 + 1, p, lambda*pow(n,2.0));
//        double *beta = cdLasso(x, y, n, p, lambda);
//        for (int i = 0; i < p; i++)
//        {
//            myfile << beta[i] << ",";
//        }
//        myfile << "\n";
//
//        delete[] beta;

        test00.delete_new();
        test00.delete_new_beta();
        test0.delete_old();
        delete[] x;
        delete[] y;
        delete[] status;
        delete[] y_pro_lasso;
        delete[] y_pro_apml0;
//        time_t result = std::time(nullptr);
//        myfile << asctime(localtime(&result)) << "\n";
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

#ifdef TRAIN_TEST
        int n0 = int(floor(n * TRAIN_TEST) + 0.5);
#else
        int n0 = int(floor(n / 2.0) + 0.5);
#endif // TRAIN_TEST
        int n1 = n - n0;
        ALPath pre_lasso(p);
        ALPath pre_apml0(p);
        XY_old test0, test1;

#ifdef MULTI_SIMULATION
        double *beta = new double[p];
        for (int i = 0; i < 15; i++) {
            beta[i] = pow(-1.0, i)* 2 * std::exp(-i / 15.0);
        }
        for (int i = 15; i < p; i++) {
            beta[i] = 0;
        }
        for (int i = 0; i < 20; i++) {
            XY_old test(n, p, i + 1, beta);
            for (int i = 0; i < n; i++) {
                y[i] = test.y[i];
                status[i] = test.status[i];
                for (int j = 0; j < p; j++) {
                    x[i][j] = test.x[i][j];
                }
            }
            test.delete_old();
#endif // MULTI_SIMULATION

        for (int z = 0; z < REPLI; z++) {
            std::shuffle(&key[0], &key[n - 1], g);
            double *y0 = new double[n0];
            double **x0 = new double*[n0];
            double *status0 = new double[n0];
            double *y1 = new double[n1];
            double **x1 = new double*[n1];
            double *status1 = new double[n1];

            if (TRAIN_TEST != 1) {

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
                test1.y = y1;
                test1.x = x1;
                test1.status = status1;
                test1.n = n1;
                test1.p = p;
            }
            else {
                for (int i = 0; i < n0; i++) {
                    y0[i] = y[i];
                    x0[i] = x[i];
                    status0[i] = status[i];
                }
            }

            test0.y = y0;
            test0.x = x0;
            test0.status = status0;
            test0.n = n0;
            test0.p = p;

            XY_new test10(test0);
#ifdef CINDEX
            if (TRAIN_TEST != 1) {
                XY_new test11(test1);
                cv_path(test10, test11, 50, pre_lasso, pre_apml0);
                test11.delete_new();
            }
#else
            cv_path(test10, 50, pre_lasso, pre_apml0);
#endif // CINDEX
            std::cout << "-------------------z1:" << z << '\n';
            test10.delete_new();
            std::cout << "-------------------z2:" << z << '\n';
            delete[] y0;
            delete[] x0;
            delete[] status0;
            delete[] y1;
            delete[] x1;
            delete[] status1;
            pre_lasso.ALprint("_lasso");
            pre_apml0.ALprint("_apml0");
        }
#ifdef MULTI_SIMULATION
        }
        delete[] beta;
#endif // MULTI_SIMULATION

//        time_t result = std::time(nullptr);
//        myfile << asctime(localtime(&result)) << "\n";
//        lasso << asctime(localtime(&result)) << "\n";
//        out << asctime(localtime(&result)) << "\n";
        myfile.close();
        lasso.close();
        out.close();

        for (int i = 0; i < n; i++) {
            delete[] x[i];
        }
        delete[] y;
        delete[] status;
        delete[] x;
        delete[] key;
    }
#endif // TEST_MODE

#ifdef TEST_MODE
    XY_old test0;
    test0.x = x;
    test0.y = y;
    test0.status = status;
    test0.n = n;
    test0.p = p;
    XY_new test00(test0);

    clock_t begin = clock();
    double* beta = cdLasso(test00.x, test00.y, test00.n1 + 1, p, 0.000001);
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << elapsed_secs << 's' << std::endl;

    test00.delete_new();
    test0.delete_old();
    delete[] beta;
#endif // TEST_MODE
    return 0;
}



