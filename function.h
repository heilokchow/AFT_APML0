#include <iomanip>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <pthread.h>
#include <memory>
#include <string>
#include <cstring>

#define USE_IDENTICAL 0
#define THREAD_VALUE 6

using namespace std;

struct c_range
{
	double low;
	double high;
};

class XY_old
{
public:
	XY_old();
	XY_old(int,int,int,double*);
	~XY_old();
	void delete_old();
	void print_old();

	double *y;
	double **x;
	double *status;
	int n;
	int p;
};

class XY_new
{
public:
	XY_new();
	XY_new(XY_old const &A);
	~XY_new();
	void delete_new();
	void delete_new_beta();
	void cross_validation(ofstream &myfile, ofstream &lasso, int maxit);
	double *get_beta();
	double *get_beta_LASSO();
	double get_lasso_lambda() const;
	double get_apml0_lambda() const;
	int get_apml0_k() const;
	double c_index(double *e_beta) const;
	double LG_all(double *beta, double lambda) const;
	double *prognostic_index(double lambda, int k);
	friend void *cv_lasso_run(void *arg);

	double *y;
	int *y1;
	int *y2;
	int *status2;
	double **x;
	int n;
	int n1;
	int p;

private:
	double *y_ori;
	double **x_ori;
	int *status_ori;

	double* best_beta;
	double* best_beta_LASSO;
	double lasso_lambda;
	double apml0_lambda;
	int apml0_k;

	void delete_new_cv();
	void seperate(int *target, int k);
	double LG(double *beta);

	double *y_train;
	int *y1_train;
	int *y2_train;
	double **x_train;
	int n_train;
	int n1_train;

	double *y_test;
	int *y1_test;
	int *y2_test;
	double **x_test;
	int n_test;
	int n1_test;
};

struct lasso_path
{
	double lambda;
	XY_new model;
	double *beta = NULL;
	int k;
	~lasso_path()
	{
//      delete[] beta; // Not used
        cout << "Thread_beta deleted" << endl;
	}
};

struct cv_lasso_path
{
	double lambda;
	XY_new model;
	double *beta = NULL;
	double **beta_all = NULL;
	double *LG;
	int *k;
	int sum;
	~cv_lasso_path()
	{
//      delete[] beta; // Not used
        for (int i = 0; i < sum; i++) {
            delete[] beta_all[i];
        }
        delete[] beta_all;
        cout << "Thread_beta_cv deleted" << endl;
	}
};

struct ALPath
{
    int path[43] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                       20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,
                       180,200,250,300,350,400,450,500};
    double lambda[43] = {0.0};
    double* beta[43];
    double cv[43] = {0.0};
    int ct[43] = {0};
    int p;

#ifdef TRUE_PARAMETER
    double TP[43] = {0.0};
    double FP[43] = {0.0};
#endif // TRUE_PARAMETER

#ifdef CINDEX
    double c_index[43] = {0.0};
#ifdef MSE
    double mse_LG[43] ={0.0};
#endif // MSE
#endif // CINDEX

    ALPath(int model_p)
    {
        p = model_p;
        for (int i = 0; i < 43; i++) {
            beta[i] = new double[p];
            for (int j = 0; j < p; j++) {
                beta[i][j] = 0;
            }
        }
    }
    ~ALPath()
    {
        for (int i = 0; i < 43; i++) {
            delete[] beta[i];
        }
        std::cout << "ALPath deleted\n";
    }

    void ALprint(std::string b)
    {
        std::ofstream AL_path;
        std::string x = "AL_path" + b + ".txt";
        char cstr[x.size() + 1];
        strcpy(cstr, x.c_str());
        AL_path.open(cstr);

        for (int i = 0; i < 43; i++) {
            AL_path << path[i] << ',' << cv[i] << ',' << lambda[i] << ',';

#ifdef TRUE_PARAMETER
            AL_path << TP[i] << ',' << FP[i] << ',';
#endif // TRUE_PARAMETER

#ifdef CINDEX
            AL_path << c_index[i] << ',';
#ifdef MSE
            AL_path << mse_LG[i] << ',';
#endif // MSE
#endif

            for (int j = 0; j < p; j++) {
                AL_path << beta[i][j] << ',';
            }
            AL_path << '\n';
        }
        AL_path.close();
    }

    void ALadd(int* k, double* lam, int sum, double** beta_path, double* LG_path);
    void ALadd(int* k, double* lam, int sum, double** beta_path, double* LG_path, const XY_new& test_model);

};

void ALPath::ALadd(int* k, double* lam, int sum, double** beta_path, double* LG_path
#ifdef CINDEX
                   , const XY_new& test_model
#endif // CINDEX
)
{
    int ini = 0;
    for (int i = 0; i < 43; i++) {
        if (k[ini] <= path[i]) {
            int flag = 0;
            for (int j = ini; j < sum; j++) {
                if(k[j] > path[i]) {
                    j--;
                    ini = j;
                    flag++;
                    ct[i]++;
                    lambda[i] = (lam[j] + lambda[i] * (ct[i] - 1))/ ct[i];
                    cv[i] = (LG_path[j] + cv[i] * (ct[i] - 1))/ ct[i];
#ifdef TRUE_PARAMETER
                    int TR = 0;
                    for (int z = 0; z < TRUE_PARAMETER; z++) {
                        if (std::abs(beta_path[j][z]) > 1e-6)
                            TR++;
                    }
                    TP[i] = (TR + TP[i] * (ct[i] - 1))/ ct[i];
                    FP[i] = (k[j] - TR + FP[i] * (ct[i] - 1))/ ct[i];
#endif // TRUE_PARAMETER

                    for (int z = 0; z < p; z++) {
                        beta[i][z] = (beta_path[j][z] + beta[i][z] * (ct[i] - 1))/ ct[i];
                    }
#ifdef CINDEX
                    double c = test_model.c_index(beta_path[j]);
                    c_index[i] = (c + c_index[i] * (ct[i] - 1))/ ct[i];
#ifdef MSE
                    double ms = test_model.LG_all(beta_path[j], 0)/(2.0 * test_model.n);
                    mse_LG[i] = (ms + mse_LG[i] * (ct[i] - 1))/ ct[i];
#endif // MSE
#endif // CINDEX
                    break;
                }
            }

            if (flag == 0) {
                int j = sum - 1;
                ini = j;
                ct[i]++;
                lambda[i] = (lam[j] + lambda[i] * (ct[i] - 1))/ ct[i];
                cv[i] = (LG_path[j] + cv[i] * (ct[i] - 1))/ ct[i];
#ifdef TRUE_PARAMETER
                int TR = 0;
                for (int z = 0; z < TRUE_PARAMETER; z++) {
                    if (std::abs(beta_path[j][z]) > 1e-6)
                        TR++;
                }
                TP[i] = (TR + TP[i] * (ct[i] - 1))/ ct[i];
                FP[i] = (k[j] - TR + FP[i] * (ct[i] - 1))/ ct[i];
#endif // TRUE_PARAMETER
                for (int z = 0; z < p; z++) {
                    beta[i][z] = (beta_path[j][z] + beta[i][z] * (ct[i] - 1))/ ct[i];
                }
#ifdef CINDEX
                double c = test_model.c_index(beta_path[j]);
                c_index[i] = (c + c_index[i] * (ct[i] - 1))/ ct[i];
#ifdef MSE
                double ms = test_model.LG_all(beta_path[j], 0)/(2.0 * test_model.n);
                mse_LG[i] = (ms + mse_LG[i] * (ct[i] - 1))/ ct[i];
#endif // MSE
#endif // CINDEX
            }
        }
    }
}

double **cholesky(double **A, int n) {
	double **L = new double *[n];
	for (int i = 0; i < n; i++) {
		L[i] = new double[n];
		for (int j = 0; j < n; j++)
			L[i][j] = 0;
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < (i + 1); j++) {
			double s = 0;
			for (int k = 0; k < j; k++)
				s += L[i][k] * L[j][k];
			L[i][j] = (i == j) ?
				sqrt(A[i][i] - s) :
				(1.0 / L[j][j] * (A[i][j] - s));
		}

	return L;
}

void swap0(double &x, double &y)
{
	double temp = x;
	x = y;
	y = temp;
}

int partition0(double *z, int low, int high)
{
	double pivot = z[high];
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++) {
		if (z[j] <= pivot) {
			i++;
			swap0(z[i], z[j]);
		}
	}
	swap0(z[i + 1], z[high]);

	return (i + 1);
}

void quickSort0(double *z, int low, int high)
{
	if (low < high) {
		int pi = partition0(z, low, high);

		quickSort0(z, low, pi - 1);
		quickSort0(z, pi + 1, high);
	}
}

c_range new_range(double *y, double low, double high, int n)
{
	double *z = new double[n];
	for (int i = 0; i < n; i++) {
		z[i] = y[i];
	}

	quickSort0(z, 0, (n - 1));

	int n1 = int (floor(low * (n - 1)) + 0.5);
	int n2 = int (floor(high* (n - 1)) + 0.5);

	c_range output;
	output.low = z[n1];
	output.high = z[n2];

	delete[] z;
	return output;
}

XY_old::XY_old()
{
	cout << "XY_old Constructed\n";
}

double XY_new::get_lasso_lambda() const
{
    return(this->lasso_lambda);
}

double XY_new::get_apml0_lambda() const
{
    return(this->apml0_lambda);
}

int XY_new::get_apml0_k() const
{
    return(this->apml0_k);
}

XY_new::XY_new()
{
	cout << "XY_new Constructed\n";
}

XY_old::XY_old(int n0, int p0, int seed0, double *beta)
{
	n = n0;
	p = p0;
	y = new double[n];
	x = new double*[n];
	status = new double[n];
	double c = 0;

	seed_seq seed{seed0};
	mt19937 e(seed);
	normal_distribution<double> normal_dist(0.0, 1.0);
	uniform_real_distribution<double> uni_dist(0.0, 1.0);

#if USE_IDENTICAL == 1
	for (int i = 0; i < n; i++) {
		x[i] = new double[p];
		y[i] = normal_dist(e);
		for (int j = 0; j < p; j++) {
			x[i][j] = normal_dist(e);
			y[i] += x[i][j] * beta[j];
		}
	}
#else
	double **CV = new double *[50];
	for (int i = 0; i < 50; i++) {
		CV[i] = new double[50];
		for (int j = 0; j < 50; j++) {
			if (i == j)
				CV[i][j] = 1.0;
			else
				CV[i][j] = 0.5;
		}
	}

	double **KCV = cholesky(CV, 50);
	double *xx = new double[p];
	for (int i = 0; i < n; i++) {
		x[i] = new double[p];
		y[i] = normal_dist(e);

		for (int j = 0; j < p; j++) {
			xx[j] = normal_dist(e);
			x[i][j] = 0;
		}

		for (int j = 0; j < 50; j++) {
			for (int k = 0; k < 50; k++)
				x[i][j] += KCV[j][k] * xx[k];
		}

		for (int j = 50; j < p; j++)
			x[i][j] = xx[j];

		for (int j = 0; j < p; j++)
			y[i] += x[i][j] * beta[j];
	}

	delete[] xx;
	for (int i = 0; i < 50; i++) {
		delete[] KCV[i];
		delete[] CV[i];
	}
	delete[] KCV;
	delete[] CV;
#endif // USE_IDENTICAL == 1

	//OLD CENSORE DETERMINATION METHOD

	//c_range c0 = new_range(y, 0.50, 0.85, n);
	//uniform_real_distribution<double> uni_dist0(c0.low, c0.high);

	for (int i = 0; i < n; i++) {
		c = uni_dist(e);
		if (0.5 < c) {
			//y[i] = c;
			status[i] = 0;
		}
		else {
			status[i] = 1;
		}
	}
}

void XY_old::print_old()
{
	std::ofstream check;
	check.open("check_file.txt");

	for (int i = 0; i < n; i++) {
		check << y[i] << "," << status[i];
		for (int j = 0; j < p; j++) {
			check << "," << x[i][j];
		}
		check << '\n';
	}
	check.close();
}

void XY_old::delete_old()
{
	for (int i = 0; i < n; i++) {
		delete[] x[i];
	}
	delete[] x;
	delete[] y;
	delete[] status;
	cout << "Manual Delete_old is called\n";
}

XY_old::~XY_old() {
	cout << "Destructor_old is called\n";
}

XY_new::XY_new(XY_old const &A)
{
	int k = 0;
	n = A.n;
	p = A.p;

	y_ori = new double[n];
	status_ori = new int[n];
	x_ori = new double*[n];
	for (int i = 0; i < n; i++) {
		x_ori[i] = new double[p];
		y_ori[i] = A.y[i];
		status_ori[i] = A.status[i];
		for (int j = 0; j < p; j++) {
			x_ori[i][j] = A.x[i][j];
		}
	}

	for (int i = 0; i < n; i++) {
		if (A.status[i] == 1) k++;
	}

	n1 = k * (n - 1);

	x = new double *[n1 + 1];
	y = new double[n1 + 1];
	y1 = new int[n1 + 1];
	y2 = new int[n1 + 1];
	status2 = new int[n1 + 1];
	int k1 = 0;

	x[n1] = new double[p];
	for (int i = 0; i < p; i++) {
		x[n1][i] = 0;
	}

	for (int i = 0; i < n; i++) {
		if (A.status[i] == 1) {
			for (int j = 0; j < n; j++) {
				if (j != i) {
					x[k1] = new double[p];
					for (int z = 0; z < p; z++) {
						x[k1][z] = A.x[i][z] - A.x[j][z];
						x[n1][z] += -x[k1][z];
					}
					y[k1] = A.y[i] - A.y[j];
					y1[k1] = i;
					y2[k1] = j;
					status2[k1] = A.status[j];
					k1++;
				}
			}
		}
	}
	cout << "New Sample Formulated\n";

	y1[n1] = 0;
	y2[n1] = 0;
	status2[n1] = 0;
	y[n1] = 1e8;
}

XY_new::~XY_new() {
	cout << "Destructor_new is called\n";
}

double XY_new::LG_all(double *beta, double lambda) const
{
	double s = 0, r = 0;
	for (int i = 0; i < n1; i++) {
		r = y[i];
		for (int j = 0; j < p; j++) {
			r += -x[i][j] * beta[j];
		}
		if (r < 0)
		{
			s += -2 * r;
		}
	}

	for (int i = 0; i < p; i++)
	{
		s += lambda * abs(beta[i]);
	}
	return s;
}

void XY_new::seperate(int *target, int k) {
	int s = 0;
	for (int i = 0; i < n1; i++) {
		int flag = 0;
		for (int j = 0; j < k; j++) {
			if (target[j] == y1[i] || target[j] == y2[i]) {
				flag = 1;
			}
		}

		if (flag == 0) {
			s++;
		}
	}

	y_train = new double[s + 1];
	y_test = new double[n1 - s + 1];
	y1_train = new int[s + 1];
	y1_test = new int[n1 - s + 1];
	y2_train = new int[s + 1];
	y2_test = new int[n1 - s + 1];

	x_train = new double *[s + 1];
	x_test = new double *[n1 - s + 1];

	int s1 = 0, s2 = 0;
	for (int i = 0; i < n1; i++) {
		int flag = 0;
		for (int j = 0; j < k; j++) {
			if (target[j] == y1[i] || target[j] == y2[i]) {
				flag = 1;
			}
		}

		if (flag == 0) {
			y_train[s1] = y[i];
			y1_train[s1] = y1[i];
			y2_train[s1] = y2[i];
			x_train[s1] = x[i];
			s1++;
		}
		else {
			y_test[s2] = y[i];
			y1_test[s2] = y1[i];
			y2_test[s2] = y2[i];
			x_test[s2] = x[i];
			s2++;
		}

	}

	x_test[s2] = new double[p];
	for (int j = 0; j < p; j++) {
		x_test[s2][j] = 0;
	}

	for (int i = 0; i < s2; i++) {
		for (int j = 0; j < p; j++) {
			x_test[s2][j] += -x_test[i][j];
		}
	}

	x_train[s1] = new double[p];
	for (int j = 0; j < p; j++) {
		x_train[s1][j] = x[n1][j] - x_test[s2][j];
	}

	y_train[s1] = 1e8;
	y_test[s2] = 0;
	y1_train[s1] = 0;
	y1_test[s2] = 0;
	y2_train[s1] = 0;
	y2_test[s2] = 0;

	n_test = k;
	n1_test = s2;
	n_train = n - k;
	n1_train = s1;
}

void XY_new::delete_new()
{
	for (int i = 0; i < (n1 + 1); i++) {
		delete[] x[i];
	}

	for (int i = 0; i < n; i++)
	{
		delete[] x_ori[i];
	}

	delete[] x;
	delete[] y;
	delete[] y1;
	delete[] y2;
	delete[] status2;
	delete[] x_ori;
	delete[] y_ori;
	delete[] status_ori;
	cout << "Manual Delete_new is called\n";
}

void XY_new::delete_new_beta()
{
	delete[] best_beta;
	delete[] best_beta_LASSO;
}

void XY_new::delete_new_cv()
{
	delete[] x_train[n1_train];
	delete[] x_train;
	delete[] y_train;
	delete[] y1_train;
	delete[] y2_train;

	delete[] x_test[n1_test];
	delete[] x_test;
	delete[] y_test;
	delete[] y1_test;
	delete[] y2_test;
	cout << "Manual Delete_new_cv is called\n";
}

double XY_new::LG(double *beta)
{
	double s = 0, r = 0;
	for (int i = 0; i < n1; i++) {
		r = y[i];
		for (int j = 0; j < p; j++) {
			r += -x[i][j]*beta[j];
		}
		if (r < 0)
		{
			s += -r;
		}
	}

	double s_train = 0;
	for (int i = 0; i < n1_train; i++) {
		r = y_train[i];
		for (int j = 0; j < p; j++) {
			r += -x_train[i][j]*beta[j];
		}
		if (r < 0)
		{
			s_train += -r;
		}
	}

	double ss = s/n - s_train/n_train;
	return ss;
}

int find_rank(double *abs_beta, double v, int low, int high)
{
	int rank;
	if (low == high) {
		rank = low;
		return rank;
	}

	if (high - low == 1) {
		if (v > abs_beta[high]) {
			rank = low;
			return rank;
		}
		else {
			rank = high;
			return high;
		}
	}

	int k = int (floor((low + high)/2.0) + 0.5);

	if (v > abs_beta[k]) {
		high = k - 1;
		rank = find_rank(abs_beta, v, low, high);
	}
	else {
		low = k;
		rank = find_rank(abs_beta, v, low, high);
	}
	return rank;
}

int *beta_rank(double *beta, int p)
{
	int *rank = new int[p];
	double *beta_copy = new double[p];

	for (int i = 0; i < p; i++) {
		beta_copy[i] = abs(beta[i]);
	}

	quickSort0(beta_copy, 0, p - 1);

	int i = 0, j = p - 1;
	while (i < j) {
		swap0(beta_copy[i], beta_copy[j]);
		i++;
		j--;
	}

	for (int i = 0; i < p; i++) {
		rank[i] = find_rank(beta_copy, abs(beta[i]), 0, p - 1);
	}

	delete[] beta_copy;
	return rank;
}

int number_nzero(double *beta, int p)
{
	int s = 0;
	for (int i = 0; i < p; i++) {
		if (abs(beta[i]) > 1e-6) {
			s++;
		}
	}
	return s;
}

void XY_new::cross_validation(ofstream &myfile, ofstream &lasso, int maxit)
{
	double lambda[7];
	lambda[0] = 0.02;
	lambda[1] = 0.06;
	lambda[2] = 0.11;
	lambda[3] = 0.18;
	lambda[4] = 0.25;
	lambda[5] = 0.30;
	lambda[6] = 0.35;

	double *stage2_beta = new double[p];
	int fold = int (floor(n / 10) + 0.5);

	double min_LG = 1e8, min_lambda = 1;
	double min_LG_LASSO = 1e8, min_lamda_LASSO = 1;
	double *LG1 = new double[maxit];
	int min_k = maxit;

	for (int m = 0; m < 7; m++) {
		for (int i = 0; i < maxit; i++) {
			LG1[i] = 0;
		}

		double LG_LASSO = 0;

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

			seperate(target, d);
			double *rep_beta = cdLasso(x_train, y_train, (n1_train + 1),
										p, lambda[m] * (n_train + 1));
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
				LG1[j] += LG(stage2_beta);
			}
			LG_LASSO += LG(rep_beta);

			delete[] rank;
			delete[] rep_beta;
			delete[] target;
			delete_new_cv();
		}

		for (int j = 0; j < maxit; j++) {
			if (min_LG > LG1[j]) {
				min_LG = LG1[j];
				min_lambda = lambda[m];
				min_k = j + 1;
			}
		}

		if (min_LG_LASSO > LG_LASSO) {
			min_LG_LASSO = LG_LASSO;
			min_lamda_LASSO = lambda[m];
		}
		cout << m << "\n";
	}
	delete[] LG1;

	apml0_k = min_k;
	lasso_lambda = min_lamda_LASSO;
	apml0_lambda = min_lambda;

	myfile << min_LG << "," << min_lambda << "," << min_k << ",";
	lasso << min_LG_LASSO << "," << min_lamda_LASSO << ",";
	best_beta = cdLasso(x, y, (n1 + 1), p, min_lambda*(n1 + 1));
	best_beta_LASSO = cdLasso(x, y, (n1 + 1), p, min_lamda_LASSO*(n1 + 1));
	int *rank = beta_rank(best_beta, p);
	int nlasso = number_nzero(best_beta_LASSO, p);
	lasso << nlasso << ",";

	for (int j = 0; j < p; j++) {
		if (rank[j] < min_k) {
			stage2_beta[j] = best_beta[j];
			myfile << stage2_beta[j] << ",";
		}
		else {
			stage2_beta[j] = 0;
			myfile << stage2_beta[j] << ",";
		}
		lasso << best_beta_LASSO[j] << ",";
	}

	best_beta = stage2_beta;
	myfile << "\n";
	lasso << "\n";
	delete[] rank;
}

double*XY_new::prognostic_index(double lambda, int k)
{
	double *stage2_beta = new double[p];
	double *y_out = new double[n];
	int k1 = 0;

	if (k == 0) {
		k1 = p;
	}
	else
	{
		k1 = k;
	}

	int fold = int(floor(n / 10) + 0.5);
	std::random_device rd;
	std::mt19937 g(rd());
	int *key = new int[n];
	for (int i = 0; i < n; i++) {
		key[i] = i;
	}
	std::shuffle(&key[0], &key[n - 1], g);

	for (int i = 0; i < 10; i++) {
		int low = fold * i;
		int high = fold * (i + 1) - 1;
		if (i == 9) {
			high = n - 1;
		}

		int d = high - low + 1;
		int *target = new int[d];
		for (int j = 0; j < d; j++) {
			target[j] = key[low + j];
		}
		seperate(target, d);
		double *rep_beta = cdLasso(x_train, y_train, (n1_train + 1),
			p, lambda * (n_train + 1));
		int *rank = beta_rank(rep_beta, p);

		for (int j_ = 0; j_ < p; j_++) {
			if (rank[j_] <= k1) {
				stage2_beta[j_] = rep_beta[j_];
			}
			else {
				stage2_beta[j_] = 0;
			}
		}

		for (int j = 0; j < d; j++) {
			y_out[target[j]] = 0;
			for (int j_ = 0; j_ < p; j_++) {
				y_out[target[j]] += x_ori[target[j]][j_] * stage2_beta[j_];
			}
		}
		delete[] rank;
		delete[] rep_beta;
		delete[] target;
		delete_new_cv();
	}

	delete[] key;
	delete[] stage2_beta;
	return y_out;
}

double*XY_new::get_beta()
{
	return best_beta;
}

double*XY_new::get_beta_LASSO()
{
	return best_beta_LASSO;
}

double XY_new::c_index(double *e_beta) const
{
	int c0 = 0;
	double c2 = 0;
	double y0;
	for (int i = 0; i < n1; i++) {
		y0 = 0;

		for (int j = 0; j < p; j++)
			y0 += x[i][j] * e_beta[j];

		if (y[i] < 0) {
			if (y0 < 0)
				c0++;
			c2++;
		}
	}

	double c1 = c0 / c2;
	return c1;
}

void *lasso_run(void *arg)
{
	lasso_path *m = (lasso_path*)arg;
	m->k = 0;
	m->beta = cdLasso(m->model.x, m->model.y, (m->model.n1 + 1), m->model.p, m->lambda*(m->model.n1 + 1));
	for (int i = 0; i < m->model.p; i++) {
		if (abs(m->beta[i]) > 1e-6)
			m->k++;
	}
	cout << "beta: " << m->beta[0] << endl;
	pthread_exit(m);
}

void *cv_lasso_run(void *arg)
{
#ifdef FOLDS
    int folds = FOLDS;
#else
    int folds = 10;
#endif
    cv_lasso_path *m = (cv_lasso_path*)arg;

    m->LG = new double[m->sum];
    for (int i = 0; i < m->sum; i++) {
        m->LG[i] = 0;
    }
    int fold = int (floor(m->model.n / folds) + 0.5);

    double *stage2_beta = new double[m->model.p];
    for (int i = 0; i < folds; i++) {
        int low = fold * i;
        int high = fold * (i + 1) - 1;
        if (i == (folds - 1)) {
            high = m->model.n - 1;
        }

        int d = high - low + 1;
        int *target = new int[d];
        for (int j = 0; j < d; j++) {
            target[j] = low + j;
        }

        m->model.seperate(target, d);
        double *rep_beta = cdLasso(m->model.x_train, m->model.y_train, (m->model.n1_train + 1),
										m->model.p, m->lambda* (m->model.n_train + 1));
			int *rank = beta_rank(rep_beta, m->model.p);
			for (int j = 0; j < m->sum; j++) {
				for (int j_ = 0; j_ < m->model.p; j_++) {
					if (rank[j_] < m->k[j]) {
						stage2_beta[j_] = rep_beta[j_];
					}
					else {
						stage2_beta[j_] = 0;
					}
				}
				m->LG[j] += m->model.LG(stage2_beta);
			}

			delete[] rank;
			delete[] rep_beta;
			delete[] target;
			m->model.delete_new_cv();
		}

        int *rank = beta_rank(m->beta, m->model.p);
        m->beta_all = new double*[m->sum];
        for (int j = 0; j < m->sum; j++) {
            m->beta_all[j] = new double[m->model.p];
            for (int j_ = 0; j_ < m->model.p; j_++) {
                if (rank[j_] < m->k[j]) {
                    m->beta_all[j][j_] = m->beta[j_];
                }
                else {
                    m->beta_all[j][j_] = 0;
                }
            }
        }
        delete[] rank;
        pthread_exit(m);
}

void cv_path(const XY_new &new_class,
#ifdef CINDEX
const XY_new &test_class,
#endif
int k, ALPath& pre_lasso, ALPath& pre_apml0)
{
	double lambda[50];
	for (int i = 0; i < 50; i++)
		lambda[49 - i] = (i + 1)*0.01;

	lasso_path *path = new lasso_path[k];
	for (int i = 0; i < k; i++) {
		path[i].lambda = lambda[i];
		path[i].model = new_class;
		cout << path[i].lambda << endl;
	}

	pthread_t tid[k];
	pthread_attr_t attr[k];

    int tv0 = (int) ceil(k/(double)THREAD_VALUE);
    lasso_path* temp;
    for (int j = 0; j < tv0; j++) {
        for (int i = j*THREAD_VALUE; i < std::min((j + 1)*THREAD_VALUE, k); i++) {
            cout << path[i].lambda << endl;
            pthread_attr_init(&attr[i]);
            pthread_create(&tid[i], &attr[i], lasso_run, &path[i]);
        }

        for (int i = j*THREAD_VALUE; i < std::min((j + 1)*THREAD_VALUE, k); i++) {
            temp = &path[i];
            pthread_join(tid[i], (void**) &temp);
        }
    }

//	for (int i = 0; i < k; i++)
//	{
//		cout << "Thread_beta: " << sizeof(path[i]) << endl;
//		delete[] path[i].beta;
//	}

//  C++14
//  std::unique_ptr<int[]> k_size = std::make_unique<int[]>(k);
//  std::unique_ptr<int[]> k_flag = std::make_unique<int[]>(k);

//  C++11
    std::unique_ptr<int[]> k_size(new int[k]);
    std::unique_ptr<int[]> k_flag(new int[k]);

    for (int i = 0; i < k; i++) {
        k_size[i] = number_nzero(path[i].beta, new_class.p);
        std::cout << k_size[i] << std::endl;
    }

    int sum = 0;
    int c = 0;
    for (int i = 0; i < k - 1; i++) {
        c = 0;
        for (int j = i + 1; j < k; j++) {
            if (k_size[j] == k_size[i])
            {
                c++;
                break;
            }
        }
        if (c == 0) {
            sum++;
            k_flag[i] = 1;
        }
    }
    sum++;
    k_flag[k - 1] = 1;

//  C++14
//  std::unique_ptr<int[]> k_sum = std::make_unique<int[]>(sum);
//  std::unique_ptr<double[]> lam_sum = std::make_unique<double[]>(sum);
//  std::unique_ptr<double*[]> beta_sum = std::make_unique<double*[]>(sum);

//  C++11
    std::unique_ptr<int[]> k_sum(new int[sum]);
    std::unique_ptr<double[]> lam_sum(new double[sum]);
    std::unique_ptr<double*[]> beta_sum(new double*[sum]);

    int q = 0;
    for (int i = 0; i < k; i++) {
        if (k_flag[i] == 1) {
            k_sum[q] = k_size[i];
            lam_sum[q] = lambda[i];
            beta_sum[q] = path[i].beta;
            std::cout << "k: " << k_size[i] << " lambda: " << lambda[i] << " beta: " << path[i].beta[0] << std::endl;
            std::cout << "k: " << k_sum[q] << " lambda: " << lam_sum[q] << " beta: " << beta_sum[q][0] << std::endl;
            q++;
        }
        else {
            delete[] path[i].beta;
        }
    }

	delete[] path;
    std::cout << "-------------------4\n";
//  CROSS VALIDDATION FOR APML0 UNDER DIFFERENT PENALTY LEVEL
    cv_lasso_path *cv_path = new cv_lasso_path[sum];
    int *k_sum_ = new int[sum];
    for (int i = 0; i < sum; i++) {
        k_sum_[i] = k_sum[i];
    }

	for (int i = 0; i < sum; i++) {
		cv_path[i].lambda = lam_sum[i];
		cv_path[i].beta = beta_sum[i];
        cv_path[i].k = k_sum_;
		cv_path[i].model = new_class;
		cv_path[i].sum = sum;
		std::cout << "cv_lambda: " << cv_path[i].lambda << std::endl;
	}

	pthread_t tid1[sum];
	pthread_attr_t attr1[sum];

	int tv = ceil(sum/(double)THREAD_VALUE);
    cv_lasso_path* temp1;
    std::cout << "sum: " << sum << "T_V: " << THREAD_VALUE << "tv: " << tv << std::endl;
	for (int j = 0; j < tv; j++) {
        for (int i = j*THREAD_VALUE; i < std::min((j + 1)*THREAD_VALUE, sum); i++) {
            std::cout << "cv_lambda: " << cv_path[i].lambda << std::endl;
            pthread_attr_init(&attr1[i]);
            pthread_create(&tid1[i], &attr1[i], cv_lasso_run, &cv_path[i]);
        }

        for (int i = j*THREAD_VALUE; i < std::min((j + 1)*THREAD_VALUE, sum); i++) {
            temp1 = &cv_path[i];
            pthread_join(tid1[i], (void**) &temp1);
        }
	}

    for (int i = 0; i < sum; i++) {
        for (int j = 0; j < sum; j++) {
            std::cout << cv_path[i].LG[j] << ' ';
//            std::cout << i << "------" << j << std::endl;
        }
//        std::cout << std::endl;
    }
    std::cout << sum << std::endl;
//  OUTPUT RESULT TO FILE
//    std::ofstream wp_lasso;
//    std::ofstream wp_apml0;

//    wp_lasso.open("wp_lasso.txt");
//    wp_apml0.open("wp_apml0.txt");

    double** lasso_t = new double*[sum];
    double** apml0_t = new double*[sum];
    double* lasso_LG = new double[sum];
    double* apml0_LG = new double[sum];
    std::cout << "-------------------0\n";
    for (int i = 0; i < sum; i++) {
        lasso_t[i] = new double[new_class.p];
        apml0_t[i] = new double[new_class.p];
        lasso_LG[i] = cv_path[i].LG[i];
        apml0_LG[i] = cv_path[i].LG[i];

        for (int z = 0; z < new_class.p; z++) {
            lasso_t[i][z] = cv_path[i].beta_all[i][z];
            apml0_t[i][z] = cv_path[i].beta_all[i][z];
        }

        for (int j = i; j < sum; j++) {
            if(cv_path[j].LG[i] < apml0_LG[i]) {
                apml0_LG[i] = cv_path[j].LG[i];
                for(int z = 0; z < new_class.p; z++) {
                    apml0_t[i][z] = cv_path[j].beta_all[i][z];
                }
            }
        }
    }

    int* k_sum__ = new int[sum];
    double* lam_sum__ = new double[sum];

    for (int i = 0; i < sum; i++) {
        k_sum__[i] = k_sum[i];
        lam_sum__[i] = lam_sum[i];
    }
    std::cout << "-------------------1\n";
#ifdef CINDEX
    pre_lasso.ALadd(k_sum__, lam_sum__, sum, lasso_t, lasso_LG, test_class);
    pre_apml0.ALadd(k_sum__, lam_sum__, sum, apml0_t, apml0_LG, test_class);
#else
    pre_lasso.ALadd(k_sum__, lam_sum__, sum, lasso_t, lasso_LG);
    pre_apml0.ALadd(k_sum__, lam_sum__, sum, apml0_t, apml0_LG);
#endif // CINDEX
    std::cout << "-------------------2\n";
//    for (int i = 0; i < sum; i++) {
//        wp_lasso << lasso_LG[i] << ',' << k_sum[i] << ',' << lam_sum[i] << ',';
//        wp_apml0<< apml0_LG[i] << ',' << k_sum[i] << ',' << lam_sum[i] << ',';
//        for (int j = 0; j < new_class.p; j++) {
//            wp_lasso << lasso_t[i][j] << ',';
//            wp_apml0<< apml0_t[i][j] << ',';
//        }
//        wp_lasso << std::endl;
//        wp_apml0 << std::endl;
//    }
    std::cout << "-------------------3\n";
    for (int i = 0; i < sum; i++) {
        delete[] lasso_t[i];
        delete[] apml0_t[i];
        delete[] beta_sum[i];
    }
    std::cout << "-------------------3.5\n";
    delete[] k_sum_;
    std::cout << "-------------------4\n";
    delete[] k_sum__;
    delete[] lam_sum__;
    delete[] lasso_t;
    delete[] apml0_t;
    delete[] lasso_LG;
    delete[] apml0_LG;
    std::cout << "-------------------5\n";
    delete[] cv_path;
//    wp_lasso.close();
//    wp_apml0.close();
    std::cout << "-------------------6\n";
}
