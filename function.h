#include <iomanip>
#include <random>
#include <cstdlib>

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

	double *y;
	double **x;
	double *status;
	int n;
	int p;
};

class XY_new 
{
public:
	XY_new(XY_old const &A);
	~XY_new();
	void delete_new();
	void delete_new_beta();
	void cross_validation(ofstream &myfile, ofstream &lasso, int maxit);
	double *get_beta();
	double *get_beta_LASSO();
	double c_index(double *e_beta);
	double LG_all(double *beta, double lambda);

	double *y;
	int *y1;
	int *y2;
	int *status2;
	double **x;
	int n;
	int n1;
	int p;

	double* best_beta;
	double* best_beta_LASSO;

private:
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

	for (int i = 0; i < n; i++) {
		x[i] = new double[p];
		y[i] = normal_dist(e);
		for (int j = 0; j < p; j++) {
			x[i][j] = normal_dist(e);
			y[i] += x[i][j] * beta[j];
		}
	}

	c_range c0 = new_range(y, 0.50, 0.85, n);
	uniform_real_distribution<double> uni_dist0(c0.low, c0.high);
	
	for (int i = 0; i < n; i++) {
		c = uni_dist0(e);
		if (y[i] > c) {
			y[i] = c;
			status[i] = 0;
		}
		else {
			status[i] = 1;
		}
	}	
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

double XY_new::LG_all(double *beta, double lambda)
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
	delete[] x;
	delete[] y;
	delete[] y1;
	delete[] y2;
	delete[] status2;
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

double*XY_new::get_beta()
{
	return best_beta;
}

double*XY_new::get_beta_LASSO()
{
	return best_beta_LASSO;
}

double XY_new::c_index(double *e_beta)
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