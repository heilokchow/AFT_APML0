#include <iomanip>
#include <random>

using namespace std;

void swap0(double &x, double &y) {
	double temp = x;
	x = y;
	y = temp;
}

int partition0(double *z, int low, int high)
{
	double pivot = z[high];
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++)
	{
		if (z[j] <= pivot)
		{
			i++;
			swap0(z[i], z[j]);
		}
	}
	swap0(z[i + 1], z[high]);

	return (i + 1);
}

void quickSort0(double *z, int low, int high)
{
	if (low < high)
	{
		int pi = partition0(z, low, high);

		quickSort0(z, low, pi - 1);
		quickSort0(z, pi + 1, high);
	}
}

struct c_range
{
	double low;
	double high;
};

class XY_old
{
public:
	XY_old(int,int,int,double*);
	~XY_old();

	void delete_old()
	{
		for (int i = 0; i < n; i++)
		{
			delete[] x[i];
		}
		delete[] x;
		delete[] y;
		delete[] status;
		cout << "Manual Delete_old is called\n";
	}

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

	void delete_new()
	{
		for (int i = 0; i < (n1 + 1); i++)
		{
			delete[] x[i];
		}
		delete[] x;
		delete[] y;
		delete[] y1;
		delete[] y2;
		cout << "Manual Delete_new is called\n";
	}

	void delete_new_cv()
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

	void seperate(int *target, int k);

	double *y;
	int *y1;
	int *y2;
	double **x;
	int n;
	int n1;
	int p;

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



c_range new_range(double *y, double low, double high, int n)
{
	double *z = new double[n];
	for (int i = 0; i < n; i++)
	{
		z[i] = y[i];
	}

	quickSort0(z, 0, (n - 1));

	int n1 = floor(low * (n - 1));
	int n2 = floor(high* (n - 1));

	c_range output;
	output.low = z[n1];
	output.high = z[n2];

	delete[] z;
	return output;
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

	for (int i = 0; i < n; i++)
	{
		x[i] = new double[p];
		y[i] = normal_dist(e);
		for (int j = 0; j < p; j++)
		{
			x[i][j] = normal_dist(e);
			y[i] += x[i][j] * beta[j];
		}
	}

	c_range c0 = new_range(y, 0.50, 0.85, n);
	uniform_real_distribution<double> uni_dist0(c0.low, c0.high);
	
	for (int i = 0; i < n; i++)
	{
		c = uni_dist0(e);
		if (y[i] > c)
		{
			y[i] = c;
			status[i] = 0;
		}
		else
		{
			status[i] = 1;
		}
	}	
}

XY_old::~XY_old()
{
	cout << "Destructor_old is called\n";
}

XY_new::XY_new(XY_old const &A)
{
	int k = 0;
	n = A.n;
	p = A.p;

	for (int i = 0; i < n; i++)
	{
		if (A.status[i] == 1) k++;
	}

	n1 = k * (n - 1);

	x = new double *[n1 + 1];
	y = new double[n1 + 1];
	y1 = new int[n1 + 1];
	y2 = new int[n1 + 1];
	int k1 = 0;

	x[n1] = new double[p];
	for (int i = 0; i < p; i++)
	{
		x[n1][i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		if (A.status[i] == 1)
		{
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				{
					x[k1] = new double[p];
					for (int z = 0; z < p; z++)
					{
						x[k1][z] = A.x[i][z] - A.x[j][z];
						x[n1][z] += -x[k1][z];
					}
					y[k1] = A.y[i] - A.y[j];
					y1[k1] = i;
					y2[k1] = j;
					k1++;
				}
			}
		}
		cout << n - i << '\n';
	}

	y1[n1] = 0;
	y2[n1] = 0;
	y[n1] = 1e6;
}

XY_new::~XY_new()
{
	cout << "Destructor_new is called\n";
}

void XY_new::seperate(int *target, int k)
{
	int s = 0;
	for (int i = 0; i < n1; i++)
	{
		int flag = 0;
		for (int j = 0; j < k; j++)
		{
			if (target[j] == y1[i] || target[j] == y2[i])
			{
				flag = 1;
			}
		}

		if (flag == 0)
		{
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
	for (int i = 0; i < n1; i++)
	{
		int flag = 0;
		for (int j = 0; j < k; j++)
		{
			if (target[j] == y1[i] || target[j] == y2[i])
			{
				flag = 1;
			}
		}

		if (flag == 0)
		{
			y_train[s1] = y[i];
			y1_train[s1] = y1[i];
			y2_train[s1] = y2[i];
			x_train[s1] = x[i];
			s1++;
		}
		else
		{
			y_test[s2] = y[i];
			y1_test[s2] = y1[i];
			y2_test[s2] = y2[i];
			x_test[s2] = x[i];
			s2++;
		}
	}

	x_test[s2] = new double[p];
	for (int i = 0; i < s2; i++)
	{
		for (int j = 0; j < p; j++)
		{
			x_test[s2][j] += -x_test[i][j];
		}
	}

	x_train[s1] = new double[p];
	for (int j = 0; j < p; j++)
	{
		x_train[s1][j] = x[n1][j] - x_test[s2][j];
	}

	y_train[s1] = 1e6;
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