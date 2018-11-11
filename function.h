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

struct XY_new
{
	double *y;
	double **x;
};

class XY_old
{
public:
	XY_old(int,int,int,double*);
	~XY_old();
	double *y;
	double **x;
	double *status;
	int n;
	int p;
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
	for (int i = 0; i < n; i++)
	{
		delete[] x[i];
	}
	delete[] x;
	delete[] y;
	delete[] status;
}


XY_new new_sample(double *log_t, double **X, double *status, int n, int p)
{
	int k = 0;
	int n1 = 0;
	XY_new output;

	for (int i = 0; i < n; i++)
	{
		if (status[i] == 1) k++;
	}

	n1 = k * (n - 1);

	output.x = new double *[n1 + 1];
	output.y = new double [n1 + 1];
	int k1 = 0;

	for (int i = 0; i < p; i++)
	{
		output.x[n1][i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		if (status[i] == 1)
		{
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				{
					output.x[k1] = new double[p];
					for (int z = 0; z < p; z++)
					{
						output.x[k1][z] = X[i][p] - X[j][p];
					}
					output.y[k1] = log_t[i] - log_t[j];
					k1++;
				}
			}
		}

		for (int j = 0; j < p; j++)
		{
			output.x[n1][j] += -output.x[i][j];
		}
	}
	output.y[n1] = 1e8;

	return output;
}

