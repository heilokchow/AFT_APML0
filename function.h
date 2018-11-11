using namespace std;

struct XY_new
{
	double *y;
	double **x;
};

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

	output.x = new double *[n1];
	output.y = new double *[n1];
	int k1 = 0;
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
				}
			}
		}
		k1++;
	}
}
