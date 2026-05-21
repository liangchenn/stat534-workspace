#include "matrices.h"
#include <lapacke.h>

// log marginal likelihood function
double marglik(int n, int p, double **data, int lenA, int *A)
{

	int i, j;

	// 1. define variables
	double **D1 = allocmatrix(n, 1);
	double **DA = allocmatrix(n, lenA);

	// response Y
	for (i = 0; i < n; i++)
	{
		D1[i][0] = data[i][0];
	}
	// data X
	for (j = 0; j < lenA; j++)
	{
		for (int i = 0; i < n; i++)
		{
			DA[i][j] = data[i][A[j] - 1];
		}
	}

	// debugmatrix(n, 1, D1);
	// printf("\n");
	// debugmatrix(n, lenA, DA);

	// 2. compute MA = I + DA^T DA
	double **DAt = transposematrix(n, lenA, DA);
	double **MA = allocmatrix(lenA, lenA);

	matrixproduct(lenA, n, lenA, DAt, DA, MA);

	for (i = 0; i < lenA; i++)
	{
		MA[i][i] += 1.0;
	}

	double log_det_MA = logdet(lenA, MA);

	// 3. compute D1^T D1
	double d1td1 = 0.0;
	for (i = 0; i < n; i++)
	{
		d1td1 += D1[i][0] * D1[i][0];
	}

	// 4. compute D1^T DA MA^-1 DA^T D1
	double **b = allocmatrix(lenA, 1);
	matrixproduct(lenA, n, 1, DAt, D1, b);
	// get MA^-1
	inverse(lenA, MA);

	double **temp = allocmatrix(lenA, 1);
	matrixproduct(lenA, lenA, 1, MA, b, temp);

	double prod = 0.0;
	for (i = 0; i < lenA; i++)
	{
		prod += b[i][0] * temp[i][0];
	}

	// 5. calc. log marignal lik
	double prior = (lenA + 2.0) / 2.0;
	double post = (n + lenA + 2.0) / 2.0;

	double q = 1.0 + d1td1 - prod;

	double log_marlik = (lgamma(post) - lgamma(prior) - 0.5 * log_det_MA - post * log(q));

	// free mem
	freematrix(n, D1);
	freematrix(n, DA);
	freematrix(lenA, DAt);
	freematrix(lenA, MA);
	freematrix(lenA, b);
	freematrix(lenA, temp);

	return (log_marlik);
}

// allocates the memory for a matrix with
// n rows and p columns
double **allocmatrix(int n, int p)
{
	int i;
	double **m;

	m = new double *[n];
	for (i = 0; i < n; i++)
	{
		m[i] = new double[p];
		memset(m[i], 0, p * sizeof(double));
	}
	return (m);
}

// frees the memory for a matrix with n rows
void freematrix(int n, double **m)
{
	int i;

	for (i = 0; i < n; i++)
	{
		delete[] m[i];
		m[i] = NULL;
	}
	delete[] m;
	m = NULL;
	return;
}

// creates the copy of a matrix with n rows and p columns
void copymatrix(int n, int p, double **source, double **dest)
{
	int i, j;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			dest[i][j] = source[i][j];
		}
	}
	return;
}

// reads from a file a matrix with n rows and p columns
void readmatrix(char *filename, int n, int p, double *m[])
{
	int i, j;
	double s;
	FILE *in = fopen(filename, "r");

	if (NULL == in)
	{
		printf("Cannot open input file [%s]\n", filename);
		exit(1);
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < p; j++)
		{
			fscanf(in, "%lf", &s);
			m[i][j] = s;
		}
	}
	fclose(in);
	return;
}

// prints the elements of a matrix in a file
void printmatrix(char *filename, int n, int p, double **m)
{
	int i, j;
	double s;
	FILE *out = fopen(filename, "w");

	if (NULL == out)
	{
		printf("Cannot open output file [%s]\n", filename);
		exit(1);
	}
	for (i = 0; i < n; i++)
	{
		fprintf(out, "%.3lf", m[i][0]);
		for (j = 1; j < p; j++)
		{
			fprintf(out, "\t%.3lf", m[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
	return;
}

// creates the transpose of the matrix m
double **transposematrix(int n, int p, double **m)
{
	int i, j;

	double **tm = allocmatrix(p, n);

	for (i = 0; i < p; i++)
	{
		for (j = 0; j < n; j++)
		{
			tm[i][j] = m[j][i];
		}
	}

	return (tm);
}

// calculates the dot (element by element) product of two matrices m1 and m2
// with n rows and p columns; the result is saved in m
void dotmatrixproduct(int n, int p, double **m1, double **m2, double **m)
{
	int i, j;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < p; j++)
		{
			m[i][j] = m1[i][j] * m2[i][j];
		}
	}

	return;
}

// calculates the product of a nxp matrix m1 with a pxl matrix m2
// returns a nxl matrix m
void matrixproduct(int n, int p, int l, double **m1, double **m2, double **m)
{
	int i, j, k;
	double s;

	for (i = 0; i < n; i++)
	{
		for (k = 0; k < l; k++)
		{
			s = 0;
			for (j = 0; j < p; j++)
			{
				s += m1[i][j] * m2[j][k];
			}
			m[i][k] = s;
		}
	}
	return;
}

void set_mat_identity(int p, double *A)
{
	int i;

	for (i = 0; i < p * p; i++)
		A[i] = 0;
	for (i = 0; i < p; i++)
		A[i * p + i] = 1;
	return;
}

// computes the inverse of a symmetric positive definite matrix

void inverse(int p, double **m)
{
	int i, j, k, info;
	double *m_copy = (double *)malloc((p * p) * sizeof(double));
	double *m_inv = (double *)malloc((p * p) * sizeof(double));

	// Flatten the matrix 'm' to 'm_copy' in column-major order
	k = 0;
	for (j = 0; j < p; j++)
	{
		for (i = 0; i < p; i++)
		{
			m_copy[j * p + i] = m[i][j]; // Note the change to column-major
		}
	}

	set_mat_identity(p, m_inv);

	// Use LAPACKE_dposv to compute the inverse
	info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', p, p, m_copy, p, m_inv, p);

	if (info != 0)
	{
		fprintf(stderr, "Something was wrong with LAPACKE_dposv [%d]\n", info);
		exit(1);
	}

	// Un-flatten the matrix 'm_inv' back to 'm'
	k = 0;
	for (j = 0; j < p; j++)
	{
		for (i = 0; i < p; i++)
		{
			m[i][j] = m_inv[j * p + i]; // Change back to row-major
		}
	}

	free(m_copy);
	free(m_inv);
}

// computes the log of the determinant of a symmetric positive definite matrix
double logdet(int p, double **m)
{
	int i, j, info;
	double *a = (double *)malloc(p * p * sizeof(double));
	double *wr = (double *)malloc(p * sizeof(double));
	double *wi = (double *)malloc(p * sizeof(double));

	// Prepare matrix in column-major order
	for (i = 0; i < p; i++)
	{
		for (j = 0; j < p; j++)
		{
			a[j * p + i] = m[i][j];
		}
	}

	// Compute eigenvalues
	info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', p, a, p, wr, wi, NULL, 1, NULL, 1);

	if (info != 0)
	{
		printf("Error in eigenvalue computation [info = %d]\n", info);
		exit(1);
	}

	double logdet = 0.0;
	for (i = 0; i < p; i++)
	{
		logdet += log(fabs(wr[i])); // Compute log of absolute values of real parts
	}

	free(a);
	free(wr);
	free(wi);

	return logdet;
}

// debug matrix
void debugmatrix(int n, int p, double **m)
{
	int i, j;
	if (n >= 5)
	{
		n = 5;
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < p; j++)
		{
			printf("%10.5lf ", m[i][j]);
		}
		printf("\n");
	}
}