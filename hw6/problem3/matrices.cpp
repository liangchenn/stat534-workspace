#include "matrices.h"

// log marginal likelihood w/ GSL
double marglik(gsl_matrix *data, int lenA, int *A)
{
	int i, j;
	int n = data->size1;

	// 1. define variables
	gsl_matrix *D1 = gsl_matrix_alloc(n, 1);
	gsl_matrix *DA = gsl_matrix_alloc(n, lenA);

	// response Y
	for (i = 0; i < n; i++)
	{
		gsl_matrix_set(D1, i, 0, gsl_matrix_get(data, i, 0));
	}
	// data X
	for (j = 0; j < lenA; j++)
	{
		for (i = 0; i < n; i++)
		{
			gsl_matrix_set(DA, i, j, gsl_matrix_get(data, i, A[j] - 1));
		}
	}

	// 2. compute MA = I + DA^T DA
	gsl_matrix *DAt = transposematrix(DA);
	gsl_matrix *MA = gsl_matrix_alloc(lenA, lenA);

	matrixproduct(DAt, DA, MA);

	for (i = 0; i < lenA; i++)
	{
		gsl_matrix_set(MA, i, i, gsl_matrix_get(MA, i, i) + 1.0);
	}

	double log_det_MA = logdet(MA);

	// 3. compute D1^T D1
	double d1td1 = 0.0;
	for (i = 0; i < n; i++)
	{
		d1td1 += gsl_matrix_get(D1, i, 0) * gsl_matrix_get(D1, i, 0);
	}

	// 4. compute D1^T DA MA^-1 DA^T D1
	gsl_matrix *b = gsl_matrix_alloc(lenA, 1);
	matrixproduct(DAt, D1, b);

	gsl_matrix *MA_inv = inverse(MA);
	gsl_matrix *temp = gsl_matrix_alloc(lenA, 1);
	matrixproduct(MA_inv, b, temp);

	double prod = 0.0;
	for (i = 0; i < lenA; i++)
	{
		prod += gsl_matrix_get(b, i, 0) * gsl_matrix_get(temp, i, 0);
	}

	// 5. calc. log marignal lik
	double prior = (lenA + 2.0) / 2.0;
	double post = (n + lenA + 2.0) / 2.0;
	double q = 1.0 + d1td1 - prod;

	double log_marlik = (lgamma(post) - lgamma(prior) - 0.5 * log_det_MA - post * log(q));

	gsl_matrix_free(D1);
	gsl_matrix_free(DA);
	gsl_matrix_free(DAt);
	gsl_matrix_free(MA);
	gsl_matrix_free(b);
	gsl_matrix_free(MA_inv);
	gsl_matrix_free(temp);

	return (log_marlik);
}

// prints the elements of a matrix in a file
void printmatrix(gsl_matrix *m)
{
	int i, j;
	for (i = 0; i < m->size1; i++)
	{
		printf("%.3lf", gsl_matrix_get(m, i, 0));
		for (j = 1; j < m->size2; j++)
		{
			printf("\t%.3lf",
				   gsl_matrix_get(m, i, j));
		}
		printf("\n");
	}
	return;
}

// creates the transpose of the matrix m
gsl_matrix *transposematrix(gsl_matrix *m)
{
	int i, j;

	gsl_matrix *tm = gsl_matrix_alloc(m->size2, m->size1);

	for (i = 0; i < tm->size1; i++)
	{
		for (j = 0; j < tm->size2; j++)
		{
			gsl_matrix_set(tm, i, j, gsl_matrix_get(m, j, i));
		}
	}

	return (tm);
}

// calculates the product of a nxp matrix m1 with a pxl matrix m2
// returns a nxl matrix m
void matrixproduct(gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m)
{
	int i, j, k;
	double s;

	for (i = 0; i < m->size1; i++)
	{
		for (k = 0; k < m->size2; k++)
		{
			s = 0;
			for (j = 0; j < m1->size2; j++)
			{
				s += gsl_matrix_get(m1, i, j) * gsl_matrix_get(m2, j, k);
			}
			gsl_matrix_set(m, i, k, s);
		}
	}
	return;
}

// computes the inverse of a positive definite matrix
// the function returns a new matrix which contains the inverse
// the matrix that gets inverted is not modified
gsl_matrix *inverse(gsl_matrix *K)
{
	int j;

	gsl_matrix *copyK = gsl_matrix_alloc(K->size1, K->size1);
	if (0 != gsl_matrix_memcpy(copyK, K))
	{
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}

	gsl_matrix *inverse = gsl_matrix_alloc(K->size1, K->size1);
	gsl_permutation *myperm = gsl_permutation_alloc(K->size1);

	if (0 != gsl_linalg_LU_decomp(copyK, myperm, &j))
	{
		printf("GSL failed LU decomposition.\n");
		exit(1);
	}
	if (0 != gsl_linalg_LU_invert(copyK, myperm, inverse))
	{
		printf("GSL failed matrix inversion.\n");
		exit(1);
	}
	gsl_permutation_free(myperm);
	gsl_matrix_free(copyK);

	return (inverse);
}

// creates a submatrix of matrix M
// the indices of the rows and columns to be selected are
// specified in the last four arguments of this function
gsl_matrix *MakeSubmatrix(gsl_matrix *M,
						  int *IndRow, int lenIndRow,
						  int *IndColumn, int lenIndColumn)
{
	int i, j;
	gsl_matrix *subM = gsl_matrix_alloc(lenIndRow, lenIndColumn);

	for (i = 0; i < lenIndRow; i++)
	{
		for (j = 0; j < lenIndColumn; j++)
		{
			gsl_matrix_set(subM, i, j,
						   gsl_matrix_get(M, IndRow[i], IndColumn[j]));
		}
	}

	return (subM);
}

// computes the log of the determinant of a symmetric positive definite matrix
double logdet(gsl_matrix *K)
{
	int i;

	gsl_matrix *CopyOfK = gsl_matrix_alloc(K->size1, K->size2);
	gsl_matrix_memcpy(CopyOfK, K);
	gsl_permutation *myperm = gsl_permutation_alloc(K->size1);
	if (0 != gsl_linalg_LU_decomp(CopyOfK, myperm, &i))
	{
		printf("GSL failed LU decomposition.\n");
		exit(1);
	}
	double logdet = gsl_linalg_LU_lndet(CopyOfK);
	gsl_permutation_free(myperm);
	gsl_matrix_free(CopyOfK);
	return (logdet);
}
