#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

void printmatrix(char *filename, gsl_matrix *m);
gsl_matrix *transposematrix(gsl_matrix *m);
void matrixproduct(gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m);
gsl_matrix *inverse(gsl_matrix *K);
gsl_matrix *MakeSubmatrix(gsl_matrix *M,
						  int *IndRow, int lenIndRow,
						  int *IndColumn, int lenIndColumn);
double logdet(gsl_matrix *K);

// log marginal lik func.
double marglik(gsl_matrix *data, int lenA, int *A);
