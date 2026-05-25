#ifndef MATRICES_H
#define MATRICES_H

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>

double getDeterminant(gsl_matrix *matrix);
gsl_matrix *getSubmatrix(gsl_matrix *matrix, int row_remove_index, int col_remove_index);

#endif
