#ifndef AI_MATRICES_H
#define AI_MATRICES_H

#include <gsl/gsl_matrix.h>

double ai_getDeterminant(gsl_matrix *matrix);
gsl_matrix *ai_getSubmatrix(gsl_matrix *matrix, int removeRow, int removeCol);

#endif
