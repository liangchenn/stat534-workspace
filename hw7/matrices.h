// matrices.h
#ifndef MATRICES_H
#define MATRICES_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

int load_data(char *filename, gsl_matrix *data);
void makeCovariance(gsl_matrix *covX, gsl_matrix *X);
gsl_matrix *makeCholesky(gsl_matrix *K);

void draw_multivariate_standard_normal(gsl_rng *rng, gsl_matrix *z);
void draw_multivariate_normal(gsl_matrix *chol,
                              gsl_rng *rng,
                              gsl_matrix *z,
                              gsl_matrix *x);
void randomMVN(gsl_rng *mystream, gsl_matrix *samples, gsl_matrix *sigma);

int simulate_draws(gsl_matrix *chol, gsl_rng *rng, gsl_matrix *draws);

int save_matrix(char *filename, gsl_matrix *matrix);

#endif
