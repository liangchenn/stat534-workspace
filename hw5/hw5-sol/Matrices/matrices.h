#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// Include the LAPACKE header for LAPACK support
#include <lapacke.h>

int LAPACKE_dgeev(int matrix_layout, char jobvl, char jobvr, int n,
                  double* a, int lda, double* wr, double* wi, double* vl,
                  int ldvl, double* vr, int ldvr);


// Declarations for matrix utility functions
double ** allocmatrix(int n,int p);
void freematrix(int n,double** m);
void copymatrix(int n,int p,double** source,double** dest);
void readmatrix(char* filename,int n,int p,double* m[]);
void printmatrix(char* filename,int n,int p,double** m);
double** transposematrix(int n,int p,double** m);
void dotmatrixproduct(int n,int p,double** m1,double** m2,double** m);
void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m);
void inverse(int p,double** m);
double logdet(int p,double** m);

void debugmatrix(int n, int p, double** m);

// log marginal likelihood
double marglik(int n, int p, double** data, int lenA, int* A);
