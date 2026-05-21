#define marglik reference_marglik
#include "../Matrices/matrices.cpp"
#undef marglik

#include "problem1_ai_marglik.cpp"

int main(int argc, char *argv[])
{
    int n = 158; // sample size
    int p = 51;  // number of variables
    int A[] = {2, 5, 10};
    int lenA = 3;
    int i;

    char default_datafilename[] = "../Matrices/erdata.txt";
    char *datafilename = default_datafilename;
    if (argc > 1)
    {
        datafilename = argv[1];
    }

    double **data = allocmatrix(n, p);
    readmatrix(datafilename, n, p, data);

    printf("Using data file: %s\n", datafilename);
    printf("AI marginal likelihood of regression [1|%d", A[0]);
    for (i = 1; i < lenA; i++)
    {
        printf(",%d", A[i]);
    }
    printf("] = %.5lf\n", marglik(n, p, data, lenA, A));

    freematrix(n, data);
    return 0;
}
