#ifndef _REGMODELS
#define _REGMODELS

// alias
typedef struct Regression *RegressionPointer;

// regression class
struct Regression
{
    int lenA;               // number of regressors
    double logmarglik;      // log marginal likelihood of regression
    int *A;                 // regressor index array
    RegressionPointer Next; // link to next regression model
};

// utils
void AddRegression(int nMaxRegs, RegressionPointer reg, int lenA, int *A, double logmarglikA, int &counter);
void SaveRegressions(const char *filename, RegressionPointer regressions);
void SaveRegressions(const char *filename, RegressionPointer regressions, const char **columnNames);
void DeleteAllRegressions(RegressionPointer regressions, int &counter);
#endif
