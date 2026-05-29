#ifndef AI_REGMODELS_H
#define AI_REGMODELS_H

#include <gsl/gsl_matrix.h>

typedef struct AIRegression *LPAIRegression;

struct AIRegression
{
    int lenA;
    int *A;
    double logmarglik;
    LPAIRegression Next;
};

int ai_sameRegression(int lenA, int *A, int lenB, int *B);
int ai_countRegressions(LPAIRegression regressions);
void ai_deleteLastRegression(LPAIRegression regressions);
void ai_addRegression(int nMaxRegs, LPAIRegression regressions, int lenA, int *A, double logmarglikA);
void ai_findBestTwoPredictorModels(int nMaxRegs, gsl_matrix *data, int p, LPAIRegression regressions);

#endif
