#include "sparse.h"

void mostraResposta(int n, double* r);
int GradConjSparse(int n, Sparse A, double* b, double* x, double tol);
int GradConj (int n, Sparse A, double* b, double* x, double tol);
//int GradConj(int n, Sparse A, double* b, double* x, double tol);
//void PreCond(int n, Sparse A, Sparse *M, double w);
//int GradConjPreCond (int n, Sparse A, double* b, double* x, double tol);
