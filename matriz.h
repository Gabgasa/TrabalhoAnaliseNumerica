#ifndef MATRIZ_H
#define MATRIZ_H

// cria (aloca) um vetor de dimens�o n, retornando seu ponteiro
double* criavet(int n);

// libera (a mem�ria) de um vetor previamente criado
void liberavet(double* v);

// calcula e retorna o produto escalar entre dois vetores de dimens�o n
double prodescalar(int n, const double* v, const double* w);

// calcula e retorna a norma-2 de um vetor de dimens�o n
double norma2(int n, const double* v);

// calcula a produto de um vetor v por um escalar s;
// o resultado deve ser armazenado no vetor w, previamente criado
void multvs(int n, const double* v, double s, double* w);


// cria (aloca) uma matriz de dimens�o m x n, retornando seu ponteiro;
// a matriz � representado por vetor de vetores linha
double** criamat(int m, int n);

// cria (aloca) uma matriz triangular inferior de dimens�o n x n;
// a matriz � representado por vetor de vetores linha:
// o primeiro vetor linha tem dimens�o 1, o segundo 2, e assim por diante
double** criamattri(int n);

// libera (a mem�ria) de uma matriz previamente criada
void liberamat(int m, double** A);

// preenche a matriz transposta de A em T, previamente criada;
// A tem dimens�o m x n; T tem dimens�o n x m
void transposta(int m, int n, const double** A, double** T);

// calcula o produto de uma matriz A (m x n) por um vetor v (m),
// resultando no vetor w (m), previamente criado
void multmv(int m, int n, const double** A, const double* v, double* w);

// calcula o produto de uma matriz A (m x n) por uma matriz B (n x q),
// armazenando o resultado na matriz C (m x q), previamente criada
void multmm(int m, int n, int q, const double** A, const double** B, double** C);

#endif
