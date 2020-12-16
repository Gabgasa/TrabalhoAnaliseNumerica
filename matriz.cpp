#include "matriz.h"
#include <math.h>
#include <stdlib.h>


double* criavet(int n) {
	double* v = (double*)malloc(n * sizeof(double));

	return v;
};


void liberavet(double* v) {
	free(v);
};


double prodescalar(int n, const double* v, const double* w) {
	double escalar = 0;
	for (int i = 0; i < n; i++) {
		escalar += v[i] * w[i];
	}

	return escalar;
};

double norma2(int n, const double* v) {
	double norma_2 = 0;
	for (int i = 0; i < n; i++) {
		norma_2 += v[i] * v[i];
	}
	norma_2 = sqrt(norma_2);

	return norma_2;
}

void multvs(int n, const double* v, double s, double* w) {
	for (int i = 0; i < n; i++)
	{
		w[i] = v[i] * s;
	}
}

double** criamat(int m, int n) {
	double** mat = (double**)malloc(m * sizeof(double*));

	for (int i = 0; i < m; i++) {
		mat[i] = (double*)malloc(n * sizeof(double));
	}

	return mat;
}

double** criamattri(int n) {
	double** mat = (double**)malloc(n * sizeof(double*));

	for (int i = 0; i < n; i++) {
		mat[i] = (double*)malloc((i + 1) * sizeof(double*));
	}

	return mat;
}

void liberamat(int m, double** A) {

	for (int i = 0; i < m; i++)
	{
		free(A[i]);
	}

	free(A);
}

// preenche a matriz transposta de A em T, previamente criada;
// A tem dimens�o m x n; T tem dimens�o n x m
void transposta(int m, int n, const double** A, double** T) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			T[i][j] = A[j][i];
		}
	}
}

// calcula o produto de uma matriz A (m x n) por um vetor v (m),
// resultando no vetor w (m), previamente criado
void multmv(int m, int n, const double** A, const double* v, double* w) {

	for (int i = 0; i < m; i++) {
		w[i] = 0;
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			w[i] += v[j] * A[i][j];
		}
	}
}

// calcula o produto de uma matriz A (m x n) por uma matriz B (n x q),
// armazenando o resultado na matriz C (m x q), previamente criada
void multmm(int m, int n, int q, const double** A, const double** B, double** C) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < q; j++) {
			C[i][j] = 0;
		}
	}

	for (int i = 0; i < m; i++) {
		for (int k = 0; k < q; k++) {
			for (int j = 0; j < n; j++) {
				C[i][k] += A[i][j] * B[j][k];
			}
		}

	}
}








