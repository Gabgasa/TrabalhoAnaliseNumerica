#include "gradconj.h"
#include "sparse.h"
#include "matriz.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void constroiMatriz(int n, Sparse *A, double* b){
	for (int i = 0; i < n; i++)
	{
		b[i] = 0;
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				A->setValue(i, j, i + 1);
				b[i] += A->getValue(i, j);
			}
			else if (j == i + 1 || j == i - 1)
			{
				A->setValue(i, j, 0.5);
				b[i] += A->getValue(i, j);
			}
			else if (j == i + 2 || j == i - 2)
			{
				A->setValue(i, j, 0.5);
				b[i] += A->getValue(i, j);
			}
			else if (j == i * 2 || i == j * 2)
			{
				A->setValue(i, j, 0.5);
				b[i] += A->getValue(i, j);
			}
			
		}
	}
}

void zeraResposta(int n, double* x){
	for(int i = 0; i<n; i++){
		x[i] = 0;
	}
}

int main(void)
{
	int n1 = 100;
	int n2 = 1000;
	int n3 = 10000;
	int n4 = 100000;
	Sparse *Mat1 = new Sparse(n1);
	Sparse *Mat2 = new Sparse(n2);
	Sparse *Mat3 = new Sparse(n3);
	Sparse *Mat4 = new Sparse(n4);

	double* b1 = criavet(n1);
	double* b2 = criavet(n2);
	double* b3 = criavet(n3);
	double* b4 = criavet(n4);

	double* x1 = criavet(n1);
	double* x2 = criavet(n2);
	double* x3 = criavet(n3);
	double* x4 = criavet(n4);
    
	zeraResposta(n1, x1);
	zeraResposta(n2, x2);
	zeraResposta(n3, x3);
	zeraResposta(n4, x4);

	constroiMatriz(n1, Mat1, b1);
	constroiMatriz(n2, Mat2, b2);
	constroiMatriz(n3, Mat3, b3);
	constroiMatriz(n4, Mat4, b4);
    
	printf("----------- N = 100 -----------\n");
    printf("Gradientes conjugados sem pre condicionador \n");
    printf("\nVezes: %d\n",GradConj(n1, *Mat1, b1, x1, pow(10, -7)));
    //mostraResposta(n1, x1);
	zeraResposta(n1, x1);
	    
    printf("Gradientes conjugados com pre condicionador \n");
	printf("\nVezes: %d\n",GradConjSparse(n1, *Mat1, b1, x1, pow(10, -7)));
	//imprimevec(n, x, "%f");
    //mostraResposta(n1, x1);

	printf("----------- N = 1000 -----------\n");
    printf("Gradientes conjugados sem pre condicionador \n");
    printf("\nVezes: %d\n",GradConj(n2, *Mat2, b2, x2, pow(10, -7)));
    //mostraResposta(n2, x2);
	zeraResposta(n2, x2);
	    
    printf("Gradientes conjugados com pre condicionador \n");
	printf("\nVezes: %d\n",GradConjSparse(n2, *Mat2, b2, x2, pow(10, -7)));
    //mostraResposta(n2, x2);

	printf("----------- N = 10000 -----------\n");
    printf("Gradientes conjugados sem pre condicionador \n");
    printf("\nVezes: %d\n",GradConj(n3, *Mat3, b3, x3, pow(10, -7)));
    //mostraResposta(n3, x3);
	zeraResposta(n3, x3);
	    
    printf("Gradientes conjugados com pre condicionador \n");
	printf("\nVezes: %d\n",GradConjSparse(n3, *Mat3, b3, x3, pow(10, -7)));

    //mostraResposta(n3, x3);

	printf("----------- N = 100000 -----------\n");
    printf("Gradientes conjugados sem pre condicionador \n");
    printf("\nVezes: %d\n",GradConj(n4, *Mat4, b4, x4, pow(10, -7)));
    //mostraResposta(n1, x1);
	zeraResposta(n4, x4);
	    
    printf("Gradientes conjugados com pre condicionador \n");
	printf("\nVezes: %d\n",GradConjSparse(n4, *Mat4, b4, x4, pow(10, -7)));
	//imprimevec(n, x, "%f");
    //mostraResposta(n4, x1);
    
    
    
    
    

	return 0;
}
