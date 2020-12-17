#include "matriz.h"
#include <math.h>
#include "sparse.h"

void mostraResposta(int n, double* r) {
    printf(" [ ");
    for (int i = 0; i < n; i++) {
        printf(" %g ", r[i]);

    }
    printf(" ]\n");
}

int GradConj (int n, Sparse A, double* b, double* x, double tol)
{	
	int count;
	double aux = 0, alpha = 0, beta;
	double *d, *Ax, *rk, *r, *auxVet;
	//double **auxMat;
	r=criavet(n);
	rk=criavet(n);
	d=criavet(n);
    auxVet=criavet(n);
    Sparse *auxMat = new Sparse(n);
//	auxMat=criamat(n,n);
	Ax=criavet(n);

    A.multSMV(n, x, Ax);
	//multmv (n,n,A, x, Ax);
	for(int i=0;i<n;i++)
	{
		r[i]=b[i]-Ax[i];
		d[i]=b[i]-Ax[i];
	}

	count = 0;
	for(int k=0;k<n;k++)
	{
		if (tol > norma2(n, r))
		{
			break;
		}
		alpha = prodescalar(n, r, r);
        A.multSMV(n, d, auxVet);
		//multmv(n, n, A, d, auxVet);
		aux = prodescalar(n, d, auxVet);
		alpha =alpha/aux;
		for (int i = 0; i < n; i++)
		{
			x[i] = x[i] + alpha*d[i];
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
                auxMat->setValue(i,j, A.getValue(i,j) * alpha);
				//auxMat[i][j] = A[i][j] * alpha;
			}
		}
        auxMat->multSMV(n, d, auxVet);
		//multmv(n, n, auxMat, d, auxVet);
		for (int i = 0; i < n; i++)
		{
			rk[i] = r[i] - auxVet[i];
		}
		beta = (prodescalar(n, rk, rk)) / (prodescalar(n, r, r));
		for (int i = 0; i < n; i++)
		{
			d[i] = rk[i] + beta*d[i];
			r[i] = rk[i];
		}
		count++;
	}
    
    
	liberavet(r);
	liberavet(rk);
	liberavet(d);
	liberavet(auxVet);
	liberavet(Ax);
	//liberamat(n, auxMat);

	return count;
}

void PreCond(int n, Sparse A, Sparse *M, double w){
	
	Sparse D(n);
    Sparse Dinv(n);
	Sparse L(n);
	Sparse U(n);
	Sparse *Aux = new Sparse(n);


	//Matrizes L*w e U*w
	for(int i = 0; i<n; i++){
		D.setValue(i,i, A.getValue(i,i));
		Dinv.setValue(i,i, 1.0/A.getValue(i,i));
		U.setValue(i,i, D.getValue(i,i));
		for(int j = 0; j < i; j++){
			if(A.getValue(i,j) == 0){
				continue;
			}
			L.setValue(i,j, A.getValue(i,j) * w);
			U.setValue(j,i, A.getValue(j,i) * w);
		}
	}   

	//L.print();
	L.multSMM(n, Dinv, Aux);
	//I + wL*Dinv
	for(int i = 0; i<n; i++){
        Aux->setValue(i,i, Aux->getValue(i,i) + 1.0);
	}
	//(I + wL*Dinv)*(D + wU)
    Aux->multSMM(n,U,M);
}

void substituicoes(int n, Sparse *A, double* b, double* z, double*d) {

	double* y;
	double s;
	y = criavet(n);
	for (int i = 0; i < n; i++)
	{
		y[i] = 0;
	}

	for (int i = 0; i < n; i++) {
		s = b[i];
		for (int j = 0; j < i; j++)
		{
            s -= A->getValue(i,j) * y[j];
		}
        y[i] = s / A->getValue(i,i);
	}

	for (int i = n - 1; i >= 0; i--)
	{
		s = y[i];
		for (int j = i + 1; j < n; j++)
		{
            s -= A->getValue(j,i) * z[j];
		}
        d[i] = z[i] = s / A->getValue(i,i);
	}
}

int GradConjSparse (int n, Sparse A, double* b, double* x, double tol, double w)
{
	int count;
	double aux = 0, alpha = 0, beta;
	double *d,*z,*zk, *Ax, *rk, *r, *r_, *b_, *auxVet;
    Sparse *M = new Sparse(n);
    Sparse *auxMat = new Sparse(n);

    
	z = criavet(n);
	zk = criavet(n);
	r_ = criavet(n);
    r = criavet(n);
	rk = criavet(n);
	d = criavet(n);
	Ax = criavet(n);
	b_ = criavet(n);
	auxVet = criavet(n);
	
	A.multSMV(n, x, Ax);
	
	count = 0;
	for (int i = 0; i<n; i++)
	{
		r[i] = b[i] - Ax[i];
	}
	PreCond(n, A, M, w);
	
	//cholesky(n, M);
	//M->print();
	substituicoes(n, M, r, d, z);
	//printf("A ");
	//mostraResposta(n, r);
	

	for (int k = 0;k < n;k++)
	{
		if (tol > norma2(n, r))
		{
			break;
		}
		alpha = prodescalar(n, r, z);
		//multmv(n, n, A, d, auxVet);
		A.multSMV(n, d, auxVet);
		
		//mostraResposta(n, auxVet);
		aux = prodescalar(n, d, auxVet);
		
		alpha = alpha/aux;
		
		for (int i = 0; i < n; i++)
		{
			x[i] = x[i] + alpha*d[i];
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
                auxMat->setValue(i,j, A.getValue(i, j) * alpha);
			}
		}
        auxMat->multSMV(n, d, auxVet);
		//multmv(n, n, auxMat, d, auxVet);
		for (int i = 0; i < n; i++)
		{
			rk[i] = r[i] - auxVet[i];
		}
        M->multSMV(n, rk, zk);
		//multmv(n, n, M, rk, zk);
		substituicoes(n, M, rk, zk, zk);

		beta = prodescalar(n, rk, zk) / prodescalar(n, r, z);
		
        for (int j = 0; j < n; j++)
		{
			auxVet[j] = beta*d[j];
		}

		for (int i = 0; i < n; i++)
		{
			d[i] = zk[i] + beta*d[i];;
			r[i] = rk[i];
			z[i] = zk[i];
		}
		count++;
	}
	
	liberavet(z);
	liberavet(zk);
	liberavet(r);
	liberavet(rk);
	liberavet(d);
	liberavet(auxVet);
	liberavet(Ax);
//	liberamat(n, auxMat);
//	liberamat(n, M);
	return count;
}

