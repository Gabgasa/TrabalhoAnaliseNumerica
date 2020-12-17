#include "sparse.h"
#include <sstream>
#include <iostream>
#include <iomanip>

Sparse::Sparse()
{
}



Sparse::Sparse(int l)
{
	matrix.resize(l);
}



Sparse::~Sparse()
{
}

std::vector<std::vector<std::pair<int, double>>> Sparse::getMatrix(){
	return matrix;
}


void Sparse::setValue(int l, int c, double val)
{
    for (int col = 0; col < matrix[l].size(); col++)
    {
		
        if (c == matrix[l][col].first)
        {
            if (val == 0)
            {
                matrix[l].erase(matrix[l].begin() + col);
                return;
            }
            matrix[l][col].second = val;
            return;
        }
    }
    if (val == 0)
    {
        return;
    }
    std::pair<int, double> aux;
    aux.first = c;
    aux.second = val;
    matrix[l].push_back(aux);
}



double Sparse::getValue(int l, int c)
{
	for (int col = 0; col < matrix[l].size(); col++)
	{
		if (c == matrix[l][col].first)
		{
			return matrix[l][col].second;
		}
	}
	return 0;
}



void Sparse::multSMV(int n, double* vec, double* w)
{
	int i, j, k;
	for (k = 0; k < n; k++)
	{
		w[k] = 0;
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < matrix[i].size(); j++)
		{
			int coluna = matrix[i][j].first;
			w[i] += matrix[i][j].second * vec[coluna];
		}
	}
}

// calcula o produto de uma matriz A (m x n) por uma matriz B (n x q),
// armazenando o resultado na matriz C (m x q), previamente criada
void Sparse::multSMM(int q, Sparse B, Sparse *C) {
	int m = matrix.size();
	int n = B.size();
	
	std::vector<std::vector<std::pair<int, double>>> matrixB = B.getMatrix();

	// for (int i = 0; i < m; i++) {
	// 	for (int k = 0; k < q; k++) {
	// 		for (int j = 0; j < n; j++) {
	// 			C->setValue(i, k, C->getValue(i, k) + getValue(i, j) * B.getValue(j, k));
	// 		}
	// 	}

	// }

	//Selecionando i linha de A
	for (int i = 0; i < m; i++) {

		//Coluna de A
		for(int j = 0; j < matrix[i].size(); j++){
			int colunaA = matrix[i][j].first;
			
			//Coluna de B
			for(int k = 0; k < matrixB[i].size(); k++){
				int colunaB = matrixB[colunaA][k].first;
				//printf("col %d ", colunaB);
				C->setValue(i, colunaB, C->getValue(i, colunaB) + getValue(i, colunaA) * B.getValue(colunaA, colunaB));

			}
		}

	}
}
int Sparse::size(){
	return matrix.size();
}



void Sparse::print()
{
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.size(); j++)
		{
			std::cout << std::setw(10) << getValue(i, j) << " ";
		}
		std::cout << std::endl;
	}
}

