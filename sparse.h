#pragma once

#include <iostream>
#include <vector>

class Sparse
{
public:

	Sparse();

	// creates a matrix with l lines
	Sparse(int l);

	~Sparse();

	// sets value for matrix at M[l][c]
	// please dont set a value to zero, why would you do this? What's the point of using this structure?
	void setValue(int l, int c, double val);

	// gets value for matrix at M[l][c], if its not in matrix, returns zero
	double getValue(int l, int c);

	// multiplies a vector with the sparce matrix
	// n is the size of the vector
	// vec is the vector
	void multSMV(int n, double* vec, double* w);

	//Devolve o número de linhas da matriz
	int size();

	// calcula o produto de uma matriz A (m x n) por uma matriz B (n x q),
	// armazenando o resultado na matriz C (m x q), previamente criada
	void multSMM(int q, Sparse B, Sparse *C);

	void print();
private:

	/* matrix:
	* vector for each line of the matrix
	* pair where int means the column of the matrix and the double means the value
	* 
	*/
	std::vector<std::vector<std::pair<int, double>>> matrix;
	
};