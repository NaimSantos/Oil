
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using Vec1D = std::vector<double>;

const int MAX_ITER = 50;
const double eps = 0.0001;

void gauss_solver(std::vector<std::vector<double>>& A, std::vector<double>& B, std::vector<double>& X);
void printarray(const std::vector<std::vector<double>>& M);
bool is_diagonal_dom(const std::vector<std::vector<double>>& M);

int main(int argc, char *argv[] ){
	
	std::vector<std::vector<double>> A {{2, 1, -1},
							            {1, 3, 2},
							            {3, -2, 7}
						};
	std::vector<double> B {3, 4, -5};
	std::vector<double> X0 {0, 0, 0};

	printarray(A);
	std::cout << is_diagonal_dom(A) << std::endl;
	//gauss_solver(A, B, X0);
	return 0;
}


void gauss_solver(std::vector<std::vector<double>>& A, std::vector<double>& B, std::vector<double>& X){
	auto m = B.size();  // numero de linhas
	auto n = m;         // numero de colunas
	
	auto Y = X;          // matriz auxiliar
	auto E = X;          // necessária para estimar o erro de uma iteração a outra
	int counter = 1;     // Contar iterações apenas pro caso da tolerancia nao ser atingida.
	//Se a matriz é diagonal dominante, a convergência é garantida
	bool teste = false;
	
	while(!teste && counter<MAX_ITER){
		teste = true;
		//std::cout << "Iteracao " << std::setprecision(10) << counter << '\n';
		for (int i = 0; i < m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j==i)
					continue;
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * X[j]);
				X[i] = Y[i]; //escreve em X a estimativa encontrada
			}
			auto res = std::fabs(((X[i] - E[i]) / X[i])) <= eps;
			teste = teste & res;
			//std::cout<< "x" << i + 1 << " = " << Y[i] << '\n';
			E[i] = X[i];
		}
		counter++;
		//std::cout << '\n';
	}
}

bool is_diagonal_dom(const std::vector<std::vector<double>>& M){
	auto m=M[0].size();
	auto n=M.size();
	if (m != n)
		std::cout << "As dimensoes nao sao compativeis" << std::endl;

	for (int i=0; i < m; i++){
		double diag = M[i][i];
		double sum = 0.0;
		for (int j = 0; j < n; j++){
			if (i==j)
				continue;
			sum = sum + std::fabs(M[i][j]);
		}
		if (sum > diag){
			return false;
		}
	}
	return true;
}

void printarray(const std::vector<std::vector<double>>& M){ 
	auto nrow=M[0].size();
	auto ncol=M.size();
	for(int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			std::cout << M[i][j] << ' '; 
		}
		std::cout << '\n'; 
	}
}


/*
	Uma implementação do método de Gauss-Siedel para resolver um sistema de equações lineares


#include <iostream>
#include <iomanip>
#include <cmath>

//#include "utilitarios.h"

void GS_Solver(double** A, const int m, const int n, double* B, const unsigned int n_eq, const float eps, double* X);

int main(){

	//Alocação: A[dim1][dim2], B[dim1], C[dim1];
	const int dim1 = 3;
	const int dim2 = 3;

	double** A{new double*[dim1] {}};
	for (int i = 0; i < dim1; ++i)
		A[i] = new double[dim2] {};

	double* B = new double [dim1];
	double* X = new double [dim1];

	//Exemplo:
	A[0][0] = 4.0; A[0][1] = 1.0; A[0][2] = -1.0;
	A[1][0] = 2.0; A[1][1] = 7.0; A[1][2] = 1.0;
	A[2][0] = 1.0; A[2][1] = -3.0; A[2][2] = 12.0;
	B[0]= 3.0; B[1]=19.0; B[2] = 31.0;
	X[0]= 0.0; X[1]=0.0; X[2]=0.0;

	const double eps = 0.00001; //define o erro permitido, usado como criterio de parada

	GS_Solver(A, dim1, dim2, B, dim1, eps, X);

	//Desalocar as matrizes:
	for (int i = 0; i < dim1; ++i)
		delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;
}

void GS_Solver(double** A, const int m, const int n, double* B, const unsigned int n_eq, const float eps, double* X){

	double* Y = new double [n_eq];	//matriz auxiliar
	double* E = new double [n_eq];	//necessária para estimar o erro de uma iteração a outra
	for (int i = 0; i < n_eq; i++)
		E[i] = X[i];
	
	unsigned int counter = 1; //Contar iterações apenas pro caso da tolerancia nao ser atingida. Se a matriz é diagonal dominante, a convergência é garantida
	bool teste = false;

	while(!teste && counter<20){
		teste = true;
		std::cout << "Iteracao " << std::setprecision(10) << counter << '\n';
		for (int i = 0; i < m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j==i)
					continue;
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * X[j]);
				X[i] = Y[i]; //escreve em X a estimativa encontrada
			}
			auto res = std::fabs(((X[i] - E[i]) / X[i])) <= eps;
			teste = teste & res;
			std::cout<< "x" << i + 1 << " = " << Y[i] << '\n';
			E[i] = X[i];
		}
		counter++;
		std::cout << '\n';
	}
	delete[] E;
	delete[] Y;
}
*/