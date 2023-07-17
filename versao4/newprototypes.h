#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <sstream>    //std::istringstream
#include <fstream>    // std::ifstream
#include <cmath>      // std::fabs, std::pow
#include <algorithm>  // std::min, std::max, std::max_element
#include <iomanip>    // std::setw
#include <chrono>     // registro de tempo

//Define o tipo "Vec1D"
using Vec1D = std::vector<double>;

//Struct usado para armazenar as propriedades do meio e do dominio da simulacao
struct DadosTotais{
	double tempo_max {5000.0}; // Tempo máximo de simulação (dias)
	double dt {0.25};           // tamanho do passo de tempo (dias)
	double Lx {1000};          // Dimensão do reservatório (m)
	int nx {12};                // número de células (adimensional)
	double tolerancia {1e-8};  // tolerância
	double dx {0.1};           // dimensão do discretização (m)
	double p_0 {15000};        // pressão inicial (Pa)
	double p_W {25000};        // pressão no contorno esquerdo (Pa)
	double p_E {15000};        // pressão no contorno direito (Pa)
	double Sg_0 {0.0};         // Saturação inicial (adimensional)
	double Sg_W {0.85};      // Saturação no contorno esquerdo (adimensional)
	double phi {0.50};         // porosidade da rocha (adimensional)
	double k {0.030};          // permeabilidade da rocha (unidade)
	double siw {0.15};         // saturação irredutivel agua
	double sor {0.15};         // saturação irredutivel oleo
	double rhowsc {1000.0};    // massa específica (kg/m^3)
	double rhow {990};         // massa específica (kg/m^3)
	double mug {1000.0};       // viscosidade da agua
	double Boref {1.4};        // Fator Volume Formação (?)
	double co {1e-6};          // compressibilidade
	double muo {50};           // Viscosidade do oleo
	double pref {5000};          //
	double Bg {50};            //
	double Alf{4.17};
	double Z{1.0};
	double T{10.0};
	double temat {10.0};

};
// Constantes do processo iterativo do Guass Siedel
const int MAX_ITER = 50;
const double eps = 1e-06;

//Funcoes relacionadas ao registro do tempo:
using milisegundos = std::chrono::milliseconds;
using SteadyTimePoint = std::chrono::time_point<std::chrono::steady_clock>;
class CustomTimer{
	private:
		SteadyTimePoint startpoint;
	public:
		CustomTimer(){
			startpoint = std::chrono::steady_clock::now();
		}
		~CustomTimer(){
			auto endpoint = std::chrono::steady_clock::now();
			auto start = std::chrono::time_point_cast<milisegundos>(startpoint).time_since_epoch().count();
			auto end = std::chrono::time_point_cast<milisegundos>(endpoint).time_since_epoch().count();
			auto total = end - start;
			std::cout << "Tempo decorrido: " << total << " millisegundos" << std::endl;
		}
};

//Uma forma de capturar precisamente o tempo:
SteadyTimePoint capture_time(){
	return std::chrono::steady_clock::now();;
}
//Retorna, em milissegundos, o tempo decorrido entre duas capturas de tempo:
double get_elapsed_time(SteadyTimePoint ti, SteadyTimePoint tf){
	return static_cast<double>(std::chrono::duration_cast<milisegundos> (tf - ti).count());
}


/*
//Funções:
DadosTotais read_input_data(std::string filename);
void gauss_solver(std::vector<Vec1D>& A, Vec1D& B, Vec1D& X);
bool is_diagonal_dom(const std::vector<Vec1D>& M);
void solver();
void propriedades_oleo(Vec1D& B, const Vec1D& P);
void propriedades_gas(Vec1D& B, const Vec1D& P);
double calc_Bo(const double P);
double calc_Bg(double P);
void TransmissibilidadeGas();
void TransmissibilidadeOleo();
double DerivadaPressaoCapilar(double p, double S);
double CalcularPcgo(double Sg);


void save_full_data(const Vec1D& X, std::string variavel);
void print_simulation_properties();
void print_array_2D(const std::vector<Vec1D>& M);
void print_vector(const Vec1D& M);








void resize_if_needed(int n);
void rel_perm();

void update(Vec1D& V0, const Vec1D& V);
double evaluate_dB(double bo, double bo0, double P, double P0);
void oil_prop(Vec1D& B, const Vec1D& P);
double calc_Bo(const double P);
int upwind(double pesq, double pdir);
double dB(double bo, double bo0, double P, double P0);
void divx (Vec1D& X, Vec1D& A, Vec1D& B);
void atualvec (Vec1D& x0, Vec1D& X);
double calc_Bg(double P);
double ca_kro (double S);
double ca_krg (double S);
double CalcularPcgo(double Sg);
void gas_prop(Vec1D& B, const Vec1D& P);
void criarvx0();
////////////////////////////////
void ccoefp();
void ccoefS();
void cvecDp();
void cvecDS();
void solucao_sist(std::vector<Vec1D>& A,std::vector<Vec1D>& C, Vec1D& B,Vec1D& D, Vec1D& X, Vec1D& Y);
void solucao_sist2(std::vector<Vec1D>& A,std::vector<Vec1D>& C, Vec1D& B,Vec1D& D, Vec1D& X, Vec1D& Y);
void termofonte();
void termofonte();
void calcularpo(std::vector<Vec1D>& A, Vec1D& B, Vec1D& X0, Vec1D& Z, Vec1D& X);
void calcularSg(std::vector<Vec1D>& A, Vec1D& B, Vec1D& Z0, Vec1D& X, Vec1D& Z);
void resolverSistemas(std::vector<Vec1D>& A,std::vector<Vec1D>& C, Vec1D& B,Vec1D& H, Vec1D& X, Vec1D& Z);

////////////////////////////////////////////
void funcionamento();
void TransmissibilidadeGas();
void TransmissibilidadeOleo();
double DerivadaPressaoCapilar(double p, double S);
void CalcularCgg();
void CalcularCgp();
void CalcularCog();
void CalcularCop();
double derivada_B_gas(double b,double p);
double derivada_B_oleo(double b,double p);
//area para novas funçoes:
void mmll(std::vector<Vec1D>& A,Vec1D& B, double i);
void vetpc(Vec1D& P,const Vec1D& S);
void vecdPc ();
//bool Sassenfeld(const std::vector<Vec1D>& M);
bool Sassenfeld(const std::vector<std::vector<double>>& matrix);

*/


