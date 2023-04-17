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

//Struct usado para armazenar as propriedades do meio e do dominio da simulacao
struct DadosTotais{
	double tempo_max {5000.0}; // Tempo máximo de simulação (dias)
	double dt {0.5};           // tamanho do passo de tempo (dias)
	double Lx {1000};          // Dimensão do reservatório (m)
	int nx {400};              // número de células (adimensional)
	double tolerancia {1e-8};  // tolerância
	double dx {0.1};           // dimensão do discretização (m)
	double p_0 {15000};        // pressão inicial (Pa)
	double p_W {25000};        // pressão no contorno esquerdo (Pa)
	double p_E {15000};        // pressão no contorno direito (Pa)
	double Sw_0 {0.2};         // Saturação inicial (adimensional)
	double Sw_W {1.0};         // Saturação no contorno esquerdo (adimensional)
	double phi {0.20};         // porosidade da rocha (adimensional)
	double k {0.030};          // permeabilidade da rocha (unidade)
	double siw {0.15};         // saturação irredutivel agua
	double sor {0.10};         // saturação irredutivel oleo
	double rhowsc {1000.0};    // massa específica (kg/m^3)
	double rhow {990};         // massa específica (kg/m^3)
	double muw {1000.0};       // viscosidade da agua
	double Boref {1.4};        // Fator Volume Formação (?)
	double co {1e-6};          // compressibilidade
	double muo {50};           // Viscosidade do oleo
	double pref {50};          //
	double Bw {50};            //
};

//Define o tipo "Vec1D"
using Vec1D = std::vector<double>;

//Atualiza um vetor V0 com os valores the um outro vetor V
void update(Vec1D& V0, const Vec1D& V);

//Atualiza um vetor "diff" com a diferença relativa entre os pares de dois vetores V e V0
void rel_diff(Vec1D& dif, const Vec1D& V, const Vec1D& V0);

//Atualiza um vetor "diff" com a diferença absoluta entre os pares de dois vetores V e V0
void abs_diff(Vec1D& dif, const Vec1D& V, double *V0);
double max_vetor(Vec1D& V);
void acumulo();
double evaluate_dB(double bo, double bo0, double P, double P0);

//Preenche um vetor "B" com os valores do fator volume formação para cada pressão em P
void oil_prop(Vec1D& B, const Vec1D& P);
double calc_Bo(const double P);
void rel_perm();

// Modelo de Corey para permeabilidade relativa óleo/água
double calc_kro(double sw);
double calc_krw(double sw);
int upwind(double pesq, double pdir);
void Transmissibilidade();
double dB(double bo, double bo0, double P, double P0);
void print_field_x(Vec1D&, double DX, int n);
void solucao_explicita();
void simulador_exp();
void solucao_implicita();
void fill_matrix();
void solver_system();
void resize_if_needed(int n);
void save_full_data(const Vec1D& X, std::string variavel);
void solver_s();
void print_simulation_properties();
DadosTotais read_input_data(std::string filename);
//Registro do tempo:
using milisegundos = std::chrono::milliseconds;
using SteadyTimePoint = std::chrono::time_point<std::chrono::steady_clock>;

class CustomTimer
{
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

			std::cout << "\nTempo decorrido: " << total << " millisegundos" << std::endl;
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