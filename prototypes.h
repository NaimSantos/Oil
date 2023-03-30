#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <sstream>  //std::istringstream
#include <fstream>    // std::ifstream
#include <cmath>      // std::fabs, std::pow
#include <algorithm>  // std::min, std::max, std::max_element
#include <chrono>     // registro de tempo

//Structs usados para armazenar as propriedades do meio e do dominio da simulacao
struct RockFluidProperties{
	//Propriedades do meio
	double porosidade_rocha {0.1};
	double permeabilidade_rocha {0.1};
	double saturacao_agua {0.1};
	double saturacao_oleo {0.1};
	//Propriedades da água:
	double densidade_superficie_sc {1000};
	double densidade_reservatorio {900};
	double viscosidade_muw {0.1};
	//Propriedades do óleo:
	double b_referencia {1.0};
	double compressibilidade {1E-6};
	double viscosidade {2.0};
};
struct SimulationProperties{
	//Características do domínio:
	int tempo_simulado {100};       // dias
	double passo_de_tempo {0.1};    // dias
	double comprimento_x {1000};    // m
	int numero_de_celulas {400};    // adimensional
	double tolerancia_tol {1E-8};
	//Condições iniciais e de contorno:
	double pressao_inicial {15000};
	double pressao_esquerda {25000};
	double pressao_direita {15000};
	double saturacaow_inicial {1.0};
	double saturacaow_esquerda {0.0};
	//Metodo de linearização:
	int linearizacao {0};
};

struct DadosTotais{
	double tempo_max {100.0};  // Tempo máximo de simulação (dias)
	double dt {0.5};           // tamanho do passo de tempo (dias)
	double Lx {150};           // Dimensão do reservatório (m)
	int nx {100};              // número de células (adimensional)
	int sim_type {1};          // tipo de simulação (adimensional)
	double tolerancia {0.001}; // tolerância
	double dx {0.1};           // dimensão do discretização (m)
	double p_0 {0.1};          // pressão inicial (Pa)
	double p_W {0.1};          // pressão no contorno esquerdo (Pa)
	double p_E {0.1};          // pressão no contorno direito (Pa)
	double Sw_0 {0.1};         // Saturação inicial (adimensional)
	double Sw_W {0.1};         // Saturação no contorno esquerdo (adimensional)
	double phi {0.1};          // porosidade da rocha (adimensional)
	double k {1000};           // permeabilidade da rocha (unidade)
	double siw {0.1};          // 
	double sor {0.1};          // 
	double rhowsc {1000.0};    // massa específica (kg/m^3)
	double rhow {1000.0};      // massa específica (kg/m^3)
	double muw {1000.0};       // 
	double Boref {1000};       // Fator Volume Formação (?)
	double co {10};            // 
	double muo {50};           // 
	double pref {50};          //
	double Bw {50};            //
};

//Define o tipo "Vec1D"
using Vec1D = std::vector<double>;

//Prototipos das funções:
RockFluidProperties read_input_fluid_rock(std::string filename);
SimulationProperties read_input_domain(std::string filename);

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
void resize_vectors();

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