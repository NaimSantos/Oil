#include <iostream>
#include "prototypes.h"

//Váriáveis da simulação

DadosTotais dados;
//Inicialização dos vetores utilizados no processo de resolução numérica:
Vec1D W(dados.nx, 0.0);                 // n elementos, inicializados em 0
Vec1D M(dados.nx, 0.0);                 // n elementos, inicializados em 0
Vec1D E(dados.nx, 0.0);                 // n elementos, inicializados em 0
Vec1D D(dados.nx, 0.0);                 // n elementos, inicializados em 0
Vec1D RESp(dados.nx, 0.0);              // n elementos, inicializados em 0
Vec1D RESs(dados.nx, 0.0);              // n elementos, inicializados em 0
//Inicialização de vetores utilizados na simulação:
Vec1D p(dados.nx, dados.p_0);           // n elementos, inicializados com p_0
Vec1D p_old(dados.nx, dados.p_0);       // n elementos, inicializados com p_0
Vec1D p_it(dados.nx, dados.p_0);        // n elementos, inicializados com p_0
Vec1D Sw(dados.nx, dados.Sw_0);         // n elementos, inicializados com Sw_0
Vec1D Sw_old(dados.nx, dados.Sw_0);     // n elementos, inicializados com Sw_0
Vec1D Sw_it(dados.nx, dados.Sw_0);      // n elementos, inicializados com Sw_0
Vec1D Cop(dados.nx, 0.0);
Vec1D Cow(dados.nx, 0.0);
Vec1D Cww(dados.nx, 0.0);
Vec1D Kro(dados.nx, 0.0);
Vec1D Krw(dados.nx, 0.0);
Vec1D Bo(dados.nx, 0.0);
Vec1D Bo_old(dados.nx, 0.0);
Vec1D To(dados.nx+1, 0.0);               // Transmissibilidade
Vec1D Tw(dados.nx+1, 0.0);               // Transmissibilidade

int main(int argc, char *argv[] ){
	std::cout << "Inicio da simulacao" << std::endl;
	//Inicia o registro do tempo:
	CustomTimer Cronometro;
	
	//Leitura das váriaveis a partir de um arquivo de entrada:
	auto input_meio=read_input_fluid_rock("entrada_rocha_fluido.txt");
	auto input_simulacao=read_input_domain("entrada_simulacao.txt");
	auto n=dados.nx;
	//Atualização de variáveis a partir dos dados de entrada:
	resize_if_needed(n);
	dados.Bw = dados.rhowsc/dados.rhow;
	dados.dx = dados.Lx/n;
	dados.pref = dados.p_0;


	// Atualizar o fator volume formação:
	oil_prop(Bo_old,p);

	//Invoca o simulador na configuração desejada:
	if (dados.sim_type==0)
		solucao_explicita();
	/*
		if (dados.sim_type==1)
		solucao_implicita();

	*/

	save_full_data(p,"Pressao");
	save_full_data(Sw,"Saturacao");
	std::cout << "Fim da simulacao" << std::endl;
	return 0;
}
void update(Vec1D& V0, const Vec1D& V){
	const auto n = V0.size();
	for(int i = 0; i < n; i++){
		V0[i] = V[i];
	}
}

void rel_diff(Vec1D& dif, const Vec1D& V, const Vec1D& V0){
	const auto n = dif.size();
	for (int i = 0; i < n; i++){
		dif[i] = std::fabs( (V[i] - V0[i])/V0[i]);
	}
}

void abs_diff(Vec1D& dif, const Vec1D& V, double *V0){
	const auto n = dif.size();
	for (int i = 0; i < n; i++){
		dif[i] = std::fabs(V[i] - V0[i]);
	}
}
double max_vetor(Vec1D& V){
	return *std::max_element(V.begin(), V.end());
}
void acumulo(){
	const auto v = dados.phi/dados.dt;
	double bl {1.0};
	for(int i = 0; i < dados.nx; i++){
		bl = dB(Bo[i],Bo_old[i],p[i],p_old[i]);
		Cop[i] = v*(1 - Sw_old[i])*bl;
		Cow[i] = -v/Bo[i];
		Cww[i] = v/dados.Bw;
	}
}
double evaluate_dB(double bo, double bo0, double P, double P0){
	double ret {0.0};
	if(std::fabs(P - P0) < dados.tolerancia)
		ret = 0.0;
	else
		ret = ( (1/bo) - (1/bo0) )/(P - P0);
	return ret;
}
void oil_prop(Vec1D& B, const Vec1D& P){
	for(int i = 0; i < dados.nx; i++){
		B[i] = calc_Bo(P[i]);
	}
}
double calc_Bo(const double P){
	return dados.Boref/( 1.0 + dados.co*(P - dados.pref));
}
void rel_perm(){
	for(int i = 0; i < dados.nx; i++){
		Kro[i] = calc_kro(Sw[i]);
		Krw[i] = calc_krw(Sw[i]);
	}
}
// Modelo de Corey para permeabilidade relativa óleo/água
double calc_kro(double sw){
	double so = 1 - sw;
	double aux {0.0};
	double ret {0.0};
	aux = (sw - dados.siw)/(1 - dados.siw);
	ret = (1 - aux*aux) * (1 - aux)*(1 - aux);
	ret = std::max(0.0,ret);
	ret = std::min(ret,1.0);
	return ret;
}
double calc_krw(double sw){
	double aux {0.0};
	double ret {0.0};
	aux = (sw - dados.siw)/(1 - dados.siw);
	ret = std::pow(aux,4);
	ret = std::max(0.0,ret);
	ret = std::min(ret,1.0);
	return ret;
}
int upwind(double pesq, double pdir){
	return (pesq > pdir) ? 1 : 0;
}
void Transmissibilidade(){
	double conv = 8.64E-2;
	const auto Gw = (conv*dados.k)/(dados.Bw*dados.muw*dados.dx*dados.dx);
	const auto Go = (conv*dados.k)/(dados.muo*dados.dx*dados.dx);

	auto Bomed = calc_Bo(dados.p_W);
	auto ud = upwind(dados.p_W, p[0]);
	auto kroud = ud*calc_kro(dados.Sw_W) + (1-ud)*Kro[0];
	auto krwud = ud*calc_krw(dados.Sw_W) + (1-ud)*Krw[0];
	To[0] = 2*Go*(kroud/Bomed);
	Tw[0] = 2*Gw*krwud;
	for(int i = 1; i < dados.nx; i++){
		Bomed = (Bo[i-1] + Bo[i])/2;
		ud    = upwind(p[i-1],p[i]);
		kroud = ud*Kro[i-1] + (1-ud)*Kro[i];
		krwud = ud*Krw[i-1] + (1-ud)*Krw[i];
		To[i] = Go*(kroud/Bomed);
		Tw[i] = Gw*krwud;
	}
	Bomed = calc_Bo(dados.p_E);
	/// S_E = Sw[nx-1] (condicao de derivada nula)
	kroud = Kro[dados.nx-1];
	krwud = Krw[dados.nx-1];
	To[dados.nx] = 2*Go*(kroud/Bomed);
	Tw[dados.nx] = 2*Gw*krwud;
}
double dB(double bo, double bo0, double P, double P0){
	double ret {0.0};
	if(std::fabs(P - P0) < dados.tolerancia)
		ret = 0.0;
	else
		ret = ((1/bo)-(1/bo0))/(P-P0);
	return ret;
}
void solucao_explicita(){
	std::cout << "Solucao via Linearizacao Explicita" << std::endl;
	auto nx=dados.nx;
	double time {0.0}, res {0.0}, resp {0.0};
	int iter = 0;
	update(Bo, Bo_old);
	while (time < dados.tempo_max){
		res = 1.0;
		iter = 0;
		rel_perm();
		Transmissibilidade();
		while (res > dados.tolerancia){
			acumulo();
			fill_matrix();
			solver_system();
			oil_prop(Bo, p);
			rel_diff(RESp, p, p_it);
			res = max_vetor(RESp);
			update(p_it ,p);
			iter++;
		}
		solver_s();
		time += dados.dt;
		update(Bo_old,Bo);
		update(p_old,p);
		update(Sw_old,Sw);
	}
}

void fill_matrix(){
	auto nx=dados.nx;
	auto bw=dados.Bw;
	M[0] = Bo[0]*To[0] + bw*Tw[0] + Bo[0]*To[1] + bw*Tw[1] + Bo[0]*Cop[0];
	E[0] = - (Bo[0]*To[1] + bw*Tw[1]);
	D[0] = Bo[0]*Cop[0]*p_old[0] + (Bo[0]*To[0] + bw*Tw[0])*dados.p_W;

	for(int i = 1; i < dados.nx-1; i++){
		W[i] = - (Bo[i]*To[i] + bw*Tw[i]);
		E[i] = - (Bo[i]*To[i+1] + bw*Tw[i+1]);
		M[i] = Bo[i]*Cop[i] - W[i] - E[i];
		D[i] = Bo[i]*Cop[i]*p_old[i];
	}
	W[nx-1] = - (Bo[nx-1]*To[nx-1] + bw*Tw[nx-1]);
	M[nx-1] = Bo[nx-1]*To[nx-1] + bw*Tw[nx-1] + Bo[nx-1]*To[nx] + bw*Tw[nx] + Bo[nx-1]*Cop[nx-1];
	D[nx-1] = Bo[nx-1]*Cop[nx-1]*p_old[nx-1] + (Bo[nx-1]*To[nx] + bw*Tw[nx])*dados.p_E;
}
void solver_system(){
	auto nx=dados.nx;
	double omega {0.0};
	for(int i = 1; i < nx; i++){
		omega = W[i]/M[i-1];
		M[i] = M[i] - omega*E[i-1];
		D[i] = D[i] - omega*D[i-1];
	}
	p[nx-1] = D[nx-1]/M[nx-1];
	for(int i = nx-2; i >= 0; i--){
		p[i] = (D[i] - E[i]*p[i+1])/M[i];
	}
}
/// Solucao explicita da saturacao da agua
void solver_s(){
	auto nx=dados.nx;
	Sw[0] = Sw_old[0] + (1.0/Cww[0]) * ( Tw[1]*(p[1] - p[0]) + Tw[0]*(dados.p_W - p[0]));
	for(int i = 1; i < nx-1; i++){
		Sw[i] = Sw_old[i] + (1.0/Cww[i]) * ( Tw[i+1]*(p[i+1] - p[i]) + Tw[i]*(p[i-1] - p[i]) );
	}
	Sw[nx-1] = Sw_old[nx-1] + (1.0/Cww[nx-1]) * ( Tw[nx]*(dados.p_E - p[nx-1]) + Tw[nx-1]*(p[nx-2] - p[nx-1]));
}

RockFluidProperties read_input_fluid_rock(std::string filename){
	RockFluidProperties  dados_entrada;
	//Nome do arquivo de entrada
	std::ifstream input_file {filename};
	if (input_file.is_open()){
		double temp_var {0.0};
		std::vector<double> Dados;
		for (std::string linha; std::getline(input_file, linha, '\n');){
			std::istringstream current(linha);
			std::string item;
			// Leitura do stream, separados por espacos, um a um:
			while(current >> item){
				double x;
				// Isto é um double?
				if(std::istringstream(item) >> x){
					// Use x e coloque no vetor:
					Dados.push_back(x);
				}
			}
		}
		input_file.close();
		//Preenche o objeto dados_entrada com as informações armazenadas em Dados:
		if (Dados.size() < 10)
			std::cerr << "Arquivo de entrada com as com propriedade de fluido e rocha contem menos que 10 informacoes";
			dados_entrada.porosidade_rocha=Dados[0];
			dados_entrada.permeabilidade_rocha=Dados[1];
			dados_entrada.saturacao_agua=Dados[2];
			dados_entrada.saturacao_oleo=Dados[3];
			dados_entrada.densidade_superficie_sc=Dados[4];
			dados_entrada.densidade_reservatorio=Dados[5];
			dados_entrada.viscosidade_muw=Dados[6];
			dados_entrada.b_referencia=Dados[7];
			dados_entrada.compressibilidade=Dados[8];
			dados_entrada.viscosidade=Dados[9];
	}
	else{
		std::cerr << "Nao foi possivel encontrar o arquivo " << filename;
		std::exit(1);
	}
	return dados_entrada;
}
SimulationProperties read_input_domain(std::string filename){
	SimulationProperties dados_entrada;
	//Nome do arquivo de entrada
	std::ifstream input_file {filename};
	if (input_file.is_open()){
		double temp_var {0.0};
		std::vector<double> Dados;
		for (std::string linha; std::getline(input_file, linha, '\n');){
			std::istringstream current(linha);
			std::string item;
			// Leitura do stream, separados por espacos, um a um:
			while(current >> item){
				double x;
				// Isto é um double?
				if(std::istringstream(item) >> x){
					// Use x e coloque no vetor:
					Dados.push_back(x);
				}
			}
		}
		input_file.close();
		//Preenche o objeto dados_entrada com as informações armazenadas em Dados:
		if (Dados.size() < 10)
			std::cerr << "Arquivo de com as informacoes do dominio contem menos que 11 informacoes";
			dados_entrada.tempo_simulado=static_cast<int>(Dados[0]);
			dados_entrada.passo_de_tempo=Dados[1];
			dados_entrada.comprimento_x=Dados[2];
			dados_entrada.numero_de_celulas=static_cast<int>(Dados[3]);
			dados_entrada.tolerancia_tol=Dados[4];
			dados_entrada.pressao_inicial=Dados[5];
			dados_entrada.pressao_esquerda=Dados[6];
			dados_entrada.pressao_direita=Dados[7];
			dados_entrada.saturacaow_inicial=Dados[8];
			dados_entrada.saturacaow_esquerda=Dados[9];
			dados_entrada.linearizacao=static_cast<int>(Dados[10]);
	}
	else{
		std::cerr << "Nao foi possivel encontrar o arquivo " << filename;
		std::exit(1);
	}
	return dados_entrada;
}

void save_full_data(const Vec1D& Vetor, std::string variavel){
	static int a {1};
	std::string init {"output_"};
	std::string filename = init + variavel + std::to_string(a) + ".txt";
	std::fstream saver{filename, std::ios::out|std::ios::trunc};

	saver << std::setw(5) << "x (m) " << std::setw(8) << variavel << std::endl;
	const auto N = Vetor.size();
	auto dx = dados.dx;
	for (int i = 0; i < N; i++){
		saver << std::fixed << std::setprecision(3) << (dx+dx*i) << " " << Vetor[i] << std::endl;
	}
	a++;
}

void resize_if_needed(int n){
	W.resize(n,0.0);
	M.resize(n,0.0);
	E.resize(n,0.0);
	D.resize(n,0.0);
	RESp.resize(n,0.0);
	RESs.resize(n,0.0);
	p.resize(n,dados.p_0);
	p_old.resize(n,dados.p_0);
	p_it.resize(n,dados.p_0);
	Sw.resize(n,dados.Sw_0);
	Sw_old.resize(n,dados.Sw_0);
	Sw_it.resize(n,dados.Sw_0);
	Cop.resize(n,0.0);
	Cow.resize(n,0.0);
	Cww.resize(n,0.0);
	Kro.resize(n,0.0);
	Krw.resize(n,0.0);
	Bo.resize(n,0.0);
	Bo_old.resize(n,0.0);
	To.resize(n+1,0.0);
	Tw.resize(n+1,0.0);
}