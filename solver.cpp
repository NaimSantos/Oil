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
	dados=read_input_data("input_data.txt");
	auto n=dados.nx;

	//Definicão de variáveis calculadas a partir dos dados de entrada:
	resize_if_needed(n);
	dados.Bw = dados.rhowsc/dados.rhow;
	dados.dx = dados.Lx/n;
	dados.pref = dados.p_0;

	//Exibição das informações do problema:
	print_simulation_properties();

	// Atualizar o fator volume formação:
	oil_prop(Bo_old,p);
	//Invoca o simulador:
	solucao_explicita();

	//Salva os dados:
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
	auto nx=dados.nx;
	double current_time {0.0}, res {0.0}, resp {0.0};
	int iter = 0;
	update(Bo, Bo_old);
	while (current_time < dados.tempo_max){
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
		current_time += dados.dt;
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

void solver_s(){
	auto nx=dados.nx;
	Sw[0] = Sw_old[0] + (1.0/Cww[0]) * ( Tw[1]*(p[1] - p[0]) + Tw[0]*(dados.p_W - p[0]));
	for(int i = 1; i < nx-1; i++){
		Sw[i] = Sw_old[i] + (1.0/Cww[i]) * ( Tw[i+1]*(p[i+1] - p[i]) + Tw[i]*(p[i-1] - p[i]) );
	}
	Sw[nx-1] = Sw_old[nx-1] + (1.0/Cww[nx-1]) * ( Tw[nx]*(dados.p_E - p[nx-1]) + Tw[nx-1]*(p[nx-2] - p[nx-1]));
}

DadosTotais read_input_data(std::string filename){
	DadosTotais dados_entrada;
	//Nome do arquivo de entrada
	std::ifstream input_file {filename};
	if (input_file.is_open()){
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
		//Preenche o objeto dados_entrada com as informações armazenadas no vetor Dados:
		if (Dados.size() < 20){
			std::cerr << "Arquivo de entrada nao contem todas as informacoes";
		}
		dados_entrada.phi=Dados[0];
		dados_entrada.k=Dados[1];
		dados_entrada.siw=Dados[2];
		dados_entrada.sor=Dados[3];
		dados_entrada.rhowsc=Dados[4];
		dados_entrada.rhow=Dados[5];
		dados_entrada.muw=Dados[6];
		dados_entrada.Boref=Dados[7];
		dados_entrada.co=Dados[8];
		dados_entrada.muo=Dados[9];
		dados_entrada.tempo_max=static_cast<int>(Dados[10]);
		dados_entrada.dt=Dados[11];
		dados_entrada.Lx=Dados[12];
		dados_entrada.nx=static_cast<int>(Dados[13]);
		dados_entrada.tolerancia=Dados[14];
		dados_entrada.p_0=Dados[15];
		dados_entrada.p_W=Dados[16];
		dados_entrada.p_E=Dados[17];
		dados_entrada.Sw_0=Dados[18];
		dados_entrada.Sw_W=Dados[19];
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
	//Se os vetores criados não sao suficientes, redimensionar para o novo tamanho:
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


void print_simulation_properties(){
	std::cout << "\nParametros e propriedades utilizadas na simulacao: \n";
	std::cout << "Comprimento do dominio: " << dados.Lx << " (m)\n";
	std::cout << "Numero de celulas: " << dados.nx << "\n";
	std::cout << "Dimensao da celulas: " << dados.dx << " (m)\n";
	std::cout << "Tempo total simulado: " << dados.tempo_max << " (dias)\n";
	std::cout << "Passo de tempo: " << dados.dt << " (dias)\n";
	std::cout << "Viscosidade do oleo: " << dados.muo << "\n";
	std::cout << "Compressibiliade do oleo: " << dados.co << "\n";
	std::cout << "B_ref: " << dados.Boref << "\n";
	std::cout << "Viscosidade da agua: " << dados.muw << "\n";
	std::cout << "Densidade (sc): " << dados.rhowsc << "\n";
	std::cout << "Densidade (w): " << dados.rhow << "\n";
	std::cout << "Porosidade do meio " << dados.phi << "\n";
	std::cout << "Permeabilidade do meio " << dados.k << "\n";
	std::cout << "Saturacao irredutivel (siw) " << dados.siw << "\n";
	std::cout << "Saturacao irredutivel (sor) " << dados.sor << "\n\n";
	std::cout << "Condicoes iniciais e de contorno: \n";
	std::cout << "Pressao inicial: " << dados.p_0 << "\n";
	std::cout << "Pressao no contorno esquerdo: " << dados.p_W << "\n";
	std::cout << "Pressao no contorno direito: " << dados.p_E << "\n";
	std::cout << "Saturacao inicial: " << dados.Sw_0 << "\n";
	std::cout << "Saturacao no contorno esquerdo: " << dados.Sw_W << "\n";
}