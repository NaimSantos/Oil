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
Vec1D Sg(dados.nx, dados.Sg_0);         // n elementos, inicializados com Sg_0
Vec1D Sg_old(dados.nx, dados.Sg_0);     // n elementos, inicializados com Sg_0
Vec1D Sg_it(dados.nx, dados.Sg_0);      // n elementos, inicializados com Sg_0
Vec1D Cop(dados.nx, 0.0);
Vec1D Cog(dados.nx, 0.0);
Vec1D Cgg(dados.nx, 0.0);
Vec1D Cgp(dados.nx, 0.0);
Vec1D Kro(dados.nx, 0.0);
Vec1D Krg(dados.nx, 0.0);
Vec1D Bo(dados.nx, 0.0);
Vec1D Bo_old(dados.nx, 0.0);
Vec1D To(dados.nx+1, 0.0);               // Transmissibilidade
Vec1D Tg(dados.nx+1, 0.0);               // Transmissibilidade
Vec1D Pc(dados.nx+1,0.0);
Vec1D X0(2*dados.nx,0.0);

std::vector<Vec1D> coef(2*dados.nx, Vec1D(2*dados.nx, 0.0));

int main(int argc, char *argv[] ){
	std::cout << "Inicio da simulacao" << std::endl;
	//Inicia o registro do tempo:
	CustomTimer Cronometro;
	
	//Leitura das váriaveis a partir de um arquivo de entrada:
	dados=read_input_data("input_data.txt");
	auto n=dados.nx;

	//Definicão de variáveis calculadas a partir dos dados de entrada:
	resize_if_needed(n);
	dados.Bg = dados.rhowsc/dados.rhow;
	dados.dx = dados.Lx/n;
	dados.pref = dados.p_0;

	//Exibição das informações do problema:
	print_simulation_properties();

	//std::cout << "Antes do preenchimento: \n";
	//print_array_2D(coef);
	//ccoef();
	//std::cout << "Depois do preenchimento: \n";
	//print_array_2D(coef);


	// Atualizar o fator volume formação:
	oil_prop(Bo_old,p);
	//Invoca o simulador:
	//solucao_explicita();

	//Salva os dados:
	//save_full_data(p,"Pressao");
	//save_full_data(Sg,"Saturacao");
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
		Cop[i] = v*(1 - Sg_old[i])*bl;
		Cog[i] = -v/Bo[i];
		Cgg[i] = v/dados.Bg;
		Cgp[i] = v*Sg_old[i]*bl;
		
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
//Vamos usar
double calc_Bo(const double P){
	return dados.Boref/( 1.0 + dados.co*(P - dados.pref));
}
void rel_perm(){
	for(int i = 0; i < dados.nx; i++){
		Kro[i] = calc_kro(Sg[i]);
		Krg[i] = calc_krg(Sg[i]);
	}
}
// Modelo de Corey para permeabilidade relativa óleo/água
double calc_kro(double Sg){
	double so = 1 - Sg;
	double aux {0.0};
	double ret {0.0};
	aux = (Sg - dados.siw)/(1 - dados.siw);
	ret = (1 - aux*aux) * (1 - aux)*(1 - aux);
	ret = std::max(0.0,ret);
	ret = std::min(ret,1.0);
	return ret;
}
double calc_krg(double Sg){
	double aux {0.0};
	double ret {0.0};
	aux = (Sg - dados.siw)/(1 - dados.siw);
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
	const auto Gw = (conv*dados.k)/(dados.Bg*dados.mug*dados.dx*dados.dx);
	const auto Go = (conv*dados.k)/(dados.muo*dados.dx*dados.dx);

	auto Bomed = calc_Bo(dados.p_W);
	auto ud = upwind(dados.p_W, p[0]);
	auto kroud = ud*calc_kro(dados.Sg_W) + (1-ud)*Kro[0];
	auto krgud = ud*calc_krg(dados.Sg_W) + (1-ud)*Krg[0];
	To[0] = 2*Go*(kroud/Bomed);
	Tg[0] = 2*Gw*krgud;
	for(int i = 1; i < dados.nx; i++){
		Bomed = (Bo[i-1] + Bo[i])/2;
		ud    = upwind(p[i-1],p[i]);
		kroud = ud*Kro[i-1] + (1-ud)*Kro[i];
		krgud = ud*Krg[i-1] + (1-ud)*Krg[i];
		To[i] = Go*(kroud/Bomed);
		Tg[i] = Gw*krgud;
	}
	Bomed = calc_Bo(dados.p_E);
	/// S_E = Sg[nx-1] (condicao de derivada nula)
	kroud = Kro[dados.nx-1];
	krgud = Krg[dados.nx-1];
	To[dados.nx] = 2*Go*(kroud/Bomed);
	Tg[dados.nx] = 2*Gw*krgud;
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
	while (current_time < 1.0){
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
		update(Sg_old,Sg);
	}
}

void funcionamento(){
	// preencher a matriz de coeficientes
	// preencher a matriz de termos independentes
	// invocar o Gauss-Siedel, utilizando as pressao/saturação inicial como estimativa
	// obtem novas estimativa para pressao/saturação
	
	// while (tempo_atual < tempo_max)
	/* 
		corrige os termos independentes
		atualiza a matriz de coeficientes
		invocar o Gauss-Siedel, utilizando as pressao/saturação obtidas
		obtem novas estimativa para pressao/saturação
	*/
}
void TransmissibilidadeGas(){
	//Tg
	double conv = 8.64E-2;
	const auto Go = (conv*dados.k)/(dados.muo*dados.dx*dados.dx);
	auto Bo_med = calc_Bo(dados.p_W);
	auto ud = upwind(dados.p_W, p[0]);
	auto kroud = ud*calc_kro(dados.Sg_W) + (1-ud)*Kro[0];
	Tg[0] = 2*Go*(kroud/Bo_med);
	for(int i = 1; i < dados.nx; i++){
		Bo_med = (Bo[i-1] + Bo[i])/2;
		ud    = upwind(p[i-1],p[i]);
		kroud = ud*Kro[i-1] + (1-ud)*Kro[i];
		Tg[i] = Go*(kroud/Bo_med);
	}
	Bo_med = calc_Bo(dados.p_E);
	kroud = Kro[dados.nx-1];
	Tg[dados.nx] = 2*Go*(kroud/Bo_med);
}
void TransmissibilidadeOleo(){
	double conv = 8.64E-2;
	const auto Go = (conv*dados.k)/(dados.muo*dados.dx*dados.dx);
	auto Bomed = calc_Bo(dados.p_W);
	auto ud = upwind(dados.p_W, p[0]);
	auto kroud = ud*calc_kro(dados.Sg_W) + (1-ud)*Kro[0];
	To[0] = 2*Go*(kroud/Bomed);
	for(int i = 1; i < dados.nx; i++){
		Bomed = (Bo[i-1] + Bo[i])/2;
		ud    = upwind(p[i-1],p[i]);
		kroud = ud*Kro[i-1] + (1-ud)*Kro[i];
		To[i] = Go*(kroud/Bomed);
	}
	Bomed = calc_Bo(dados.p_E);
	/// S_E = Sg[nx-1] (condicao de derivada nula)
	kroud = Kro[dados.nx-1];
	To[dados.nx] = 2*Go*(kroud/Bomed);
}
void DerivadaPressaoCapilar(){
	//para implementar
}
void CalcularCgg(){
	auto constante = dados.phi/dados.dt;
	for(int i=0; i < dados.nx; i++){
		Cgg[i]= constante*1/(calc_Bo(P[i])); // B do gas
	}
	
}
void CalcularCgp(){
	auto constante = dados.phi/dados.dt;
	for(int i=0; i < dados.nx; i++){
		Cgp[i] = derivada_B_gas()*Sg[i];
	}
}
void CalcularCog(){
	auto constante = - dados.phi/dados.dt;
	for(int i=0; i < dados.nx; i++){
		Cog[i] = constante*1/(calc_Bo(P[i])); // B do oleo
	}
}
void CalcularCop(){
	auto constante = dados.phi/dados.dt;
	for(int i=0; i < dados.nx; i++){
		Cop[i] = derivada_B_oleo()*(1-Sg[i]);
	}
}
void derivada_B_gas(){
	//para implementar
	
}
void derivada_B_oleo(){
	//para implementar
	
}

//montar vetor de incogintas inicial com 2*nx	
void criarvx0 (){
	for (int w=0; w<=2*dados.nx;w=w+2){
		X0[w]=dados.Sg_0;
		X0[w+1]=dados.p_0;
	}
}
//matriz dos coeficientes
void ccoef(){
	// analizar que pode ter coef errados 
	int m=2*dados.nx;
	int n=2*dados.nx;
	
	//Preencher os vetores: Tg, To, Cgg, Cgp, Cog, Cop, Pc
	TransmissibilidadeGas(); //Tg
	DerivadaPressaoCapilar(); //Pc
	CalcularCgg(); //Cgg
	CalcularCgp(); //Cgp
	TransmissibilidadeOleo(); //To
	CalcularCog(); //Cog
	CalcularCop(); //Cop

	// Condiçoes de contorno esquerda
	// saturação
	coef[0][0]= 1; //Cgg[0]+Tg[0-1]*Pc[0-1]+Tg[0]Pc[0];
	coef[0][1]= 1; //Tg[0]+Tg[0-1]+Cgp[0];
	coef[0][2]= 1; // -(Tg[0-1]*Pc[0];
	coef[0][3]= - Tg[0];
	// pressão
	coef[1][0]= Cog[0];
	coef[1][1]= 1; // To[0]+To[0-1]+Cop[0];
	coef[1][3]= -To[0];

	for(int i=2; i<m-2; i=i+2){
		//preenchimento saturação 
		int k = i/2; //os vetores têm nx de tamanho, mas a matrix tem 2nx
		int j = i-2;
		coef[i][j]= - (Tg[k-1]*Pc[k-1]);
		coef[i][j+1]= -Tg[k-1];
		coef[i][j+2]= Cgg[k]+Tg[k-1]*Pc[k-1]+Tg[k]*Pc[k];
		coef[i][j+3]= Tg[k]+Tg[k-1]+Cgp[k];
		coef[i][j+4]= -(Tg[k-1]*Pc[k]);
		coef[i][j+5]= - Tg[k];
		// pressão
		j=i-1;
		coef[i+1][j]=-To[k-1];
		coef[i+1][j+1]= Cog[k];
		coef[i+1][j+2]= To[k]+To[k-1]+Cop[k];
		coef[i+1][j+4]= -To[k];
	}
	// saturação
	coef[m-2][m-4] = -(Tg[nx-2]*Pc[nx-2];
	coef[m-2][m-3] = -Tg[nx-2];
	coef[m-2][m-2] = Cgg[nx-1]+Tg[nx-2]*Pc[nx-2]+Tg[nx-1]Pc[nx-1];
	coef[m-2][m-1] = Tg[nx-1]+Tg[(nx-2]+Cgp[nx-1];
	// pressão
	coef[m-1][m-3] = -To[nx-2];
	coef[m-1][m-2] = Cog[nx-1]; 
	coef[m-1][m-1] = To[nx-1]+To[nx-2]+Cop[nx-1];
	
}
//vetor de termos constantes
void cvecD(){
	int m=2*dados.nx;
	D[0]= 0; // Cgp[0]*p_old[0]+Cgg[0]*Sg_old[0] Tg[0-1]*P_E+ td[0-1]*Pc[0-1]*S_E;
	D[1]= 0; // Cop[0]*p_old[0]+Cog[0]*Sg_old[0]+To[0-1]*P_E;
	for (int i=2;i<m-2;i=i+2){
		D[i]=Cgp[i/2]*p_old[i/2]+Cgg[i/2]*Sg_old[i/2];
		D[i+1]=Cop[i/2]*p_old[i/2]+Cog[i/2]*Sg_old[i/2];
	}
	D[m-2]=0; //Cgp[nx-1]*p_old[nx-1]+Cgg[nx-1]*Sg_old[nx-1]+Tg[nx-1]*P_D+ Tg[nx-2]*Pc[nx-2]*S_D;
	D[m-1]= 0; //Cop[nx-1]*p_old[nx-1]+Cog[nx-1]*Sg_old[nx-1]+Cog[nx-1]*Sg_old[mx-1]+To[nx-1]*P_D;
}
void fill_matrix(){
	auto nx=dados.nx;
	auto bg=dados.Bg;
	M[0] = Bo[0]*To[0] + bg*Tg[0] + Bo[0]*To[1] + bg*Tg[1] + Bo[0]*Cop[0];
	E[0] = - (Bo[0]*To[1] + bg*Tg[1]);
	D[0] = Bo[0]*Cop[0]*p_old[0] + (Bo[0]*To[0] + bg*Tg[0])*dados.p_W;

	for(int i = 1; i < dados.nx-1; i++){
		W[i] = - (Bo[i]*To[i] + bg*Tg[i]);
		E[i] = - (Bo[i]*To[i+1] + bg*Tg[i+1]);
		M[i] = Bo[i]*Cop[i] - W[i] - E[i];
		D[i] = Bo[i]*Cop[i]*p_old[i];
	}
	W[nx-1] = - (Bo[nx-1]*To[nx-1] + bg*Tg[nx-1]);
	M[nx-1] = Bo[nx-1]*To[nx-1] + bg*Tg[nx-1] + Bo[nx-1]*To[nx] + bg*Tg[nx] + Bo[nx-1]*Cop[nx-1];
	D[nx-1] = Bo[nx-1]*Cop[nx-1]*p_old[nx-1] + (Bo[nx-1]*To[nx] + bg*Tg[nx])*dados.p_E;
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
// não é necessaria 
void solver_s(){
	auto nx=dados.nx;
	Sg[0] = Sg_old[0] + (1.0/Cgg[0]) * ( Tg[1]*(p[1] - p[0]) + Tg[0]*(dados.p_W - p[0]));
	for(int i = 1; i < nx-1; i++){
		Sg[i] = Sg_old[i] + (1.0/Cgg[i]) * ( Tg[i+1]*(p[i+1] - p[i]) + Tg[i]*(p[i-1] - p[i]) );
	}
	Sg[nx-1] = Sg_old[nx-1] + (1.0/Cgg[nx-1]) * ( Tg[nx]*(dados.p_E - p[nx-1]) + Tg[nx-1]*(p[nx-2] - p[nx-1]));
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
		dados_entrada.mug=Dados[6];
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
		dados_entrada.Sg_0=Dados[18];
		dados_entrada.Sg_W=Dados[19];
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
	Sg.resize(n,dados.Sg_0);
	Sg_old.resize(n,dados.Sg_0);
	Sg_it.resize(n,dados.Sg_0);
	Cop.resize(n,0.0);
	Cog.resize(n,0.0);
	Cgg.resize(n,0.0);
	Cgp.resize(n,0.0);
	Kro.resize(n,0.0);
	Krg.resize(n,0.0);
	Bo.resize(n,0.0);
	Bo_old.resize(n,0.0);
	To.resize(n+1,0.0);
	Tg.resize(n+1,0.0);
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
	std::cout << "Viscosidade da agua: " << dados.mug << "\n";
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
	std::cout << "Saturacao inicial: " << dados.Sg_0 << "\n";
	std::cout << "Saturacao no contorno esquerdo: " << dados.Sg_W << "\n";
}

void print_array_2D(const std::vector<Vec1D>& M){
	
	std::cout << "Entered the print array 2D function" << std::endl;
	
	auto nrow=M.size();
	auto ncol=M[0].size();
	std::cout << "Numero de linhas eh: " << nrow << std::endl;
	std::cout << "Numero de colunas eh: " << ncol << std::endl;
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			std::cout << M[i][j] << ' ';
		}
		std::cout << '\n';
	}
	std::cout << "Exit the print array 2D function" << std::endl;
}

void print_vector(const Vec1D& M){
	for (int i=0; i < M.size(); i++){
		std::cout << M[i] << ' ' ;
	}
	std::cout << '\n';
}