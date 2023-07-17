#include "newprototypes.h"

//Váriáveis da simulação
DadosTotais dados;

//Vetores para armazenamento:
Vec1D p(dados.nx, dados.p_0);           // n elementos, inicializados com p_0
Vec1D p_old(dados.nx, dados.p_0);       // n elementos, inicializados com p_0
Vec1D p_it(dados.nx, dados.p_0);        // n elementos, inicializados com p_0
Vec1D Sg(dados.nx, dados.Sg_0);         // n elementos, inicializados com Sg_0
Vec1D Sg_old(dados.nx, dados.Sg_0);     // n elementos, inicializados com Sg_0
Vec1D Sg_it(dados.nx, dados.Sg_0);      // n elementos, inicializados com Sg_0
//Vetores para as propriedas:
Vec1D Cop(dados.nx, 0.0);
Vec1D Cog(dados.nx, 0.0);
Vec1D Cgg(dados.nx, 0.0);
Vec1D Cgp(dados.nx, 0.0);
Vec1D Kro(dados.nx, 0.0);
Vec1D Krg(dados.nx, 0.0);
Vec1D Bo(dados.nx, 0.0);
Vec1D Bo_old(dados.nx, 0.0);
Vec1D Bg(dados.nx, 0.0);
Vec1D Bg_old(dados.nx, 0.0);
Vec1D To(dados.nx+1, 0.0);               // Transmissibilidade oleo
Vec1D Tg(dados.nx+1, 0.0);               // Transmissibilidade gas
Vec1D Pc(dados.nx+1,0.0);
Vec1D Pc_old(dados.nx+1,0.0);
Vec1D dPc(dados.nx,0.0);
Vec1D dBo(dados.nx,0.0);
//Matrizes de coeficientes:
std::vector<Vec1D> coefp(dados.nx, Vec1D(dados.nx, 0.0));
std::vector<Vec1D> coefS(dados.nx, Vec1D(dados.nx, 0.0));
//Termos independentes:
Vec1D indepen_p(dados.nx, 0.0);
Vec1D indepen_S(dados.nx, 0.0);



void preencher_coef_saturacao(){
	int m=dados.nx;
	
	coefS[0][0]=Cgg[0]+Tg[0]*dPc[0]+Tg[0]*dPc[0];
	coefS[0][1]= -Tg[0]*dPc[0]; 
	for(int i=1; i<m-1; i=i+1){
		coefS[i][i-1]= - (Tg[i-1]*dPc[i-1]);
		coefS[i][i]  = Cgg[i]+Tg[i-1]*dPc[i-1]+Tg[i]*dPc[i];
		coefS[i][i+1]= -(Tg[i]*dPc[i]);
	}
	coefS[m-1][m-2] = -Tg[m-2]*dPc[m-2];
	coefS[m-1][m-1] =Cgg[m-1]+Tg[m-2]*dPc[m-2]+Tg[m-1]*dPc[m-1];
}
void preencher_indep_saturacao(){
	int m=dados.nx;

	indepen_S[0]=-Cgp[0]*(p[0]-p_old[0])+Cgg[0]*Sg_old[0]+ Tg[0]*(dados.p_W-p[0])+Tg[0]*(p[1]-p[0]) + Tg[0]*dPc[0]*dados.Sg_W;
	for (int i=1;i<m-1;i=i+1){
		indepen_S[i]=-Cgp[i]*(p[i]-p_old[i])+Cgg[i]*Sg_old[i]+ Tg[i-1]*(p[i-1]-p[i])+Tg[i]*(p[i+1]-p[i]);
	}
	indepen_S[m-1]=-Cgp[m-1]*(p[m-1]-p_old[m-1])+Cgg[m-1]*Sg_old[m-1]+ Tg[m-2]*(p[m-2]-p[m-1])+Tg[m-1]*(dados.p_E-p[m-1])+ Tg[m-2]*dPc[m-2]*dados.Sg_0;
}


void preencher_coef_pressao(){
	int m=dados.nx;

	coefp[0][0] =  To[0]+To[0]+Cop[0];
	coefp[0][1] = -To[0];
	for(int i=1; i<m-1; i=i+1){
		coefp[i][i-1] = -To[i-1]; 
		coefp[i][i]   =  To[i]+To[i-1]+Cop[i];                        
		coefp[i][i+1] = -To[i];                           
	}
	coefp[m-1][m-2] = -To[m-2];             
	coefp[m-1][m-1] = To[m-1]+To[m-2]+Cop[m-1];
}
void preencher_indep_pressao(){
	int m=dados.nx;

	indepen_p[0]= Cop[0]*p_old[0]+Cog[0]*Sg_old[0]-Cog[0]*Sg[0] +To[0]*dados.p_W;
	for (int i=1;i<m-1;i++){
		indepen_p[i] = Cop[i]*p_old[i]+Cog[i]*Sg_old[i]- Cog[i]*Sg[i];
	}
	indepen_p[m-1]= Cop[m-1]*p_old[m-1]+Cog[m-1]*Sg_old[m-1] - Cog[m-1]*Sg[m-1]+To[m-1]*dados.p_E;
}

int upwind(double pesq, double pdir){
	return (pesq > pdir) ? 1 : 0;
}
double ca_kro (double S){
	double frac= ((1.0 - S)*dados.sor)/(1.0 - dados.sor);
	double exp= (2.0 + 3.0*dados.Alf)/dados.Alf;
	return std::pow(frac,exp);
}
double ca_krg (double S){
	double frac=(1.0 - S - dados.sor)/(1.0 - dados.sor);
	double frac2=(S)/(1.0-dados.sor);
	double exp=(2.0 + dados.Alf)/dados.Alf;
	return std::pow(frac2,2)*(1-std::pow(frac,exp));
}
void permeabilidades_relativas(){
	for(int i = 0; i < dados.nx; i++){
		Kro[i] = ca_kro(Sg[i]);
		Krg[i] = ca_krg(Sg[i]);
	}
}

double calc_Bo(const double P){
	return dados.Boref/( 1.0 + dados.co*(P - dados.pref));
}
void propriedades_oleo(Vec1D& B, const Vec1D& P){
	for(int i = 0; i < dados.nx; i++){
		B[i] = calc_Bo(P[i]);
	}
}
double calc_Bg(double P){
	return 0.0282*((dados.T*dados.Z)/P);
}
void propriedades_gas(Vec1D& B, const Vec1D& P){
	for(int i = 0; i < dados.nx; i++){
		B[i] = calc_Bg(P[i]); 
	}
}
void TransmissibilidadeGas(){
	double conv =8.64E-2;
	const auto Gg = (conv*dados.k)/(dados.muo*dados.dx*dados.dx);
	auto Bg_med = calc_Bg(dados.p_W);
	auto ud = upwind(dados.p_W, p[0]);
	auto krgud =Krg[0];
	Tg[0] =Gg*(krgud/Bg_med); 
	for(int i = 1; i < dados.nx; i++){
		Bg_med = (Bg[i-1] + Bg[i])/2;
		ud    = upwind(p[i-1],p[i]);
		krgud = ud*Krg[i-1] + (1-ud)*Krg[i];
		Tg[i] = Gg*(krgud/Bg_med);
		std::cout << "Tg= " << krgud << std::endl;
	}
	Bg_med = calc_Bg(dados.p_E);
	krgud = Krg[dados.nx-1];
	Tg[dados.nx] = Gg*(krgud/Bg_med);
}
void TransmissibilidadeOleo(){
	double conv =8.64E-2;
	const auto Go = (conv*dados.k)/(dados.muo*dados.dx*dados.dx);
	auto Bomed = calc_Bo(dados.p_W);
	auto ud = upwind(dados.p_W, p[0]);
	auto kroud = Kro[0];
	To[0] = Go*(kroud/Bomed);
	for(int i = 1; i < dados.nx; i++){
		Bomed = (Bo[i-1] + Bo[i])/2;
		ud    = upwind(p[i-1],p[i]);
		kroud = ud*Kro[i-1] + (1-ud)*Kro[i];
		To[i] = Go*(kroud/Bomed);
	}
	Bomed = calc_Bo(dados.p_E);
	kroud = Kro[dados.nx-1];
	To[dados.nx] =Go*(kroud/Bomed);// tirei o 2
}
double CalcularPcgo(double Sg){
	double Pct = 1.0; 
	double Sorg = dados.sor;
	double lambda = 2.0; 
	return Pct*(std::pow(((1.0-Sorg)/((1.0-Sg)-Sorg)),(1.0/lambda)));
}
double DerivadaPressaoCapilar(double p, double S){
	return ((1/p)-(1/CalcularPcgo(S-0.05)))/(0.05);
}
double derivada_B_gas(double b, double p){
	return ((1/b)-(1/calc_Bg(p-10.0)))/(10.0);
}
double derivada_B_oleo(double b, double p){
	double ret {0.0};
	double cte = 1.5;
	ret=(b-calc_Bo(p-cte))/cte;
	double bot=b;
	double botmenos=calc_Bo(p-cte);
	return ret;
}
void CalcularCgg(){
	auto constante = dados.phi/dados.dt;
	for(int i=0; i < dados.nx; i++){
		Cgg[i]= constante/(calc_Bg(p[i]));
	}
	
}
void CalcularCgp(){
	auto constante = dados.phi/dados.dt;
	for(int i=0; i < dados.nx; i++){
		Cgp[i] = constante*derivada_B_gas(Bg[i],p[i])*Sg[i];
	}
}
void CalcularCog(){
	auto constante = - dados.phi/dados.dt;
	for(int i=0; i < dados.nx; i++){
		Cog[i] = constante/(calc_Bo(p[i]));
	}
}
void CalcularCop(){
	auto constante = - dados.phi/dados.dt;
	for(int i=0; i < dados.nx; i++){
		dBo[i] = derivada_B_oleo(calc_Bo(p[i]),p[i]);
		Cop[i] = constante*dBo[i]*(1.0-Sg[i]);
	}	
}
void vetor_pressao_capilar(Vec1D& p,const Vec1D& S){
	for (int i=0;i<dados.nx;i++){
		p[i]= CalcularPcgo(S[i]);
	}
}
void vetor_derivada_pressao_capilar (){
	for(int i=0;i<dados.nx;i++){
		dPc[i]=DerivadaPressaoCapilar(Pc[i],Sg[i]);
	}
}



void update_new_vector_with_old(Vec1D& V0, const Vec1D& V){
	const auto n = V0.size();
	for(int i = 0; i < n; i++){
		V0[i] = V[i];
	}
}

void resize_if_needed(int n){
	//Se os vetores criados não sao suficientes, redimensionar para o novo tamanho:
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
	Bg.resize(n,0.0);
	Bg_old.resize(n,0.0);
	To.resize(n+1,0.0);
	Tg.resize(n+1,0.0);
	Pc.resize(n+1,0.0);	
	Pc_old.resize(n+1,0.0);	
	dPc.resize(n,0.0);	
	dBo.resize(n,0.0);
	indepen_p.resize(n,0.0);
	indepen_S.resize(n,0.0);


	//Resizing: coefp
	auto nlines = coefp.size();
	auto ncol= coefp[0].size();
	if (nlines != ncol)
		std::cout<< "A matriz de coeficientes coefp nao e quadrada" << std::endl;
	else
		{
		if (nlines != n){
			std::cout << "Redimensionando a matriz de coeficientes coefp para " << n << "x" << n << std::endl;
			coefp.resize(n);
				for (int i = 0; i < n; ++i)
					coefp[i].resize(n,0.0);
		}
	}
	//Resizing: coefS
	nlines = coefS.size();
	ncol= coefS[0].size();
	if (nlines != ncol)
		std::cout<< "A matriz de coeficientes coefS nao e quadrada" << std::endl;
	else
		{
		if (nlines != n){
			std::cout << "Redimensionando a matriz de coeficientes coefS para " << n << "x" << n << std::endl;
			coefS.resize(n);
				for (int i = 0; i < n; ++i)
					coefS[i].resize(n,0.0);
		}
	}
}



DadosTotais read_input_data(std::string filename){
	DadosTotais dados_entrada;
	//Nome do arquivo de entrada
	std::ifstream input_file {filename};
	if (input_file.is_open()){
		Vec1D Dados;
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
		dados_entrada.muo=Dados[9]; //muo
		dados_entrada.tempo_max=static_cast<int>(Dados[10]);
		dados_entrada.dt=Dados[11];
		dados_entrada.Lx=Dados[12];
		dados_entrada.nx=static_cast<int>(Dados[13]);
		dados_entrada.tolerancia=Dados[14];
		dados_entrada.p_0=Dados[15];//inicail
		dados_entrada.p_W=Dados[16];//esquerda
		dados_entrada.p_E=Dados[17];//direita
		dados_entrada.Sg_0=Dados[18];//inicial
		dados_entrada.Sg_W=Dados[19];//esquerda
		dados_entrada.Alf=Dados[20];//
		dados_entrada.Z=Dados[21];
		dados_entrada.T=Dados[22];
	}
	else{
		std::cerr << "Nao foi possivel encontrar o arquivo " << filename;
		std::exit(1);
	}
	return dados_entrada;
}
void gauss_solver(std::vector<Vec1D>& A, Vec1D& B, Vec1D& X){
	//A = coeficientes, B = termos independentes, X = estimativas iniciais, sera sobreescrito
	auto m = B.size();   // numero de linhas
	auto n = m;          // numero de colunas
	Vec1D Y = X;          // matriz auxiliar
	Vec1D E = X;         // necessária para estimar o erro de uma iteração a outra
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
			teste = teste && res;
			E[i] = X[i];
		}
		counter++;
	}
}
bool is_diagonal_dom(const std::vector<Vec1D>& M){
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
void print_array_2D(const std::vector<Vec1D>& M){
	auto nrow=M.size();
	auto ncol=M[0].size();
	std::cout << "Exibindo matriz 2D (" << nrow << "x" << ncol << ")"<< std::endl;
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			std::cout << std::fixed << std::setprecision(6) << M[i][j] << ' ';
		}
		std::cout << '\n';
	}
}
void print_vector(const Vec1D& M){
	for (int i=0; i < M.size(); i++){
		std::cout << M[i] << ' ' ;
	}
	std::cout << '\n';
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



void preencher_vetores_aux(){
	permeabilidades_relativas();
	TransmissibilidadeGas();          //vetor Tg
	TransmissibilidadeOleo();         //vetor To
	CalcularCgg();                    //vetor Cgg
	CalcularCog();                    //vetor Cog
	vetor_pressao_capilar(Pc_old,Sg); //vetor Pc
	vetor_derivada_pressao_capilar ();//vetor dPc
}


void solver(){
	propriedades_oleo(Bo,p);
	propriedades_gas(Bg,p);
	
	double tempo_atual = 0.0;
	while (tempo_atual < 1.0){
		//Preencher os vetores: Tg, To, Cgg, Cgp, Cog, Cop, Pc:
		std::cout << "Exibindo Tg:" << std::endl;
		print_vector(Tg);
		std::cout << "Exibindo Bo_old:" << std::endl;
		print_vector(Bo_old);
		preencher_vetores_aux();
		
		std::cout << "Exibindo Tg:" << std::endl;
		print_vector(Tg);
		
		
		//Prencher os vetores: coefp, indepen_p, coefS, indepen_S
		preencher_coef_pressao();
		preencher_indep_pressao();
		preencher_coef_saturacao();
		preencher_indep_saturacao();

		//Exibição:
		std::cout << "Tempo atual: " << tempo_atual << std::endl;
		std::cout << "Matriz de p: " << std::endl;
		print_array_2D(coefp);
		std::cout << "Matriz de Sg: " << std::endl;
		print_array_2D(coefS);
		


		//Armazenamento das informações, antes da resolução dos sistemas
		update_new_vector_with_old(Bg, Bg_old);
		update_new_vector_with_old(Bo, Bo_old);
		update_new_vector_with_old(Sg_old,Sg);
		update_new_vector_with_old(p_old,p);
		update_new_vector_with_old(Pc,Pc_old);
		
		//Executa a soluçao dos sistemas:
		gauss_solver(coefp,indepen_p,p);  //escreve em p as novas estimativas
		gauss_solver(coefS,indepen_S,Sg); //escreve em Sg as novas estimativas
		
		//Corrige:

		
		tempo_atual += dados.dt;
	}
}


int main(int argc, char *argv[] ){
	std::cout << "Inicio da simulacao..." << std::endl;
	//Inicia o registro do tempo:
	CustomTimer Cronometro;

	//Leitura das váriaveis a partir de um arquivo de entrada:
	dados=read_input_data("input_data.txt");
	int n=dados.nx;
	dados.dx = dados.Lx/n;
	dados.pref = dados.p_0;
	//Exibição das informações adquiridas:
	print_simulation_properties();
	
	//Verificação das dimensoes dos vetores:
	resize_if_needed(n);

	//Execução do solver
	solver();

	//Salva os dados:
	save_full_data(p,"Pressao");
	save_full_data(Sg,"Saturacao");
	std::cout << "\nFim da simulacao" << std::endl;
	return 0;
}
