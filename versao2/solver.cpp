#include <iostream>
#include "prototypes.h"

//Váriáveis da simulação
DadosTotais dados;
//Inicialização de vetores utilizados na simulação:
Vec1D Dp(dados.nx, 0.0);                // n elementos, inicializados em 0
Vec1D DS(dados.nx, 0.0);                // n elementos, inicializados em 0
Vec1D RESp(dados.nx, 0.0);              // n elementos, inicializados em 0
Vec1D RESs(dados.nx, 0.0);              // n elementos, inicializados em 0
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
Vec1D Bg(dados.nx, 0.0);
Vec1D Bg_old(dados.nx, 0.0);
Vec1D To(dados.nx, 0.0);                 // Transmissibilidade
Vec1D Tg(dados.nx+1, 0.0);               // Transmissibilidade
Vec1D Pc(dados.nx+1,0.0);
Vec1D X0(2*dados.nx,0.0);
Vec1D X(2*dados.nx,0.0);
Vec1D dPc(dados.nx,0.0);
Vec1D dBo(dados.nx,0.0);
Vec1D soma(2*dados.nx,0.0);
Vec1D beta(2*dados.nx,0.0);
Vec1D diagonal(2*dados.nx,0.0);
Vec1D diferenca(2*dados.nx,0.0);
Vec1D Pc_old(dados.nx+1,0.0);
std::vector<Vec1D> coefp(dados.nx, Vec1D(dados.nx, 0.0));
std::vector<Vec1D> coefS(dados.nx, Vec1D(dados.nx, 0.0));
std::vector<Vec1D> tnot(dados.nx+1, Vec1D(dados.tempo_max/100, 0.0));
Vec1D digcoef(2*dados.nx,0.0);
Vec1D q(dados.nx, 0.0);
Vec1D qp(dados.nx, 0.0);


int main(int argc, char *argv[] ){
	std::cout << "Inicio da simulacao..." << std::endl;
	//Inicia o registro do tempo:
	CustomTimer Cronometro;
	//Leitura das váriaveis a partir de um arquivo de entrada:
	dados=read_input_data("input_data.txt");
	int n=dados.nx; //numero de celulas, lido do arquivo de entrada
	dados.dx = dados.Lx/n;
	dados.pref = dados.p_0;
	 //solucao_sist(coefp,coefS,Dp,DS,p,Sg);
	 
	 
	oil_prop(Bo_old,p);
	gas_prop(Bg_old,p);
	update(Bo, Bo_old);
	update(Bg, Bg_old);
	update(p_old,p);
	update(Sg_old,Sg);
	 
	funcionamento();
	//Salva os dados:
	save_full_data(p,"Pressao");
	save_full_data(Sg,"Saturacao");
	std::cout << "\nFim da simulacao" << std::endl;
	return 0;
}


void funcionamento(){
	double temat {0.0};
	ccoefS();
	ccoefp();
	cvecDS();
	cvecDp();

	while (temat < 1.0){
	
		std::cout << "\n\ntemat " << temat << std::endl;
		std::cout << "coeficientes de p:" << std::endl;
		print_array_2D(coefp);
		std::cout << "coeficientes de S:" << std::endl;
		print_array_2D(coefS);
		oil_prop(Bo,p);
		gas_prop(Bg,p);
		ccoefS();
		ccoefp();
		resolverSistemas(coefp,coefS,Dp,DS,p,Sg);
		//Verificar a convergência de p e S:
		
		update(p_old,p);
		update(Sg_old,Sg);
		temat=temat+dados.dt;
		
		update(Bo, Bo_old);
		update(Bg, Bg_old);
	}
}


/*
Eq para saturação
-Tg[i-1]*dPc[i-1]*Sg[i-1]+(Cgg[i]+Tg[i-1]*dPc[i-1]+Tg[i]*dPc[i])*Sg[i]-Tg[i]*dPc[i]*Sg[i+1]=


Eq para saturação
-Tg[i-1]*Pc[i-1]*Sg[i-1]+(Cgg[i]+Tg[i-1]*Pc[i-1]+Tg[i]*Pc[i])*Sg[i]-Tg[i1]*Pc[i]*Sg[i+1]=
-Cgp[i]*(po[i]-po_old[i])+Cgg[i]*Sg_old[i]+ Tg[i-1]*(p[i-1]-p[i])+Tg[i]*(p[i+1]-p[i])
*/


// separando
void ccoefp(){
	// analizar que pode ter coef errados 
	int m=dados.nx;
	int n=dados.nx;
	double c;
	//preencher vetores kr
	oil_prop(Bo_old,p);
	update(Bo, Bo_old);
	rel_perm();
	//Preencher os vetores: Tg, To, Cgg, Cgp, Cog, Cop, Pc
	TransmissibilidadeOleo(); //To
	CalcularCog(); //Cog
	CalcularCop(); //Cop
	mmll(tnot,Tg,c);
	// pressão
	coefp[0][0]=To[0]+To[0]+Cop[0];
	coefp[0][1]=-To[0];
	for(int i=1; i<m-1; i=i+1){
		// pressão
		int j=i-1;
		int k=i;
		coefp[i][j]=-To[k-1]; 
		coefp[i][j+1]= To[k]+To[k-1]+Cop[k];                        
		coefp[i][j+2]= -To[k];                           
	}
	// pressão
	coefp[m-1][m-2] =-To[dados.nx-2];             
	coefp[m-1][m-1] =To[dados.nx-1]+To[dados.nx-2]+Cop[dados.nx-1];
}


void ccoefS(){
	// analizar que pode ter coef errados 
	int m=dados.nx;
	int n=dados.nx;
	double c;
	//preencher vetores kr
	gas_prop(Bg_old,p);
	update(Bg, Bg_old);
	vetpc(Pc_old,Sg);
	update(Pc,Pc_old);
	rel_perm();
	//Preencher os vetores: Tg, To, Cgg, Cgp, Cog, Cop, Pc
	TransmissibilidadeGas(); //Tg
	vetpc(Pc_old,Sg);
	update(Pc,Pc_old);
	vecdPc ();//dPc
	CalcularCgg(); //Cgg
	CalcularCog(); //Cog
	mmll(tnot,Tg,c);
	// Condiçoes de contorno na esquerda (linhas 0 e 1 da matrix de coeficientes)
	// saturação
	coefS[0][0]=Cgg[0]+Tg[0]*dPc[0]+Tg[0]*dPc[0];
	coefS[0][1]= -Tg[0]*dPc[0]; 

	for(int i=1; i<m-1; i=i+1){
		//preenchimento saturação 
		int j = i-1;
		int k=i;
		coefS[i][j]= - (Tg[k-1]*dPc[k-1]);
		coefS[i][j+1]= Cgg[k]+Tg[k-1]*dPc[k-1]+Tg[k]*dPc[k];// diagonal princ
		coefS[i][j+2]= -(Tg[k]*dPc[k]);
	}
	// saturação
	coefS[m-1][m-2] = -Tg[dados.nx-2]*dPc[dados.nx-2];
	coefS[m-1][m-1] =Cgg[dados.nx-1]+Tg[dados.nx-2]*dPc[dados.nx-2]+Tg[dados.nx-1]*dPc[dados.nx-1];
}


void cvecDS(){
	int m=dados.nx;
	/////////////////////////////////////////////// esse era i-1 ////
	DS[0]=-Cgp[0]*(p[0]-p_old[0])+Cgg[0]*Sg_old[0]+ Tg[0]*(dados.p_W-p[0])+Tg[0]*(p[1]-p[0]) + Tg[0]*dPc[0]*dados.Sg_W;
	//-Cgp[0]*p_old[0]+Cgg[0]*Sg_old[0]+ Tg[0]*dPc[0]*dados.Sg_W -(-Tg[0]*dados.p_W+(Tg[0]+Tg[0]+Cgp[0])*p[0] - Tg[0]*p[1]);
	for (int i=1;i<m-1;i=i+1){
		//int k=i;//2;
		DS[i]=-Cgp[i]*(p[i]-p_old[i])+Cgg[i]*Sg_old[i]+ Tg[i-1]*(p[i-1]-p[i])+Tg[i]*(p[i+1]-p[i]);
		//Cgp[i]*p_old[i]+Cgg[i]*Sg_old[i]-(-Tg[i-1]*p[i-1]+(Tg[i]+Tg[i-1]+Cgp[i])*p[i] - Tg[i]*p[i+1]);
	}
	DS[m-1]=-Cgp[m-1]*(p[m-1]-p_old[m-1])+Cgg[m-1]*Sg_old[m-1]+ Tg[m-2]*(p[m-2]-p[m-1])+Tg[m-1]*(dados.p_E-p[m-1])+ Tg[dados.nx-2]*dPc[dados.nx-2]*dados.Sg_0;
	//Cgp[dados.nx-1]*p_old[dados.nx-1]+Cgg[dados.nx-1]*Sg_old[dados.nx-1] -(-Tg[m-2]*p[m-2]+(Tg[m-1]+Tg[m-2]+Cgp[m-1])*p[m-1] - Tg[m-1]*dados.p_E)+ Tg[dados.nx-2]*dPc[dados.nx-2]*dados.Sg_0;
}


void termofonte(){// não ta em uso
	int m=dados.nx;
	q[0]= -Tg[0]*dados.p_W+(Tg[0]+Tg[0]+Cgp[0])*p[0] - Tg[0]*p[1];
	// -(Tg[0]*dados.p_W -(Cgp[0]+Tg[0]+Tg[1])*X[0]+ Tg[1]*X[1]);
	q[1]= -Tg[0]*p[0]+(Tg[1]+Tg[0]+Cgp[1])*p[1] - Tg[1]*p[2];
	for (int i=2;i<m-1;i++){
		q[i]= -Tg[i-1]*p[i-1]+(Tg[i]+Tg[i-1]+Cgp[i])*p[i] - Tg[i]*p[i+1];
		//-(-Tg[i-1]*X[i-1] +(Cgp[i]+Tg[i-1]-Tg[i+1])*X[i]+ Tg[i+1]*X[i+1]);
	}
	q[m-1]=-Tg[m-2]*p[m-2]+(Tg[m-1]+Tg[m-2]+Cgp[m-1])*p[m-1] - Tg[m-1]*dados.p_E;
	//-Tg[m-2]*X[m-2] +(Cgp[m-1]+Tg[m-2]+Tg[m-1])*X[m-1]- Tg[m-1]*dados.p_E;
	
}

// Vetor de termos constantes
void termofontep(){// não ta em uso
	int m=dados.nx;
	for(int i=1;i<m;i++){
		qp[i]=Cog[i]*Sg[i];
	}
}
void cvecDp(){
	int m=dados.nx;
	Dp[0]= Cop[0]*p_old[0]+Cog[0]*Sg_old[0]-Cog[0]*Sg[0] +To[0]*dados.p_W; 
	for (int i=1;i<m-1;i++){
		Dp[i]=Cop[i]*p_old[i]+Cog[i]*Sg_old[i]- Cog[i]*Sg[i];
	}
	Dp[m-1]= Cop[dados.nx-1]*p_old[dados.nx-1]+Cog[dados.nx-1]*Sg_old[dados.nx-1] - Cog[dados.nx-1]*Sg[dados.nx-1]+To[dados.nx-1]*dados.p_E;
}


void update(Vec1D& V0, const Vec1D& V){
	const auto n = V0.size();
	for(int i = 0; i < n; i++){
		V0[i] = V[i];
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

double dB(double bo, double bo0, double P, double P0){
	double ret {0.0};
	if(std::fabs(P - P0) < dados.tolerancia)
		ret = 0.0;
	else
		ret = ((1/bo)-(1/bo0))/(P-P0);
	return ret;
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
void rel_perm(){
	
	for(int i = 0; i < dados.nx; i++){
		Kro[i] = ca_kro(Sg[i]);
		Krg[i] = ca_krg(Sg[i]);
	}
	//Krg[0]=ca_krg(dados.Sg_W);
	//Kro[0]=ca_kro(dados.Sg_W);
}



double calc_Bg(double P){
	return 0.0282*((dados.T*dados.Z)/P);
}

void gas_prop(Vec1D& B, const Vec1D& P){
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
	auto kroud =Kro[0];
	To[0] =Go*(kroud/Bomed);
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

double DerivadaPressaoCapilar(double p, double S){
	return ((1/p)-(1/CalcularPcgo(S-0.05)))/(0.05);
}

//Pressao capilar
double CalcularPcgo(double Sg){
	double Pct = 1.0; 
	double Sorg = dados.sor;
	double lambda = 2.0; 
	return Pct*(std::pow(((1.0-Sorg)/((1.0-Sg)-Sorg)),(1.0/lambda)));
}
void vetpc(Vec1D& p,const Vec1D& S){
	for (int i=0;i<dados.nx;i++){
		p[i]= CalcularPcgo(S[i]);
	}
}
void vecdPc (){
	for(int i=0;i<dados.nx;i++){
		dPc[i]=DerivadaPressaoCapilar(Pc[i],Sg[i]);
	}
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
		//dBo[i]=1/calc_Bo(p[i]);
		dBo[i]=derivada_B_oleo(calc_Bo(p[i]),p[i]);

		Cop[i] =constante*dBo[i]*(1.0-Sg[i]);
	}		
}

double derivada_B_gas(double b, double p){
	return ((1/b)-(1/calc_Bg(p-10.0)))/(10.0);
}

double derivada_B_oleo(double b, double p){

	double ret;
	double cte = 1.5;
	ret=(b-calc_Bo(p-cte))/cte;
	double bot=b;
	double botmenos=calc_Bo(p-cte);

	return ret;
	
}

//montar vetor de incogintas inicial com 2*nx
void criarvx0 (){
	for (int w=0; w<=2*dados.nx;w=w+2){
		X0[w]=dados.Sg_0;
		X0[w+1]=dados.p_0;
	}
}

void divx (Vec1D& X, Vec1D& A, Vec1D& B){
	int k=0;
	for (int w=0;w<=2*dados.nx;w=w+2){
		A[k]=X[w];
		B[k]=X[w+1];
		k=k+1;
	}

}


void atualvec(Vec1D& x0, Vec1D& X){
	auto w = X.size();
	for (int i=0; i<w; i++){
		x0[i]=X[i];
	}
}

void solucao_sist(std::vector<Vec1D>& A,std::vector<Vec1D>& C, Vec1D& B,Vec1D& D, Vec1D& X, Vec1D& Y){// esse aqui ta funcioan
	ccoefp();
	ccoefS();
	cvecDp();
	//cvecDS();
	auto W = X;          // necessária para estimar o erro de uma iteração a outra
	auto Z = Y;
	int counter = 1; 
	int m=dados.nx;
	bool testesol = false;
	/*gauss_solver(A,B,X);
	W = X;
	termofonte(X);
	gauss_solver(C,D,Y);
	*/
	cvecDp();
	gauss_solver(A,B,X);
	cvecDS();
	gauss_solver(C,D,Y);
	W = X;
	Z = Y;
	while(!testesol && counter<MAX_ITER){
		testesol = true;
		/*gauss_solver(A,B,X);
		termofonte(X);
		gauss_solver(C,D,Y);
		*/	
		cvecDp();
		gauss_solver(A,B,X);
		cvecDS();
		gauss_solver(C,D,Y);
		
		for (int i = 0; i < m; i++){
			auto res = std::fabs(((X[i] - W[i]) / X[i])) <= eps;
			auto ress = std::fabs(((Y[i] - Z[i]) / Y[i])) <= eps;
			testesol = testesol & res & ress;
			W[i] = X[i];
			Z[i] = Y[i];
		}
		counter++;
	}
}

void solucao_sist2(std::vector<Vec1D>& A,std::vector<Vec1D>& C, Vec1D& B,Vec1D& D, Vec1D& X, Vec1D& Y){
	ccoefp();
	ccoefS();
	cvecDp();
	termofonte();
	cvecDS();
	gauss_solver(A,B,X);
	termofonte();
	cvecDS();
	gauss_solver(C,D,Y);
}
void gauss_solver(std::vector<Vec1D>& A, Vec1D& B, Vec1D& X){
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

void calcularpo(std::vector<Vec1D>& A, Vec1D& B, Vec1D& X0, Vec1D& Z, Vec1D& X) { //matrz coefp,vetor D, vetor p entrada, vetor sg, vetor p saida
	auto m = B.size();  // numero de linhas
	auto n = m;         // numero de colunas
	auto Y = X0;          // matriz auxiliar
	auto E = X0;
	
	for (int i = 0; i < m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j==i)
					continue;
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * X[j]);
				X[i] = Y[i]; //escreve em X a estimativa encontrada
			}
	}
}

// Função para calcular o segundo sistema usando informações do primeiro sistema
void calcularSg(std::vector<Vec1D>& A, Vec1D& B, Vec1D& Z0, Vec1D& X, Vec1D& Z) { // matriz coefs,vetor D, vetor sg, vetor p
	auto m = B.size();  // numero de linhas
	auto n = m;         // numero de colunas
	auto Y = Z0;          // matriz auxiliar
	auto E = Z0;

	for (int i = 0; i < m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j==i)
					continue;
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * Z[j]);
				Z[i] = Y[i]; //escreve em X a estimativa encontrada
			}
	}
}

void resolverSistemas(std::vector<Vec1D>& A,std::vector<Vec1D>& C, Vec1D& B,Vec1D& H, Vec1D& X, Vec1D& Z) {
	gauss_solver(A, B, X);
	gauss_solver(C, H, Z);
}


bool is_diagonal_dom(const std::vector<Vec1D>& M){
	auto m=M[0].size(); auto n=M.size();
	if (m != n)
		std::cout << "A matriz teste nao e quadrada" << std::endl;

	for (int i=0; i < m; i++){
		double diag = std::fabs(M[i][i]);
		diagonal[i]=std::fabs(diag);
		double sum = 0.0;
		for (int j = 0; j < n; j++){
			if (i==j)
				continue;
			sum = sum + std::fabs(M[i][j]);
			
		}
		soma[i]=sum;
		diferenca[i]=diagonal[i]-soma[i];
		if (sum > diag){
			std::cout << "A matriz nao e estritamente diagonal dominante" << std::endl;
				std::cout << "linha : " << i << " \n";
			std::cout << "soma : " << sum << " \n";
			std::cout << "diagonal : " << M[i][i] << " \n";
			return false;
		}
	}
	std::cout << "A matriz e estritamente diagonal dominante";
	return true;
}

/*bool Sassenfeld(const std::vector<Vec1D>& M){
	auto m=M[0].size(); auto n=M.size();
	if (m != n)
		std::cout << "A matriz teste nao e quadrada" << std::endl;
	//for (int i=0;i<m;i++){
		
	//}
	int i=0;
	beta[i]=1;
	double diag = std::fabs(M[i][i]);
	double sum = 0.0;
	
	for (int j = 0; j < n; j++){
			if (i==j)
				continue;
			sum = sum + std::fabs(M[i][j]/M[i][i]);
			
		}
	beta[i]=sum;
	for (int i=1; i < m; i++){
		double diag = std::fabs(M[i][i]);
		double sum = 0.0;
		
	for (int j = 0; j < n; j++){
			if (i==j)
				continue;
			sum = sum + std::fabs(M[i][j]/diag)*beta[j];
			
		}
	beta[i]=sum;
	}
	double maior;
	maior =beta[0];
	for (int i=1;i<m;i++){
		if(beta[i]>maior)
			maior =soma[i];
	}
	if (maior<1)
		std::cout << "\n \n gaus converge\n ";
	
	else std::cout << "\n gaus NAO converge\n ";
	return true;
		
}*/
bool Sassenfeld(const std::vector<std::vector<double>>& matrix) {
	int size = matrix.size();
	std::vector<double> beta(size); // Vetor para armazenar os valores de beta

	// Inicialização do vetor beta com zeros
	for (int i = 0; i < size; i++) {
		beta[i] = 0.0;
	}

	// Cálculo dos valores de beta
	for (int i = 0; i < size; i++) {
		double sum = 0.0;
		for (int j = 0; j < size; j++) {
		    if (j != i) {
		        sum += std::abs(matrix[i][j]) * beta[j];
		    }
		}
		beta[i] = sum / std::abs(matrix[i][i]);
	}

	// Verificação da convergência
	double maxBeta = beta[0];
	for (int i = 1; i < size; i++) {
		if (beta[i] >= 1) {
			std::cout << "A matriz nao e estritamente diagonal dominante" << std::endl;
		    return false; // Não converge
		}
		if (beta[i] > maxBeta) {
		    maxBeta = beta[i];
		}
	}
	std::cout << "A matriz e estritamente diagonal dominante";
	return true; // Converge
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
	//W.resize(n,0.0); M.resize(n,0.0);
	//E.resize(n,0.0);
	Dp.resize(n,0.0);
	DS.resize(n,0.0);
	RESp.resize(n,0.0); RESs.resize(n,0.0);
	p.resize(n,dados.p_0); p_old.resize(n,dados.p_0);
	p_it.resize(n,dados.p_0); Sg.resize(n,dados.Sg_0);
	Sg_old.resize(n,dados.Sg_0); Sg_it.resize(n,dados.Sg_0);
	Cop.resize(n,0.0); Cog.resize(n,0.0);
	Cgg.resize(n,0.0); Cgp.resize(n,0.0);
	Kro.resize(n,0.0); Krg.resize(n,0.0);
	Bo.resize(n,0.0); Bo_old.resize(n,0.0);
	To.resize(n+1,0.0); Tg.resize(n+1,0.0);
	Pc.resize(n+1,0.0);
	X0.resize(2*n,0.0); X.resize(2*n,0.0);

	int nlines= coefp.size();
	int ncol= coefp[0].size();
	int newsize=2*n;
	if (nlines != ncol)
		std::cout<< "A matriz de coeficientes nao e quadrada" << std::endl;
	else
		{
		if (nlines != newsize){
			std::cout << "Redimensionando a matriz de coeficientes para " << newsize << "x" << newsize << "("<< n <<" celulas lidas no arquivo de entrada)" << std::endl;
			coefp.resize(newsize);
				for (int i = 0; i < newsize; ++i)
					coefp[i].resize(newsize,0.0);
		}
	}
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
///montar matriz linha por linha
void mmll(std::vector<Vec1D>& A, Vec1D& B, double i){
	auto m = B.size(); // 
	A[0][i]=i;
	for(int j= 1; j < m; j++){
	A[j][i]=B[j-1];
	}
}

//Funções utilitarias:
void print_array_2D(const std::vector<Vec1D>& M){
	auto nrow=M.size();
	auto ncol=M[0].size();
	std::cout << "Exibindo matriz 2D (" << nrow << "x" << ncol << ")"<< std::endl;
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			std::cout << M[i][j] << ' ';
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


//As funcões a seguir são do código referência e não serão mais usadas:
/*
//Inicialização dos vetores utilizados no processo de resolução numérica:
Vec1D W(dados.nx, 0.0);                 // n elementos, inicializados em 0
Vec1D M(dados.nx, 0.0);                 // n elementos, inicializados em 0
Vec1D E(dados.nx, 0.0);                 // n elementos, inicializados em 0
Vec1D D(dados.nx, 0.0);                 // n elementos, inicializados em 0

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
void solver_s(){
	auto nx=dados.nx;
	Sg[0] = Sg_old[0] + (1.0/Cgg[0]) * ( Tg[1]*(p[1] - p[0]) + Tg[0]*(dados.p_W - p[0]));
	for(int i = 1; i < nx-1; i++){
		Sg[i] = Sg_old[i] + (1.0/Cgg[i]) * ( Tg[i+1]*(p[i+1] - p[i]) + Tg[i]*(p[i-1] - p[i]) );
	}
	Sg[nx-1] = Sg_old[nx-1] + (1.0/Cgg[nx-1]) * ( Tg[nx]*(dados.p_E - p[nx-1]) + Tg[nx-1]*(p[nx-2] - p[nx-1]));
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

	*/
