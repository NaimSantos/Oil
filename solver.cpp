#include <iostream>
#include <vector>
#include <cmath>      //std::fabs, std::pow
#include <algorithm>  //std::max, std::max_element, std::min

using Vec1D = std::vector<double>;

//Protótipos:
void update(Vec1D& V0, const Vec1D& V);
void rel_diff(Vec1D& dif, const Vec1D& V, const Vec1D& V0);
void abs_diff(Vec1D& dif, const Vec1D& V, double *V0);
double max_vetor(Vec1D& V);
void acumulo(Vec1D& Cop, Vec1D& Cow, Vec1D& Cww, const Vec1D& Bo, const Vec1D& Bo_old, const Vec1D& p, const Vec1D& p_old);
double evaluate_dB(double bo, double bo0, double P, double P0);
void oil_prop(Vec1D& B, const Vec1D& P);
double calc_Bo(const double P);
void rel_perm(const Vec1D& Sw, Vec1D& Kro, Vec1D& Krw);
double calc_kro(double sw);
double calc_krw(double sw);
int upwind(double pesq, double pdir);
void Transmissibilidade(const Vec1D& Kro, const Vec1D& Krw, const Vec1D& To, const Vec1D& Tw);
double dB(double bo, double bo0, double P, double P0)
void acumulo(const Vec1D& Bo, const Vec1D& Bo_old, const Vec1D& p, const Vec1D& p_old, const Vec1D& Sw_old, Vec1D& Cop, Vec1D& Cow, Vec1D& Cww);
void print_field_x(Vec1D&, double DX, int n);
//Váriáveis da simulação
double tempo_max {0.0};
double dt {0.5};
double Lx {0.1};
int nx {100};
int sim_type {1};
double tol {0.0001};
double dx {0.1};
double p_0 {0.1};
double p_W {0.1};
double p_E {0.1};
double Sw_0 {0.1};
double Sw_W {0.1};
double phi {0.1};
double k {1000};
double siw {0.1};
double sor {0.1};
double rhowsc {1000.0};
double rhow {1000.0};
double muw {1000.0};
double Boref {1000};
double co {10};
double muo {50};
double pref {50};
double Bw {50};

int main(int argc, char *argv[] ){
	//Leitura das váriaveis a partir de um arquivo de entrada:

	//Inicialização de vetores e variáveis:
	Vec1D p(nx, p_0);                  // nx elementos, inicializados com p_0
	Vec1D p_old(nx, p_0);              // nx elementos, inicializados com p_0
	Vec1D p_it(nx, p_0);               // nx elementos, inicializados com p_0s
	Vec1D Sw(nx, Sw_0);                // nx elementos, inicializados com Sw_0
	Vec1D Sw_old(nx, Sw_0);            // nx elementos, inicializados com Sw_0
	Vec1D Sw_it(nx, Sw_0);             // nx elementos, inicializados com Sw_0
	Vec1D Cop(nx, 0.0);
	Vec1D Cow(nx, 0.0);
	Vec1D Cww(nx, 0.0);
	Vec1D Kro(nx, 0.0);
	Vec1D Krw(nx, 0.0);
	Vec1D Bo(nx, 0.0);
	Vec1D Bo_old(nx, 0.0);
	Vec1D To(nx+1, 0.0);               //Transmissibilidade
	Vec1D Tw(nx+1, 0.0);               //Transmissibilidade

	Bw = rhowsc/rhow;
	dx = Lx/nx;
	pref = p_0;
	

	return 0;
}

void update(Vec1D& V0, const Vec1D& V){
	constexpr auto n = dif.size();
	for(int i = 0; i < n; i++){
		V0[i] = V[i];
	}
}
void rel_diff(Vec1D& dif, const Vec1D& V, const Vec1D& V0){
	constexpr auto n = dif.size();
	for (int i = 0; i < n; i++){
		dif[i] = std::fabs( (V[i] - V0[i])/V0[i]);
	}
}
void abs_diff(Vec1D& dif, const Vec1D& V, double *V0){
	constexpr auto n = dif.size();
	for (int i = 0; i < n; i++){
		dif[i] = std::fabs(V[i] - V0[i]);
	}
}
double max_vetor(Vec1D& V){
	return std::max_element(V.begin(), V.end());
}
void acumulo(Vec1D& Cop, Vec1D& Cow, Vec1D& Cww, const Vec1D& Bo, const Vec1D& Bo_old, const Vec1D& p, const Vec1D& p_old){
	double v = phi/dt;
	double bl;
	for(int i = 0; i < nx; i++){
		bl = evaluate_dB(Bo[i],Bo_old[i],p[i],p_old[i]);
		Cop[i] = v*(1 - Sw_old[i])*bl;
		Cow[i] = -v/Bo[i];
		Cww[i] = v/Bw;
	}
}
double evaluate_dB(double bo, double bo0, double P, double P0){
	double ret = 0.0;
	if(std::fabs(P - P0) < tol){
		ret = 0.0;
	}
	else{
		ret = ( (1/bo) - (1/bo0) )/(P - P0);
	}
	return ret;
}
void oil_prop(Vec1D& B, const Vec1D& P){
	for(int i = 0; i < nx; i++){
		B[i] = calc_Bo(P[i]);
	}
}
double calc_Bo(const double P){
	return Boref/( 1.0 + co*(P - pref));
}
void rel_perm(const Vec1D& Sw, Vec1D& Kro, Vec1D& Krw){
	for(int i = 0; i < nx; i++){
		Kro[i] = calc_kro(Sw[i]);
		Krw[i] = calc_krw(Sw[i]);
	}
}
/// Modelo de Corey para permeabilidade relativa óleo/água
double calc_kro(double sw){
	double so = 1 - sw;
	double aux {0.0};
	double ret {0.0};
	aux = (sw - siw)/(1 - siw);
	ret = (1 - aux*aux) * (1 - aux)*(1 - aux);
	ret = std::max(0,ret);
	ret = std::min(ret,1.0);
	return ret;
}
double calc_krw(double sw){
	double aux {0.0};
	double ret {0.0}
	aux = (sw - siw)/(1 - siw);
	ret = std::pow(aux,4);
	ret = std::max(0,ret);
	ret = std::min(ret,1.0);
	return ret;
}
int upwind(double pesq, double pdir){
	return (pesq > pdir) ? 1 : 0;
}
void Transmissibilidade(const Vec1D& Kro, const Vec1D& Krw, const Vec1D& To, const Vec1D& Tw){
	double conv = 8.64E-2;
	const auto Gw = (conv*k)/(Bw*muw*dx*dx);
	const auto Go = (conv*k)/(muo*dx*dx);

	auto Bomed = calc_Bo(p_W);
	auto ud = upwind(p_W, p[0]);
	auto kroud = ud*calc_kro(Sw_W) + (1-ud)*Kro[0];
	auto krwud = ud*calc_krw(Sw_W) + (1-ud)*Krw[0];
	To[0] = 2*Go*(kroud/Bomed);
	Tw[0] = 2*Gw*krwud;
	for(int i = 1; i < nx; i++){
		Bomed = (Bo[i-1] + Bo[i])/2;
		ud    = upwind(p[i-1],p[i]);
		kroud = ud*Kro[i-1] + (1-ud)*Kro[i];
		krwud = ud*Krw[i-1] + (1-ud)*Krw[i];
		To[i] = Go*(kroud/Bomed);
		Tw[i] = Gw*krwud;
	}
	Bomed = calc_Bo(p_E);
	/// S_E = Sw[nx-1] (condicao de derivada nula)
	kroud = Kro[nx-1];
	krwud = Krw[nx-1];
	To[nx] = 2*Go*(kroud/Bomed);
	Tw[nx] = 2*Gw*krwud;
}
double dB(double bo, double bo0, double P, double P0){
	double ret {0.0};
	if(std::fabs(P - P0) < tol)
		ret = 0.0;
	else
		ret = ((1/bo)-(1/bo0))/(P-P0);
	return ret;
}
void acumulo(const Vec1D& Bo, const Vec1D& Bo_old, const Vec1D& p, const Vec1D& p_old, const Vec1D& Sw_old, Vec1D& Cop, Vec1D& Cow, Vec1D& Cww){
	const auto v = phi/dt;
	double bl {1.0};
	for(int i = 0; i < nx; i++){
		bl = dB(Bo[i],Bo_old[i],p[i],p_old[i]);
		Cop[i] = v*(1 - Sw_old[i])*bl;
		Cow[i] = -v/Bo[i];
		Cww[i] = v/Bw;
	}
}
void print_field_x(Vec1D&, double DX, int n){

}