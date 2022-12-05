#include <iostream>
#include <cmath>
#include <vector>
#include <string>

using VecDouble = std::vector<double>;
//Function prototypes:
double calc_krw(double Sw);
double calc_kro(double Sw);
double calc_phi(double P);
double calc_Bg(double P);
double calc_Bw(double P);
double calc_Bo(double P);
void calc_C(int n, VecDouble& Cop, VecDouble& Cow, VecDouble& Cwp, VecDouble& Cww, const VecDouble& Bo1, const VecDouble& Bo2, const VecDouble& Bw1, const VecDouble& Bw2, const VecDouble& phi2, const VecDouble& Sw2);

//Variáveis:
int nx = 20;              //descrição
int n;                    //descrição
int i;                    //descrição
int j;                    //descrição
int k;                    //descrição
int total_iter1;          //descrição
int total_iter2;          //descrição
double tol = 0.0001;      //descrição
double kx = 0.25;         //descrição
double viscg = 0.00002;   //viscosidade do gás
double visco = 0.007;     //viscosidade do óleo
double viscw = 0.005;     //viscosidade da água
double kro;               //permeabilidade relativa (óleo)
double krg;               //permeabilidade relativa (gás)
double krw;               //descrição
double phi = 0.20;        //descrição
double cg;                //descrição
double tmax = 100.0;      //descrição
double tempo = 0.0;       //descrição
double deltat = 0.05;     //descrição
double alphac = 1.0;      //descrição
double betac = 0.0000864; //fator de conversão
double Lx = 4000.0;       //Dimensão no eixo x
double Ly = 20.0;         //Dimensão no eixo y
double Lz = 10.0;         //Dimensão no eixo z
constexpr double Ax {Ly*Lz};                //descrição
constexpr double deltax {Lx/nx};            //descrição
constexpr double Vb {deltax*Ax};                //descrição
double p_init = 25000.0;  //descrição
double S_init = 0.2;      //descrição
double erro1;             //descrição
double erro2;             //descrição
double Bo = 1.1;          //descrição
double p_esq = 45000.0;   //pressão no lado esquerdo
double p_dir = 25000.0;   //pressão no lado direito
double S_esq = 0.85;      //descrição
double p_aux;             //descrição
double T_sc = 293.15;     //descrição
double p_sc = 101.325;    //descrição
double T = 400.0;         //
constexpr double der_1_Bg {T_sc/(p_sc*T)};  //descrição
double p0 = 25000.0;      //descrição
double Bo0 = 1.2;         //
double Bw0 = 1.1;         //descrição
double co = 0.00000005;   //descrição
double cw = 0.000000001;  //descrição
double phi0 = 0.2;        //descrição
double cphi = 0.00000001; //descrição
constexpr double um_Bo_linha {Bo0/co};       //descrição
constexpr double um_Bw_linha {Bo0/cw};       //descrição
constexpr double phi_linha {phi0*cphi};      //descrição
//Outras variáveis
double Sor=0.15;          //Saturação do óleo ?
double Sgc=0.1;           //Saturação do gás ?
double nw=2.0;            //descrição
double no=2.0;            //descrição
double krw_max=0.80;      //descrição
double kro_max=0.60;      //descrição
double Swi=0.2;           //descrição
double Swr=0.2;           //descrição
double alpha=-0.01;       //descrição

int main(int argc,char *argv[]){

	//Vetores utilizados
	VecDouble Cog (nx+2, 0.0);
	VecDouble Cgp (nx+2, 0.0);
	VecDouble Cop (nx+2, 0.0);
	VecDouble Cwp (nx+2, 0.0);
	VecDouble Cow (nx+2, 0.0);
	VecDouble Cww (nx+2, 0.0);
	VecDouble Sg1 (nx+2, 0.0);
	VecDouble Sg2 (nx+2, 0.0);
	VecDouble Sw1 (nx+2, S_init);
	VecDouble Sw2 (nx+2, S_init);
	VecDouble po1 (nx+2, p_init);
	VecDouble po2 (nx+2, p_init);
	VecDouble po3 (nx+2, p_init);
	VecDouble Bg1 (nx+2, 0.0);
	VecDouble Bg2 (nx+2, 0.0);
	VecDouble Bw1 (nx+2, 0.0);
	VecDouble Bw2 (nx+2, 0.0);
	VecDouble Bo1 (nx+2, 0.0);
	VecDouble Bo2 (nx+2, 0.0);
	VecDouble phi1 (nx+2, 0.0);
	VecDouble phi2 (nx+2, 0.0);
	VecDouble Tg (nx+2, 0.0);
	VecDouble To (nx+2, 0.0);
	VecDouble Tw (nx+2, 0.0);
	VecDouble qo (nx+2, 0.0);
	VecDouble qw (nx+2, 0.0);
	//Verificar as saturações em 0, 1 e nx
	//Verificar po's em 0


	//Time iteration
	while(tempo<=tmax){
		tempo+=deltat;
		printf("%lf\n\n",tempo);

		for(i=1;i<=nx;i++){
			Bo1[i]=calc_Bo(po1[i]);
			Bw1[i]=calc_Bw(po1[i]);
			phi1[i]=calc_phi(po1[i]);
		}
		//printf("here 1\n");
		total_iter1=0;

		while(total_iter1<=1000){
			//printf("here 2\n");
			total_iter1++;
			erro1=-1.0;
			for(n=1;n<=nx;n++){
				Bo2[n]=calc_Bo(po2[n]);
				Bw2[n]=calc_Bw(po2[n]);
				phi2[n]=calc_phi(po2[n]);
				calc_C(n,Cop,Cow,Cwp,Cww,Bo1,Bo2,Bw1,Bw2,phi2,Sw2);
			}
			//printf("here 3\n");
			i=0;
			krw=calc_krw(S_esq);
			kro=calc_kro(S_esq);
			To[i]=(2.0*betac*Ax*kx*kro)/(visco*deltax*calc_Bo(p_esq));
			Tw[i]=(2.0*betac*Ax*kx*krw)/(viscw*deltax*calc_Bw(p_esq));
			//printf("here 4\n");
			for(i=1;i<nx;i++){
				krw=calc_krw(Sw2[i]);
				kro=calc_kro(Sw2[i]);
				To[i]=(betac*Ax*kx*kro)/(visco*deltax*0.5*(Bo2[i]+Bo2[i+1]));
				Tw[i]=(betac*Ax*kx*krw)/(viscw*deltax*0.5*(Bw2[i]+Bw2[i+1]));
			}
			//printf("here 5\n");
			i=nx;
			krw=calc_krw(Sw2[i]);
			kro=calc_kro(Sw2[i]);
			To[i]=(2.0*betac*Ax*kx*kro)/(visco*deltax*calc_Bo(p_dir));
			Tw[i]=(2.0*betac*Ax*kx*krw)/(viscw*deltax*calc_Bw(p_dir));

			total_iter2=0;

			//printf("here 6\n");
			while(total_iter2<=1000){
				erro2=-1.0;
				total_iter2++;
				//printf("here 7\n");
				i=1;

				p_aux=(((Bo2[i]*Cop[i])+(Bw2[i]*Cwp[i]))*po1[i])+
				(((Bo2[i]*To[i])+(Bw2[i]*Tw[i]))*po3[i+1])+(((Bo2[i]*To[i-1])+(Bw2[i]*Tw[i-1]))*p_esq);
				p_aux/=((((Bo2[i]*Cop[i])+(Bw2[i]*Cwp[i])))+(((Bo2[i]*To[i])+(Bw2[i]*Tw[i])))+(((Bo2[i]*To[i-1])+(Bw2[i]*Tw[i-1]))));
				erro2=maximo(erro2,fabs(p_aux-po3[i]));
				po3[i]=p_aux;

				printf("%.12lf\n",po3[i]);

				for(i=2;i<nx;i++){
					p_aux=(((Bo2[i]*Cop[i])+(Bw2[i]*Cwp[i]))*po1[i])+
					(((Bo2[i]*To[i])+(Bw2[i]*Tw[i]))*po3[i+1])+(((Bo2[i]*To[i-1])+(Bw2[i]*Tw[i-1]))*po3[i-1]);
					p_aux/=((((Bo2[i]*Cop[i])+(Bw2[i]*Cwp[i])))+(((Bo2[i]*To[i])+(Bw2[i]*Tw[i])))+(((Bo2[i]*To[i-1])+(Bw2[i]*Tw[i-1]))));
					erro2=maximo(erro2,fabs(p_aux-po3[i]));
					po3[i]=p_aux;

					//printf("%lf\t%lf\n",po3[i],erro2);
					printf("%.12lf\n",po3[i]);
				}
				i=nx;

				p_aux=(((Bo2[i]*Cop[i])+(Bw2[i]*Cwp[i]))*po1[i])+
				(((Bo2[i]*To[i])+(Bw2[i]*Tw[i]))*p_dir)+(((Bo2[i]*To[i-1])+(Bw2[i]*Tw[i-1]))*po3[i-1]);
				p_aux/=((((Bo2[i]*Cop[i])+(Bw2[i]*Cwp[i])))+(((Bo2[i]*To[i])+(Bw2[i]*Tw[i])))+(((Bo2[i]*To[i-1])+(Bw2[i]*Tw[i-1]))));
				erro2=maximo(erro2,fabs(p_aux-po3[i]));
				po3[i]=p_aux;

				printf("%.12lf\n",po3[i]);
				printf("\n ERRO %lf \n",erro2);
				getchar();

				if(erro2<tol)
					break;
			}
			i=1
			Sw2[i]=Sw1[i]+(((Tw[i]*(po3[i+1]-po3[i]))+(Tw[i-1]*(p_esq-po3[i]))-(Cwp[i]*(po3[i]-po1[i])))/Cww[i]);
			printf("%lf\n",Sw2[i]);

			for(i=2;i<nx;i++){
				Sw2[i]=Sw1[i]+(((Tw[i]*(po3[i+1]-po3[i]))+(Tw[i-1]*(po3[i-1]-po3[i]))-(Cwp[i]*(po3[i]-po1[i])))/Cww[i]);
				printf("%lf\n",Sw2[i]);
			}
			i=nx;
			Sw2[i]=Sw1[i]+(((Tw[i]*(p_dir-po3[i]))+(Tw[i-1]*(po3[i-1]-po3[i]))-(Cwp[i]*(po3[i]-po1[i])))/Cww[i]);
			printf("%lf\n\n",Sw2[i]);
			erro1=-1.0;

			for(i=1;i<=nx;i++){
				erro1=maximo(erro1,fabs(po2[i]-po3[i]));
				po2[i]=po3[i];
			}

			if(erro1<tol){
				printf("\n ERRO!!!\n");
				getchar();
				break;
				}
		}
		for(i=1;i<=nx;i++){
			Sw1[i]=Sw2[i];
			po1[i]=po2[i]=po3[i];
		}
	}

	arq1=fopen("result.dat","w+");
	for(i=1;i<=nx;i++)
		fprintf(arq1,"%lf\t %lf\t %lf\n",(i-0.5)*deltax,po3[i],Sw2[i]);
	fclose(arq1);

	printf("\nSimulacao pronta\n");
}

double calc_Bg(double P){
	return (p_sc*T)/(T_sc*P);
}
double calc_krw(double Sw){
	return krw_max*std::pow((Sw-Swi)/(1.0-Swi-Sor),nw);
}
double calc_kro(double Sw){
	return kro_max*std::pow((1.0-Sw-Swi)/(1.0-Swi-Sor),no);
}
/*double calc_krg(double Sg){
	return krg_max*std::pow((Sg-Sgc)/(1.0-Sgc-Sorg),ng);
}
double calc_kro(double Sg){
	return kro_max*std::pow((1.0-Sg-Sorg)/(1.0-Sorg),nog);
}*/
double calc_Bw(double P){
	return Bw0/(1.0+(cw*(P-p0)));
}
double calc_Bo(double P){
	return Bo0/(1.0+(co*(P-p0)));
}
double calc_phi(double P){
	return phi0*(1.0+(cphi*(P-p0)));
}
void calc_C(int n, VecDouble& Cop, VecDouble& Cow, VecDouble& Cwp, VecDouble& Cww, const VecDouble& Bo1, const VecDouble& Bo2, const VecDouble& Bw1, const VecDouble& Bw2, const VecDouble& phi2, const VecDouble& Sw2){
	Cop[n]=(Vb/(alphac*deltat))*((phi_linha/Bo1[n])+(phi2[n]*um_Bo_linha))*(1.0-Sw2[n]);
	Cow[n]=-(Vb/(alphac*deltat))*(phi2[n]/Bo2[n]); 
	Cwp[n]=(Vb/(alphac*deltat))*((phi_linha/Bw1[n])+(phi2[n]*um_Bw_linha))*(Sw2[n]); 
	Cww[n]=(Vb/(alphac*deltat))*(phi2[n]/Bw2[n]);
}
/*
double calc_krw(double Sw){
	auto S=(Sw-Swi)/(1.0-Swi-Swr);
	return 0.2*S*S*S;
}
double calc_kro(double Sw){
	auto S=(Sw-Swi)/(1.0-Swi-Swr);
	return (1.0-S)*(1.0-S)*(1.0-S);
}
*/