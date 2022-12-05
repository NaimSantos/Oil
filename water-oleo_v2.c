#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>

int nx,n,i,j,k,total_iter1,total_iter2;
double tol, kx, viscg, visco, viscw,kro, krg,  krw,Ax, deltax, Lx, phi, cg, tmax, tempo, deltax,deltat,alphac,betac,Vb,Ly,Lz,p_init, S_init,erro1,erro2;
double Bo, p_esq,p_dir, S_esq,p_aux,T_sc,p_sc,der_1_Bg;
double *Cog, *Cgp, *Cop, *Cwp, *Cow, *Cww, *Sg1, *Sg2, *Sg3, *po1, *po2, *po3, *Bg1, *Bg2, *Tg, *To, *Tw, T, *Sw1, *Sw2, *Bw1, *Bw2, *Bo1, *Bo2;
double *qw, *qo, *phi1, *phi2;
double phi_linha,um_Bo_linha,um_Bw_linha;
double p0, Bo0, co, Bw0, cw, phi0, cphi;

FILE *arq1;

double maximo(double a,double b);

double calc_krw(double Sw);
double calc_kro(double Sw);
double calc_phi(double P);

double calc_Bg(double P);
double calc_Bw(double P);
double calc_Bo(double P);
void calc_C();

double Sor=0.15,Sgc=0.1,nw=2.0,krw_max=0.80,no=2.00,kro_max=0.60; 
double Swi=0.2, Swr=0.2;

int main(int argc,char *argv[])
{
  Bo0=1.2;
  co=0.00000005;
  Bw0=1.1;
  cw=0.000000001;
  phi0=0.2;
  cphi=0.00000001;
  T=400.0;
  p_sc=101.325;
  T_sc=293.15;
  nx=20;
  tol=0.0001; 
  kx=0.25; 
  viscg=0.00002; 
  visco=0.007; 
  viscw=0.005; 
  Lx=4000.0; 
  phi=0.20;
  tmax=100.0; 
  deltat=0.05;
  alphac=1.0;
  betac=0.0000864;
  tempo=0.0;
  Lz=10.0;
  Ly=20.0;
  Ax=Ly*Lz;
  deltax=Lx/nx;
  Vb=deltax*Ax;
  der_1_Bg=T_sc/(p_sc*T);
  Bo=1.1;
  double alpha=-0.01;

  p0=p_init=25000.0;
  S_init=0.2;
  p_esq=45000.0;
  p_dir=25000.0;
  S_esq=0.85;

  Cog=calloc(nx+2,sizeof(double));
  Cgp=calloc(nx+2,sizeof(double));
  Cop=calloc(nx+2,sizeof(double));
  Cwp=calloc(nx+2,sizeof(double));
  Cow=calloc(nx+2,sizeof(double));
  Cww=calloc(nx+2,sizeof(double));
  Sg1=calloc(nx+2,sizeof(double));
  Sg2=calloc(nx+2,sizeof(double));
  Sw1=calloc(nx+2,sizeof(double));
  Sw2=calloc(nx+2,sizeof(double));
  po1=calloc(nx+2,sizeof(double));
  po2=calloc(nx+2,sizeof(double));
  po3=calloc(nx+2,sizeof(double));
  Bg1=calloc(nx+2,sizeof(double));
  Bg2=calloc(nx+2,sizeof(double));
  Bw1=calloc(nx+2,sizeof(double));
  Bw2=calloc(nx+2,sizeof(double));
  Bo1=calloc(nx+2,sizeof(double));
  Bo2=calloc(nx+2,sizeof(double));
  phi1=calloc(nx+2,sizeof(double));
  phi2=calloc(nx+2,sizeof(double));
  Tg=calloc(nx+2,sizeof(double));
  To=calloc(nx+2,sizeof(double));
  Tw=calloc(nx+2,sizeof(double));
  qo=calloc(nx+2,sizeof(double));
  qw=calloc(nx+2,sizeof(double));

  for(i=1;i<=nx;i++)
  {
   po1[i]=po2[i]=po3[i]=p_init;
   Sw1[i]=Sw2[i]=S_init;
  }

  Sw1[i]=Sw2[i]=S_init;

  while(tempo<=tmax)
  {
    tempo+=deltat;

    printf("%lf\n\n",tempo);

    for(i=1;i<=nx;i++)
    {//printf("here 0\n");
      Bo1[i]=calc_Bo(po1[i]);//printf("here 0 a\n");
      Bw1[i]=calc_Bw(po1[i]);//printf("here 0 b\n");
      phi1[i]=calc_phi(po1[i]); //printf("here 0 c\n");
    }
//printf("here 1\n");
    total_iter1=0;

    while(total_iter1<=1000)
    {//printf("here 2\n");
      total_iter1++;
      erro1=-1.0;      

      for(n=1;n<=nx;n++)
      {
        Bo2[n]=calc_Bo(po2[n]);
        Bw2[n]=calc_Bw(po2[n]);
        phi2[n]=calc_phi(po2[n]);
        calc_C();
      }
//printf("here 3\n");
      i=0;
      krw=calc_krw(S_esq);
      kro=calc_kro(S_esq);
      To[i]=(2.0*betac*Ax*kx*kro)/(visco*deltax*calc_Bo(p_esq));
      Tw[i]=(2.0*betac*Ax*kx*krw)/(viscw*deltax*calc_Bw(p_esq));
 //printf("here 4\n");
      for(i=1;i<nx;i++)
      {
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
      while(total_iter2<=1000)
      {
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

        for(i=2;i<nx;i++)
        {
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

        if(erro2<tol) break;
      }

      i=1;

      Sw2[i]=Sw1[i]+(((Tw[i]*(po3[i+1]-po3[i]))+(Tw[i-1]*(p_esq-po3[i]))-(Cwp[i]*(po3[i]-po1[i])))/Cww[i]);
 
 printf("%lf\n",Sw2[i]);
 
      for(i=2;i<nx;i++)
      {
        Sw2[i]=Sw1[i]+(((Tw[i]*(po3[i+1]-po3[i]))+(Tw[i-1]*(po3[i-1]-po3[i]))-(Cwp[i]*(po3[i]-po1[i])))/Cww[i]);
 printf("%lf\n",Sw2[i]);
      }

      i=nx;

      Sw2[i]=Sw1[i]+(((Tw[i]*(p_dir-po3[i]))+(Tw[i-1]*(po3[i-1]-po3[i]))-(Cwp[i]*(po3[i]-po1[i])))/Cww[i]);
      
       printf("%lf\n\n",Sw2[i]);
      
      erro1=-1.0;

      for(i=1;i<=nx;i++)
      {
       erro1=maximo(erro1,fabs(po2[i]-po3[i]));
       po2[i]=po3[i];
      }

      if(erro1<tol)
      {
        printf("\n ERRO!!!\n");
        getchar();
        break;
      } 
    }
    
    for(i=1;i<=nx;i++)
    {
      Sw1[i]=Sw2[i];
      po1[i]=po2[i]=po3[i];
    }
  }

  arq1=fopen("result.dat","w+");
  for(i=1;i<=nx;i++) 	  fprintf(arq1,"%lf\t %lf\t %lf\n",(i-0.5)*deltax,po3[i],Sw2[i]);
  fclose(arq1);

  free(Cog);
  free(Cgp);
  free(Cop);
  free(Cwp);
  free(Cow);
  free(Cww);
  free(Sg1);
  free(Sg2);
  free(Sw1);
  free(Sw2);
  free(po1);
  free(po2);
  free(po3);
  free(Bg1);
  free(Bg2);
  free(Bo1);
  free(Bo2);
  free(Bw1);
  free(Bw2);
  free(phi1);
  free(phi2);
  free(Tg);
  free(To);
  free(Tw);
  free(qo);
  free(qw);
 
  printf("\nSimulacao pronta\n");
}

double maximo(double a,double b)
{
  return (a>b)? a : b;
} 

double calc_Bg(double P)
{
  return (p_sc*T)/(T_sc*P);
}  

/*double calc_krg(double Sg)
{
  return krg_max*pow((Sg-Sgc)/(1.0-Sgc-Sorg),ng);                   
}

double calc_kro(double Sg)
{
  return kro_max*pow((1.0-Sg-Sorg)/(1.0-Sorg),nog);                   
}*/

double calc_krw(double Sw)
{
  return krw_max*pow((Sw-Swi)/(1.0-Swi-Sor),nw);                   
}

double calc_kro(double Sw)
{
  return kro_max*pow((1.0-Sw-Swi)/(1.0-Swi-Sor),no);                   
}

double calc_Bw(double P)
{
  return Bw0/(1.0+(cw*(P-p0)));			
}

double calc_Bo(double P)
{
  return Bo0/(1.0+(co*(P-p0)));		
}

double calc_phi(double P)
{
  return phi0*(1.0+(cphi*(P-p0)));	
}

void calc_C()
{
  phi_linha=phi0*cphi;
  um_Bo_linha=Bo0/co;
  um_Bw_linha=Bo0/cw;
  Cop[n]=(Vb/(alphac*deltat))*((phi_linha/Bo1[n])+(phi2[n]*um_Bo_linha))*(1.0-Sw2[n]);
  Cow[n]=-(Vb/(alphac*deltat))*(phi2[n]/Bo2[n]); 
  Cwp[n]=(Vb/(alphac*deltat))*((phi_linha/Bw1[n])+(phi2[n]*um_Bw_linha))*(Sw2[n]); 
  Cww[n]=(Vb/(alphac*deltat))*(phi2[n]/Bw2[n]);
}

//double calc_krw(double Sw)
//{
//  double S;
//  S=(Sw-Swi)/(1.0-Swi-Swr);
//  return 0.2*S*S*S;                   
//}
//
//double calc_kro(double Sw)
//{
//  double S;
//  S=(Sw-Swi)/(1.0-Swi-Swr);
//  return (1.0-S)*(1.0-S)*(1.0-S);                   
//}
