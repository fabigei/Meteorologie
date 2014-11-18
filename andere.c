

// Lucas.Hoeppler@physik.uni-muenchen.de 11.11.2014 8:55
/*
pda_uebung_10112014_b.c
*/
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gnuplot_i.h"
 
#define c_p 1005.0           // [J kg^-1 K^-1]
#define g 9.80665            // [m s^-2]
#define R_a 287.0
#define sigma 5.67E-8
#define c 2.99E8         // [m s^-1]
#define h 6.626E-34      // [J s]
#define kB 1.3806E-23    // [J K^-1]
 
// l2-norm as cancellation condition----------------------------------------------
double norm (double* T_old, double* T, int n) {
  int i=0;
  double accum_old = 0.;
  double accum_new=0.;
  for (i = 0; i < n; ++i) {
    accum_old += T_old[i] * T_old[i];
    accum_new += T[i]*T[i];	
  }
  return sqrt(fabs(accum_old-accum_new));
}
 
// Planck function ---------------------------------------------------------------
double B (double T, double lambda_1, double lambda_2)  {
 
  int N_lambda = 1000;
  double B = 0.0;
  double lambda=0.0;
  double dlambda = 0.0;
 
  dlambda = (lambda_2 - lambda_1) / N_lambda;  
 
  for (lambda = lambda_1 + dlambda/2.0; lambda<lambda_2; lambda+=dlambda)
    B +=  2 * h * pow(c,2) / pow(lambda,5) * 1.0 / (exp (h * c / (lambda * kB * T)) - 1.0) * dlambda;
 
  return B;
}
 
// Radiative transfere function --------------------------------------------------
int rt (int nlyr, double B_surf, double *B, double *dtau, double *Eup, double *Edn)  {
 
  int i = 0;
  int Nmu = 10;
  double mu = 1.0;
  double dmu = 0.0;
  double *Lup;
  double *Ldn;
 
  Lup = calloc(nlyr+1, sizeof(double));
  Ldn = calloc(nlyr+1, sizeof(double));
 
  for (i=0; i<=nlyr; i++)  {
    Eup[i] = 0.0;
    Edn[i] = 0.0;
  }
 
  dmu = mu / Nmu;	// angle segments
 
// downward ...
  for (mu=dmu/2.0; mu<1; mu+=dmu)  {
    for (i=0; i<nlyr; i++)  {
      Ldn[i+1] = Ldn[i] * exp(- dtau[i] / mu) + B[i] * (1 - exp(- dtau[i] / mu));
      Edn[i+1] = Edn[i+1] + 2 * M_PI * Ldn[i+1] * mu * dmu;
    }
  }
 
// upward ...
  for (mu=dmu/2.0; mu<1; mu+=dmu)  {
    Lup[nlyr] = B_surf;
    for (i=nlyr; i>0; i--)  {
      Lup[i-1] = Lup[i] * exp(- dtau[i-1] / mu) + B[i-1] * (1 - exp(- dtau[i-1] / mu));
      Eup[i-1] = Eup[i-1] + 2 * M_PI * Lup[i-1] * mu * dmu;
    }
    Eup[nlyr] = Eup[nlyr] + 2 * M_PI * Lup[nlyr] * mu * dmu;
  }
 
  
 
  free(Lup);
  free(Ldn);
 
  return 0; 
}
 
 
// main function -----------------------------------------------------------------
int main()  {
 
  int nlyr = 100;
  int i = 0; int j = 0;
  int bands = 3;
  int plotcounter=0;
 
  double T_start = 100.0;          	// [K]
  double kappa = 0.2857;           	// [1]
  double deltap = 0.0;             	// [Pa]
  double p_surf = 100000.0;        	// [Pa]
  double p_TOA = 10.0;             	// [Pa]
  double deltat = 600.0;           	// [s]
  double E_s = 235.0;              	// [W m^-2]
  double z_s = 8000.0;			// [m]
  double time = 0.0;			// [s]
  double heating_max=0;			// [K/s]
 
  double residium=2; double value_limit=0.01; // value_limit [K] whole atmosphere
 
  gnuplot_ctrl *g1;
  g1 = gnuplot_init();
 
  double **B_bands = calloc(bands, sizeof(double*));
  for (i=0; i < bands;i++)
    B_bands[i] = calloc (nlyr, sizeof(double));
 
  double **beta = calloc(bands, sizeof(double*));
  for (i=0; i < bands;i++)
    beta[i] = calloc (nlyr, sizeof(double));
 
  double **dtau = calloc(bands, sizeof(double*));
  for (i=0; i < bands;i++)
    dtau[i] = calloc (nlyr, sizeof(double));
 
  double *z_mean  = calloc(nlyr,    sizeof(double));
  double *z       = calloc(nlyr+1,  sizeof(double));
  double *p_mean  = calloc(nlyr,    sizeof(double));
  double *T       = calloc(nlyr,    sizeof(double));
  double *p       = calloc(nlyr+1,  sizeof(double));
  double *T_a     = calloc(nlyr,    sizeof(double));
  double *deltaz  = calloc(nlyr,    sizeof(double));
  double *B_surf  = calloc(bands,   sizeof(double));
  double *Eup     = calloc(nlyr+1,  sizeof(double));
  double *Edn     = calloc(nlyr+1,  sizeof(double));
  double *lambda  = calloc(bands+1, sizeof(double));
  double *tau     = calloc(bands,   sizeof(double));
  double *Eup_tot = calloc(nlyr+1,  sizeof(double));
  double *Edn_tot = calloc(nlyr+1,  sizeof(double));
  double *beta_0  = calloc(bands,   sizeof(double));
  double *E_net   = calloc(nlyr+1,  sizeof(double));
  double *deltaE  = calloc(nlyr,    sizeof(double));
  double *deltaT  = calloc(nlyr,    sizeof(double));
  double *T_old	  = calloc(nlyr,    sizeof(double));
 
// spectrum range. Maximal absorption of CO2: 15 mu m
  lambda[0] = 1.0E-6;
  lambda[1] = 8.0E-6;
  lambda[2] = 12.0E-6;
  lambda[3] = 500.0E-6;
 
  tau[0] = 10.0;
  tau[1] = 0.5;
  tau[2] = 1.0;
 
 
  deltap = (p_surf - p_TOA) / nlyr;
 
  for (i=0; i<=nlyr; i++)
    p[i] = p_TOA + i * deltap;
 
  for (i=0; i<nlyr; i++)  {
    T[i] = T_start;			// temperature each level
    T_old[i]=-1; // 
    p_mean[i] = (p[i] + p[i+1]) / 2.0; 	// pressure each layer
  }
 
// cancellation condition: residium > value_limit
// value_limit: 0.01 -> effective error 0.01/100 [K] each level
  while (residium > value_limit){
 
    time += deltat; 
 
    for (i=0; i<nlyr; i++)
      deltaz[i] = deltap * R_a * T[i] / p_mean[i] / g;
    for (i=nlyr-1; i>0; i--)
      z[i-1] = z[i] + deltaz[i];
 
// integration of all wavelength bands
    for (i=0; i<bands; i++)  {
      beta_0[i] = tau[i] / -z_s / (exp(-z[0] / z_s) - 1.0);
      B_surf[i] = B(T[nlyr-1], lambda[i], lambda[i+1]);
      for (j=0; j<nlyr; j++) {
        	B_bands[i][j] = B(T[j], lambda[i], lambda[i+1]);
        	beta[i][j] = beta_0[i] * exp(- (z[j] + z[j+1]) / 2.0 / z_s);
        	dtau[i][j] = beta[i][j] * deltaz[j];
      }
      rt(nlyr, B_surf[i], B_bands[i], dtau[i], Eup, Edn);
      for (j=0; j<=nlyr; j++){
	       Edn_tot[j] += Edn[j];
	       Eup_tot[j] += Eup[j];
      }
    }
 
    for (i=0; i<=nlyr; i++)  {
      E_net[i] = Edn_tot[i] - Eup_tot[i];
  //  printf(" i%d Edn_tot%f Eup_tot%f E_net%f\n", i, Edn_tot[i], Eup_tot[i], E_net[i]);
    }
  sleep(1);
 // Heating rate each layer
    for (i=0; i<nlyr; i++)  {
      deltaE[i] = E_net[i] - E_net[i+1];
      
      deltaT[i] = g / c_p * deltaE[i] / deltap; //where delta T [K/s]
    }

// Heating rate at the surface
    deltaT[nlyr-1] += (E_s + Edn_tot[nlyr] - Eup_tot[nlyr]) * g / deltap / c_p;  
    
    for (i=0; i<nlyr; i++)  {
      printf("i %d deltaT %f \n", i,deltaT[i]);
    }


// How strong is the strongest heating?
    heating_max=0; // where heating_max [K/s]
    for (i=0; i<nlyr; i++){
      if (fabs(deltaT[i])>heating_max){
	       heating_max=fabs(deltaT[i]);
      }       
    }
 
    deltat=2./heating_max; //Effective duration that layer heats 2K up.
 
// Convert heating rate -> temperature
    for (i=0; i<nlyr; i++){
      T[i]+=deltaT[i]*deltat; //[K]
    }
 
 
//Convection
    for (i=0; i<nlyr-1; i++)  {
      T_a[i] = T[i+1] * pow(p_mean[i] / p_mean[i+1], kappa);
 
      if (T_a[i] > T[i])  {
	T[i] = (T[i] + T[i+1])/2.0 - (T[i+1] - T_a[i]) / 2.0;
	T[i+1] = (T[i] + T[i+1])/2.0 + (T[i+1] - T_a[i]) / 2.0;
      }
    }
 
    for (i=0; i<=nlyr; i++)  {
      Eup_tot[i] = 0.0;
      Edn_tot[i] = 0.0;
    }
 
 
// gnuplot plot
    plotcounter+=1;
 
    if (plotcounter > 10) {
    gnuplot_resetplot  (g1);
    gnuplot_plot_xy   (g1, T, z, nlyr, "Temperature") ;
      plotcounter=0;
    }
 
// calculate residium via l2 norm
    residium=norm(T_old,T,nlyr);
    printf("t %.0f [s]		residium %f [K]	\n", time, residium);
 
// Update -----------------------------------   
    for (i=0; i<nlyr; i++){ 
      T_old[i]=T[i];
    }
 
  }
 
 
// Write values into data "file"
  FILE* datei;
  datei = fopen("file","w");
 
  for (i=0;i<nlyr;i++)   {
    z_mean[i] = (z[i] + z[i+1]) / 2.0;
    fprintf(datei,"%f	%f\n",T[i],z_mean[i]);
  }
 
  gnuplot_resetplot  (g1);
  gnuplot_plot_xy   (g1, T, z, nlyr, "Temperature") ;
 
  return 0;
 
} // main


