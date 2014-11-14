#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>     /* atoi */


using namespace std;


const double pressure0=100000; //in Pa
const double gasConst=8.3144621; //J/(K*mol)
const double gravitationalConst=9.81; //m/s^2
const double kappaConst=2./7.; 
const double cpConst=1005;
const double molarMassDryAir= 28.97; // g/mol
const double z_scale=8000;


double calculatePotTemp(double temp, double pressure);
void print_vec(vector<double> vec);
double calcDeltaTempSolarSurfaceHeating(double deltaTime,double eSolarSurface,double deltaPressure);
void calcDeltaTempVec(double deltaTime,vector<double> Edelta,double deltaPressure,vector<double> deltaTempVec);
void radiativeTransfer(int nlayer,int nmu, vector<double> opt_thick_vec,double LSurface,vector<double> sb_vec,vector<double>& Edelta );
double planck(double temp, double lambda_dn, double lambda_up );
void planckVec(vector<double> tempVec, double lambda_dn, double lambda_up, vector <double> & bbVec);
void getHeightfromPressure(vector<double> pressureVec, vector<double> heightVec);
void calcOpticalThickFromHeight(double beta0,vector<double> he_vec,vector <double> &opt_thick_vec);
void calcOpticalThickfromPressure( double beta0,double deltaPressure, vector<double> &opt_thick_vec);
void calcOpticalThickfromPressure( double beta0,double deltaPressure, vector<double> &opt_thick_vec);

int main(int argc, char** args){
	
	int nLayer,nMu;
	double e0=1362; // W/m^2
	double albedo=0.3;
	double eSolarSurface=e0/4 *(1-albedo);
	double deltaTime=60; //in s
	double beta0=0.0019;
	vector<double> mid_pressure_vec;
	vector<double> mid_temp_vec;
	vector<double> height_vec;
	vector<double> opt_thick_vec;
	vector<double> plankRadiation_vec;
	vector<double> Edelta;
	double deltaPressure;
	double LSurface; 
	if (argc!=3){ 
		cout<<"Try again -- Usage: nLayer nMu"<<endl; 
		return 0;
	}	
	else {
		nLayer= atoi(args[1]);
		nMu=atoi(args[2]);
		}
	
	
	
	
	mid_pressure_vec.resize(nLayer-1);
	mid_temp_vec.resize(nLayer-1);
	opt_thick_vec.resize(nLayer-1);
	plankRadiation_vec.resize(nLayer-1);


	fill(mid_temp_vec.begin(),mid_temp_vec.end(),0);
	
	
	deltaPressure=pressure0/(nLayer-1);
	for(int i=0;i<mid_pressure_vec.size();i++){
		mid_pressure_vec[i]=i*deltaPressure +deltaPressure/2.;
		//cout<<mid_pressure_vec[i]<<endl;
	}
	
	getHeightfromPressure(mid_pressure_vec,height_vec);
		
	calcOpticalThickfromPressure(beta0,deltaPressure,opt_thick_vec);


	//mid_temp_vec.end()+=calcDeltaTempSolarSurfaceHeating(deltaTime,eSolarSurface,deltaPressure);
	

	//The surface Temp needs to be adapted. (At least I think so...) Maybe its temp_vec[nlayer-1] ?
	LSurface=planck(288.15,12,100)/M_PI;
	
	
	//Still have to create the temp_vec. 
	//Probably a loop over the  lambda seqments.
	planckVec(mid_temp_vec, 12, 100,plankRadiation_vec);
	
	radiativeTransfer(nLayer, nMu, opt_thick_vec,LSurface,plankRadiation_vec, Edelta );
  	

	
	
	
	return 0;
}


void getHeightfromPressure(vector<double> pressureVec, vector<double> heightVec){
	heightVec.resize(pressureVec.size());
	for(int i=0;i<pressureVec.size();i++){
		heightVec[i]= -8000*log(pressureVec[i]/100000);  // 100.000 Pa Surface pressure
		
	//if (heightVec[i]<1e-10) heightVec[i]=0;
	//cout<<i<<"   "<<heightVec[i]<<"           "<< pressureVec[i]<<endl;
	}


}

double getPotTempFromTemp(double temp, double pressure){

	return temp*pow(pressure/pressure0,kappaConst);

}

double getTempFromPotTemp(double potTemp, double pressure){
	return potTemp*pow(pressure/pressure0,(-1)*kappaConst);
}


double calcDeltaTempSolarSurfaceHeating(double deltaTime,double eSolarSurface,double deltaPressure){
	
	return deltaTime*eSolarSurface*gravitationalConst/deltaPressure/cpConst;

}





void print_vec(vector<double> vec){
  for(int i=0;i<vec.size();i++)	cout<<vec[i]<<endl;
}

void radiativeTransfer(int nlayer,int nmu, vector<double> opt_thick_vec,double LSurface,vector<double> sb_vec,vector<double>& Edelta ){

  vector<double> Edn, Eup, Enet;
  vector<double>L;
  
  L.resize(nlayer+1);
  Edn.resize(nlayer+1);
  Eup.resize(nlayer+1);
  Enet.resize(nlayer+1);
  Edelta.resize(nlayer);
  

  if(Edn.size()!=nlayer+1){
	  Edn.resize(nlayer+1);
	  Eup.resize(nlayer+1);
	  fill(Edn.begin(),Edn.end(),0);
	  fill(Eup.begin(),Eup.end(),0);
	}
 

    for(double mu=1./nmu/2;mu<=1;mu+=1./nmu){
    
		L[0]=0;
		for(int i=0;i<nlayer;i++){
		  L[i+1]=L[i]*exp(-opt_thick_vec[i]/mu)+sb_vec[i]*(1-exp(-opt_thick_vec[i]/mu));
		  Edn[i+1]=Edn[i+1]+L[i+1]*2*M_PI*mu*1./nmu;
		}
		fill(L.begin(),L.end(),sqrt(-1));
		L[nlayer]= LSurface;//planck(288.15,lambda_dn,lambda_up)/M_PI;        //sigma*pow(288.15,4)/M_PI;
	    
		Eup[nlayer]=Eup[nlayer]+L[nlayer]*2*M_PI*mu*1./nmu;
	  
		for(int i=nlayer;i>0;i--){
		  L[i-1]=L[i]*exp(-opt_thick_vec[i-1]/mu)+sb_vec[i-1]*(1-exp(-opt_thick_vec[i-1]/mu));
		  Eup[i-1]=Eup[i-1]+L[i-1]*2*M_PI*mu*1./nmu;
		}
    
    }

    for (int i=0;i<nlayer+1;i++){
      Enet[i]= Edn[i]-Eup[i];
		}

	for (int i=0;i<nlayer;i++){
		Edelta[i]= Enet[i]-Enet[i+1];
	}
	
}

double planck(double temp, double lambda_dn, double lambda_up ){
	double hcc=5.955214e-17;//Wm*2
 	double hc=1.986446e-25;//Jm
  	double k_b=1.38065e-23;// J/K
  	double sum=0;
  	double lambda;
 	for(int i=lambda_dn;i<=lambda_up;i++){
   		lambda=i*1e-6;
    	sum+=2*hcc/pow(lambda,5)/(exp(hc/k_b/temp/lambda)-1)*1e-6;
	}
 	return M_PI* sum; 
}

void planckVec(vector<double> tempVec, double lambda_dn, double lambda_up, vector <double> & bbVec){
  for (int i=0;i<tempVec.size();i++){
  	bbVec[i]=planck(tempVec[i],lambda_dn,lambda_up);
	}
}



void calcOpticalThickfromPressure( double beta0,double deltaPressure, vector<double> &opt_thick_vec){
	for(int i=0;i<opt_thick_vec.size();i++){
	  opt_thick_vec[i]=beta0*z_scale/pressure0*deltaPressure;
	}
}

void calcOpticalThickFromHeight(double beta0,vector<double> he_vec,vector <double> &opt_thick_vec){
	//Skalenhoehe in m
	//opt_thick_vec.resize(nlayer);
	for(int i=0;i<opt_thick_vec.size();i++){
	  opt_thick_vec[i]=beta0* exp(-(he_vec[i]/z_scale))*(he_vec[i]-he_vec[i+1]);
	}
}


void calcDeltaTempVec(double deltaTime,vector<double> Edelta,double deltaPressure,vector<double> deltaTempVec){
		
	for(int i=0;i<Edelta.size();i++){

		deltaTempVec[i]=deltaTime*Edelta[i]*gravitationalConst/deltaPressure/cpConst;


	}

}