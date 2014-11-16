#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>     /* atoi */
#include <algorithm>
#include <unistd.h>
 	


//#include "./gnuplot_i-2.10/src/gnuplot_i.h"


using namespace std;


const double pressure0=100000; //in Pa
const double gasConst=8.3144621; //J/(K*mol)
const double gravitationalConst=9.81; //m/s^2
const double kappaConst=7/2.; 
const double cpConst=1005;
const double molarMassDryAir= 28.97; // g/mol
const double z_scale=8000;


double calculatePotTemp(double temp, double pressure);
void print_vec(vector<double> vec);
double calcDeltaTempSolarSurfaceHeating(double deltaTime,double eSolarSurface,double deltaPressure);
void calcDeltaTempVec(double deltaTime,vector<double> Edelta,double deltaPressure,vector<double> &deltaTempVec);
void radiativeTransfer(int nlayer,int nmu, vector<double> opt_thick_vec,double LSurface,vector<double> sb_vec,vector<double>& Edelta );
double planck(double temp, double lambda_dn, double lambda_up );
void planckVec(vector<double> tempVec, double lambda_dn, double lambda_up, vector <double> & bbVec);
void getHeightfromPressure(vector<double> mid_pressure_vec,vector<double> mid_temp_vec, vector<double>& heightVec);
void calcOpticalThickFromHeight(double beta0,vector<double> he_vec,vector <double> &opt_thick_vec);
void calcOpticalThickfromPressure( double beta0,double deltaPressure, vector<double> &opt_thick_vec);
void calcOpticalThickfromPressure( double beta0,double deltaPressure, vector<double> &opt_thick_vec);
void calcMiddleVec(std::vector<double> inputVec,std::vector<double> & outputVec);
void getPotTempVecFromTempVec(vector<double> TempVec,vector<double> pressureVec,vector<double>&PotTempVec);
void getTempVecFromPotTempVec(vector<double>PotTempVec,vector<double> pressureVec,vector<double> & TempVec);

template <typename T>
void vectorAdd(vector<T> inputVecA,vector<T> inputVecB, vector<T> & outputVec){
	if(inputVecA.size()!=inputVecB.size()) {
		cout<<"Vectors must be of equal length."<<endl;
	}

	else{
		outputVec.resize(inputVecA.size());
		for(int i=0;i<inputVecA.size();i++)
			outputVec[i]=inputVecA[i]+inputVecB[i];
	}
}





int main(int argc, char** args){
	


	 //gnuplot_ctrl *g1;

	 //g1 = gnuplot_init();



	int nLayer,nMu;
	double e0=1362; // W/m^2
	double albedo=0.3;
	double eSolarSurface=e0/4 *(1-albedo);
	double deltaTime=8640; //in s
	double timeStep;
	double MaxTime=8640*2;//3.15569e7; //3.15569e7=1a in seconds
	double beta0=0.0019;
	vector<double> mid_pressure_vec;
	
	vector<double> mid_height_vec;
	vector<double> opt_thick_vec;
	vector<double> plankRadiation_vec;
	vector<double> Edelt_vec;
	vector<double> mid_temp_vec;
	vector<double> middleDeltaTemp_vec;
	vector<double> mid_pot_temp_vec;
	
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
	
	
	
	
	mid_pressure_vec.resize(nLayer);
	mid_temp_vec.resize(nLayer);
	opt_thick_vec.resize(nLayer);
	plankRadiation_vec.resize(nLayer);
	mid_height_vec.resize(nLayer);

	fill(mid_temp_vec.begin(),mid_temp_vec.end(),100);
	
	
	deltaPressure=pressure0/(nLayer-1);
	for(int i=0;i<mid_pressure_vec.size();i++){
		mid_pressure_vec[i]=i*deltaPressure +deltaPressure/2.;
		//cout<<mid_pressure_vec[i]<<endl;
	}
	
	
		
	calcOpticalThickfromPressure(beta0,deltaPressure,opt_thick_vec);

	//mid_temp_vec.end()+=calcDeltaTempSolarSurfaceHeating(deltaTime,eSolarSurface,deltaPressure);
	for(int timeStep=0;timeStep<=MaxTime;timeStep+=deltaTime){

	//The surface Temp needs to be adapted. (At least I think so...) Maybe its temp_vec[nlayer] ?
	LSurface=planck(mid_temp_vec[nLayer-1],4,100)/M_PI;
		

	//Still have to create the temp_vec. 
	//Probably a loop over the  lambda seqments.
	planckVec(mid_temp_vec, 4 , 100,plankRadiation_vec);
	
	radiativeTransfer(nLayer, nMu, opt_thick_vec,LSurface,plankRadiation_vec, Edelt_vec);


	Edelt_vec[nLayer-1]+=eSolarSurface;

	print_vec(Edelt_vec);

	

	calcDeltaTempVec(deltaTime,Edelt_vec,deltaPressure,middleDeltaTemp_vec);
	
	//TODO: *size(Edelt_vec),size(deltaTem_vec),size(mTemp_vec) must fit together
	//calcMiddleVec(deltaTemp_vec,MiddleDeltaTemp_vec);



	vectorAdd(middleDeltaTemp_vec,mid_temp_vec,mid_temp_vec);

	cout<<endl;

	print_vec(mid_temp_vec);

	getPotTempVecFromTempVec(mid_temp_vec, mid_pressure_vec,mid_pot_temp_vec);


	sort(mid_pot_temp_vec.begin(), mid_pot_temp_vec.end(), greater<double>());


	getTempVecFromPotTempVec(mid_pot_temp_vec,mid_pressure_vec,mid_temp_vec);
	
 	cout<<endl<<endl<<endl;
	/*if (int(timeStep/deltaTime)%10 == 0) {
      
		getHeightfromPressure(mid_pressure_vec,mid_temp_vec,mid_height_vec);

      gnuplot_resetplot  (g1);  start with new plot rather than plotting into exisiting one 
      gnuplot_setstyle   (g1, "linespoints");      // draw lines and points 
      gnuplot_set_xlabel (g1, "temperature [K]");   //xaxis label 
      gnuplot_set_ylabel (g1, "pressure Pa");    // yaxis label 
      
     
    	
    	
		

    	double *mid_temp_array=&mid_temp_vec[0];
    	double *mid_pot_temp_array=&mid_pot_temp_vec[0];

    	double *mid_height_array=&mid_height_vec[0];
    	double *mid_pressure_array=&mid_pressure_vec[0];
     	gnuplot_plot_xy(g1,mid_temp_array, mid_height_array, nLayer, "Temperature");
     	sleep(1); wait a second 
    }
 	 







	
	gnuplot_close (g1) ; 	
	 */
	
} 
	//print_vec(mid_pot_temp_vec);
	//print_vec(mid_temp_vec);

	//print_vec(mid_temp_vec);

	

	//TODO:	*loop over time steps
	
	
	 
	
	
	
	return 0;
}


void getHeightfromPressure(vector<double> mid_pressure_vec,vector<double> mid_temp_vec, vector<double>& heightVec){
	
	double deltaPressure=mid_pressure_vec[3]-mid_pressure_vec[2];
	heightVec.resize(mid_pressure_vec.size());
	heightVec[0]=0;
	for(int i=1;i<mid_pressure_vec.size();i++){
		heightVec[i]= heightVec[i-1] +deltaPressure* 287.0*mid_temp_vec[i]/mid_pressure_vec[i] /gravitationalConst; // 100.000 Pa Surface pressure
		
	//if (heightVec[i]<1e-10) heightVec[i]=0;
	//cout<<i<<"   "<<heightVec[i]<<"           "<< pressureVec[i]<<endl;
	}


}

double getPotTempFromTemp(double temp, double pressure){

	return temp*pow(pressure/pressure0,(-1)*kappaConst);

}

double getTempFromPotTemp(double potTemp, double pressure){
	return potTemp*pow(pressure/pressure0,kappaConst);
}

//TODO: Evt. klasse/template um allgemeine converter zu haben.
void getPotTempVecFromTempVec(vector<double> TempVec,vector<double> pressureVec,vector<double>&PotTempVec){
	PotTempVec.resize(TempVec.size());
	for(int i=0;i<TempVec.size();i++){
		PotTempVec[i]=getPotTempFromTemp(TempVec[i],pressureVec[i]);
	}
}

void getTempVecFromPotTempVec(vector<double>PotTempVec,vector<double> pressureVec,vector<double> & TempVec){
	TempVec.resize(PotTempVec.size());
	for(int i=0;i<TempVec.size();i++){
		TempVec[i]=getTempFromPotTemp(PotTempVec[i],pressureVec[i]);
	}
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


void calcDeltaTempVec(double deltaTime,vector<double> Edelta,double deltaPressure,vector <double>  &deltaTempVec){

	deltaTempVec.resize(Edelta.size());
	
	for(int i=0;i<Edelta.size();i++){

		deltaTempVec[i]=deltaTime*Edelta[i]*gravitationalConst/deltaPressure/cpConst;

	}

}

void calcMiddleVec(std::vector<double> inputVec,std::vector<double> & outputVec){
	outputVec.resize(inputVec.size()-1);
	for(int i =0;i<inputVec.size()-1;i++){
		outputVec[i]=(inputVec[i]+inputVec[i+1])/2;
	}
}

