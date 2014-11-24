#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>     /* atoi */
#include <algorithm>
#include <unistd.h>
 	


#include "./gnuplot_i-2.10/src/gnuplot_i.h"


using namespace std;


const double pressure0=100000; //in Pa
const double gasConst=8.3144621; //J/(K*mol)
const double gravitationalConst=9.81; //m/s^2
const double kappaConst=2/7.; 
const double spezGasConstAir=287.; //in J/(kg*K)
const double cpConst=1005;
const double molarMassDryAir= 28.97; // g/mol
const double z_scale=8000;
const double c=2.99e8;         // [m s^-1]
const double h=6.626e-34;      // [J s]
const double k_B=1.3806e-23;    // [J K^-1]


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

void getOpticalThickfromPressure( double tau,int nLayer, vector<double> &opt_thick_vec);





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

	double e0=1362; // W/m^2
	double albedo=0.3;
	double eSolarSurface=e0/4 *(1-albedo);

	gnuplot_ctrl *g1;


int main(int argc, char** args){
	
	bool do_plots=false;


	//if(do_plots==true){
		 g1 = gnuplot_init();
	//}


	int nLayer,nMu;

	double deltaTime=6000; //in s
	double timeStep;
	double MaxTime=3.15569e7; //3.15569e7=1a in seconds
	double beta0=0.0019;
	vector<double> mid_pressure_vec;
	
	vector<double> mid_height_vec;
	vector<double> opt_thick_vec;
	vector<double> planckRadiation_vec;
	vector<double> Edelt_vec;
	vector<double> mid_temp_vec;
	vector<double> middleDeltaTemp_vec;
	vector<double> mid_pot_temp_vec;
	
	double tau;
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
	planckRadiation_vec.resize(nLayer);
	mid_height_vec.resize(nLayer);

	fill(mid_temp_vec.begin(),mid_temp_vec.end(),250.);
	

	deltaPressure=pressure0/(nLayer); //(nLayer-1)
	for(int i=0;i<mid_pressure_vec.size();i++){
		mid_pressure_vec[i]=i*deltaPressure +deltaPressure/2.;
		//cout<<mid_pressure_vec[i]<<endl;
	}
	


	vector<double> bands_vec;
	bands_vec.resize(4); // Consider 3 bands: H2O absorption, atmospheric window, CO2 abs
	// lower boarders in mu m  @ToDo: describe that planck functions convert mycro_meter in meter
	bands_vec[0]= 4.0e-6;
	bands_vec[1]= 8.0e-6;
	bands_vec[2]= 12.0e-6;
	bands_vec[3]= 101.0e-6;
	// optical thickness for each band
	
		
	vector <double> Edelt_vec_3Bands;
	Edelt_vec_3Bands.resize(nLayer);
		
for (tau =0; tau <= 10; tau++){
		//mid_temp_vec.end()+=calcDeltaTempSolarSurfaceHeating(deltaTime,eSolarSurface,deltaPressure);
	for(int timeStep=0;timeStep<=MaxTime;timeStep+=deltaTime){

			

		fill(Edelt_vec_3Bands.begin(),Edelt_vec_3Bands.end(),0.);
		for (int i=0;i<bands_vec.size()-1;i++){
			
			/*
			if (i==0) {
				beta0= 0.0019;
				tau=10.;
			}
			if (i==1){
				beta0= 5e-5;
				tau=0.5;
			}
			if (i==2) {
				beta0= 6.8e-4;
				tau=1.;
			}
			*/

			//cout << i <<"  "<< beta0 << endl;
			//cout << "band lower "<< bands_vec[i]<<" band upper" << bands_vec[i+1]-1<< endl;
	
			//calcOpticalThickfromPressure(beta0,deltaPressure,opt_thick_vec);
			getOpticalThickfromPressure(  tau, nLayer, opt_thick_vec);

			LSurface=planck(mid_temp_vec[nLayer-1],bands_vec[i],bands_vec[i+1]);// /M_PI;
			
			// cout<<std::scientific;
		//  cout<<"L surf  "<< LSurface<<endl;
			//Still have to create the temp_vec. 
			//Probably a loop over the  lambda seqments.
			planckVec(mid_temp_vec, bands_vec[i] , bands_vec[i+1]-1,planckRadiation_vec);
	
			radiativeTransfer(nLayer, nMu, opt_thick_vec,LSurface,planckRadiation_vec, Edelt_vec);
			
			for (int j=0; j<Edelt_vec.size();j++){
				Edelt_vec_3Bands[j]+=Edelt_vec[j];
			}
			

		}
		Edelt_vec_3Bands[nLayer-1]= Edelt_vec_3Bands[nLayer-1]+eSolarSurface;

			
			//print_vec(Edelt_vec);

		

			calcDeltaTempVec(deltaTime,Edelt_vec_3Bands,deltaPressure,middleDeltaTemp_vec);
		
			//middleDeltaTemp_vec[nLayer-1]+=eSolarSurface * gravitationalConst/deltaPressure/cpConst;
			//cout<<"ssssss"<<endl;
			//print_vec(middleDeltaTemp_vec);
			//deltaT[nlyr-1] += (E_s + Edn_tot[nlyr] - Eup_tot[nlyr]) * g / deltap / c_p;  
			//TODO: *size(Edelt_vec),size(deltaTem_vec),size(mTemp_vec) must fit together
			//calcMiddleVec(deltaTemp_vec,MiddleDeltaTemp_vec);

		
			vectorAdd(middleDeltaTemp_vec,mid_temp_vec,mid_temp_vec);
			//for (int i=0; i<nLayer; i++){
			//	mid_temp_vec[i]=mid_temp_vec[i]+middleDeltaTemp_vec[i];
			//}


			//print_vec(middleDeltaTemp_vec);
			//cout<<endl;
		
			getPotTempVecFromTempVec(mid_temp_vec, mid_pressure_vec,mid_pot_temp_vec);


			
			sort(mid_pot_temp_vec.begin(), mid_pot_temp_vec.end(), greater<double>());


			getTempVecFromPotTempVec(mid_pot_temp_vec,mid_pressure_vec,mid_temp_vec);
	 		
	 		


			

	 		//live ploting
	 		if (do_plots== true){
				if (int(timeStep/deltaTime)%10 == 0) {
		      
					getHeightfromPressure(mid_pressure_vec,mid_temp_vec,mid_height_vec);

		      		gnuplot_resetplot  (g1);  //start with new plot rather than plotting into exisiting one 
		     		gnuplot_setstyle   (g1, "linespoints");      // draw lines and points 
		      		gnuplot_set_xlabel (g1, "temperature [K]");   //xaxis label 
		      		gnuplot_set_ylabel (g1, "height [m]");    // yaxis label 
		      
		     
		    	
		    	
				

		    		double *mid_temp_array=&mid_temp_vec[0];
		    		double *mid_pot_temp_array=&mid_pot_temp_vec[0];

		    		double *mid_height_array=&mid_height_vec[0];
		    		double *mid_pressure_array=&mid_pressure_vec[0];
		     		gnuplot_plot_xy(g1,mid_temp_array, mid_height_array, nLayer, "Temperature");
		     		//sleep(1); // wait a second 
		    		double *Edelt_array_3Bands=&Edelt_vec_3Bands[0];
		    
		     		//gnuplot_plot_xy(g1,Edelt_array_3Bands,mid_height_array,nLayer,"Edelta");
	    			}

	 		} 

	

			
	 
		}
		cout<<"tau " <<tau<<"        Bodentemp "<<mid_temp_vec[nLayer-1]<<endl;
}
gnuplot_close (g1) ; 
	//print_vec(mid_pot_temp_vec);
	//print_vec(mid_temp_vec);

	//print_vec(mid_temp_vec);
	
	 
	//TODO: 3Band and Optical Thickness
	
	
	return 0;
}


void getHeightfromPressure(vector<double> mid_pressure_vec,vector<double> mid_temp_vec, vector<double>& heightVec){
	
	double deltaPressure=mid_pressure_vec[3]-mid_pressure_vec[2];
	heightVec.resize(mid_pressure_vec.size());
	heightVec[mid_pressure_vec.size()-1]=0;
	for(int i=mid_pressure_vec.size()-2;i>=0;i--){
		heightVec[i]= heightVec[i+1] +deltaPressure* spezGasConstAir*mid_temp_vec[i]/mid_pressure_vec[i] /gravitationalConst; // 100.000 Pa Surface pressure
		
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
  
	fill(Edn.begin(),Edn.end(),0);
	fill(Eup.begin(),Eup.end(),0);
/*
  if(Edn.size()!=nlayer+1){
	  Edn.resize(nlayer+1);
	  Eup.resize(nlayer+1);
	  fill(Edn.begin(),Edn.end(),0);
	  fill(Eup.begin(),Eup.end(),0);
	}
*/
 
 
    for(double mu=1./nmu/2;mu<=1;mu+=1./nmu){
     
		L[0]=0;
		for(int i=0;i<nlayer;i++){
		  L[i+1]=L[i]*exp(-opt_thick_vec[i]/mu)+sb_vec[i]*(1-exp(-opt_thick_vec[i]/mu));
		  Edn[i+1]=Edn[i+1]+L[i+1]*2*M_PI*mu*1./nmu;
		}
		fill(L.begin(),L.end(),sqrt(-1));
		L[nlayer]= LSurface;//planck(288.15,lambda_dn,lambda_up)/M_PI;        //sigma*pow(288.15,4)/M_PI;
	     //cout<<std::scientific;
		  //cout<<"Lsuf  "<< L[nlayer]<<endl;
		Eup[nlayer]=Eup[nlayer]+L[nlayer]*2*M_PI*mu*1./nmu;
	  	
		for(int i=nlayer;i>0;i--){
		  L[i-1]=L[i]*exp(-opt_thick_vec[i-1]/mu)+sb_vec[i-1]*(1-exp(-opt_thick_vec[i-1]/mu));
		 
		  Eup[i-1]=Eup[i-1]+L[i-1]*2*M_PI*mu*1./nmu;
		}
    
    }


	

			
	 
		
//gnuplot_close (g1) ; 
    //cout<<"Edn"<<endl;
    //print_vec(Edn);
    //cout<<"Eup"<<endl;
    //print_vec(Eup);
    for (int i=0;i<nlayer+1;i++){
      Enet[i]= Edn[i]-Eup[i];
	}
	
	for (int i=0;i<nlayer;i++){
		Edelta[i]= Enet[i]-Enet[i+1];
	//cout<<i<<"  "<<Edelta[i]<<endl;
	}
	Edelta[nlayer-1]=Edelta[nlayer-1]+Edn[nlayer]-Eup[nlayer];
	

	
	

}
/*
double planck(double temp, double lambda_dn, double lambda_up ){
	double hcc=5.955214e-17;//Wm*2
 	double hc=1.986446e-25;//Jm
  	double k_b=1.38065e-23;// J/K
  	double sum=0;
  	double lambda;
 	for(double i=lambda_dn+0.125;i<=lambda_up;i+=0.25){
   		lambda=i*1e-6;
    	sum+=2*hcc/pow(lambda,5)/(exp(hc/k_b/temp/lambda)-1)*0.25e-6; // TODO: dlambda einführen!
    		}
 	return M_PI*sum; 
}

*/



double planck(double temp, double lambda_dn, double lambda_up ){
	double hcc=5.955214e-17;//Wm*2
 	double hc=1.986446e-25;//Jm
  	double k_b=1.38065e-23;// J/K
  	double sum=0;
  	double deltaLambda;
  	double nLambda=100;

  	deltaLambda= (lambda_up-lambda_dn)/nLambda;

 	for(double i=lambda_dn+deltaLambda/2;i<lambda_up;i+=deltaLambda){
   		//lambda=i*1e-6;
    	sum+=2*h*c*c/pow(i,5)/(exp(h*c/k_b/temp/i)-1.)*deltaLambda; // TODO: dlambda einführen!
    		}
 	return sum; 
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


void getOpticalThickfromPressure( double tau,int nLayer, vector<double> &opt_thick_vec){
	for(int i=0;i<opt_thick_vec.size();i++){
	  opt_thick_vec[i]=tau/nLayer;
	}
}