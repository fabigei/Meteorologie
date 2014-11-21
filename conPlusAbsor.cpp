#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>     /* atoi */
#include <algorithm>
#include <unistd.h>
#include <fstream>
#include <string>
#include <sstream>


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


double e0=1362; // W/m^2
double albedo=0.3;
double eSolarSurface=e0/4 *(1-albedo);




double calculatePotTemp(double temp, double pressure);
void print_vec(vector<double> vec);
double calcDeltaTempSolarSurfaceHeating(double deltaTime,double eSolarSurface,double deltaPressure);
void calcDeltaTempVec(double deltaTime,vector<double> Edelta,double deltaPressure,vector<double> &deltaTempVec);
void radiativeTransfer(int nlayer,int nmu, vector<double> opt_thick_vec,double LSurface,vector<double> sb_vec,vector<double>& Edelta );
double planck(double temp, double lambda_dn, double lambda_up );
void planckVec(vector<double> tempVec, double lambda_dn, double lambda_up, vector <double> & bbVec);
void getHeightfromPressure(vector<double> PressureVecM,vector<double> mid_temp_vec, vector<double>& heightVec);


void calcMiddleVec(std::vector<double> inputVec,std::vector<double> & outputVec);
void getPotTempVecFromTempVec(vector<double> TempVec,vector<double> pressureVec,vector<double>&PotTempVec);
void getTempVecFromPotTempVec(vector<double>PotTempVec,vector<double> pressureVec,vector<double> & TempVec);

void getOpticalThickfromPressure( double tau,int nLayer, vector<double> &opt_thick_vec);

void readAbsorptionFile(string filename,vector<double> &waveNumVec,vector<vector<double> > &optThickMat);
void reverseOptThickMat(vector<vector<double> > & optThickMat);
void getWaveLengthFromWaveNumVec(vector<double> waveNumVec,vector<double> & waveLengthVec);
void createPressureVecB(int nLayer, vector<double> & borderPressureVec);
void createPressureVecM(int nLayer, vector<double> & borderPressureVec);


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




gnuplot_ctrl *g1;




int main(int argc, char** args){
	
	bool do_plots=true;


	//if(do_plots==true){
		 g1 = gnuplot_init();
	//}


	int nLayer,nMu;

	double deltaTime=600; //in s
	double timeStep;
	double MaxTime=3.15569e7; //3.15569e7=1a in seconds
	
	
	vector<double> mid_height_vec;
	vector<double> opt_thick_vec;
	vector<double> planckRadiation_vec;
	vector<double> Edelt_vec;
	vector<double> mid_temp_vec;
	vector<double> middleDeltaTemp_vec;
	vector<double> mid_pot_temp_vec;
	vector<double> waveNumVec;
	vector<vector<double> >	optThickMat;
	vector<double> waveLengthVec;
	vector<double> PressureVecB;
	vector<double> PressureVecM;


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
		nLayer=20;
		cout<<"Setting number of layers to 20 for testing. To change it, one first hast interpolate for the optical thickness."<<endl;
		}
	

	cout<<"Reading in the absorption file..."<<endl;


	readAbsorptionFile("co2.dtau",waveNumVec,optThickMat);
	//cout<<"optThickMatSize:  "<<optThickMat[2].size()<<endl;
	getWaveLengthFromWaveNumVec(waveNumVec,waveLengthVec); // wavelengthVec [0]= 2e-5; waveLengthVec[waveLengthVec.size()]=1.17647e-05


	reverseOptThickMat(optThickMat);
	createPressureVecB(nLayer,  PressureVecB);
	
	cout<<"				.....Done"<<endl;

	//cout<< waveLengthVec[0] << "     "<< waveLengthVec[waveLengthVec.size()-1]<<endl;
	//print_vec(optThickMat[1]);
	PressureVecM.resize(nLayer);
	mid_temp_vec.resize(nLayer);
	opt_thick_vec.resize(nLayer);
	planckRadiation_vec.resize(nLayer);
	mid_height_vec.resize(nLayer);

	fill(mid_temp_vec.begin(),mid_temp_vec.end(),280.);
	

	deltaPressure=pressure0/(nLayer); //(nLayer-1)
	

	createPressureVecM(nLayer,PressureVecM);
	/*
	for(int i=0;i<PressureVecM.size();i++){
		PressureVecM[i]=i*deltaPressure +deltaPressure/2.;
	}
*/


	/*
	vector<double> bands_vec;
	bands_vec.resize(4); // Consider 3 bands: H2O absorption, atmospheric window, CO2 abs
	// lower boarders in mu m  @ToDo: describe that planck functions convert mycro_meter in meter
	bands_vec[0]= 4;
	bands_vec[1]= 8;
	bands_vec[2]= 12;
	bands_vec[3]= 101;
	// optical thickness for each band
	*/
	
	vector <double> Edelt_vec_3Bands;
	Edelt_vec_3Bands.resize(nLayer);
	



		//mid_temp_vec.end()+=calcDeltaTempSolarSurfaceHeating(deltaTime,eSolarSurface,deltaPressure);
	for(int timeStep=0;timeStep<=MaxTime;timeStep+=deltaTime){

			//for(int i=0;i<planckRadiation_vec.size();i++){
			//	cout<<planckRadiation_vec[i]<<"      "<<mid_temp_vec[i]<<endl;

			//}

		fill(Edelt_vec_3Bands.begin(),Edelt_vec_3Bands.end(),0.);
		for (int i=waveLengthVec.size()-1;i>=0;i--){
			//cout << i << endl;
			opt_thick_vec= optThickMat[i];
			//print_vec(opt_thick_vec);
			//cout<<"waveLengthVec "<< waveLengthVec[i]<<endl; // achtung!! planck funktion arbeitet eiegntlich mit mycron!!!!!
			LSurface=planck(mid_temp_vec[nLayer-1],waveLengthVec[i],waveLengthVec[i])/M_PI;
			
			//Still have to create the temp_vec. 
			//Probably a loop over the  lambda seqments.
			planckVec(mid_temp_vec, waveLengthVec[i], waveLengthVec[i],planckRadiation_vec);
			


			//print_vec(planckRadiation_vec);
			radiativeTransfer(nLayer, nMu, opt_thick_vec,LSurface,planckRadiation_vec, Edelt_vec);	
			//print_vec(Edelt_vec);
			for (int j=0; j<Edelt_vec.size();j++){
				Edelt_vec_3Bands[j]+=Edelt_vec[j];
			}

		}






		Edelt_vec_3Bands[nLayer-1]= Edelt_vec_3Bands[nLayer-1]+eSolarSurface;

			
			//print_vec(Edelt_vec);

			

			calcDeltaTempVec(deltaTime,Edelt_vec_3Bands,deltaPressure,middleDeltaTemp_vec);
			//print_vec(Edelt_vec_3Bands);
			//print_vec(mid_temp_vec);

			vectorAdd(middleDeltaTemp_vec,mid_temp_vec,mid_temp_vec);
			


			getPotTempVecFromTempVec(mid_temp_vec, PressureVecM,mid_pot_temp_vec);

			sort(mid_pot_temp_vec.begin(), mid_pot_temp_vec.end(), greater<double>());

			getTempVecFromPotTempVec(mid_pot_temp_vec,PressureVecM,mid_temp_vec);
	 		
	 		


			

	 		//live ploting
	 		if (do_plots== true){
				if (int(timeStep/deltaTime)%10 == 0) {
		      
					getHeightfromPressure(PressureVecM,mid_temp_vec,mid_height_vec);

		      		gnuplot_resetplot  (g1);  //start with new plot rather than plotting into exisiting one 
		     		gnuplot_setstyle   (g1, "linespoints");      // draw lines and points 
		      		gnuplot_set_xlabel (g1, "temperature [K]");   //xaxis label 
		      		gnuplot_set_ylabel (g1, "height [m]");    // yaxis label 
		      
		     
		    	
		    	
				

		    		double *mid_temp_array=&mid_temp_vec[0];
		    		double *mid_pot_temp_array=&mid_pot_temp_vec[0];

		    		double *mid_height_array=&mid_height_vec[0];
		    		double *mid_pressure_array=&PressureVecM[0];
		    	
		     		gnuplot_plot_xy(g1,mid_temp_array, mid_height_array, nLayer, "Temperature");
		     		//sleep(1); // wait a second 
		    		double *Edelt_array_3Bands=&Edelt_vec_3Bands[0];
		    
		     		//gnuplot_plot_xy(g1,Edelt_array_3Bands,mid_height_array,nLayer,"Edelta");
	    			}

	 		} 

	

			
	 
		}
gnuplot_close (g1) ; 
	//print_vec(mid_pot_temp_vec);
	//print_vec(mid_temp_vec);

	//print_vec(mid_temp_vec);
	
	 
	//TODO: 3Band and Optical Thickness
	
	
	return 0;
}


void getHeightfromPressure(vector<double> PressureVecM,vector<double> mid_temp_vec, vector<double>& heightVec){
	
	double deltaPressure=PressureVecM[3]-PressureVecM[2];
	heightVec.resize(PressureVecM.size());
	heightVec[PressureVecM.size()-1]=0;
	for(int i=PressureVecM.size()-2;i>=0;i--){
		heightVec[i]= heightVec[i+1] +deltaPressure* spezGasConstAir*mid_temp_vec[i]/PressureVecM[i] /gravitationalConst; // 100.000 Pa Surface pressure
		
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

	Edelta[nlayer-1]=Edelta[nlayer-1]+Edn[nlayer]-Eup[nlayer];
	
}



double planck(double temp, double lambda_dn, double lambda_up ){
	double hcc=5.955214e-17;//Wm*2
 	double hc=1.986446e-25;//Jm
  	double k_b=1.38065e-23;// J/K
  	double sum=0;
  	double lambda;
 	for(double i=lambda_dn;i<=lambda_up;i++){
   		lambda=i;//*1e-6;
    	sum+=2*hcc/pow(lambda,5)/(exp(hc/k_b/temp/lambda)-1)*1e-6;

	}
 	return M_PI* sum; 
}

void planckVec(vector<double> tempVec, double lambda_dn, double lambda_up, vector <double> & bbVec){
  for (int i=0;i<tempVec.size();i++){
  	bbVec[i]=planck(tempVec[i],lambda_dn,lambda_up);
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


void readAbsorptionFile(string filename,vector<double> & waveNumVec,vector<vector<double> > &optThickMat){
	int neededVecLength=164135; //Just for making the code faster; Not obligatroy. 
	ifstream file(filename);
	string str; 
	string temporary;
	
	optThickMat.reserve(neededVecLength);
	vector<double> temporary_vec;
	   
	waveNumVec.reserve(neededVecLength);
		
	while (getline(file, str)){
	   	
	   	stringstream s (str);
	    s>>temporary;
	    waveNumVec.push_back(stod(temporary));

	    while(s >> temporary){
	    	temporary_vec.push_back(stod(temporary));
	    }
	    
		optThickMat.push_back(temporary_vec);
		//cout<<"SSSSS"<<endl;
		temporary_vec.resize(0);
		temporary_vec.reserve(22);
	}
	
}


//In the original File the second coulum (first in optThickMat) corresponds to the highest pressure, 
//since we start counting from the top we need the lowest pressure in the first matrix column.

void reverseOptThickMat(vector<vector<double> > & optThickMat){
	for(int i=0;i<optThickMat.size();i++){
		reverse(optThickMat[i].begin(),optThickMat[i].end());		
	}
}


void getWaveLengthFromWaveNumVec(vector<double> waveNumVec,vector<double> & waveLengthVec){
	waveLengthVec.resize(waveNumVec.size());
	for(int i =0;i<waveNumVec.size();i++){
		waveLengthVec[i]=1e-2/waveNumVec[i]; //1e-2 because the wavenumber is given in 1/cm
	}
}



void createPressureVecB(int nLayer, vector<double> & borderPressureVec){
	double pressureStepSize=pressure0/nLayer;
	borderPressureVec.resize(nLayer+1);
	for(int i=0;i<nLayer+1;i++){
		borderPressureVec[i]=pressureStepSize*i;
	}
}


void createPressureVecM(int nLayer, vector<double> & borderPressureVec){
	double pressureStepSize=pressure0/nLayer;
	borderPressureVec.resize(nLayer);
	for(int i=0;i<nLayer;i++){
		borderPressureVec[i]=pressureStepSize*(i+0.5);
	}


}