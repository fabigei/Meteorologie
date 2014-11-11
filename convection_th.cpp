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
const double zScale= 8000; // Scaleheight in m

double calculatePotTemp(double temp, double pressure);
void print_vec(vector<double> vec);
double calcDeltaTempSolarSurfaceHeating(double deltaTime,double eSolarSurface,double deltaPressure);
void rt(int nlayer,int nmu, vector<double> opt_thick_vec,double LSurface,vector<double> sb_vec,vector<double>& Edelta );
double calcPlanck(double temp, double lambda_dn, double lambda_up );
void calcPlanckVec(vector<double> tempVec, double lambda_dn, double lambda_up, vector <double> & bbVec);
void getHeightfromPressure(vector<double> pressureVec, vector<double> heightVec);
void calcOpticalThickfromPressure(double beta0,double deltaPressure, vector<double> &opt_thick_vec);



int main(int argc, char** args){
	
	int nLayer;
	double e0=1362; // W/m^2
	double albedo=0.3;
	double eSolarSurface=e0/4 *(1-albedo);
	double deltaTime=60; //in s
	double deltaPressure;

	double beta0=0.001;
	double lambda_dn=1;
	double lambda_up=100;

	vector<double> mid_pressure_vec;
	vector<double> temp_vec;
	vector<double> height_vec;
	vector<double> opt_thick_vec; 
	vector<double> bb_vec;

	if (argc==1){ 
		cout<<"Try again -- Enter the number of layers you want to use."<<endl; 
		return 0;
	}	
	else {
		nLayer= atoi(args[1]);
		}
	
	
	
	
	mid_pressure_vec.resize(nLayer);
	temp_vec.resize(nLayer);
	fill(temp_vec.begin(),temp_vec.end(),0);
	opt_thick_vec.resize(nLayer-1);
	fill(opt_thick_vec.begin(),opt_thick_vec.end(),0);
	bb_vec.resize(nLayer-1);

	
	deltaPressure=pressure0/(nLayer-1);
	for(int i=0;i<nLayer-1;i++){
		mid_pressure_vec[i]=i*deltaPressure +deltaPressure/2.;
	}
	

	calcOpticalThickfromPressure(beta0,deltaPressure,opt_thick_vec);

	print_vec(opt_thick_vec);
	double sum=0;
	for (int i=0;i<nLayer;i++){
		sum+=opt_thick_vec[i];

	}
	//cout<< "Sum "<< sum<< endl;
	
	calcPlanckVec(temp_vec,lambda_dn,lambda_up,bb_vec);
	//print_vec(bb_vec);

	//////LSurface extra berechnen!!////
	//Lsurface=planck(288.15,lambda_dn,lambda_up)/M_PI; 







	temp_vec[nLayer-1]+=calcDeltaTempSolarSurfaceHeating(deltaTime,eSolarSurface,deltaPressure);
	
	
	
	//cout<<temp_vec[nLayer-1]<<endl;
	
	
	
	
	
	return 0;
}

/* Noch nicht ganz richtig....Probleme mit der Hoehe
void getHeightfromPressure(vector<double> pressureVec, vector<double> heightVec){
	heightVec.resize(pressureVec.size());
	for(int i=0;i<pressureVec.size();i++){
		heightVec[i]= -8000*log(pressureVec[i]/100000);  // 100.000 Pa Surface pressure
		
	if (heightVec[i]<1e-10) heightVec[i]=0;
	cout<<heightVec[i]<<"           "<< pressureVec[i]<<endl;
	}


}
*/


double calculatePotTemp(double temp, double pressure){

	return temp*pow(pressure/pressure0,kappaConst);

}


double calcDeltaTempSolarSurfaceHeating(double deltaTime,double eSolarSurface,double deltaPressure){
	
	return deltaTime*eSolarSurface*gravitationalConst/deltaPressure/cpConst;

}


void print_vec(vector<double> vec){
  for(int i=0;i<vec.size();i++){
    cout<<vec[i]<<endl;
  }
}

void calcOpticalThickfromPressure( double beta0,double deltaPressure, vector<double> &opt_thick_vec){
	for(int i=0;i<opt_thick_vec.size();i++){
	  opt_thick_vec[i]=beta0*zScale/pressure0*deltaPressure;
	}
}


void rt(int nlayer,int nmu, vector<double> opt_thick_vec,double LSurface,vector<double> sb_vec,vector<double>& Edelta ){
  	
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

  //cout<<"Edn(TOA) " << Edn[0] <<"    Eup(TOA) "<< Eup[0] <<endl;
  //cout<<"Edn(SURF)"<< Edn[nlayer] << "    Eup(SURF) "<< Eup[nlayer]<<endl;
  // cout<<nmu<<"  "<<Eup[nlayer]<<endl;

     for (int i=0;i<nlayer+1;i++){
      Enet[i]= Edn[i]-Eup[i];
		}

	for (int i=0;i<nlayer;i++){
		Edelta[i]= Enet[i]-Enet[i+1];
	}
	
}

double calcPlanck(double temp, double lambda_dn, double lambda_up ){
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

void calcPlanckVec(vector<double> tempVec, double lambda_dn, double lambda_up, vector <double> & bbVec){
  for (int i=0;i<tempVec.size();i++){
  	bbVec[i]=calcPlanck(tempVec[i],lambda_dn,lambda_up);

	}
  
}
