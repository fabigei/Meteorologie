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

double calculatePotTemp(double temp, double pressure);
void print_vec(vector<double> vec);
double calcDeltaTempSolarSurfaceHeating(double deltaTime,double eSolarSurface,double deltaPressure);
int main(int argc, char** args){
	
	int nLayer;
	double e0=1362; // W/m^2
	double albedo=0.3;
	double eSolarSurface=e0/4 *(1-albedo);
	double deltaTime=60; //in s
	vector<double> pressure_vec;
	vector<double> temp_vec;
	vector<double> 
	double deltaPressure;
	
	if (argc==1){ 
		cout<<"Try again -- Enter the number of layers you want to use."<<endl; 
		return 0;
	}	
	else {
		nLayer= atoi(args[1]);
		}
	
	
	
	
	pressure_vec.resize(nLayer);
	temp_vec.resize(nLayer);
	fill(temp_vec.begin(),temp_vec.end(),0);
	
	
	deltaPressure=pressure0/(nLayer-1);
	for(int i=0;i<nLayer;i++){
		pressure_vec[i]=i*deltaPressure;
	}
	
	
	//print_vec(pressure_vec);
	cout<<temp_vec[nLayer-1]<<endl;
	
	
	
	
	temp_vec[nLayer-1]+=calcDeltaTempSolarSurfaceHeating(deltaTime,eSolarSurface,deltaPressure);
	
	
	
	cout<<temp_vec[nLayer-1]<<endl;
	
	
	
	
	
	return 0;
}




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
