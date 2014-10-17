#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>     /* atoi */



using namespace std;

void cre_temp(int nlayer, vector<double> & temp_vec); 
void cre_midtemp(int nlayer ,vector<double> & midtemp_vec);
double SB_law(double temp);
double Rad_Tran(double b,double thickness,double l,double mu);
void create_optical_thickness(vector<double> & op_thick, vector<double> & height, double beta0);
void create_height(vector<double> & height, int nlayer,double height_max);

int main(int argc, char** args){
	
	int Nlayer; //Number of  layers 
	int Nmu;
	
	
	vector<double> midtemp_vec; //Temp vector 
	vector<double> bol_vec_up;	//Boltzmann radiation surface to top
	vector<double> height;		//Height vector
	vector<double> op_thick;	//Optical Thickness Vector 
	vector<double> l_up;
	
	double l_surf;
	double height_max=11000;
	double beta0=0.001;
	double sum;
	
	
	//Reading in the number of layers, as an additional argument
	if (argc==1){ 
		cout<<"Try again -- Enter the number of layers you want to use."<<endl; 
		return 0;
	}	
	else {
		Nlayer= atoi(args[1]);
		Nmu=atoi(args[2]);
		}
	
	bol_vec_up.resize(Nlayer);
	l_up.resize(Nlayer);
	
	//Creating the height vector
	create_height( height, Nlayer, height_max);
	
	//Creating an vector holding the middle temp of each layer=> midtem_vec
	cre_midtemp(Nlayer, midtemp_vec);
	
	//Createing the vec witht he optical thickness
	create_optical_thickness(  op_thick, height, beta0);

	
	
	
	//Radiance at the surface, only depending on the temp.
	l_surf=SB_law(288.);
	

	
	for(int i=0;i<Nlayer;i++) {	
		bol_vec_up[i]=SB_law(midtemp_vec[i]);
		//cout<<"Bolvec "<<bol_vec_up[i]<<endl;
		}
	
	sum=0;
	//Rad_Tran(double b,double thickness,double l)
	for(double mu=0;mu<=1;mu+=1./(double)Nmu){
		sum+=mu*Rad_Tran(bol_vec_up[Nlayer-1],op_thick[Nlayer-1],l_surf,mu)/Nmu;
	}
	l_up[Nlayer-1]=2*M_PI*sum/Nmu;


	for(int i=Nlayer-2;i>=0;i--){
		sum=0;
		for(double mu=0;mu<=1;mu+=1./(double)Nmu){
			sum+=mu*Rad_Tran(bol_vec_up[i],op_thick[i],l_up[i+1],mu);
		}
		l_up[i]=2*M_PI*sum/Nmu;
	}
	//for(int i=0;i<Nlayer;i++){
	//		cout<<"L_up "<<l_up[i]<<endl;
	//}
	cout<<"L_up "<<l_up[0]<<"  "<<l_up[Nlayer-1]<<endl;
	double z_scal=8000;
	
	
	return 0;
}	//End of main function






	
// Stefan Boltzmann law; Return the balck body radiation dependent on the temp
double SB_law(double temp){	
	return 5.6704e-8 * pow(temp,4.0);
}


//Schwarzschild equation	
double Rad_Tran(double b,double thickness,double l,double mu){
	return l*exp(-thickness/mu)+b*(1-exp(-thickness/mu)); 
}

void create_optical_thickness(  vector<double> & op_thick, vector<double> & height, double beta0){
	
	int nlayer=height.size()-1;
	double z_scal=8000;
	op_thick.resize(nlayer);
	for(int i=0;i<nlayer;i++){
		op_thick[i]=beta0*(height[i+1]+height[i])*exp(-(height[i+1]+height[i])/(2*z_scal));
		//cout<<"Optical thickness "<<op_thick[i]<<endl;
	}
}
	
	
void create_height(vector<double> & height, int nlayer,double height_max){
	height.resize(nlayer+1);
	for(int i=0;i<nlayer;i++){
		height[i]= height_max/nlayer*(nlayer-i);
	}
	//for(int i=0;i<=nlayer;i++){
	//	cout<<"height "<<height[i]<<endl;
	//}
	
}
