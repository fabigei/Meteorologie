#include <iostream>
#include <vector>
#include <stdlib.h>     /* atoi */
#include <math.h>

using namespace std;



void create_temp(int nlayer, vector<double> & temp_vec); 
double SB_law(double temp);
double Rad_Tran(double b,double thickness,double l);
void create_optical_thickness(vector<double> & op_thick, vector<double> & height, double beta0);
void create_height(vector<double> & height, int nlayer,double height_max);


int main(int argc, char** args){
	
	double height_max=11000;
	int Nlayer; //Number of Layers
	//double optical_thickness;
	vector <double> temp_vec; 		//Will be filled with international standard atmospheric temp.
	vector <double> bol_vec_up;	//Will be filled with BB radiation values	
	vector <double> rad_vec_up;		//The actuall result, 
	vector <double> height_vec;
	vector <double> op_thick_vec;
	//Reading in the number of layers as an argument when executing the program
	if (argc==1){ 
		cout<<"Try again -- Enter the number of layers you want to use."<<endl; 
		return 0;
	}	
	else {
		Nlayer= atoi(args[1]);
		}
	
	
	
	//optical_thickness=1/(double)Nlayer; 
	
	
	create_temp(Nlayer,temp_vec);
	
	create_height(height_vec, Nlayer, height_max);
	create_optical_thickness(op_thick_vec,height_vec,0.01);
	
	
	bol_vec_up.resize(temp_vec.size());
	rad_vec_up.resize(Nlayer);
	
	//Calculating the boltzmann radiation, surface to top. 
	for(int i=0;i<Nlayer+1;i++){
		bol_vec_up[i]=SB_law(temp_vec[i]);
		cout<<"bolvecup "<<bol_vec_up[i]<<endl;
	}
	
	rad_vec_up[Nlayer-1]=bol_vec_up[Nlayer];
	for(int i=Nlayer-1;i>0;i--){
		rad_vec_up[i-1]=Rad_Tran(bol_vec_up[i],op_thick_vec[i],rad_vec_up[i]);
		cout<< bol_vec_up[i]<<"  "<<op_thick_vec[i]<<"  "<<rad_vec_up[i+1]<<endl;
	}
	
	for(int i=0;i<=Nlayer;i++){
		cout<<"Radian "<<i<<"	"<<rad_vec_up[i]<<endl;
	}
	//cout<<"up "<<rad_vec_up[0]<<endl;
	//cout<<"down "<<rad_vec_up[Nlayer]<<endl;
	
	
	return 0;
}


//Fill the temp vec with temperatures according to the std ISA

void create_temp(int nlayer ,vector<double> & temp_vec){
	
	double temp_down=288.15;
	double temp_up= 273.15-56.5 ;
	//temp_vec.resize(nlayer)
	double tempdiff=(temp_down-temp_up)/nlayer;
	
	for(int  i=0;i<nlayer+1;i++){
		temp_vec.push_back(temp_up+tempdiff*i);
		cout<<"Tempvec "<<temp_vec[i]<<endl;
	}

	
	
}
	
	
// Stefan Boltzmann law; Return the balck body radiation dependent on the temp
double SB_law(double temp){	
	return 5.6704e-8 * pow(temp,4.0);
}


//Schwarzschild equation	
double Rad_Tran(double b,double thickness,double l){
	return l*exp(-thickness)+b*(1-exp(-thickness)); 
}

void create_optical_thickness(  vector<double> & op_thick, vector<double> & height, double beta0){
	
	int nlayer=height.size()-1;
	double z_scal=8000;
	op_thick.resize(nlayer);
	for(int i=0;i<nlayer;i++){
		op_thick[i]=beta0*exp(-(height[i+1]+height[i])/2 /z_scal);
	}
}
	
	
void create_height(vector<double> & height, int nlayer,double height_max){
	height.resize(nlayer+1);
	for(int i=0;i<nlayer+1;i++){
		height[i]= height_max/nlayer*i;
	}
}
