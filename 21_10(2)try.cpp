#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>     /* atoi */
#include <algorithm>

const double sigma=5.6704e-8; // Stefan Boltzmann Konstante

using namespace std;

void init_sys(int nlayer, vector<double> & h_vec, vector <double> & t_vec, vector<double> & mt_vec,vector<double> & sb_vec); // Function to initialize height,temperature vector, also a vector with the average temp of each layer
void print_vec(vector<double> vec);
void rt(int nlayer,int nmu,vector <double>  temp_vec, vector <double> opt_thick_vec,vector<double> sb_vec,vector<double>& Edn, vector<double>&Eup );
double planck(double temp );

int main(int argc, char** args){
	
	
	int Nlayer; //Number of  layers 
	int Nmu;
	double beta0=0.001;
	double z_scale=8000;
	vector<double> he_vec;
	vector<double> te_vec;
	vector<double> mt_vec;
	vector<double> sb_vec;
	vector<double> bb_vec;
	vector<double> opt_thick_vec;
	vector<double> Edn;
	vector<double> Eup;
		//Reading in the number of layers, as an additional argument
	if (argc==1){ 
		cout<<"Try again -- Enter the number of layers you want to use."<<endl; 
		return 0;
	}	
	else {
		Nlayer= atoi(args[1]);
		Nmu= atoi(args[2]);
		}
		
	
	init_sys(Nlayer,he_vec,te_vec,mt_vec,sb_vec);
        
	opt_thick_vec.resize(Nlayer);
	for(int i=0;i<Nlayer;i++){
	  opt_thick_vec[i]=beta0* exp(-(he_vec[i]+he_vec[i+1])/2/z_scale)*(he_vec[i]-he_vec[i+1]);
	}
	  
	//rt(Nlayer, Nmu, te_vec,  opt_thick_vec, sb_vec, Edn, Eup );
	double p=planck(288.15);
	bb_vec.resize(Nlayer);
	for(int i=0;i<Nlayer;i++){
		bb_vec[i]=planck(mt_vec[i]);
	}
	print_vec(bb_vec);
	//print_vec(mt_vec);
}

void init_sys(int nlayer, vector<double> & h_vec, vector <double> & t_vec,vector <double> & mt_vec,vector<double> & sb_vec){
	


	double h_mp= 20000.; // Mesopausenhoehe [km]
	double h_0=0.;
	double dh= (h_mp-h_0)/nlayer;
       
	double t_down= 288.15; // Temp in K, from wikipedia
	double t_trop11= 273.15-56.5;	//Var namen beinhalten hoehe in km
	double t_stra20= 273.15-56.5;
	double t_stra32= 228.65;
	double t_stra47= 270.65;
	double t_meso51= 270.65;
 	double t_meso71= 214.65;
	double t_meso85= 186.95;
	double m ;
	double t;
	h_vec.resize(nlayer+1);
	t_vec.resize(nlayer+1);
	mt_vec.resize(nlayer);
	sb_vec.resize(nlayer);
	

	for(int i=0; i<=nlayer;i++){
		h_vec[i]=(i*dh);
		
		if(h_vec[i] <=11000){
		m=(t_trop11-t_down)/(11000-0);
		t=t_down;
		 }
		else if(h_vec[i] <=20000){
		m=(t_stra20-t_trop11)/(20000-11000);
		t=t_stra20-m*20000;	
 		 }
		else if(h_vec[i] <=32000){
		m=(t_stra32-t_stra20)/(32000-20000);
		t=t_stra32-m*32000;
    		 }
		else if(h_vec[i] <=47000){
		m=(t_stra47-t_stra32)/(47000-32000);
		t=t_stra47-m*47000;
    		 }
		else if(h_vec[i] <=51000){
		m=(t_meso51-t_stra47)/(51000-47000);
		t=t_meso51-m*51000;
    		 }
		else if(h_vec[i] <=71000){
		m=(t_meso71-t_meso51)/(71000-51000);
		t=t_meso71-m*71000;
    		 }
		else if(h_vec[i] <=85000){
		m=(t_meso85-t_meso71)/(85000-71000);
		t= t_meso85-m*85000;
		}
		else {
		cout<<"Height "<< h_vec[i] << " not considered in this model." <<endl;
		exit(1); 
		}
	
	t_vec[i]=(m*h_vec[i]+t);

	}
	for(int i=0;i<nlayer;i++){
		mt_vec[i]=(t_vec[i]+t_vec[i+1])/2;
		sb_vec[i]=sigma*pow(mt_vec[i],4)/M_PI;
	}
	reverse(t_vec.begin(),t_vec.end());
	reverse(mt_vec.begin(),mt_vec.end());
	reverse(h_vec.begin(),h_vec.end());
	reverse(sb_vec.begin(),sb_vec.end());
}

void print_vec(vector<double> vec){
  for(int i=0;i<vec.size();i++){
    cout<<vec[i]<<endl;
  }

}

void rt(int nlayer,int nmu,vector <double>  temp_vec, vector <double> opt_thick_vec,vector<double> sb_vec,vector<double>& Edn, vector<double>&Eup ){
  



   vector<double>L;
  L.resize(nlayer+1);
  Edn.resize(nlayer+1);
  Eup.resize(nlayer+1);
  fill(Edn.begin(),Edn.end(),0);
  fill(Eup.begin(),Eup.end(),0);
  
 

      for(double mu=1./nmu/2;mu<=1;mu+=1./nmu){
    
	L[0]=0;
	for(int i=0;i<nlayer;i++){
	  L[i+1]=L[i]*exp(-opt_thick_vec[i]/mu)+sb_vec[i]*(1-exp(-opt_thick_vec[i]/mu));
	  Edn[i+1]=Edn[i+1]+L[i+1]*2*M_PI*mu*1./nmu;
	}
	fill(L.begin(),L.end(),sqrt(-1));
	L[nlayer]=sigma*pow(288.15,4)/M_PI;
    
	Eup[nlayer]=Eup[nlayer]+L[nlayer]*2*M_PI*mu*1./nmu;
  
	for(int i=nlayer;i>0;i--){
	  L[i-1]=L[i]*exp(-opt_thick_vec[i-1]/mu)+sb_vec[i-1]*(1-exp(-opt_thick_vec[i-1]/mu));
	  Eup[i-1]=Eup[i-1]+L[i-1]*2*M_PI*mu*1./nmu;
	}
    
      }

  cout<<"Edn(TOA) " << Edn[0] <<"    Eup(TOA) "<< Eup[0] <<endl;
  cout<<"Edn(SURF)"<< Edn[nlayer] << "    Eup(SURF) "<< Eup[nlayer]<<endl;
  // cout<<nmu<<"  "<<Eup[nlayer]<<endl;

}

double planck(double temp,double lambdadn, double lambdaup ){
  double hcc=5.955214e-17;//Wm*2
  double hc=1.986446e-25;//Jm
  double k_b=1.38065e-23;// J/K
  double sum=0;	
  double lambda;
 for(double i=4;i<=100;i++){
   lambda=i*1e-6;
     sum+=2*hcc/pow(lambda,5)/(exp(hc/k_b/temp/lambda)-1)*1e-6;
  }

 return M_PI* sum;
  
}
