#include<iostream>
#include<vector>

using namespace std;

class raditaveTransfer{
private:
	

	std::vector<double> opticalThickness;
	std::vector<double> temp_vec;
	

	int nLayer;
	int nMu;

	//Intern 
	std::vector<double>	Eup,Edn;


public:
	void calcRadiativeTransfer();
	void showEdn();



};





void raditaveTransfer::calcRadiativeTransfer(){

	

}

void raditaveTransfer::showEdn(){

	for(int i=0;i<Edn.size();i++){
		cout<<Edn[i]<<endl;
	}
}


int main(){

	raditaveTransfer air;
	air.calcRadiativeTransfer();

	air.showEdn();

	return 0;
}








void rt(int nlayer,int nmu, double beta0, double lambda_dn, double lambda_up, vector <double> he_vec,  vector <double>  temp_vec, vector <double> opt_thick_vec,vector<double> sb_vec,vector<double>& Edn, vector<double>&Eup ){
  
	opt_thick_vec.resize(nlayer);
	for(int i=0;i<nlayer;i++){
	  opt_thick_vec[i]=beta0* exp(-(he_vec[i]+he_vec[i+1])/2/z_scale)*(he_vec[i]-he_vec[i+1]);
	}


   vector<double>L;
  L.resize(nlayer+1);
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
	L[nlayer]= planck(288.15,lambda_dn,lambda_up)/M_PI;        //sigma*pow(288.15,4)/M_PI;
    
	Eup[nlayer]=Eup[nlayer]+L[nlayer]*2*M_PI*mu*1./nmu;
  
	for(int i=nlayer;i>0;i--){
	  L[i-1]=L[i]*exp(-opt_thick_vec[i-1]/mu)+sb_vec[i-1]*(1-exp(-opt_thick_vec[i-1]/mu));
	  Eup[i-1]=Eup[i-1]+L[i-1]*2*M_PI*mu*1./nmu;
	}
    
      }

  //cout<<"Edn(TOA) " << Edn[0] <<"    Eup(TOA) "<< Eup[0] <<endl;
  //cout<<"Edn(SURF)"<< Edn[nlayer] << "    Eup(SURF) "<< Eup[nlayer]<<endl;
  // cout<<nmu<<"  "<<Eup[nlayer]<<endl;


	cout<<beta0<<"  "<< Edn[nlyer]<<endl;
}