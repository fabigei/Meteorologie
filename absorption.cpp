#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>




using namespace std;


void print_vec(vector<double> vec);
void readAbsorptionFile(string filename,vector<double> &waveNumVec,vector<vector<double> > &optThickMat);
void reverseOptThickMat(vector<vector<double> > & optThickMat);
void getWaveLenghtFromWaveNumVec(vector<double> waveNumVec,vector<double> & waveLengthVec);
void createBorderPressureVec(int nLayer, vector<double> & borderPressureVec);


int main(){

	

	vector<vector<double> > optThickMat;
	vector<double> waveNumVec,waveLengthVec;
	vector<double> borderPressureVec;


	readAbsorptionFile("co2.dtau",waveNumVec,optThickMat);
	//cout<<"optThickMatSize:  "<<optThickMat[2].size()<<endl;
	getWaveLenghtFromWaveNumVec(waveNumVec,waveLengthVec);

	reverseOptThickMat(optThickMat);
	createBorderPressureVec(optThickMat[1].size(),  borderPressureVec);

	for(int i=0;i<optThickMat[1].size();i++){
		cout<<(borderPressureVec[i]+borderPressureVec[i+1])/2.<<"   "<<optThickMat[2][i]<<endl;
	}

	return 0;
}



void print_vec(vector<double> vec){
  for(int i=0;i<vec.size();i++)	cout<<vec[i]<<endl;
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


void getWaveLenghtFromWaveNumVec(vector<double> waveNumVec,vector<double> & waveLengthVec){
	waveLengthVec.resize(waveNumVec.size());
	for(int i =0;i<waveNumVec.size();i++){
		waveLengthVec[i]=1e-2/waveNumVec[i];
	}
}



void createBorderPressureVec(int nLayer, vector<double> & borderPressureVec){
	double pressureStepSize=1000000./nLayer;
	borderPressureVec.resize(nLayer+1);
	for(int i=0;i<nLayer+1;i++){
		borderPressureVec[i]=pressureStepSize*i;
	}
}
