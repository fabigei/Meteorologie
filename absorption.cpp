#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>




using namespace std;


void print_vec(vector<double> vec);
void readInAbsorptionFile(string filename,vector<double> &waveNumVec,vector<vector<double> > &optThickMat);



int main(){

	

	vector<vector<double> > optThickMat;
	vector<double> waveNumVec;

	readInAbsorptionFile("co2.dtau",waveNumVec,optThickMat);
	print_vec(optThickMat[0]);


	return 0;
}



void print_vec(vector<double> vec){
  for(int i=0;i<vec.size();i++)	cout<<vec[i]<<endl;
}

void readInAbsorptionFile(string filename,vector<double> & waveNumVec,vector<vector<double> > &optThickMat){

	ifstream file(filename);
	string str; 
	string temporary;
	
	optThickMat.reserve(164135);
	vector<double> temporary_vec;
	   
	waveNumVec.reserve(164135);
		
	while (getline(file, str)){
	    	
	   	stringstream s (str);
	    s>>temporary;
	    waveNumVec.push_back(stod(temporary));

	    while(s >> temporary){
	    			
	    	temporary_vec.push_back(stod(temporary));
	    
	    }
		optThickMat.push_back(temporary_vec);
		temporary_vec.resize(0);
		temporary_vec.reserve(22);
	}
	for(int i=0;i<optThickMat.size();i++){
		reverse(optThickMat[i].begin(),optThickMat[i].end());		
	}
}