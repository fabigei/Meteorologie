class opticalThickness{
	
private:
	double zscale=8000; //in meter

	double beta0; 
	std::vector<double> height_vec;
	std::vector<double> opticalThickness_vec;

public:
	void calcOpticalThickness(double betaZero,std::vector<double> );



};