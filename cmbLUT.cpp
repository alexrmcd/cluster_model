//////// cmbSpect.cpp ////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iomanip> 


#include "Constants.h"
/*
class LUT{

	double size;
	std::vector<double> table;

	LUT() : size(1000) , 
			scale(0.01),
			table(size) 
			{};

};
*/








double CMB_bbSpectrum(double nu){
	//double hplanckGeV = J2Gev*hplanck;
	//kb *= J2Gev;
	double T = 2.73;
	//nu = pow(10, nu );
	//double nu = eps/(hplanck*J2Gev);
	//std::cout << "nu = " << nu << std::endl;
	double CMB_bbSpectrum = 8*pi* pow(nu, 3.0)/pow(clight, 3.0)	/(	exp(hplanck * nu /(kb * T )) - 1	);
	//std::cout << "nu("<< eps <<")  = " << nu <<  " spectrum = " << CMB_bbSpectrum << std::endl;;//CMB_bbSpectrum  << std::endl;

	return CMB_bbSpectrum;
}




std::vector<double> createLUT(){
	// iteration timer start
	std::clock_t Gstart;
	double Greensduration;
	Gstart = std::clock();
	int Ga ; 
	///////before algorithm
	std::cout << "creating LUT..." <<std::endl; 
	double nu_max = 1e15 ;
	double nu_min =  1e1 ;
	
	double size = 10000;
	double scale =  (log10(nu_max) - log10(nu_min))/size;//nu_min*pow( 10 , ( log10(nu_max) - log10(nu_min) )/size );
	std::cout << "scale: " << scale << std::endl;


	std::vector<double> CMB_LUT( size );

		for(int i = 0; i <= size  ; ++i){
			
			CMB_LUT[i] = CMB_bbSpectrum(pow(10, log10(nu_min) + i*scale )  );
			
			//std::cout << i <<" "<< pow( 10, log10(nu_min) + i*scale )<< " " << CMB_LUT[i]<<std::endl;
		}

	////////after algorithm
	Greensduration = (std::clock()  -  Gstart)/(double) CLOCKS_PER_SEC;
	std::cout << "CMB_LUT time = " << Greensduration <<std::endl;

	return CMB_LUT;
}

main(){



	std::ostringstream makefilename;
	makefilename << "nuCMB_LUT.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());
	std::vector<double> cmb= createLUT();

	int n_nu = 1000;
	for (int i = 0 ; i <= n_nu ; ++i  ){

		double nu_max =  1e15 ;
		double nu_min =  1e1 ;
		
		double size = 10000;
		double scale =  (log10(nu_max) - log10(nu_min))/size;



		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i)) ;
		double nu_int = (int)(log10(nu/nu_min) /scale);
		//std::cout<< nu << " " << log10(nu/nu_min) /scale <<"  " <<nu_int<<std::endl;
		
		double data = cmb[nu_int] ;
		
		file << nu << "\t" <<  data <<std::endl;
		std::cout << i << "/" << n_nu  << " " << nu << "  "<< nu_int << "\t" << data << std::endl;

	};
	


}