
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
//////////////
//  g++ -o P_IC P_IC.cpp -lgsl -lgslcblas
//////////////

double G(double eps, double E_gamma, double boost){

	//double boost = E/(me * pow( clight , 2))
	double GAMMA = 4 * eps * boost/me;
	double q = E_gamma / (GAMMA * (boost * me  - E_gamma) );
	

	double G = 2 * q* log(q) + (1+2*q)*( 1 - q ) ;//+ pow(GAMMA*q , 2)*(1-q) / ( 2*(1+GAMMA*q ) );//

		//std::cout << "G = " << G <<" ";
		//std::cout << "1/4*boost^2 = " << 1/(4*boost*boost) << " q = " << q<< std::endl;
		//std::cout << "boost = " << boost <<"  " ;
		//std::cout <<"eps = "<< eps <<" GAMMA = " << GAMMA << " q = " << q <<std::endl;
	
	return G;
}




double IC_cross_Section( double eps, double E_gamma, double E){

	//eps *= 1e-9;

	double boost = E/me;
	//std::cout << "boost (sigma) = "<< boost <<std::endl;
	
	double sigma_thompson = 6.6524e-25;

	double sigma = 3 * sigma_thompson/(4* eps * pow(boost ,2))* G(eps, E_gamma, boost);

	return sigma;
}

double CMB_bbSpectrum(double eps){
	//double hplanckGeV = J2Gev*hplanck;
	//kb *= J2Gev;
	double T = 2.73;

	double nu_CMB = eps/(hplanck*J2Gev);
	double CMB_bbSpectrum = 8*pi *pow(nu_CMB, 2.0)/pow(clight, 3.0)	/(	exp(hplanck * nu_CMB /(kb * T )) - 1	);
	std::cout << "nu("<< eps <<")  = " << nu_CMB <<  " spectrum = " << CMB_bbSpectrum << std::endl;//CMB_bbSpectrum  << std::endl;

	return CMB_bbSpectrum;
}

double dpIC(double eps, void * params ){

	//double nu = 0.02;	//Ghz
	std::vector<double> pICParams = *(std::vector<double> *)params;
	double E_gamma = pICParams[0] ;
	double E = pICParams[1] ;
	

	double dpIC = CMB_bbSpectrum(eps) *  IC_cross_Section( eps, E_gamma,  E);
	return dpIC;
}


double gslInt_pIC(double nu, double E){			//int over eps

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> pICParams (2);

	double E_gamma = hplanck * nu * J2Gev;
	std::cout << "E_gamma =  "<< E_gamma<< std::endl;
	double boost = E/me;
	std::cout <<"BOOST = " << boost <<std::endl;

	double eps_max = E_gamma / (1 - E_gamma/(boost*me) );
	double eps_min = E_gamma / (4 * boost*( boost - E_gamma/me));

	std::cout << "eps_max = " << eps_max <<" eps_min = "	<< eps_min << std::endl;




	pICParams[0] = E_gamma;
	pICParams[1] = E;

	gsl_function F;
	F.function = &dpIC;
	F.params = &pICParams;
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, eps_min, eps_max, 0, 1e-3, 1000, //x?
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= clight * E_gamma;
		std::cout << "pIC( nu = " << nu <<", "<< "E_gamma = " << E_gamma << " , E = "<<E<<" ) = "<< result<<std::endl;


	return result;

}


main(){
	std::cout << gslInt_pIC( 1e7	, 40) <<"\n"<< std::endl;


}