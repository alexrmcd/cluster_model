//////// root_dv.cpp 4-26-16 ////////////


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <gsl/gsl_integration.h>
#include <iomanip> 

#include "Constants.h"


class Cluster{  
	
	public:
	std::string name;
	double z ;
	double rh; //halo radius
	int ch ;
	double B0 ; 				//microGauss
	double rcore ; 		//0.291; //core radius in Mpc for coma


	Cluster() :	
		name("Coma"),
        z(0.0232),             //redshift
        rh(0.415),            //halo radius Mpc
        ch(25),            //darksusy channel
        B0(4.7), rcore(0.404)      //Bfield Params


	{	//everything labelled or Coma, sgould work in user options
		std::cout << "creating " << name << " cluster... " <<std::endl;
	}

};


Cluster c;


double bloss(double E ){ 												//overload so that we can just have bloss(E), could also have same as b(E, r) but set r=0 ??  kinda sloppy
	double ne = 1.3e-3;
	double Bmu = 1;													//should us average or something later

	double bloss = 0.0254*pow(Bmu, 2.0)*E*E 								//bsyn bfield_model(r)
					+ 0.25 /** pow(1 + c.z, 4 )*/*E*E  					//bIC
					+ 1.51*ne*(0.36 + log(E/me/ne) )						//brem , in Emma's code  has + 1.51*n*(0.36 + log(E/me) )*E
					+ 6.13*ne*( 1 + log(E/me/ne)/75); 					//bcoul

	bloss = bloss ;//	*1e-16;	

	return bloss;
};


double D(double E){
	//	std::cout << "D" << std::endl;

	double Bmu = 1;	
	double alpha = 1.0/3.0; //close to 1/3??
	double db = pow(20.0 , 2.0/3.0); //just a scaling factor
	double D0 = 3.1e28; // cm/s

	double D = db * D0 * pow(E, alpha)/pow(Bmu , 1.0/3.0);

	return D;
}



double du(double Eu, void * params){

	double du = 1.0/bloss(Eu);

	return du;
}



double U(double E, double mx) {

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	gsl_function F;
	F.function = &du;

	gsl_integration_qags (&F, E, mx, 0, 1e-3, 5000, 
	                w, &result, &error); 

	gsl_integration_workspace_free (w);
	//std::cout << "U(E, mx) : " << result << "  E :" << E   << std::endl;
	return result;

}



double dv(double E , void * params){

	double dv = D(E)/bloss(E);

	return dv;
}




double v( double E,  double mx ){
	
	double max =   U(E, mx);
	//std::cout << "max: " << max << std::endl;
	double min = 0;


		gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	gsl_function F;
	F.function = &dv;


	gsl_integration_qags (&F, min, max, 0, 1e-3, 5000, 
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= 1e16;
	return result;
}



double root_dv(double E, double Ep, double mx ){


	double root_dv  =   sqrt( ( v( E, mx) - v(Ep, mx))   )/mpc2cm * 1000 ;

	return root_dv;
}


main(){

	double mx  = 100;
	double E = 1 ;
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm

	std::cout << root_dv(E, mx , mx) << std::endl; 
	for (int i = 0 ; i < mx + 1 ; ++i ){
		//std::cout << i << " v(Ep) =  "  << sqrt ( v( i, mx)  )/mpc2cm * 1000 << std::endl;
		
		//root_dv(E, i, mx) ;
		std::cout << i << " rdv = "   << root_dv(E, i, mx) << std::endl;

	}


	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;



}



//			std::cout << "		" << 		<< std::endl;



/*







	*/