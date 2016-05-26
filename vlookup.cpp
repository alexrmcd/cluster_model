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


double bloss(double E ){ 	
	//std::cout <<  "E in bloss: " << E <<std:: endl;											//overload so that we can just have bloss(E), could also have same as b(E, r) but set r=0 ??  kinda sloppy
	double ne = 1.3e-3;
	double Bmu = 1;													//should us average or something later

	double bloss = 0.0254*pow(Bmu, 2.0)*E*E 								//bsyn bfield_model(r)
					+ 0.25 /** pow(1 + c.z, 4 )*/*E*E  		//bIC
					+ 1.51*ne*(0.36 + log(E/me/ne) )						//brem , in Emma's code  has + 1.51*n*(0.36 + log(E/me) )*E
					+ 6.13*ne*( 1 + log(E/me/ne)/75); 				//bcoul



	return bloss;
};


double D(double E){
	//	std::cout << "D" << std::endl;

	double Bmu = 1;	
	double alpha = 0.9; //close to 1/3??
	double db = pow(20.0 , 2.0/3.0); //just a scaling factor
	double D0 = 3.1e28; // cm/s

	double D = db * D0 * pow(E, alpha)*pow(Bmu , -1.0/3.0);

	return D;
}


double dv(double E , void * params){
	//std::cout << "E in dv: " << E<< std::endl;
	double dv = D(E)/bloss(E);

	return dv;
}




double v( double E,  double mx ){
	//std::cout << "E in v: " << E<< std::endl;
	/*double max =   U(E, mx);
	std::cout << "max: " << max << std::endl;
	double min = mx;// 0;
	std::cout << "min: " << min << std::endl;
*/
		gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &dv;


	gsl_integration_qags (&F, E, mx, 0, 1e-3, 1000, 
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= 1e16;
	return result;
}



double root_dv(double E, double Ep, double mx ){


	double root_dv  =  sqrt( ( v( E, mx) - v(Ep, mx))   )/mpc2cm * 1000 ;
	root_dv *= 1e8;
	return root_dv;
}

void createLUT(double n_vLUT){
		// iteration timer start
	std::clock_t vstart;
	double vduration;
	vstart = std::clock();
	int va ; 
	///////before algorithm
	//std::cout << "creating LUT..." <<std::endl; 
	double mx = 88.7361;

	std::vector<double> vlookup( n_vLUT + 1 );
double scale = mx/n_vLUT;
	//std::cout << root_dv(E, 5 , mx) << std::endl; 
	for (int j = 0 ; j < n_vLUT +1; ++j ){

		double dE = mx/n_vLUT;
		//std::cout << j << "  " << j*dE<< std::endl;
		vlookup[j] = v(j*dE, mx);
		//std::cout << j << ", " << j*dE << " scaled = "<< (int)(j*dE/scale)<<std::endl;
		//if(vlookup[j] > vlookup[0])
		//std::cout <<j << " rootv(E = " << j*dE << ") = " << sqrt(vlookup[j])/mpc2cm*1000 << std::endl;
	}
	//std::cout << "vlookup created..." << c.vlookup[] <<std::endl; 
	////////after algorithm
	vduration = (std::clock()  -  vstart)/(double) CLOCKS_PER_SEC;
	//std::cout << "vlookup time = " << vduration <<std::endl;
	//std::cout << vlookup[0] <<" " <<vlookup[115] <<" "<<vlookup[1000]<< " "<<std::endl;

double E = (int)(me/scale) ;
double Ep = (int)(100/scale) ;
	std::cout<< n_vLUT << " V(E = "<<E<<") = "<< vlookup[me] << ", V(Ep = "<<Ep<<") = " << vlookup[Ep]<< " rdv = " << sqrt(vlookup[E]- vlookup[Ep])/mpc2cm*1000<< std::endl;;
}

main(){


	double E = me;
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm
	double Ep = 0.01265;
	std::cout << E << std::endl;
	std::cout << Ep << std::endl;
	double scale = 0.01;
	E = (int)(E/scale)*scale;
	std::cout << (int)(Ep/scale);
	//Ep = (int)(Ep/scale)*scale;
	std::cout << E << std::endl;
	std::cout << Ep << std::endl;

createLUT(10);
createLUT(100);
createLUT(1000);
createLUT(10000);



/*
	//std::cout << root_dv(E, 5 , mx) << std::endl; 
	for (int i = 0 ; i < n +1; ++i ){
		//std::cout << i << " v(Ep) =  "  << sqrt ( v( i, mx)  )/mpc2cm * 1000 << std::endl;
		double dE = mx/n;

		vlookup[i] = v(i*dE, mx);

		std::cout << "rootdv(" << i*dE << ") = " << 1e8*sqrt(vlookup[i])/mpc2cm*1000 <<std::endl; 

		//root_dv(E, i, mx) ;
		//std::cout << i << " rdv = "   << root_dv(E, i, mx) << std::endl;

	}*/

	//std::cout   << v(1,mx) - v(5,mx) << std::endl;
	//std::cout << " Umax(E = .511) : " << v(0,5,  mx) << std::endl;


	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;



}



//			std::cout << "		" << 		<< std::endl;



/*







	*/