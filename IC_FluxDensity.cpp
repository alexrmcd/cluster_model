////////////////////////////////////////////////////
////		FluxDensity.cpp 			5-19-16  ///
////////////////////////////////////////////////////

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

/* to compile and run, use command 

g++ -o IC_FluxDensity IC_FluxDensity.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
*/

//////////////////////////// routines for calling fortran/darksusy stuff /////////////////////////

extern"C" {	 								//interface with fortran code to initialize darksusy

void dsinit_();

}	

extern"C" {									// dshayield gives injection spectrum
double dshayield_(double *mwimp, double *emuthr, int *ch,  int *yieldk, int *istat);
}


//////////////////////////////////////////////////////////////////////////////////////////////////


// consider putting this in its own header?? problem is that we need to access and change member variables
class Cluster{  
	
	public:
	std::string name;
	double z ;
	double rh; //halo radius

	double B0 ; 				//microGauss
	double rcore ; 		//0.291; //core radius in Mpc for coma
	double root_dv;

	double alpha;		// D(E) ~ E^alpha

	Cluster() :	name(""),
				z(0.0232), 						//redshift
				rh(0.415),						//halo radius Mpc
				B0(4.7), rcore(0.291)	,			//Bfield Params
				root_dv(0.035)	,
				alpha(1.0/3.0)				//DIffusion parameter, not really a cluster thing but easy access is good

	{	//everything labelled or Coma, sgould work in user options
		std::cout << "creating cluster... " << std::endl;
	}

	double bfield_model (double r) {

		double beta = 0.75;
		double eta = 0.5;
		double rc = rcore*mpc2cm;
		double B_field = B0 * pow(( 1 + r*r/(rc*rc)),(-1.5*beta*eta));		// Storm et al 2013 

		return B_field;

	};	

	double bloss(double E , double r){

		double bloss = 0.0253*pow(bfield_model(r), 2)*E*E 					//bIC bfield_model(r)
						+ 0.265 * pow(1 + z, 4 )*E*E  					//bIC
						+ 1.51*n*(0.36 + log(E/me/n) )			//+ 1.51*n*(0.36 + log(E/me) )*E						//bbrem Note this is most likely incorrect, no factor of E, and should be E/me/nin log
						+ 6.13*n*( 1 + log(E/me/n)/75); 					//bcoul

		bloss = bloss *1e-16;	

		return bloss;
	};


	double bloss(double E ){ 												//overload so that we can just have bloss(E), could also have same as b(E, r) but set r=0 ??  kinda sloppy
		double ne = 1.3e-3;
		double Bmu = 1;													//should us average or something later

		double bloss = 0.0254*pow(Bmu, 2.0)*E*E 								//bIC bfield_model(r)
						+ 0.25 /** pow(1 + c.z, 4 )*/*E*E  					//bIC
						+ 1.51*ne*(0.36 + log(E/me/ne) )						//brem , in Emma's code  has + 1.51*n*(0.36 + log(E/me) )*E
						+ 6.13*ne*( 1 + log(E/me/ne)/75); 					//bcoul

		bloss = bloss ;//	*1e-16;	

		return bloss;
	};


	double DM_profile(double r){

		double rhos = 0.039974 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
		double rs = 0.404*mpc2cm;  //DM scale radius in Mpc Storm13  (Mathematica)


		double rho = rhos / ( r/rs * pow(1 + r/rs , 2 )); //x? + 1e-100

		return rho;

	};

	double D(double E){

		double Bmu = 1;	
		//double alpha = 1.0/3.0; //close to 1/3??
		double db = pow(20.0 , 2.0/3.0); //just a scaling factor
		double D0 = 3.1e28; // cm/s

		double D = db * D0 * pow(E, alpha)/pow(Bmu , 1.0/3.0);

		return D;
	};



};

class Particle {
	public:

	int ch;
	double mx;
	double sv;

	Particle() : ch(0) , mx(0), sv(4.7e-25){}

};

Cluster c;
Particle p;


///////////////////        root_dv        ////////////////////////////////


double dv(double E , void * params){

	double dv = c.D(E)/c.bloss(E);

	return dv;
}



double v( double E ){



	gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &dv;

	gsl_set_error_handler_off();
	gsl_integration_qags (&F, E, p.mx, 0, 1e-1, 1000, 
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= 1e16;
	return result;
}

double root_dv(double Ep,  double vE){

	double root_dv  =   sqrt( ( vE ) - v(Ep)   ) ;

	return root_dv;
}

////////////////////////////////////////////////

double ddist(double z  , void * params){

	double distint =  mpc2cm * clight / ( H0 * sqrt( OmegaM * pow(1 + z , 3)  + OmegaL ) ); 

	return distint;

}

//distance as a function of redshift
double Dist(){ 

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &ddist;

	gsl_integration_qags (&F, 0, c.z, 0, 1e-3, 1000,
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;
 		
}


double rconst(double rcm){

	double thetaB = 25; // beam size in arcsec
	double dist_z = Dist() / (1 + c.z);


	double rb = 0.5 * dist_z *thetaB/(3600) * (pi/180); //beam size in cm


	double rconst;

	if( rb < rcm){
		rconst = rcm;
	}

	else if(rb > rcm){
		rconst = rb;
	};
	
	return rconst;

}


double dGreens(double rp,  void * params){ 
	std::vector<double> greenParam = *(std::vector<double> *)params;
	double ri = greenParam[0] ;
	double root_dv = greenParam[1];

	double dGreens = rp/ri* pow( c.DM_profile(rp) , 2.0) * (exp( - pow( (rp-ri)/( 2*root_dv ) , 2)) - exp( - pow( ( rp + ri)/(2*root_dv) , 2)) ) ;

	return dGreens;

}
/*
double GreenSum (double rp, void * params) {  //called by ddIC

	std::vector<double> greenParam = *(std::vector<double> *)params;



	double root_dv  = greenParam[1];



	double rh = c.rh * mpc2cm ;
	int imNum = 7; //number of image pairs + 1
	double Gsum = 0 ;


	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);

		Gsum += pow(-1, i) * dGreens(rp, ri, root_dv);

	}

	Gsum *=  ;

	return Gsum;

}*/


double gslInt_Green(double ri,  double root_dv){
	/*	///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	*/	///////////

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> greenParam (2);

	double rh = c.rh*mpc2cm;
	greenParam[0] = ri ;
	greenParam[1] = root_dv ;


	gsl_function F;
	F.function = &dGreens;
	F.params = &greenParam;

	gsl_integration_qags (&F, 1e-16, rh, 0, 1e-1, 1000, //x?
	                w, &result, &error); 

	gsl_integration_workspace_free (w);


	/*		///////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "greens duration: " << duration << std::endl;
	*/	///////
	//std::cout << "greens = " << result << std::endl;
	return result;

}


double darksusy (double Ep){

	int yieldk = 151;
	int istat;

	double ds = dshayield_(&p.mx, &Ep, &p.ch, &yieldk, &istat);

	return ds;
}

double ddiffusion(double Ep, void * params){
	std::vector<double> diffusionParams = *(std::vector<double> *)params;

	double E = diffusionParams[0];
	double ri = diffusionParams[1] ;
	double vE = diffusionParams[2];
	double rootdv = diffusionParams[3];

	//double rootdv = root_dv( Ep, vE); // 0.035*mpc2cm ; //  	//		 

	double ddiffusion = darksusy(Ep) ;//* (1.0/rootdv) * gslInt_Green(ri, rootdv);

	return ddiffusion;

}

double gslInt_diffusion( double E,  double ri){			// int over Ep

	double vE = v(E);
	double Ep = (p.mx + E) /2;
	//std::cout << Ep << std::endl;
	double rootdv = sqrt(vE); //gives max value for rootdv//root_dv(Ep, vE); //
	//std::cout << rootdv << std::endl;
	//std::cout << "umax:  " << umax << std::endl;
		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> diffusionParams (4);

	diffusionParams[0] = E;
	diffusionParams[1] = ri;
	diffusionParams[2] = vE;
	diffusionParams[3] = rootdv;

	gsl_function F;
	F.function = &ddiffusion;
	F.params = &diffusionParams; 								//pass Ep to rootdv(), pass r from dndeeq as well, 
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, E, p.mx, 0, 1e-1, 1000, 
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);
		//std::cout << "result: " <<result << std::endl;
	return result;

}

double dndeeq(double E, double r ){

	/*//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	*///////before algorithm


	double rh = c.rh * mpc2cm ;
	int imNum = 4; //number of image pairs + 1, total points = 2*imNum + 1
	double diffsum = 0 ;
/*

	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);

		diffsum += pow(-1, i) * gslInt_diffusion(E, ri);

	}

	*/
	double dndeeq = (1 / c.bloss(E,r))*gslInt_diffusion(E, r);	//pow(4*pi , -1.0/2.0)* pow(c.DM_profile(r) , -2.0)*diffsum ;
	/*////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	if (duration > 0.1)*/
	//std::cout << "dndeeq(E = "<< E  << " , r = " << r/mpc2cm*1000 << " ) = " << dndeeq <<  std::endl;
	//dndeeq *= sqrt(4*pi) * pow(c.DM_profile(r) , 2.0)/ (2*imNum+1);
	return dndeeq;

}

double G(double eps, double E_gamma, double boost){

	//double boost = E/(me * pow( clight , 2))
	double GAMMA = 4 * eps * boost/me;
	double q = E_gamma / (GAMMA * (boost * me  - E_gamma) );
	

	double G = 2 * q* log(q) + (1+2*q)*( 1 - q ) + pow(GAMMA*q , 2)*(1-q) / ( 2*(1+GAMMA*q ) );//
	//if (G<0)
	//std::cout << "G( eps = " << eps <<", "<< "E_gamma = " << E_gamma << " , boost = " << boost << " ) = "<< G<<std::endl;	
	return G;
}




double IC_cross_Section( double eps, double E_gamma, double E){


	double boost = E/me;

	double sigma_thompson = 6.6524e-25;

	double sigma = 3 * sigma_thompson/(4* eps * pow(boost ,2))* G(eps, E_gamma, boost);
//if (sigma<0)
//	std::cout << "sigma( eps = " << eps <<", "<< "E_gamma = " << E_gamma << " , E = "<<E<<" ) = "<< sigma<<std::endl;
	return sigma;
}

double CMB_bbSpectrum(double eps){
	//double hplanckGeV = J2Gev*hplanck;
	//kb *= J2Gev;
	double T = 2.73;

	double nu = eps/(hplanck*J2Gev);
	//std::cout << "nu = " << nu << std::endl;
	double CMB_bbSpectrum = 8*pi* pow(nu, 2.0)/pow(clight, 3.0)	/(	exp(hplanck * nu /(kb * T )) - 1	);
	//std::cout << "nu("<< eps <<")  = " << nu <<  " spectrum = " << CMB_bbSpectrum << std::endl;;//CMB_bbSpectrum  << std::endl;

	return CMB_bbSpectrum;
}

double dpIC(double eps, void * params ){

	//double nu = 0.02;	//Ghz
	std::vector<double> pICParams = *(std::vector<double> *)params;
	double E_gamma = pICParams[0] ;
	double E = pICParams[1] ;
	

	double dpIC = CMB_bbSpectrum(eps) *  IC_cross_Section( eps, E_gamma,  E);
	//std::cout << dpIC<<std::endl;
	return dpIC;
}


double gslInt_pIC(double nu, double E){			//int over eps

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> pICParams (2);

	double E_gamma = hplanck * nu * J2Gev;
	double boost = E/me;
	pICParams[0] = E_gamma;
	pICParams[1] = E;

	double eps_max = E_gamma / (1 - E_gamma/(boost*me) );
	double eps_min = E_gamma / (4 * boost*( boost - E_gamma/me));

	gsl_function F;
	F.function = &dpIC;
	F.params = &pICParams;
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, eps_min, eps_max,  0, 1e-3, 1000, //x?
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= clight * E_gamma;
	
	//std::cout << "pIC( nu = " << nu <<", "<< "E_gamma = " << E_gamma << " , E = "<<E<<" ) = "<< result<<std::endl;
	
	return result;

}


double djIC(double E , void * params){

	std::vector<double> jICParams = *(std::vector<double> *)params;
	double nu = jICParams[0];
	double r = jICParams[1];

	double djIC = 2* gslInt_pIC(nu, E)* dndeeq(E , r);
//std::cout << "pIC = " << gslInt_pIC(nu, E) << "  dndeeq = " << dndeeq(E,r) <<std::endl;

	return djIC;
}


double gslInt_jIC(double nu, double r){ 				// int over E

			///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
		///////////


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;
	std::vector<double> jICParams (2);

	jICParams[0] = nu;
	jICParams[1] = r;

	gsl_function F;
	F.function = &djIC;
	F.params = &jICParams;

	gsl_integration_qags (&F, me, p.mx, 0, 1e-2, 1000,
	                    w, &result, &error); 


	gsl_integration_workspace_free (w);
	
	//std::cout << " . ";
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;

	//std::cout << "jIC( r = " << r/mpc2cm*1000 <<" )  = "<<result<<" duration: " << duration <<std::endl;  //      ~30s-60s	

	return result;

}


double dsIC( double r, void * params ){

	double nu = *(double *)params;


	double dist_z = Dist() / (1+c.z);

	double sICIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(c.DM_profile(r) , 2)*gslInt_jIC(nu, r);

		//std::cout << "sICIntegrand = " <<sICIntegrand<<std::endl;
	
	return sICIntegrand;
}


double gslInt_sIC( double nu, double r ){				// int over r

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;



	gsl_function F;
	F.function = &dsIC;
	F.params = &nu;

	gsl_integration_qags (&F, 1e-16, r, 0, 1e-3, 1000,
	                w, &result, &error); 






	//std::cout <<"result = " << result <<std::endl;
	result *= 0.00160218 * p.sv/(8* pi*pow( p.mx , 2.0 ));
	//	std::cout <<"new result  = " << result <<std::endl;
	return result;

}

void runFlux(double mx, double r){ //runs through values of root_dv
	
	p.mx = mx;
	int n_nu = 50;


	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm

	std::ostringstream makefilename;
	makefilename << "IC_SED_NSD_JUNK" << p.mx << "Gev_coma.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){


		double nu_min = 1e10;
		double nu_max = 1e25;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));
		double data = gslInt_sIC( nu , r)*nu ;
		file << log10(nu) << "\t" <<  data <<std::endl;
		std::cout << i<< "/" << n_nu << " " << log10(nu) << "\t" << data << std::endl;

	};
	
	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}
//NOTE that we need to set a cross section for this, use eq 1 in STorm for Snu.

main(){
	dsinit_(); //initialixe DarkSUSY
	

	double r = c.rh*mpc2cm;
	p.ch = 25;							//darksusy channel
	//std::cout << c.bfield_model(r) << std::endl;
 	runFlux( 40, r );


	p.ch = 13;							//darksusy channel
	//std::cout << c.bfield_model(r) << std::endl;
	runFlux( 81, r ) ;



}




