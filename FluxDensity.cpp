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

g++ -o FluxDensity FluxDensity.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
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

		double bloss = 0.0253*pow(bfield_model(r), 2)*E*E 					//bsyn bfield_model(r)
						+ 0.265 * pow(1 + z, 4 )*E*E  					//bIC
						+ 1.51*n*(0.36 + log(E/me/n) )			//+ 1.51*n*(0.36 + log(E/me) )*E						//bbrem Note this is most likely incorrect, no factor of E, and should be E/me/nin log
						+ 6.13*n*( 1 + log(E/me/n)/75); 					//bcoul

		bloss = bloss *1e-16;	

		return bloss;
	};


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

	Particle() : ch(0) , mx(0), sv(0){}

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
double GreenSum (double rp, void * params) {  //called by ddsyn

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

	double ddiffusion = darksusy(Ep) * (1.0/rootdv) * gslInt_Green(ri, rootdv);

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


	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);

		diffsum += pow(-1, i) * gslInt_diffusion(E, ri);

	}

	double dndeeq = pow(4*pi , -1.0/2.0)*(1 / c.bloss(E,r))* pow(c.DM_profile(r) , -2.0)*diffsum ;	//gslInt_diffusion(E, r);	
	/*////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	if (duration > 0.1)*/
	//std::cout << "dndeeq(E = "<< E  << " , r = " << r/mpc2cm*1000 << " ) = " << dndeeq <<  std::endl;
	
	return dndeeq;

}

//SYnchrotron emmission spectral function form Cola2006.
double fff(double x){

	double fff = 1.25 * pow( x , 1.0/3.0) * exp(-x) * pow((648 + x*x) , 1.0/12.0);

	return fff;
}

double dpsyn(double theta, void * params ){

	//double nu = 0.02;	//Ghz
	std::vector<double> psynParams = *(std::vector<double> *)params;
	double E = psynParams[0];
	double r = psynParams[1];
	double nu = psynParams[2];


	double psyn0 = 1.46323e-25 ; // Gev/s/Hz
	double x0 = 62.1881 ;			// dimensionless constant
	double nu_em = ( 1 + c.z )* nu/1000; // (observing freq)*(1+z)/100 convert from Ghz to MHz 


	//std::cout << E << " , " << r/mpc2cm*1000 << std::endl;

	double x = x0 *nu_em / ( c.bfield_model( r ) * pow( E , 2) );
	//std::cout << "x = " << x  << " B = " << c.bfield_model( r ) << std::endl;
	double dpsyn = psyn0 * c.bfield_model(r) * 0.5 * pow(  sin(theta) , 2) * fff( x  /sin(theta) ); 
	
	
	//std::cout << " dpsyn = "<< dpsyn << ", B = " << c.bfield_model(r) << std::endl;
	return dpsyn;

}

double gslInt_psyn(  double nu, double E, double r){			//int over theta

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> psynParams (3);

	psynParams[0] = E;
	psynParams[1] = r;
	psynParams[2] = nu;

	gsl_function F;
	F.function = &dpsyn;
	F.params = &psynParams;

	gsl_integration_qags (&F, 1e-16, pi, 0, 1e-3, 1000, //x?
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}


double djsyn(double E , void * params){

	std::vector<double> jsynParams = *(std::vector<double> *)params;
	double nu = jsynParams[0];
	double r = jsynParams[1];

	double djsyn = 2* gslInt_psyn(nu, E, r)* dndeeq(E , r);


	return djsyn;
}


double gslInt_jsyn(double nu, double r){ 				// int over E

			///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
		///////////


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;
	std::vector<double> jsynParams (2);

	jsynParams[0] = nu;
	jsynParams[1] = r;

	gsl_function F;
	F.function = &djsyn;
	F.params = &jsynParams;

	gsl_integration_qags (&F, me, p.mx, 0, 1e-2, 1000,
	                    w, &result, &error); 


	gsl_integration_workspace_free (w);

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "alpha = "<< c.alpha << ", jsyn( r = " << r/mpc2cm*1000 <<" ) duration: " << duration <<std::endl;  //      ~30s-60s
	//std::cout << " . ";


	return result;

}


double dssyn( double r, void * params ){

	double nu = *(double *)params;


	double dist_z = Dist() / (1+c.z);

	double ssynIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(c.DM_profile(r) , 2)*gslInt_jsyn(nu, r);	
	
	return ssynIntegrand;
}


double gslInt_ssyn( double nu, double r ){				// int over r

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;



	gsl_function F;
	F.function = &dssyn;
	F.params = &nu;

	gsl_integration_qags (&F, 1e-16, r, 0, 1e-2, 1000,
	                w, &result, &error); 
	 result *= GevJy*sv/(8* pi*pow( p.mx , 2.0 ));
	return result;

}

void runFlux(double mx, double r){ //runs through values of root_dv
	
	p.mx = mx;
	int n = 50;


	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm

	std::ostringstream makefilename;
	makefilename << "Flux_" << p.mx << "Gev_coma.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n + 1; ++i  ){


		double nu_min = 10;
		double nu_max = 10e5;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n * i));
		double data = gslInt_ssyn( nu , r);
		file << nu << "\t" <<  data <<std::endl;
		std::cout << "\n" << nu << "\t" << data << std::endl;

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




