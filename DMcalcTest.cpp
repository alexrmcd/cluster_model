// 4-19-16
//
// Diffusion Equation with spatial diffusion, diffusion const independent of spatial coordinate -> D(E). Uses GSL and Darksusy
//
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

g++ -o DMcalc DMcalc.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
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
	int ch ;
	double B0 ; 				//microGauss
	double rcore ; 		//0.291; //core radius in Mpc for coma
	double root_dv;

	Cluster() :	name(""),
				z(0), 						//redshift
				rh(0),						//halo radius Mpc
				ch(0), 						//darksusy channel
				B0(0), rcore(0)	,			//Bfield Params
				root_dv(0.035)					//DIffusion parameter, not really a cluster thing but easy access is good

	{	//everything labelled or Coma, sgould work in user options
		std::cout << "creating cluster... " << std::endl;
	}

	double bfield_model (double r) {

		double beta = 0.75;
		double eta = 0.5;

		double B_field = B0 * pow(( 1 + r*r/(rcore*rcore)),(-1.5*beta*eta));		// Storm et al 2013 

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
		double alpha = 1.0/3.0; //close to 1/3??
		double db = pow(20.0 , 2.0/3.0); //just a scaling factor
		double D0 = 3.1e28; // cm/s

		double D = db * D0 * pow(E, alpha)/pow(Bmu , 1.0/3.0);

		return D;
	};



};

Cluster c;


///////////////////        root_dv        ////////////////////////////////

double du(double Eu, void * params){

	double du = 1.0/c.bloss(Eu);

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

	return result;

}

double dv(double E , void * params){

	double dv = c.D(E)/c.bloss(E);

	return dv;
}


double v( double E,  double mx , double umax){
	
	//double max =  U(E, mx);
	//std::cout << "in v: "<< max <<std::endl;

	double min = 0;


		gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	gsl_function F;
	F.function = &dv;


	gsl_integration_qags (&F, min, umax, 0, 1e-3, 5000, 
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= 1e16;
	return result;
}

double v( double E,  double mx ){
	
	double max =  U(E, mx);
	//std::cout << "in v: "<< max <<std::endl;

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

double root_dv(double E, double Ep, double mx, double umax ){

	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm

	double root_dv  =   sqrt( ( v( E, mx, umax) - v(Ep, mx))   ) ;

	double root_dvcm  = root_dv/mpc2cm *1000;

	////////////
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	//std::cout << "root_dv:  " << root_dvcm <<std::endl;

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
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	gsl_function F;
	F.function = &ddist;

	gsl_integration_qags (&F, 0, c.z, 0, 1e-3, 5000,
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


double dGreens(double rp,  double ri , double r , double root_dv){ 


	double dGreens = ( 1/root_dv ) * rp/ri * (exp( - pow( (rp-ri)/( 2*root_dv ) , 2)) 
		- exp( - pow( ( rp + ri)/(2*root_dv) , 2)) ) * pow( c.DM_profile(rp),2)/pow( c.DM_profile(r),2);

	return dGreens;

}

double GreenSum (double rp, void * params) {  //called by ddsyn

	std::vector<double> greenParam = *(std::vector<double> *)params;


	double r  = greenParam[0];
	double root_dv  = greenParam[1];


	double rh = c.rh * mpc2cm ;
	int imNum = 10; //number of image pairs
	double Gsum = 0 ;

	/*	///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	*/	///////////


	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);

		Gsum += pow(-1, i) * dGreens(rp, ri, r, root_dv);

	}

	double GreenSum = pow(4*pi , -1.0/2.0)*Gsum ;

	/*	///////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "greens duration: " << duration << std::endl;
	*/	///////

	return GreenSum;

}


double gslInt_GreenSum(double r,  double root_dv){


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	std::vector<double> greenParam (3);

	double rh = c.rh*mpc2cm;
	greenParam[0] = r ;
	greenParam[1] = root_dv ;

	gsl_function F;
	F.function = &GreenSum;
	F.params = &greenParam;

	gsl_integration_qags (&F, 1e-16, rh, 0, 1e-3, 5000, //x?
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}





double darksusy (double Ep, double mx){

	int yieldk = 151;
	int istat;

	double ds = dshayield_(&mx, &Ep, &c.ch, &yieldk, &istat);

	return ds;
}

double ddiffusion(double Ep, void * params){
	std::vector<double> diffusionParams = *(std::vector<double> *)params;





	double E = diffusionParams[0];
	double mx = diffusionParams[1] ;
	double r = diffusionParams[2] ;
	double umax = diffusionParams[3];








	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm


	double rootdv = root_dv(E, Ep, mx, umax);//	0.035*mpc2cm ; // 





	double ddiffusion = darksusy(Ep, mx)* gslInt_GreenSum(r, rootdv) ;//;
		////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	//std::cout << "ddiffusion:  " << duration <<std::endl;

	return ddiffusion;

}

double gslInt_diffusion( double E, double mx , double r){

	double umax = U(E,mx);
	//std::cout << "umax:  " << umax << std::endl;
		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	std::vector<double> diffusionParams (4);

	diffusionParams[0] = E;
	diffusionParams[1] = mx;
	diffusionParams[2] = r;
	diffusionParams[3] = umax;

	gsl_function F;
	F.function = &ddiffusion;
	F.params = &diffusionParams; 								//pass Ep to rootdv(), pass r from dndeeq as well, 
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, E, mx, 0, 1e-3, 5000, 
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);
		//std::cout << "result: " <<result << std::endl;
	return result;

}

double dndeeq(double mx, double E, double r ){



	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm

	double dndeeq = 	(1 / c.bloss(E,r))* gslInt_diffusion(E, mx, r);	



		////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "dndeeq time:  " << duration <<std::endl;

	return dndeeq;

}

//SYnchrotron emmission spectral function form Cola2006.
double fff(double x){

	double fff = 1.25 * pow( x , 1.0/3.0) * exp( -x )* pow((648 + x*x) , 1.0/12.0);

	return fff;
}

double dpsyn(double theta, void * params ){

	double nu = 1.4;	//Ghz

	double psyn0 = 1.46323e-25 ; // Gev/s/Hz
	double x0 = 62.1881 ;			// dimensionless constant
	double nu_em = ( 1 + c.z )* nu; // (observing freq)*(1+z)

	std::vector<double> beta = *(std::vector<double> *)params;

	double x = x0 *nu_em / ( c.bfield_model( beta[1] ) * pow( beta[0] , 2) );
	
	double dpsyn = psyn0 * c.bfield_model(beta[1]) * 0.5 * pow(  sin(theta) , 2)* fff( x  /sin(theta) ); 
	
	
	
	return dpsyn;

}

double gslInt_psyn(  double E, double r){

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	std::vector<double> beta (2);

	beta[0] = E;
	beta[1] = r;


	gsl_function F;
	F.function = &dpsyn;
	F.params = &beta;

	gsl_integration_qags (&F, 1e-16, pi, 0, 1e-3, 5000, //x?
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}


double djsyn(double E , void * params){

	std::vector<double> alpha = *(std::vector<double> *)params;


	double djsyn = 2* gslInt_psyn(E, alpha[1])* dndeeq(alpha[0], E , alpha[1]);


	return djsyn;
}


double gslInt_jsyn(double mx, double r){

			///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
		///////////


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	std::vector<double> alpha (2);

	alpha[0] = mx;
	alpha[1] = r;


	gsl_function F;
	F.function = &djsyn;
	F.params = &alpha;

	gsl_integration_qags (&F, me, mx, 0, 1e-3, 5000,
	                    w, &result, &error); 


	gsl_integration_workspace_free (w);

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "jsyn duration: " << duration <<std::endl;  //      ~30s-60s

	std::cout << "ints: " << w -> size << std::endl;


	return result;

}


double dssyn( double r, void * params ){

	double mx = *(double *)params;

	double dist_z = Dist() / (1+c.z);

	double ssynIntegrand = 4 *pi /pow(dist_z , 2) * pow(r,2)  *  pow(c.DM_profile(r) , 2)*gslInt_jsyn(mx, r);	
	
	return ssynIntegrand;
}


double gslInt_ssyn( double r , double mx){

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;



	gsl_function F;
	F.function = &dssyn;
	F.params = &mx;

	gsl_integration_qags (&F, 1e-16, r, 0, 1e-3, 5000,
	                w, &result, &error); 

	return result;

}


double min_flux(double r){

	double thetaB = 25.0; // beam size in arcsec
	double frms  = 1e-5; //noise per beam in Jy

	double dist_z = Dist() / (1.0 + c.z);

	double thetaH = r/dist_z * 180.0/pi * 3600.0;

	double min_flux = 4.0 * log(2.0) * frms * pow(thetaH/thetaB, 2.0); 

	return min_flux;
}


double Calc_sv(double mx, double r){ // potentially add ch, z here?

	double Sin = gslInt_ssyn(r, mx) * GeVJy ; 
	double Sout = min_flux(r);

	double sv = 8*pi * mx*mx * (Sout/Sin);

	return sv ; 

}

void runComa(int ch){
	c.name = "Coma";
	c.ch = ch;	 						//darksusy channel
	c.z = 0.0232; 						//redshift
	c.rh = 0.415;						//halo radius Mpc
	c.B0 = 4.7;							//	
	c.rcore = 0.291*mpc2cm;				//


	double rcm = c.rh * mpc2cm ; 
	double rmax = rconst(rcm);

	std::string channel;

	if(c.ch == 13){
		channel = "WW";
	}
	else if(c.ch == 15){
		channel = "ee";
	}
	else if(c.ch == 17){
		channel = "mumu";
	}
	else if(c.ch == 19){
		channel = "tt";
	}
	else if(c.ch == 25){
		channel = "bb";
	};

	/*
	std::ostringstream makefilename;
	makefilename << c.name << "_" << channel << "DV.txt" ;
	std::string filename = makefilename.str();
	*/
	std::cout << "Running mx = 50" << std::endl;
	std::cout << "mx = 50" << Calc_sv(50, rmax) <<std::endl;

	int n_mx = 100 ;//number of mx values used

	//std::ofstream file(filename.c_str());
	for (int i = 0 ; i < n_mx +1 ; ++i){

		// iteration timer start
		std::clock_t start;
		double duration;
		start = std::clock();
		int a ; 
		///////before algorithm

			double mx_min = 5;
			double mx_max = 1000;

			double mx = mx_min * ( exp(    (log(mx_max) - log(mx_min))/ n_mx * i));

			//file << mx << "\t" << Calc_sv(mx,rmax) <<std::endl;
			std::cout << "sv( " << mx << " ) = " << Calc_sv(mx, rmax) << std::endl;
		////////after algorithm
		duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
		std::cout << c.ch << " time = " << i << " " << duration <<std::endl;

	}

//end runComa()
}


main(){
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm
	
		dsinit_(); //initialixe DarkSUSY

		/*runComa(13);
		runComa(15);
		runComa(17);
		runComa(19);*/
		runComa(25);

	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}









/*
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm




	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;




	*/