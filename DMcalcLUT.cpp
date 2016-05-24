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

	double B0 ; 				//microGauss
	double rcore ; 		//0.291; //core radius in Mpc for coma
	double root_dv;

	double alpha;		// D(E) ~ E^alpha

	double n;
	std::vector<double> vlookup;

	Cluster() :	name(""),
				z(0), 						//redshift
				rh(0),						//halo radius Mpc
				B0(0), rcore(0)	,			//Bfield Params
				root_dv(0.035)	,
				alpha(1.0/3.0)	,			//DIffusion parameter, not really a cluster thing but easy access is good
				n(1000000),
				vlookup(n)



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
		double ne = 1e-3;

		double bloss = 0.0253*pow(bfield_model(r), 2)*E*E 					//bsyn bfield_model(r)
						+ 0.265 * pow(1 + z, 4 )*E*E  					//bIC
						+ 1.51*ne*(0.36 + log(E/me/ne) )			//+ 1.51*n*(0.36 + log(E/me) )*E						//bbrem Note this is most likely incorrect, no factor of E, and should be E/me/nin log
						+ 6.13*ne*( 1 + log(E/me/ne)/75); 					//bcoul

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
	//std::cout << "rdv = " << root_dv/mpc2cm*1000 << std::endl;
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
	//double vE = diffusionParams[2];
	double rootdv = diffusionParams[2];

	//double rootdv = root_dv( Ep, vE); // 0.035*mpc2cm ; //  	//		 

	double ddiffusion = darksusy(Ep) * (1.0/rootdv)  * gslInt_Green(ri, rootdv);

	return ddiffusion;

}

double gslInt_diffusion( double E,  double ri){			// int over Ep

	//double vE = v(E);

	//double Ep = (p.mx + E) /2;
	//std::cout << E << std::endl;
	//double scale = 0.001;
	//double Es = (int)(E*c.n/p.mx);
	//std::cout <<  sqrt(vE)/mpc2cm*1000 << "   "<<sqrt(c.vlookup[ E*c.n/p.mx ])/mpc2cm*1000 <<std::endl;;
	//std::cout << Es << " LUT = " << sqrt(c.vlookup[ Es ])/mpc2cm*1000 << std::endl;
	double rootdv =sqrt(c.vlookup[E*c.n/p.mx])/mpc2cm ;//0.035*mpc2cm ;//sqrt(c.vlookup[ Es ]); //gives max value for rootdv//root_dv(Ep, vE); //

	//std::cout << "umax:  " << umax << std::endl;
		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> diffusionParams (3);

	diffusionParams[0] = E;
	diffusionParams[1] = ri;
	//diffusionParams[2] = vE;
	diffusionParams[2] = rootdv;

	gsl_function F;
	F.function = &ddiffusion;
	F.params = &diffusionParams; 								//pass Ep to rootdv(), pass r from dndeeq as well, 
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, E, p.mx, 0, 1e-1, 1000, 
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

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
	if (duration > 0.1)
	std::cout << "dndeeq(E = "<< E  << " , r = " << r/mpc2cm*1000 << " ) --> " << duration << std::endl;
	*/

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

	std::vector<double> psynParams = *(std::vector<double> *)params;

	double x = x0 *nu_em / ( c.bfield_model( psynParams[1] ) * pow( psynParams[0] , 2)	 );
	//std::cout << "x = " << x  << " B = " << c.bfield_model(psynParams[1]) << std::endl;
	double dpsyn = psyn0 * c.bfield_model(psynParams[1]) * 0.5 * pow(  sin(theta) , 2)* fff( x  /sin(theta) ); 
	
	

	return dpsyn;

}

double gslInt_psyn(  double E, double r){			//int over theta

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> psynParams (2);

	psynParams[0] = E;
	psynParams[1] = r;


	gsl_function F;
	F.function = &dpsyn;
	F.params = &psynParams;

	gsl_integration_qags (&F, 1e-16, pi, 0, 1e-3, 1000, //x?
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}


double djsyn(double E , void * params){

	double r = *(double *)params;


	double djsyn = 2* gslInt_psyn(E, r)* dndeeq(E , r);


	return djsyn;
}


double gslInt_jsyn(double r){ 				// int over E

		///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
		///////////


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;


	gsl_function F;
	F.function = &djsyn;
	F.params = &r;

	gsl_integration_qags (&F, me, p.mx, 0, 1e-2, 1000,
	                    w, &result, &error); 


	gsl_integration_workspace_free (w);

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;

	std::cout << "alpha = "<< c.alpha << ", jsyn( r = " << r/mpc2cm*1000 <<" ) = "<< result <<"duration: " << duration <<std::endl;  //      ~30s-60s



	return result;

}


double dssyn( double r, void * params ){

	double dist_z = Dist() / (1+c.z);

	double ssynIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(c.DM_profile(r) , 2)*gslInt_jsyn(r);	
	
	return ssynIntegrand;
}


double gslInt_ssyn( double r ){				// int over r

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;



	gsl_function F;
	F.function = &dssyn;
	//F.params = &mx;

	gsl_integration_qags (&F, 1e-16, r, 0, 1e-2, 1000,
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


double Calc_sv(double r){ // potentially add ch, z here?
	
	double Sin =  gslInt_ssyn(r) * GeVJy ; 
	double Sout = min_flux(r);

	double sv = 8*pi * pow(p.mx, 2) * (Sout/Sin);

	return sv ; 

}

void runComa(int ch){
	
	p.ch = ch;							//darksusy channel

	c.name = "Coma"; 						
	c.z = 0.0232; 						//redshift
	c.rh = 0.415;						//halo radius Mpc
	c.B0 = 4.7;			
	//c.alpha = 1.0/3.0;			//	
	c.rcore = 0.291*mpc2cm;				//

	double rcm = c.rh * mpc2cm ; 
	double rmax = rconst(rcm);

	double mx_min = 5;
	double mx_max = 1000;
	double data;

	std::string channel;

	if(p.ch == 13){
		channel = "WW";
	}
	else if(p.ch == 15){
		channel = "ee";
	}
	else if(p.ch == 17){
		channel = "mumu";
	}
	else if(p.ch == 19){
		channel = "tt";
	}
	else if(p.ch == 25){
		channel = "bb";
	};

	
	std::ostringstream makefilename;
	makefilename << c.name << "_express_" << channel << "_alpha_" <<c.alpha <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());


	int n_mx = 35 ;//number of mx values used




	for (int i = 0 ; i < n_mx + 1 ; ++i){

		// iteration timer start
		std::clock_t start;
		double duration;
		start = std::clock();
		int a ; 
		///////before algorithm

			p.mx = mx_min * ( exp(    (log(mx_max) - log(mx_min))/ n_mx * i));

				// iteration timer start
				std::clock_t vstart;
				double vduration;
				vstart = std::clock();
				int va ; 
				///////before algorithm
				std::cout << "creating LUT..." <<std::endl; 
				//std::cout << root_dv(E, 5 , mx) << std::endl; 
				for (int j = 0 ; j < c.n +1; ++j ){

					double dE = p.mx/c.n;
					//std::cout << j << "  " << j*dE<< std::endl;
					c.vlookup[j] = v(j*dE);
					//if (j % 1000 == 0)
					//	std::cout << "." ; 
					//if(c.vlookup[j] < 0)
					//std::cout <<p.mx << " rootv(E = " << j*dE << ") = " << sqrt(c.vlookup[j])/mpc2cm*1000 << std::endl;
				}
				//std::cout << "vlookup created..." << c.vlookup[] <<std::endl; 
				////////after algorithm
				vduration = (std::clock()  -  vstart)/(double) CLOCKS_PER_SEC;
				std::cout << "vlookup time = " << vduration <<std::endl;

			data = Calc_sv(rmax);
			file << p.mx << "\t" <<  data <<std::endl;
			std::cout << "sv( " << p.mx << " ) = " << data << std::endl;

		////////after algorithm
		duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
		std::cout << i << "/"<< n_mx << " ";
		std::cout << p.ch << ", alpha = " << c.alpha << ", time = " << duration <<std::endl;


	};

//end runComa()
}


void alphaPrint(){

	std::cout << c.alpha << std::endl;
}


main(){
	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	int a ; 
	///////before algorithm
	
		dsinit_(); //initialixe DarkSUSY

		//runComa(13);
		//runComa(15);
		//runComa(17);
		//runComa(19);
		//runComa(25);

		c.alpha = 1.0/3.0;
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