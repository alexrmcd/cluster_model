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

};

Cluster c;

double bfield_model (double r) {

	double beta = 0.75;
	double eta = 0.5;

	double B_field = c.B0 * pow(( 1 + r*r/(c.rcore*c.rcore)),(-1.5*beta*eta));		// Storm et al 2013 

	return B_field;

};	



double bloss(double E , double r){

	double bloss = 0.0253*pow(bfield_model(r), 2)*E*E 					//bsyn bfield_model(r)
					+ 0.265 * pow(1 + c.z, 4 )*E*E  					//bIC
					+ 1.51*n*(0.36 + log(E/me/n) )			//+ 1.51*n*(0.36 + log(E/me) )*E						//bbrem Note this is most likely incorrect, no factor of E, and should be E/me/nin log
					+ 6.13*n*( 1 + log(E/me/n)/75); 					//bcoul

	bloss = bloss *1e-16;	

	return bloss;
};




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


double DM_profile(double r){

    double rhos = 0.039974 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
    double rs = 0.404*mpc2cm;  //DM scale radius in Mpc Storm13  (Mathematica)


    double rho = rhos / ( r/rs * pow(1 + r/rs , 2 )); //x? + 1e-100

    return rho;

};



double dGreens(double rp, void * params ){ 

	std::vector<double> greenParam = *(std::vector<double> *)params;

	double root_dv = c.root_dv*mpc2cm; // 
	double dGreens = pow(root_dv , -1) * rp/greenParam[0] * (exp( - pow( (rp-greenParam[0])/(2*root_dv) , 2)) 
		- exp( - pow( ( rp + greenParam[0])/(2*root_dv) , 2)) ) * pow( DM_profile(rp),2)/pow( DM_profile(greenParam[1]),2);

	return dGreens;

}



double gslInt_Greens(double ri , double r, double rh){


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	std::vector<double> greenParam (2);


	greenParam[0] = ri;
	greenParam[1] = r;

	gsl_function F;
	F.function = &dGreens;
	F.params = &greenParam;

	gsl_integration_qags (&F, 1e-16, rh, 0, 1e-3, 5000, //x?
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}


double Greens (double r) {  //called by ddsyn

	double rh = c.rh * mpc2cm ;
	int imNum = 20; //number of image pairs
	double Gsum = 0 ;

	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);

		Gsum += pow(-1, i) * gslInt_Greens(ri, r, rh);

	}

	double Greens = pow(4*pi , -1.0/2.0)*Gsum ;
	//std::cout << Greens << std::endl;

	return Greens;

};



double ddarksusy (double Ep, void * params){

   int yieldk = 151;
   int istat;
   double mx = *(double *)params; 
   
   double ds = dshayield_(&mx, &Ep, &c.ch, &yieldk, &istat);


   return ds;
}



double gslInt_ds( double Ep, double mx){
		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error;

	gsl_function F;
	F.function = &ddarksusy;
	F.params = &mx; 								//pass Ep to rootdv(), pass r from dndeeq as well, 

	gsl_integration_qags (&F, Ep, mx, 0, 1e-3, 5000, 
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);
	return result;

}



double dndeeq(double mx, double E, double r ){

	double dndeeq = 	(1 / bloss(E,r))  * gslInt_ds(E, mx);		; 

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

	double x = x0 *nu_em / ( bfield_model( beta[1] ) * pow( beta[0] , 2) );
	
	double dpsyn = psyn0 * bfield_model(beta[1]) * 0.5 * pow(  sin(theta) , 2)* fff( x  /sin(theta) ); 
	
	// fff(x) doesnt like negative arguments, sin(pi) ~ -2e-13, take abs for now(requires cmath) 4-22-16 removed abs by using less precisepi value	
	
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

	//std::cout << "psyn(5): " << result <<std::endl; 

	return result;

}


double djsyn(double E , void * params){

	std::vector<double> alpha = *(std::vector<double> *)params;


	double djsyn = 2* gslInt_psyn(E, alpha[1])* dndeeq(alpha[0], E , alpha[1]);

   	//std::cout << "E: " << E <<std::endl;
	return djsyn;
}


double gslInt_jsyn(double mx, double r){
		//std::cout << "jsyn " <<std::endl;
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

	return result;

	//std::cout << "gslInt_jsyn: " << std::setprecision(19) << result << std::endl;
}


double dssyn( double r, void * params ){

	double mx = *(double *)params;

	double dist_z = Dist() / (1+c.z);

	double ssynIntegrand = 4 *pi /pow(dist_z , 2) * pow(r,2) *Greens(r) *  pow(DM_profile(r) , 2)*gslInt_jsyn(mx, r);	
	//double ssynIntegrand = 4 *pi /pow(dist_z , 2) * r * r * pow(DM_profile(r) , 2)*gslInt_jsyn(mx, r);	
	return ssynIntegrand;
}


double gslInt_ssyn( double r , double mx){
		//std::cout <<"ssyn" <<std::endl;
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
		std::cout << channel << std::endl;
	}
	else if(c.ch == 15){
		channel = "ee";
		std::cout << channel << std::endl;
	}
	else if(c.ch == 17){
		channel = "mumu";
		std::cout << channel << std::endl;
	}
	else if(c.ch == 19){
		
		channel = "tt";
		std::cout << channel << std::endl;
	}
	else if(c.ch == 25){
		
		channel = "bb";
		std::cout << channel << std::endl;
	};

	
	std::ostringstream makefilename;
	makefilename << c.name << "_" << channel << ".txt" ;
	std::string filename = makefilename.str();
	


	int n_mx = 100 ;//number of mx values used

	std::ofstream file(filename.c_str());
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

			file << mx << "\t" << Calc_sv(mx,rmax) <<std::endl;
			//std::cout << "sv( " << mx << " ) = " << Calc_sv(mx, rmax) << std::endl;
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

		runComa(13);
		runComa(15);
		runComa(17);
		runComa(19);
		runComa(25);

	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}

