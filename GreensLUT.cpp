/////////////////////////////////////////////////////////
///				GreensLUT.cpp					  5-25-16 ///
/////////////////////////////////////////////////////////

/*
g++ -o GreensLUT GreensLUT.cpp  -lgsl -lgslcblas

*/



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




double DM_profile(double r){

    double rhos = 0.039974 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
    double rs = 0.404*mpc2cm;  //DM scale radius in Mpc Storm13  (Mathematica)
    double rc = 0.241 *mpc2cm;
    
    // NFW
    //double rho = rhos / ( r/rs * pow(1 + r/rs , 2 )); //x? + 1e-100
	

    /*
    // Berkert
    double rho = rhos / ((1+ r/rc) * (1 + pow(r/rc , 2 ))); //shouldnt use rhos here, different rho
	*/
	double rho_N04 =0.08296;		//use Colafrancesco 2006 Eq 5 and Eq 33 
	double rs_N04 = 0.28*mpc2cm; 
	// N04
	double rho = rho_N04 * exp(-2.0/0.17 * ( pow(r/rs_N04, 0.17) - 1 ));
 




    return rho;

};



double dGreens(double rp, void * params ){ 

	std::vector<double> greenParam = *(std::vector<double> *)params;
	double ri = greenParam[0];
	double r = greenParam[1];
	double root_dv = greenParam[2]; // 


	double dGreens = pow(root_dv , -1) * rp/ri * (exp( - pow( (rp-ri)/(2*root_dv) , 2)) 
		- exp( - pow( ( rp + ri)/(2*root_dv) , 2)) ) * pow( DM_profile(rp),2)/pow( DM_profile(r),2);

	return dGreens;

}



double gslInt_Greens(double ri , double r, double root_dv, double rh){


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> greenParam (3);


	greenParam[0] = ri;
	greenParam[1] = r;
	greenParam[2] = root_dv;

	gsl_function F;
	F.function = &dGreens;
	F.params = &greenParam;
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, 1e-16, rh, 0, 1e-5, 1000, //x?
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}


double Greens (double r, double root_dv) {  //called by ddsyn

	double rh = 0.415 * mpc2cm ;
	int imNum = 7; //number of image pairs
	double Gsum = 0 ;
	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);

		Gsum += pow(-1, i) * gslInt_Greens(ri, r, root_dv, rh);

	}

	double Greens = pow(4*pi , -1.0/2.0)*Gsum ;


	return Greens;

};



std::vector< std::vector<double> >  createLUT(double n_r, double n_rootdv){
	// iteration timer start
	std::clock_t Gstart;
	double Greensduration;
	Gstart = std::clock();
	int Ga ; 
	///////before algorithm
	std::cout << "creating LUT..." <<std::endl; 

	std::vector< std::vector<double> > Greens_lookup( n_r , std::vector<double>(n_rootdv) );

	double rh = 415*mpc2cm/1000;
	double r_scale = rh/n_r;

	double rootdv_max = 40*mpc2cm/1000;
	double rootdv_scale = rootdv_max/n_rootdv;
	

	std::cout << r_scale/mpc2cm*1000 << " , rdv_max = "<<rootdv_max<<" n_rootdv = " <<n_rootdv <<" rdv_scale = " << rootdv_scale/mpc2cm*1000 << std::endl; 
	

	for (int i = 1 ; i < n_r ; ++i ){
		
		double dr = rh/n_r;
		
		for(int j = 1; j < n_rootdv ; ++j){
			double drootdv = rootdv_max/n_rootdv;

			Greens_lookup[i][j] = Greens(i*dr, j*rootdv_scale);

			
		}

		std::cout << i << "/" << n_r <<std::endl;
	};
	////////after algorithm
	Greensduration = (std::clock()  -  Gstart)/(double) CLOCKS_PER_SEC;
	std::cout << "Glookup time = " << Greensduration <<std::endl;

	return Greens_lookup;
}


void runGreens(double r){ //runs through values of root_dv
	//int n = 10000;
	int n_r = 415;
	int n_rootdv = 1000;
	double root_dv = 0.04*mpc2cm;

	double rk = mpc2cm/1000*r;
	double h = root_dv/10000;

	double rh = 415*mpc2cm/1000;
	double rootdv_max = 40*mpc2cm/1000;
	double r_scale = rh/n_r;
	double rootdv_scale = rootdv_max/n_rootdv;


	std::vector< std::vector<double> >  glookup = createLUT(n_r, n_rootdv);
	std::ostringstream makefilename;
	makefilename << "Greens_NFW_r" << r << "kpc" <<"GLUTtest.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int ix = 0 ; ix < 10000 ; ++ix  ){
		double r_int  = (int)(rk/r_scale);
		double rdv_int = (int)( h*ix/rootdv_scale);
		
		double data = glookup[r_int][rdv_int];

		file << h*ix/mpc2cm*1000 << "\t" <<  data <<std::endl;
		std::cout << h*ix/mpc2cm*1000 << "\t" << data << std::endl;


		if(ix % 1000 == 0)
			std::cout << ix/1000 << "/" << 10000.0 /100000.0<< std::endl;;		


	};
	
}


void runGreens(double rdv, int rdv_switch){
	//int n = 10000;
	int n_r = 415;
	int n_rootdv = 1000;
	double root_dv = rdv*mpc2cm/1000.0;

	//double rk = mpc2cm/1000*r;
	//double h = root_dv/10000;

	double rh = 415*mpc2cm/1000;
	double rootdv_max = 40*mpc2cm/1000;
	double r_scale = rh/n_r;
	double rootdv_scale = rootdv_max/n_rootdv;


	std::vector< std::vector<double> >  glookup = createLUT(n_r, n_rootdv);
	std::ostringstream makefilename;
	makefilename << "Greens_N04_rdv" << rdv << "kpc" <<"GLUTtest.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int ir = 0 ; ir < n_r ; ++ir  ){
		double r_int  = (int)(ir*mpc2cm/1000.0/r_scale);
		double rdv_int = (int)( root_dv/rootdv_scale);
		
		double data = glookup[r_int][rdv_int];

		file << ir << "\t" <<  data <<std::endl;
		std::cout << ir << "\t" << data << std::endl;


		if(ir % 1000 == 0)
			std::cout << ir/1000 << "/" << 10000.0 /100000.0<< std::endl;;		


	};
	

}

main(){





//createLUT(415, 1000);


	runGreens(5, 1);
	std::cout << "1/5" <<std ::endl;
	runGreens(10, 1);
	std::cout << "2/5" <<std ::endl;
	runGreens(15, 1);
	std::cout << "3/5" <<std ::endl;
	runGreens(25, 1);
	std::cout << "4/5" <<std ::endl;
	runGreens(35, 1);
	std::cout << "5/5" <<std ::endl;
	//runGreens(500);
	//std::cout << "7/7" <<std ::endl;


}