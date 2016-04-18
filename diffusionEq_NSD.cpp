// 4-13-16
//
// Diffusion Equation with spatial diffusion, diffusion const independent of spatial coordinate -> D(E).
//
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <ctime>


/* to compile and run, use command 

g++ -o diffusionEq diffusionEq.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -ldarksusy -lFH -lHB -lgfortran
*/

//////////////////////////// routines for calling fortran/darksusy stuff /////////////////////////

extern"C" {	 								//interface with fortran code to initialize darksusy

void dsinit_();

}	

extern"C" {
double dshayield_(double *mwimp, double *emuthr,int *ch,  int *yieldk, int *istat);
}


//////////////////////////////////////////////////////////////////////////////////////////////////

double bfield_model (double r) {


        float Bo = 4.7; //microGauss
        float rcore = 8.9793e23 ; //0.291; //core radius in Mpc for coma
        float beta = 0.75;
        float eta = 0.5;

        double B_field = Bo * pow(( 1 + r*r/(rcore*rcore)),(-3*beta*eta/2));// Storm et al 2013 



        return B_field;


};	

double bloss(double E , double r){

	double z = 0.0232 ; //for Coma
	double me = 0.000511; //mass electron in Gev
	double n = 0.001; //electron density (all these straight from Storm code  cm^-3, Cola has n = 1.3e-3)

	double B = 1; //from colafrancesco estimate, should use b_field from cluster or use average, 

	double bloss = 0.0254*E*E*pow(bfield_model(r) , 2) 						//bsyn
					+ 0.25 * pow(1+z, 4 )*E*E  						//bIC
					+ 1.51*n*(0.36 + log(E/me/n) )		//bbrem
					+ 6.13*( 1 + log(E/me/n)/0.75); 	//bcoul

	bloss = bloss *1e-16;	

	return bloss;
};





double U(double Emax, double E0, int ns, double r){

	double h = (Emax - E0)/ns;
	double simp1 = 1/bloss( E0 + h/2, r);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += 1/bloss( E0 + h * i + h/2 , r);
		simp2 += 1/bloss( E0 + h * i , r);

	}

	double result = h/6 * (1/bloss(E0 , r) + 1/bloss(Emax , r) + 4 * simp1 + 2 * simp2); // in 1/(10^-16 Gev/s)

	std::cout << result << std::endl;
};



double D(double E){

	double D0 = 3.1; // 10^28 cm^2/s
	double db = 20; // no reason to be int, estimate from colafranceso for Coma
	double Bmu = 1; // in mG, estimate from colafrancesco, same thing as B in bloss

	double D = D0 * pow(db, 2/3)/pow(Bmu,1/3) * pow(E,1/3);
};



double distint(double z){
	double mpc2cm = 3.085678e24; 	//Mpc to cm
	double clight = 299792.458;		// km/s
	double H0 = 70.4;				// km/s/Mpc
	double OmegaM = 0.272 ; 		
	double OmegaL = 0.728 ;

	double distint =  mpc2cm* clight / ( H0 * sqrt( OmegaM * pow(1 + z , 3)  + OmegaL ) ); 

	return distint;

}


//distance as a function of redshift
double Dist(double z){ 

	double mpc2cm = 3.085678e24; 	//Mpc to cm
	double clight = 299792.458;		// km/s
	double H0 = 70.4;				// km/s/Mpc
	double OmegaM = 0.272 ; 		
	double OmegaL = 0.728 ;

	int ns = 100; // 100 works , 1000 maybe? , 10000 nogo

	double h = z/ns;
	double simp1 = distint(h/2);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += distint(h * i + h/2);
		simp2 += distint( h * i);

	}

	double result =  h/6 * (distint(0) + distint(z) + 4 * simp1 + 2 * simp2);


	



	return result;
 		
}





double DM_profile(double r){

    double rhos = .039974 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
    double rs = 1.245518e24 ;//.404 ;//DM sacale radius in Mpc Storm13  (Mathematica)


    double rho = rhos / ( (r + 1e-100 )/rs * pow(1 + r/rs , 2 ));


    return rho;

};




double GreensIntegrand(double rp, double ri, double root_dv, double r){

	double GreensIntegrand = rp/ri * (exp( - pow( (rp-ri)/(2*root_dv) , 2)) - exp( - pow( (rp+ri)/(2*root_dv) , 2)) ) * pow( DM_profile(rp),2)/pow( DM_profile(r),2);

	return GreensIntegrand;


}




double Integrate_Green(double rf, double r0, double ri ,int ns, double root_dv, double r){

	double h = (rf - r0)/ns;
	double simp1 = GreensIntegrand( r0 + h/2, ri, root_dv, r);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += GreensIntegrand( r0 + h * i + h/2,  ri, root_dv, r);
		simp2 += GreensIntegrand( r0 + h * i, ri, root_dv, r);

	}

	double result = h/6 * (GreensIntegrand(r0, ri, root_dv, r) + GreensIntegrand(rf, ri, root_dv, r) + 4 * simp1 + 2 * simp2);

	//std::cout << "Integration result: "<< result << std::endl;

	return result;
	//std::cout << result << std::endl;

}



double Greens (double r, double root_dv, int ns) {

	double mpc2cm = 3.085678e24;
	double rh = 0.415 * mpc2cm ; //Coma Halo radius (Storm 2013)
	double rs = 0.404 * mpc2cm; //Coma scale radius (mathematica jfactor)
	int imNum = 5; //number of image pairs
	double pi = 3.1415926;

	//double root_dv = Root_dv(0.02);

	double Gsum = 0 ;

	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);


		//std::cout << "r" << i << ":  " << ri <<", Gsum: "<< Gsum << std::endl;

		Gsum += pow(-1, i) * Integrate_Green(rh, 0, ri, ns, root_dv, r);


	}

	double Greens = 1/pow(4*pi , 0.5)/root_dv*Gsum ;

	//std::cout << "Gsum " << Gsum << std::endl;


	return Greens;


};


void integrandTest(double root_dv, double r){

	for (int ix = -100; ix < 100; ++ix){


		double dx = .01;

		std::cout << "Integrand (x = " << ix*dx << "): "<< GreensIntegrand(0.1, ix*dx ,  root_dv, r) << std::endl;

	}

	
}

void runGreens(double r, int ns){ //runs through values of root_dv
	
	double root_dv = 0.025;
	int n = 100;

	double h = root_dv/n;

	for (int ix = 0 ; ix < n + 1; ++ix  ){
		//Greens(0.05, h*ix, ns )

		std::cout << h*ix << ":  " << Greens(r, h*ix, ns ) << std::endl;
	};
	
}



double IntDarksusy(double mx, double E, int ns){

	
	int ch = 25;
	int yieldk = 151;
	int istat;

	double h = (mx - E)/ns;
	double E0 = E + h/2;
	double simp1 = dshayield_(&mx, &E0, &ch,  &yieldk, &istat); //f( E + h/2);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		double E1 = E + h*i + h/2;
		double E2 = E + h*i;
		simp1 += dshayield_(&mx, &E1, &ch,  &yieldk, &istat);	//f( E + h * i + h/2);
		simp2 += dshayield_(&mx, &E2, &ch,  &yieldk, &istat);	//f( E + h * i);

	}

	double result = h/6 * (dshayield_(&mx, &E, &ch,  &yieldk, &istat) + dshayield_(&mx, &mx, &ch,  &yieldk, &istat) + 4 * simp1 + 2 * simp2);

	//std::cout << " dshayield: " << dshayield_(&mx, &E, &ch,  &yieldk, &istat) <<std::endl;	
	//std::cout << "integration result: " << result << std::endl;
	return result;

}





double dndeeq(double mx, double E, double r ){
	int ns = 100	;
	int imNum = 10;
	double mpc2cm = 3.085678e24;
	double root_dv = 0.025*mpc2cm;

	//double dndeeq = 	(1 / bloss(E, r)) * Greens(r, root_dv, imNum) * IntDarksusy(mx, E, ns)		;
	double dndeeq = 	(1 / bloss(E, r))  * IntDarksusy(mx, E, ns)		;

	return dndeeq;

}



double rundndeeq(double mx, double E, double r ){
	//std::ofstream rundndeeq("dndeeq.txt");
	double n = 100;
	double h = (mx-E)/n;

	for (int ix = 0 ; ix < n + 1; ++ix  ){
		//Greens(0.05, h*ix, ns )

		std::cout << h*ix << ":  " << dndeeq(mx, h*ix, r ) << std::endl;

	};
	//rundndeeq.close();


}




//SYnchrotron emmission spectral function form Cola2006.
double fff(double x){

	double fff = 1.25 * pow( x , 1/3) * exp( -x ) * pow((648 + x*x) , 1/12);


	return fff;
}






double psynIntegrand(double theta, double E, double r){


	double bmu = 1; 	//muGauss

	double z = 0.0232;


	double nu = 1.4;	//Ghz

	double psyn0 = 1.463232e-25 ; // Gev/s/Hz
	double x0 = 62.1881 ;			// dimensionless constant
	double nu_em = ( 1 + z )* nu; // (observing freq)*(1+z)



	double x = x0 *nu_em / (bfield_model(r) * E*E) ;
	double psynIntegrand = psyn0 * bfield_model(r) * 0.5 * pow( sin(theta) , 2) * fff( x  / std::abs(sin(theta)) ); // fff(x) doesnt like negative arguments, sin(pi) ~ -2e-13, take abs for now(requires cmath)

	return psynIntegrand;


}





double psyn(double mx, double E, double r){

	double pi = 3.14159265359;
	//double dndeeq = 1;
	int ns = 100;

	double h = (pi)/ns;
	double simp1 = psynIntegrand( h/2, E ,r);

	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += psynIntegrand( h * i + h/2 , E , r);
		simp2 += psynIntegrand(  h * i , E , r);

	}

	double result = h/6 * (psynIntegrand(0 , E , r) + psynIntegrand(pi , E , r) + 4 * simp1 + 2 * simp2);
	

	double psyn = 2* result * dndeeq(mx, E, r);
	

	return psyn;

}


double jsyn(double mx, double r){


	int ns = 100;
	double me = 0.000511; //electron mass Gev

	double h = (mx - me)/ns;
	double simp1 = psyn( mx , me + h/2, r);

	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += psyn( mx , me + h * i + h/2 ,r );
		simp2 += psyn(  mx, me + h * i, r);

		//std::cout << "psyn( n = " << i << ") "<<std::endl;

	};	

	double result = h/6 * (psyn(mx, me, r) + psyn(mx, mx, r) + 4 * simp1 + 2 * simp2);
	


	double jsyn = result;
	

	return jsyn;

}





double ssynIntegrand( double mx, double r){
	double pi = 3.14159265359;
	double z = 0.0232;

	double dist_z = Dist(z) / (1+z);


	//jsyn(r) << std::endl;
	double ssynIntegrand = 4 *pi /pow(dist_z , 2) * pow(r,2) * pow(DM_profile(r) , 2)*jsyn(mx, r);	


	


	return ssynIntegrand;
}



double ssyn(double mx, double r){

	int ns = 100;
	

	double h = r/ns;
	double simp1 = ssynIntegrand(  mx, h/2);

	double simp2 = 1e-16;

	for (int i = 0 ; i < ns ; ++i){


		simp1 += ssynIntegrand( mx, h * i + h/2 );
		simp2 += ssynIntegrand(mx,  h * i + 1e-16);

		//std::cout << mx << " :"<<"ssynIntegrand( n = " << i << ") = " << simp1 <<std::endl;

	};	


	double result = h/6 * (ssynIntegrand(mx, 1e-16) + ssynIntegrand(mx, r) + 4 * simp1 + 2 * simp2);
	//std::cout << "result "<< result << std::endl;

	double ssyn = result;

	return ssyn;


}




double min_flux(double r){
	double pi = 3.14159265359;
	double z = 0.0232;
	double thetaB = 25; // beam size in arcsec
	double frms  = 1e-5; //noise per beam in Jy

	double dist_z = Dist(z) / (1+z);

	double thetaH = r/dist_z * 180/pi * 3600;

	double min_flux = 4 * log(2) * frms * pow(thetaH/thetaB, 2); 


	return min_flux;


}





double rconst(double rcm){

	double pi = 3.14159265359;
	double z = 0.0232;
	double thetaB = 25; // beam size in arcsec
	double dist_z = Dist(z) / (1+z);


	double rb = .5 * dist_z *thetaB/(3600) * (pi/180); //beam size in cm

	double rconst;

	if( rb < rcm){
		rconst = rcm;
	}

	else if(rb > rcm){
		rconst = rb;
	};


	return rconst;

}




double Calc_sv(double mx, double r){ // potentially add ch, z here?

	double pi = 3.14159265359;
	double GeVJy=1.6e20 ;   //GeV/s/Hz/cm^2 to Jy




	double Sin = ssyn(mx, r) * GeVJy ; 
	double Sout = min_flux(r);



	double sv = 8*pi * mx*mx * (Sout/Sin);


	return sv ; 

}




main(){
	//timer start
	std::clock_t start;

	double duration;

	start = std::clock();
	int a ; 

	///////before algorithm

	dsinit_(); //initialixe DarkSUSY

	int n_mx = 100 ;//number of mx values used



	double mpc2cm = 3.085678e24;
	double r = 0.415; //ROI for cluster/object (Coma here)

	double rcm = r * mpc2cm ; 

	double rmax = rconst(rcm);


	
	std::ofstream file("coma_ch25_rdv25_NSD.txt");	
	for (int i = 0 ; i < n_mx +1 ; ++i){
			// timer start
			std::clock_t start;

			double duration;
			start = std::clock();
			int a ; 

			///////before algorithm


         double mx_min = 5;
         double mx_max = 1000;

         double mx = mx_min * ( exp(    (log(mx_max) - log(mx_min))/ n_mx * i));

         file << mx << "\t" << Calc_sv(mx,rmax) <<std::endl;
         std::cout << "sv( " << mx << " ) = " << Calc_sv(mx, rmax) << std::endl;


         	////////after algorithm
			duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;

			std::cout << "time:  " << duration <<std::endl;

	}




	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;

	std::cout << "Total time:  " << duration <<std::endl;
}

