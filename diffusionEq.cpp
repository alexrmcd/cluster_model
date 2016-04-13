// 4-13-16
//
// Diffusion Equation with spatial diffusion, diffusion const independent of spatial coordinate -> D(E).
//
#include <iostream>
#include <math.h>


double bloss(double E){

	double z = 0.0232 ; //for Coma
	double me = 0.000511; //mass electron in Gev
	double n = 0.001; //electron density (all these straight from Storm code, Cola has n = 1.3e-3)

	double B = 1; //from colafrancesco estimate, should use b_field from cluster or use average, 

	double bloss = 0.0254*E*E*B*B 						//bsyn
					+ 0.25*E*E  						//bIC
					+ 1.51*n*(0.36 + log(E/me/n) )		//bbrem
					+ 6.13*( 1 + log(E/me/n)/0.75); 	//bcoul

	return bloss;
};





double U(double Emax, double E0, int ns){

	double h = (Emax - E0)/ns;
	double simp1 = 1/bloss( E0 + h/2);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += 1/bloss( E0 + h * i + h/2);
		simp2 += 1/bloss( E0 + h * i);

	}

	double result = h/6 * (1/bloss(E0) + 1/bloss(Emax) + 4 * simp1 + 2 * simp2); // in 1/(10^-16 Gev/s)

	std::cout << result << std::endl;
};



double D(double E){

	double D0 = 3.1; // 10^28 cm^2/s
	double db = 20; // no reason to be int, estimate from colafranceso for Coma
	double Bmu = 1; // in mG, estimate from colafrancesco, same thing as B in bloss

	double D = D0 * pow(db, 2/3)/pow(Bmu,1/3) * pow(E,1/3);
};



/*double Root_dv(){      // eventually should be some function of D(E), see colafranceso 


	double root_dv = 0.01; //Mpc
	return root_dv;
};*/






double DM_profile(double r){

    double rhos = .04 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
    double rs = .404 ;//DM sacale radius in Mpc Storm13  (Mathematica)


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

	double rh = 0.415; //Coma Halo radius (Storm 2013)
	double rs = 0.404; //Coma scale radius (mathematica jfactor)
	int imNum = 20; //number of image pairs
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
	
	double root_dv = 0.035;
	int n = 100;

	double h = root_dv/n;

	for (int ix = 0 ; ix < n + 1; ++ix  ){
		//Greens(0.05, h*ix, ns )

		std::cout << h*ix << ":  " << Greens(r, h*ix, ns ) << std::endl;
	};
	
}

























main(){

	int ns = 10000;
	double r = 0.01;
	runGreens(r, ns);

	//integrandTest( 0.35 , 0.05);

	//std::cout << "Greens: " << Greens(0.05, 0.01, ns) <<std::endl;
}

