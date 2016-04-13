





#include <vector>
#include <iostream>
#include <math.h>



class Distribution {
    public:


		int nr ;
		std::vector <double> value;
        std::vector <double> s;                                                         
		
		Distribution(int nr) : nr(16), value(nr, 0) , s(nr, 0){}  //value ex Bfield.value= B0. s is energy spectrum at each point, useful for particle distributions

        void print(){
            for (int ir = 0; ir < nr ; ++ir){
                    //value[ix + iy * nx] = Bfield_model(ix, iy);
                    std::cout << "( r = " << ir << ")" << "  " << value[ir] << std::endl;;
                        //  std::cout << bfield.value[ix + iy*nx] << "  ";
            //std::cout << std::endl;
            };
        };

};

double bfield_model (double r) {


        float Bo = 5.0; //microGauss
        float rcore = 20; //core radius in kpc
        float beta = 6;
        float eta = 0.5;

        double B_field = Bo * pow(( 1 + r*r/(rcore*rcore)),(-3*beta*eta/2));// Storm et al 2013 



        return B_field;


};








class Cluster {

    public:
        int nr;
        Distribution bfield;
        Distribution  rhox;


        Cluster() : nr(16), bfield(16) , rhox(16) {
            
            std::cout << "creating cluster..." << std::endl;

        };



       
    	//Dm profile, default to NFW
        double DM_profile(double r){

            double rhos = .04 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
            double rs = 1.25 ;//DM sacale radius in cm Storm13  (Mathematica)


            double rho = rhos / ( (r + 1e-100 )/rs * pow(1 + r/rs , 2 ));


            return rho;

        };


        //Dm profile, options for optimistic or conservative NO OPTIONS ADDED YET
        double DM_profile(double r, int op){

            double rhos = .04; //DM char density in Gev /cm^3 Storm13 (Mathematica)
            double rs = 1.25; //DM sacale radius in cm Storm13  (Mathematica)
            


            double rho = rhos / ( (r + 1e-100 )/rs * pow(1 + r/rs , 2 ));

            return rho;
        };






};





void create_cluster(){

    Cluster c;
    

    std::cout << "calculating bfield... " <<std::endl;
    for (int ir = 0; ir < c.nr ; ++ir){                       

            c.bfield.value[ir] = bfield_model(ir);

    }

    std::cout << "calculating DM_profile... " <<std::endl;
    for (int ir = 0; ir < c.nr ; ++ir){                       
        
            c.rhox.value[ir] = c.DM_profile(ir);

        

    }


     c.rhox.print();


    

    std::cout<< "Cluster created!" << std::endl;

}





void Diffusion(){

  double D0 = 10; //not realistic value at all... Check colafracesco

  Distribution dnde(16); //initialize dnde

  dnde.value[15] = 0; // boundary condition, vanish at rh
  
  double dr = .1; //step size in r


  for ( int ir = 0 ; ir < 16 ; ++ir ){
    

  }




}














int main() {


    
    create_cluster();

 
};