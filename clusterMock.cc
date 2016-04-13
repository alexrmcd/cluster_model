





#include <vector>
#include <iostream>
#include <math.h>



class Distribution {
    public:


		int nx , ny ;
		std::vector <double> value;
        std::vector <double> s;                                                         // access data elements by (x, y) = (nx*y + x)
		
		Distribution(int nx, int ny) : nx(16), ny(16), value(nx*ny, 0) , s(nx*ny, 0){}  //value ex Bfield.value= B0. s is energy spectrum at each point, usefule for particle distributions

        void print(){
            for (int ix = 0; ix < nx ; ix++){
                for ( int iy=0; iy < ny ; iy ++ ){
                    //value[ix + iy * nx] = Bfield_model(ix, iy);
                    std::cout << "(" << ix << ',' << iy << ")" << "  " << value[ix + iy * nx] << std::endl;;
                        //  std::cout << bfield.value[ix + iy*nx] << "  ";
                };
            //std::cout << std::endl;
            };
        };

};

double bfield_model (double x, double y) {


        double r =  sqrt ( x*x + y*y );
        float Bo = 5.0; //microGauss
        float rcore = 20; //core radius in kpc
        float beta = 6;
        float eta = 0.5;

        double B_field = Bo * pow(( 1 + r*r/(rcore*rcore)),(-3*beta*eta/2));// Storm et al 2013 



        return B_field;


};








class Cluster {

    public:
        int nx , ny ;
        Distribution bfield;
        Distribution  rhox;


        Cluster() : nx(16), ny(16), bfield(16,16) , rhox(16,16) {
            
            std::cout << "creating cluster..." << std::endl;

        };



       
    	//Dm profile, default to NFW
        double DM_profile(double x, double y){

            double rhos = .04 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
            double rs = 1.25 ;//DM sacale radius in cm Storm13  (Mathematica)
            
            double r = sqrt( x*x + y*y);

            double rho = rhos / ( (r + 1e-100 )/rs * pow(1 + r/rs , 2 ));


            return rho;

        };


        //Dm profile, options for optimistic or conservative NO OPTIONS ADDED YET
        double DM_profile(double x, double y, int op){

            double rhos = .04; //DM char density in Gev /cm^3 Storm13 (Mathematica)
            double rs = 1.25; //DM sacale radius in cm Storm13  (Mathematica)
            
            double r = sqrt( x*x + y*y);

            double rho = rhos / ( (r + 1e-100 )/rs * pow(1 + r/rs , 2 ));

            return rho;
        };






};





void create_cluster(){

    Cluster c;
    

    std::cout << "calculating bfield... " <<std::endl;
    for (int ix = 0; ix < c.nx ; ix++){                       
        for ( int iy=0; iy < c.ny ; iy ++ ){
            c.bfield.value[ix + iy * c.nx] = bfield_model(ix, iy);

        };

    }

    std::cout << "calculating DM_profile... " <<std::endl;
    for (int ix = 0; ix < c.nx ; ix++){                       
        for ( int iy=0; iy < c.ny ; iy ++ ){
            c.rhox.value[ix + iy * c.nx] = c.DM_profile(ix, iy);

        };

    }





    c.rhox.print();

    std::cout<< "Cluster created!" << std::endl;

}




















int main() {


    
    create_cluster();

 
};