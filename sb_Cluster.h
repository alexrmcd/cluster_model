#ifndef CLUSTER_H
#define CLUSTER_H

//      #include "Distribution.h"

#include<vector>
#include<iostream>



class Distribution {
	private:
		int nx , ny ;
		std::vector <double> data;
	public:
		
		Distribution(int nx, int ny) : nx(10), ny(10), data(nx*ny, 0) {}  // access data elements by (x, y) = (nx*y + x)

};



class Cluster {
    private:
        Distribution bfield_;
        Distribution energies_;

    public:
    	Cluster();
    	Distribution bfield;
    	
};



#endif