#ifndef CLUSTER_H
#define CLUSTER_H

#include "Distribution.h"

#include<vector>

class Cluster {
    private:
        Distribution bfield_;
        Distribution energies_;

    public:
    	Cluster(): Cluster(Distribution(10,10)) {}
    	


};


#endif