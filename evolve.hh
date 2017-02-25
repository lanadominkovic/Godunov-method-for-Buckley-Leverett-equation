#ifndef EVOLVE_HH
#define EVOLVE_HH

#include <cassert>
/**
 * 
 * @param gv  = LeafGridView
 * @param current = rješenje na trenutnom vremenskom sloju
 * @param next = rješenje na sljedećem vremenskom sloju
 * @param bc = vrijednost rubnog uvjeta (na sljedećem vremenskom sloju)
 * @param lambda = dt/dx.
 */

double f(double u, double a){
    return u*u/(u*u + a*(1-u)*(1-u)); 
}

double df(double u, double a){
    return (-2*a*u*u + 2*a*u)/( ((1+a)*u*u - 2*a*u + a)*((1+a)*u*u - 2*a*u + a) ); 
}

template <typename GV, typename Vec, typename VMap>
void godunov(GV const & gv, Vec const & current, Vec & next, double dt, VMap const & vmapper, double a){

    int dim = GV::dimension;
    assert(dim == 1);
    
    // po svim elementima
    auto it = gv.template begin<0>();
    for (; it != gv.template end<0>(); ++it){
        double dx = it->geometry().volume();
        double lambda = dt/dx;
         
        int left  = vmapper.subIndex(*it,0,dim);
        int right = vmapper.subIndex(*it,1,dim);
	
        next[ left ]  +=  - lambda * f(current[ left ], a);
        next[ right ] +=    lambda * f(current[ left ], a);
        
    }
   
}
#endif /* EVOLVE_HH */
