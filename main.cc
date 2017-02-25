
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <string>
#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>  
#include <dune/common/exceptions.hh>          
 
#include <dune/grid/onedgrid.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include "driver.hh"

int main(int argc, char** argv)
{
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    // ==============  Čitanje ulaznih podataka =========================
    // defaultno ime ulazne datoteke
    std::string filename("1D.input");
    // ime ulazne datoteke možemo dati i kao argument komandne linije
    if (argc > 1) filename = argv[1];
    
    Dune::ParameterTree input_data;  // ulazni podaci/parametri
    try {
        Dune::ParameterTreeParser::readINITree(filename, input_data);
     }
     catch (...) {
           std::cerr << "The configuration file \"" << filename << "\" "
                        "could not be read. Exiting..." << std::endl;
	       std::exit(1);
      }
    
    int n     	 = input_data.get<int>("N");
    double left  = input_data.get<double>("L");
    double right = input_data.get<double>("R");
    double T     = input_data.get<double>("T");
    double dt    = input_data.get<double>("dt");
    double bc	 = input_data.get<double>("bc");
    double in    = input_data.get<double>("in");
    double a     = input_data.get<double>("a");
    // ==============  Kraj čitanje ulaznih podataka =========================
     

    // Konstrukcija mreže
    using Grid = Dune::OneDGrid; 
    Grid grid(n, left, right);
    auto gv = grid.leafGridView();  
    
    // simulacija se vrši u driver rutini
    driver(gv, T, dt, bc, in, a); 
    
    return 0;
}
