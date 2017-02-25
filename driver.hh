#ifndef DRIVER_HH
#define DRIVER_HH


#include <vector>
#include <cassert>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "evolve.hh"

/**
 * 
 * @param gv = leaf grid view
 * @param T  = krajnje vrijeme simulacije
 * @param dt = vremenski korak
 * @param vel = polje brzine (konstanta)
 */


template <typename GV>
void driver(GV & gv, double T, double dt, double bc, double in, double a)
{
  const int dim = GV::dimension; // dim = 1
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<typename GV::Grid, Dune::MCMGElementLayout> el_mapper(gv.grid());
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<typename GV::Grid, Dune::MCMGVertexLayout>  ve_mapper(gv.grid());
  
  assert( el_mapper.size() + 1 == ve_mapper.size() );
  int N = ve_mapper.size(); // broj vrhova
  
  // current = u^n, next = u^{n+1}
  std::vector<double> current, next, shock;
  current.resize(N);
  next.resize(N);
  
  
  // Inicijalizija početnim uvjetom -> s njim dobijemo current i pomoću njega računamo next
  int index = 0;
  for(auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it){ 
      current[index] = in;
      index++;
  }

  if(bc==1.0){ // za bc=1 znamo izracunati u*
        shock.resize(N);
  	double u = std::sqrt(a/(1+a));  
  	for (int i = 0; i < shock.size(); i++ ) shock[i] = u;
  }

  if(bc > in and bc == 1) std::cout << "Brzina soka=" << df( shock[0], a) << std::endl;

  next = current;
  double time = 0.0;
  int k_out = 0;

   // Vremenska petlja
  while( time <= T){

      // evolucija sheme u novi vremenski korak
	godunov(gv, current, next, dt, ve_mapper, a);
	next[0] = bc;
	int n = next.size() - 1;
	next[n] = in;

	// priprema za sljedeću iteraciju
	time += dt;
	current = next;
	
	// Ispis trenutne iteracije
	Dune::VTKWriter<GV> vtkwriter(gv);
	std::string fname("godunov-"); 
	fname += std::to_string(k_out);   
	vtkwriter.addVertexData(current, "numerical");
	if (bc == 1.0)	vtkwriter.addVertexData(shock, "u*");
	vtkwriter.write(fname, Dune::VTK::ascii );
	k_out++;
  }
  
}

#endif /* DRIVER_HH */

