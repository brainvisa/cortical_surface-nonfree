#ifndef AIMS_ANNEAL_H
#define AIMS_ANNEAL_H
#include <cortical_surface/structuralanalysis/minimization.h>
#include <cortical_surface/structuralanalysis/cliques.h>


class Anneal: public SurfaceBased_StructuralAnalysis{
  private:

  public:
    Anneal(){}
    Anneal(Graph &primal);
    void Run(int verbose=1);
    void Step( std::vector<int> &random, long double temp, uint &mod );
};




#endif

