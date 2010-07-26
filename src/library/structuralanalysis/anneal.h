#ifndef AIMS_ANNEAL_H
#define AIMS_ANNEAL_H
#include <cortical_surface/structuralanalysis/minimization.h>
#include <cortical_surface/structuralanalysis/cliques.h>

using namespace aims;
using namespace carto;
using namespace std;


class Anneal: public SurfaceBased_StructuralAnalysis{
  private:

  public:
    Anneal(){}
    Anneal(Graph &primal);
    void Run(int verbose=1);
    void Step(vector<int> &random, long double temp, uint &mod);
};




#endif

