#ifndef AIMS_ICM_H
#define AIMS_ICM_H
#include <cortical_surface/structuralanalysis/minimization.h>
#include <cortical_surface/structuralanalysis/meshdistance.h>
#include <cortical_surface/structuralanalysis/cliques.h>

using namespace aims;
using namespace carto;
using namespace std;


class ICM: public SurfaceBased_StructuralAnalysis{
  private:

  public:
    ICM(){}


    void Run();
    void Step(uint &mod);

};




#endif

