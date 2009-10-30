#ifndef AIMS_OLDANNEAL_H
#define AIMS_OLDANNEAL_H
#include <cortical_surface/structuralanalysis/minimization.h>
#include <cortical_surface/structuralanalysis/meshdistance.h>
#include <cortical_surface/structuralanalysis/cliques.h>

using namespace aims;
using namespace carto;
using namespace std;


class OldAnneal: public SurfaceBased_StructuralAnalysis{
  private:

  public:
    OldAnneal(){}
    OldAnneal(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons);
    void Run(int verbose=1);
    void Step(vector<int> &random, long double temp, uint &mod);
};




#endif

