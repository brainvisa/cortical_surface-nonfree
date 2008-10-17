#ifndef AIMS_ANNEAL_H
#define AIMS_ANNEAL_H
#include "minimization.h"
#include "meshdistance.h"
#include "cliques.h"

using namespace aims;
using namespace carto;
using namespace std;


class Anneal: public SurfaceBased_StructuralAnalysis{
  private:

  public:
    vector<uint> ipscliques;
    Anneal(){}
    Anneal(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons);
    void Run();
    void Step(vector<int> &random, long double temp, uint &mod);
};




#endif

