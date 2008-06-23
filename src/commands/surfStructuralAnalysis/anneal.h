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
    Anneal(){}
    Anneal(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon);
    void Run();
    void Step(vector<int> &random, long double temp, uint &mod);
};




#endif

