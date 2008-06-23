#ifndef AIMS_ICM_H
#define AIMS_ICM_H
#include "minimization.h"
#include "meshdistance.h"
#include "cliques.h"

using namespace aims;
using namespace carto;
using namespace std;


class ICM: public SurfaceBased_StructuralAnalysis{
  private:

  public:
    ICM(){}
    ICM(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon);

    void Run();
    void Step(uint &mod);

};




#endif

