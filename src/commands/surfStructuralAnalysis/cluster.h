#ifndef AIMS_CLUSTER_H
#define AIMS_CLUSTER_H
#include "meshdistance.h"
#include "cliques.h"
#include "minimization.h"

using namespace aims;
using namespace carto;
using namespace std;


class SWC: public SurfaceBased_StructuralAnalysis{
  private:
    void getListeTriangles();
    void Step(vector<int> &random, long double temp, uint &mod);
    vector<uint> getCliquesTurnedOn(float temp);
    vector<int> getCompConn(vector<uint> &indicesCliques);
  public:
    SWC(){}
    SWC(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon);
    void Run();
};




#endif

