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
    vector<uint> getCliquesTurnedOn(float temp, vector<uint> &indicesCliques);
    vector<int> getCompConn(vector<uint> &indicesCliques, set<uint> &listeSites);
    vector<set<uint> > getCompConnVector(vector<int> &comp);
//     void getSpecialMesh();
    long double getCompacite(set<uint> &comp, bool verb=true);


  public:
    SWC(){}
    SWC(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons);
    void Run();
    void Run2();
};




#endif

