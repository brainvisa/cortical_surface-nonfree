#ifndef AIMS_CLUSTER_H
#define AIMS_CLUSTER_H
#include <cortical_surface/structuralanalysis/meshdistance.h>
#include <cortical_surface/structuralanalysis/cliques.h>
#include <cortical_surface/structuralanalysis/minimization.h>

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

    void Run();
    void Run2();
};




#endif

