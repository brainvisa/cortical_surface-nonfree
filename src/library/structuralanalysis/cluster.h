#ifndef AIMS_CLUSTER_H
#define AIMS_CLUSTER_H
#include <cortical_surface/structuralanalysis/meshdistance.h>
#include <cortical_surface/structuralanalysis/cliques.h>
#include <cortical_surface/structuralanalysis/minimization.h>



class SWC: public SurfaceBased_StructuralAnalysis{
  private:
    void getListeTriangles();
    void Step(std::vector<int> &random, long double temp, uint &mod);
    std::vector<uint> getCliquesTurnedOn(float temp, std::vector<uint> &indicesCliques);
    std::vector<int> getCompConn(std::vector<uint> &indicesCliques, std::set<uint> &listeSites);
    std::vector<std::set<uint> > getCompConnVector(std::vector<int> &comp);
//     void getSpecialMesh();
    long double getCompacite(std::set<uint> &comp, bool verb=true);


  public:
    SWC(){}

    void Run();
    void Run2();
};




#endif

