#ifndef AIMS_MINIMIZATION_H
#define AIMS_MINIMIZATION_H
#include "cliques.h"

using namespace aims;
using namespace carto;
using namespace std;



enum typesMinim{
  ICM, ANNEAL, CUSTOM
};

class SurfaceBased_StructuralAnalysis{
  protected:
    long double energy;
    uint nbsujets;
    vector<int> labels;
    vector<Site *> sites;
    vector<Clique> cliques;
    vector<vector<int> > cliquesDuSite;

    void ShortSummaryLabels();


  public:
    void MinimizationSetup(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon);
    SurfaceBased_StructuralAnalysis(){}
    SurfaceBased_StructuralAnalysis(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon);
    void setModelParameters(float _ddweight=2.0, float _intrapsweight = 10.0, float _simweight=3.0, float _lsweight=0.002, float _ddx2 = 4.0, float _ddx1 = 2.0, float _ddh=0.0001);
    
    void RunMinimization(int type);
    double getLabelEnergy(int label, int type=UNKNOWN);
    double getTypeEnergy(int type);
    double getTotalEnergy();
    void SummaryLabels();
    
    void Initialization();
};




#endif

