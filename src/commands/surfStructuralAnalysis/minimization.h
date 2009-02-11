#ifndef AIMS_MINIMIZATION_H
#define AIMS_MINIMIZATION_H
#include "cliques.h"

using namespace aims;
using namespace carto;
using namespace std;

enum typesValid{
  PERMUT, BOOTSTRAP, RANDOM, NEIGHBORHOOD
};

enum typesMinim{
  ICM, ANNEAL, CUSTOM
};

class SurfaceBased_StructuralAnalysis{
  protected:
    float Esimil;
    vector<uint> ipscliques;
    
    long double energy;
    uint nbsujets;
    vector<int> labels;
    vector<pair<Point2df,Point2df> > labelsZones;
    vector<set<uint> > zonesListesBlobs;
    vector<set<uint> > listeZones; // attention les indices de listeZones sont décalés de 1 par rapport à labelsZones (à cause du label 0 qui recouvre tout l'espace 2D)
    vector<Site *> sites;
    vector<Clique> cliques;
    vector<vector<int> > cliquesDuSite;

    void ShortSummaryLabels();


  public:
    string energypath, recuitpath;
    void MinimizationSetup(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons);
    SurfaceBased_StructuralAnalysis(){}
    SurfaceBased_StructuralAnalysis(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons);
    void setModelParameters(float _ddweight=0.0, float _intrapsweight = 0.0, float _simweight=0.0, float _lsweight=0.0, float _ddx2 = 0.0, float _ddx1 = 0.0, float _ddh=0.0, float _ddweight2=0.0, float _dd2x1=0.0, float _dd2x2=0.0);
    
    void RunMinimization(int type);
    long double getLabelEnergy(int label, int type=UNKNOWN);
    long double getTypeEnergy(int type);
    long double getTotalEnergy();
//     double getTotalEnergyLastChance(uint site, uint newlabel);

    void SummaryLabels();
    void StoreToGraph(Graph &primal);
    void Initialization();


    long double getCompacite(set<uint> &liste, bool verb);
    vector<float> getPseudoSamplesPermut(vector<uint> &listeBlobs);
    vector<float> getPseudoSamplesBootstrap(vector<uint> &listeBlobs);
    vector<float> getPseudoSamplesFullBootstrap(vector<uint> &listeBlobs);

    vector<long unsigned int> creerHisto(vector<double> &samples, uint histosize, float *mini, float *step);
    void getRandomLabelsEnergies(long double nb, FILE *f);
    void Validation(int type=PERMUT);
};




#endif

