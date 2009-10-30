#ifndef AIMS_VALIDATION_H
#define AIMS_VALIDATION_H
#include <cortical_surface/structuralanalysis/minimization.h>

using namespace aims;
using namespace carto;
using namespace std;

enum typesValid{
  PERMUT, BOOTSTRAP, RANDOM, CLUSTERS
};

enum typeHisto{
  HORIZONTAL,VERTICAL
};

class StructuralAnalysis_Validation{
  protected:
    SurfaceBased_StructuralAnalysis *ssb;
    
  public:
    StructuralAnalysis_Validation(SurfaceBased_StructuralAnalysis *ssb2){ssb = ssb2; }
    vector<int> getCompConn(vector<uint> &indicesCliques, set<uint> &listeSites);
    vector<set<uint> > getCompConnVector(vector<int> &comp);

//     long double getCompacite(set<uint> &liste, bool verb);
//     vector<float> getPseudoSamplesPermut(vector<uint> &listeBlobs);
    vector<double> getPseudoSamplesBootstrap(vector<uint> &listeBlobs, uint type=0);
    vector<double> getPseudoSamplesFullBootstrap(vector<uint> &listeBlobs, uint type =0);
    vector<double> getPseudoSamplesFullBootstrap2(vector<uint> &listeBlobs);

    vector<int> creerHisto(vector<double> &samples, uint histosize, float *mini, float *step);
    vector<int> creerHisto2(vector<double> &samples, double step, float *mini);

    void printHisto(vector<int> &histo, float mini, float step, int type=VERTICAL, FILE *f=NULL);
    void printFile(vector<double> &samples, FILE *f);


//     void getRandomLabelsEnergies(long double nb, FILE *f);
//     void Validation(int type=PERMUT);

    void ValidTestLastChance();
    void ValidTest();
    void ValidAround();
    void followCompConn(set<uint> comp, uint last, uint card, set<uint>&graphe,set<set<uint> > &composantes);
    vector<double> getCaracSample(vector<uint> &comp);
    vector<double> getBackup(vector<uint> &comp);
    uint nbcombinaisons(set<uint> &graphe, uint card);
    double WalshTest(vector<double> &samplesdist,int r);
};






















#endif
