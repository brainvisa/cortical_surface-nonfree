#ifndef AIMS_VALIDATION_H
#define AIMS_VALIDATION_H
#include <cortical_surface/structuralanalysis/minimization.h>



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
    std::vector<int> getCompConn(std::vector<uint> &indicesCliques, std::set<uint> &listeSites);
    std::vector<std::set<uint> > getCompConnVector(std::vector<int> &comp);


//     vector<double> getPseudoSamplesBootstrap(vector<uint> &listeBlobs, uint type=0);
//     vector<double> getPseudoSamplesFullBootstrap(vector<uint> &listeBlobs, uint type =0);
//     vector<double> getPseudoSamplesFullBootstrap2(vector<uint> &listeBlobs);

    std::vector<int> creerHisto(std::vector<double> &samples, uint histosize, float *mini, float *step);
    std::vector<int> creerHisto2(std::vector<double> &samples, double step, float *mini);

    void printHisto(std::vector<int> &histo, float mini, float step, int type=VERTICAL, FILE *f=NULL);
    void printFile(std::vector<double> &samples, FILE *f);




    void ValidAround();
    void followCompConn(std::set<uint> comp, uint last, uint card, std::set<uint>&graphe, std::set<std::set<uint> > &composantes);
    std::vector<double> getCaracSample(std::vector<uint> &comp);
    std::vector<double> getBackup(std::vector<uint> &comp);
    uint nbcombinaisons(std::set<uint> &graphe, uint card);
    double WalshTest(std::vector<double> &samplesdist,int r);
};






















#endif
