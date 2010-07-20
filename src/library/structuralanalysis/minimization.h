#ifndef AIMS_MINIMIZATION_H
#define AIMS_MINIMIZATION_H
#include "cliques.h"


enum typesMinim{
  ICM, ANNEAL, CUSTOM
};



class SurfaceBased_StructuralAnalysis{
  public:
    float Esimil;
    std::vector<uint> ipscliques;
    
    long double energy;
    uint nbsujets;
    std::vector<int> labels;
    std::vector<std::pair<Point2df,Point2df> > labelsZones;
    std::vector<std::set<uint> > zonesListesBlobs;
    std::vector<std::set<uint> > listeZones; // attention les indices de listeZones sont décalés de 1 par rapport à labelsZones (à cause du label 0 qui recouvre tout l'espace 2D)
    std::vector<Site *> sites;
    std::vector<Clique> cliques;
    std::vector<std::vector<int> > cliquesDuSite;

    void ShortSummaryLabels();

    std::string energypath, recuitpath;
    uint run,save;
    void noLabelsZones();
    void regionLabelsZones(); //vector<pair<Point2df,Point2df> > &labelsZones, vector<set<uint> > &zonesListesBlobs, vector<set<uint> > &listeZones );
    void MinimizationSetup(Graph &primal);
    

    SurfaceBased_StructuralAnalysis(){}
    SurfaceBased_StructuralAnalysis(Graph &primal);

    void setModelParameters(float _ddweight=0.0, float _intrapsweight = 0.0, float _simweight=0.0, float _lsweight=0.0, float _ddx2 = 0.0, float _ddx1 = 0.0, float _ddh=0.0);
    
    void RunMinimization(int type);
    long double getLabelEnergy(int label, int type=UNKNOWN);
    long double getClusterEnergy(std::vector<uint> &composante);
    long double getTypeEnergy(int type);
    long double getTotalEnergy();
//     double getTotalEnergyLastChance(uint site, uint newlabel);

    void SummaryLabels();
    void StoreToGraph(Graph &primal);
    void Initialization();


};




#endif

