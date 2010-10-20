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
    
  public:
    std::vector<int> getCompConn ( SurfaceBased_StructuralAnalysis &ssb, 
                                   std::vector<uint> &indicesCliques, 
                                   std::set<uint> &listeSites ) ;
    std::vector<std::set<uint> > getCompConnVector ( std::vector<int> &comp );

    void printFile ( std::vector<double> &samples, FILE *f);

    void ValidAround ( SurfaceBased_StructuralAnalysis &ssb );
    
    std::vector<double> getCaracSample ( SurfaceBased_StructuralAnalysis &ssb, 
                                         std::vector<uint> &comp );
    std::vector<double> getBackup ( SurfaceBased_StructuralAnalysis &ssb,
                                    std::vector<uint> &comp );
    uint nbcombinaisons ( SurfaceBased_StructuralAnalysis &ssb, 
                          std::set<uint> &graphe, 
                          uint card );
    
    double WalshTest ( std::vector<double> &samplesdist, int r );
    
    StructuralAnalysis_Validation() {}
};






















#endif
