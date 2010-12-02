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
    double getRankPercentile ( std::vector< double > &samplesT,
                               std::set<uint> &activblobs, 
                               SurfaceBased_StructuralAnalysis &ssb ) ;
    std::vector<uint> getComponent ( SurfaceBased_StructuralAnalysis &ssb,
                                     std::set<uint> &activblobs, 
                                     std::set< uint > &forbidden,
                                     std::vector< int > &tirage ) ;
    std::vector<std::vector<uint> > ValidAround ( uint label,
                                                  SurfaceBased_StructuralAnalysis &ssb,
                                                  uint number_of_samples = 10000 );
    void saveSignificanceInfo ( uint label,
                                SurfaceBased_StructuralAnalysis &swc, 
                                std::vector<std::vector<uint> > & composantes, 
                                FILE *f );
    
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
