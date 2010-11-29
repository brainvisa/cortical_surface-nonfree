#ifndef AIMS_MINIMIZATION_H
#define AIMS_MINIMIZATION_H
#include <cortical_surface/structuralanalysis/sites.h>
#include <cortical_surface/structuralanalysis/cliques.h>
#include <cortical_surface/structuralanalysis/blobs.h>


enum typesMinim {
    ICM, ANNEAL, CUSTOM
};


enum typeDistance { 
    DISTANCE_3DEUCLIDIAN, DISTANCE_LATITUDES, DISTANCE_LONGITUDES, DISTANCE_LATLON, DISTANCE_BOUNDING_BOXES_LATLON
};

class SurfaceBased_StructuralAnalysis {
    public:

        std::vector<uint> ipscliques;
        int globalclique;

        double energy;
        uint nbsujets;
        std::vector<int> labels;
        std::vector<std::pair<Point2df,Point2df> > labelsZones;
        std::vector<std::set<uint> > zonesListesBlobs;
        std::vector<std::set<uint> > listeZones; // attention les indices de listeZones sont décalés de 1 par rapport à labelsZones (à cause du label 0 qui recouvre tout l'espace 2D)
        std::vector<Site *> sites;
        std::vector<Clique> cliques;
        std::vector<std::vector<int> > cliquesDuSite;

        void ShortSummaryLabels();

        std::string energyPath, labelsPath;
        uint save;
        void noLabelsZones( int number_of_labels = 21 );
        void regionLabelsZones();

        void PrintCliquesNumbers(){
            uint nb_cl_sim = 0, nb_cl_dd = 0, nb_cl_intraps = 0, nb_cl_lower = 0;

            for ( uint i = 0 ; i < cliques.size() ; i++ ) {
                if ( cliques[i].type == SIMILARITY ) nb_cl_sim++;
                else if ( cliques[i].type == DATADRIVEN ) nb_cl_dd++;
                else if ( cliques[i].type == BESTLOWERSCALE ) nb_cl_lower++;
                else if ( cliques[i].type == INTRAPRIMALSKETCH ) nb_cl_intraps++;
            }
            std::cout << " done (" << nb_cl_sim << " similarity cliques ; " <<
                nb_cl_dd << " datadriven cliques ; " <<
                nb_cl_lower << " lower scale cliques ; " <<
                nb_cl_intraps << " maximal order cliques ; " <<
                cliques.size() << " cliques in total)" << std::endl;
            std::cout << nbsujets << " subjects" << std::endl;
        }

        SurfaceBased_StructuralAnalysis(){}

        void setModelParameters ( float _ddweight = 0.0,
                                  float _intrapsweight = 0.0,
                                  float _simweight = 0.0,
                                  float _lsweight = 0.0,
                                  float _ddx1 = 0.0,
                                  float _ddx2 = 0.0,
                                  float _simx1 = 0.0,
                                  float _simx2 = 0.0,
                                  float _lsx1 = 0.0,
                                  float _lsx2 = 0.0,
                                  float _ddh = 0.0 );

        double getLabelEnergy ( int label, int type = UNKNOWN );
        double getClusterEnergy ( std::vector<uint> &composante );
        double getTypeEnergy ( int type );
        double getTotalEnergy ();

        void SummaryLabels();
        void StoreToGraph(Graph &primal);
        void Initialization( bool initLabels = true );


        static void ConvertSSBlobsToSites ( std::vector<surf::ScaleSpaceBlob *> &ssblobs, std::vector<Site *> &sites );

        void importGraphNodesFromBlobs ( std::vector<surf::ScaleSpaceBlob *> &ssblobs ) {
            ConvertSSBlobsToSites ( ssblobs, this->sites );
        }

        static void GetSimilarityCliquesFromSSBCliques ( std::vector<surf::Clique> &ssbcliques,
                                        std::vector<Site *> &sites,
                                        std::vector<Clique> &cliques,
                                        std::vector<std::vector<int> > &cliquesDuSite ) ;

        static std::vector<surf::Clique> BuildSimilarityCliques ( std::vector<surf::ScaleSpaceBlob *>   &ssblobs,
                                                         std::vector<std::vector<surf::GreyLevelBlob *> > &matchingblobs ) ;


        static std::vector<surf::Clique> BuildSimilarityCliques3D ( std::vector<surf::ScaleSpaceBlob *>   &ssblobs,
                                                           GroupData &data,
                                                           float threshold = 5.0,
                                                           int typeDistance = DISTANCE_3DEUCLIDIAN );

        void importGraphCliquesFromSSBCliques ( std::vector<surf::Clique> &ssbcliques ) {
            GetSimilarityCliquesFromSSBCliques( ssbcliques, this->sites, this->cliques, this->cliquesDuSite );
        }

        void importGraphNodesAndCliques ( std::vector<surf::ScaleSpaceBlob *> &ssblobs, std::vector<surf::Clique> &ssbcliques ) {
            ConvertSSBlobsToSites ( ssblobs, this->sites );

            GetSimilarityCliquesFromSSBCliques( ssbcliques, this->sites, this->cliques, this->cliquesDuSite );

            BuildGlobalClique( this->sites, this->cliquesDuSite, this->cliques );
            BuildDataDrivenCliques( this->sites, this->cliquesDuSite, this->cliques );
            BuildMaximalOrderCliques( this->sites, this->cliquesDuSite, this->cliques );
            BuildLowerScaleCliques( this->sites, this->cliquesDuSite, this->cliques );

            for ( uint i = 0 ; i < this->cliques.size() ; i++ ) {
                if ( this->cliques[i].type == INTRAPRIMALSKETCH )
                    this->ipscliques.push_back(i);
                if ( this->cliques[i].type == GLOBAL )
                    this->globalclique = i;
            }

            std::set<std::string> subjects;
            for ( uint i = 0 ; i < this->sites.size() ; i++ )
                subjects.insert(this->sites[i]->subject);
            nbsujets = subjects.size();
            assert( nbsujets == ipscliques.size() );
        }


};




#endif

