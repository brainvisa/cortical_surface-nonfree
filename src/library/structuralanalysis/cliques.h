#ifndef AIMS_CLIQUES_H
#define AIMS_CLIQUES_H

#include <math.h>
#include <cortical_surface/structuralanalysis/sites.h>




enum typesCliques {
    DATADRIVEN, BESTLOWERSCALE, INTRAPRIMALSKETCH, SIMILARITY, UNKNOWN, DATADRIVEN2, GLOBAL
};

class Clique{
    public:
        static float ddweight, intrapsweight, globalweight, simweight, lsweight, ddx1, ddx2, simx1, simx2, ddh;

        short type;
        std::vector< Site * > blobs;
        double energie;
        float similarity, distance;
        std::map<int, uint> labelscount;
        std::map<int, std::set<std::string> > subjectscount;

        float computeEnergy ( bool save, uint CLIQUESNBSUJETS ) {
            float energy = -1.0;
            switch ( type ) {
                case DATADRIVEN:
                    ASSERT( blobs.size() == 1 );
                    if ( blobs[0]->label != 0 ){
                        if ( energy > ddx1 )
                            energy = ddh;
                        else if ( energy < ddx2 )
                            energy = 1.0;
                        else
                            energy = (1.0 - ddh) / (ddx2 - ddx1) * ( blobs[0]->t - ddx2 ) + 1.0;
                    }
                    else {
                        energy = 0.0;
                    }
                    energy *= ddweight;
                    energy *= CLIQUESNBSUJETS;

                break;
                case BESTLOWERSCALE:
                    ASSERT(blobs.size()==1);
                    if (blobs[0]->label != 0)
                    energy = lsweight * blobs[0]->trep;
                    else
                    energy = 0.0;
                    energy *= CLIQUESNBSUJETS;
                    break;
                case INTRAPRIMALSKETCH:
                    energy = 0;
                    for ( uint i = 1 ; i < labelscount.size() ; i++ ) {
                        if ( labelscount[i] <= 1 )
                            energy += 0;
                        else
                            energy += intrapsweight * (labelscount[i]-1);
                    }
                    energy *= CLIQUESNBSUJETS;
                break;
                case SIMILARITY:
                    ASSERT( blobs.size() == 2 );
                    if ( blobs[0]->label == blobs[1]->label && blobs[0]->label != 0 ) {
                        // ATENTION MODE OVERLAP (INVERSE DE DIST)
                        if ( distance > simx2 )
                            energy = -1.0;
                        else if ( distance < simx1 )
                            energy = 0.0;
                        else 
                            energy = (- 1.0) / (simx1 - simx2) * ( distance - simx1 ) - 1.0;

                        energy *= simweight;
                    }
                    else {
                        energy = 0.0;
                    }

                break;
                case GLOBAL:
                    energy = 0;
                    assert( subjectscount.size() > 0 );
                    for ( uint j = 1 ; j < subjectscount.size() ; j++ ) {
                        assert( subjectscount[j].size() <= CLIQUESNBSUJETS );
//                        if ( subjectscount[j].size() == CLIQUESNBSUJETS )
//                            energy += 0;
//                        else
                            energy += - globalweight * subjectscount[j].size();
                    }
                    energy *= CLIQUESNBSUJETS;

                break;
            }
            if (save) energie = energy;
                return energy;
        }

        float updateEnergy( uint node, int old, bool save, uint CLIQUESNBSUJETS ) {

            float energy = 0.0;
            float _intrapsweight, _globalweight;

            switch ( type ) {

                case DATADRIVEN:
                    if ( old == 0 && blobs[0]->label != 0 )
                        energy = computeEnergy( false, CLIQUESNBSUJETS );
                    else if (old != 0 && blobs[0]->label == 0)
                        energy = - energie;
//                    energy = 0.0;
                break;
                case BESTLOWERSCALE:
                    if (old == 0 && blobs[0]->label != 0)
                    energy = computeEnergy(false, CLIQUESNBSUJETS);
                    else if (old != 0 && blobs[0]->label == 0)
                    energy = -energie;
                    break;
                case SIMILARITY:
                    ASSERT((uint)blobs.size()==2);
                        // ASSERT(((uint)blobs[0]->index == (uint)node || (uint)blobs[1]->index == (uint)node));
                        // if ((uint)blobs[0]->index == (uint)node) index = 0;
                        // else if ((uint)blobs[1]->index == (uint)node) index = 1;
                    if ( energie * energie < 0.0001 ) {
                        if ( (uint) blobs[0]->label == (uint) blobs[1]->label && (uint) blobs[0]->label != 0 )
                            energy = computeEnergy( false, CLIQUESNBSUJETS );
                    }
                    else {
                        if ( ( blobs[0]->label != blobs[1]->label ) || ( blobs[1]->label == 0 || blobs[0]->label == 0 ) )
                            energy = -energie;
                    }

                break;
                case INTRAPRIMALSKETCH:
                    _intrapsweight = intrapsweight;
                    uint i;
                    for ( i = 0 ; i < blobs.size() && (uint) blobs[i]->index != (uint) node ; i++ )
                    {}
                    ASSERT( i < blobs.size() );
                    if ( old == blobs[i]->label )
                        energy = 0.0;
                    else {
                        if ( old == 0 )
                            energy = 0.0;
                        else if ( labelscount[ old ] > 1 )
                            energy += -_intrapsweight;

                        if ( blobs[i]->label == 0 )
                            energy += 0.0;
                        else if ( labelscount[ blobs[i]->label ] > 0 )
                            energy += _intrapsweight;
                    }
                    energy *= CLIQUESNBSUJETS;
                    if (save){
                        labelscount[blobs[i]->label]++;
                        labelscount[old]--;
                    }
                break;
                case GLOBAL:
                    _globalweight = globalweight;
                    uint j, k;
                    for ( j = 0 ; j < blobs.size() && (uint) blobs[j]->index != (uint) node ; j++ )
                    {}
                    ASSERT( j < blobs.size() );
                    k = 0;
                    while ( k < blobs.size() && !(blobs[k]->subject == blobs[j]->subject && blobs[k]->label == old ) ){
                        k++;
                    }
                    if ( old == blobs[j]->label )
                        energy = 0.0;
                    else {
                        if ( old == 0 && blobs[j]->label == 0 )
                            energy = 0.0;
                        else if ( old != 0 ) { // && blobs[j]->label == 0 ) {
//                            if ( subjectscount[ blobs[j] ].size() < CLIQUESNBSUJETS )

                            if ( k < blobs.size() )
                                energy += 0.0;
                            else 
                                energy += _globalweight;
                        }

                        if ( blobs[j]->label != 0 ) {
                            if (subjectscount[blobs[j]->label].find(blobs[j]->subject) == subjectscount[blobs[j]->label].end()) {
                                energy += -_globalweight;
                            }
                        }
//                        else if ( old != 0 && blobs[j]->label != 0 ) {
//                            if (subjectscount[old].find(blobs[j]->subject) == subjectscount[blobs[j]->label].end()) {
//                                energy += -_globalweight;
//                            }
//                            if (subjectscount[blobs[j]->label].find(blobs[j]->subject) == subjectscount[blobs[j]->label].end()) {
//                                energy += -_globalweight;
//                            }
//                        }
                    }
                    energy *= CLIQUESNBSUJETS;
            //           energy = 0.0;
                    if (save){
                        if ( k == blobs.size() )
                            subjectscount[old].erase(blobs[j]->subject);
                        
                        subjectscount[blobs[j]->label].insert(blobs[j]->subject);
                        
                    }
                break;
            }
            if (save) energie += energy;
            return energy;
        }

        void updateLabelsCount();

        void updateSubjectsCount();

        static void setParameters ( float _ddweight, 
                                    float _intrapsweight, 
                                    float _simweight, 
                                    float _lsweight, 
                                    float _ddx2, 
                                    float _ddx1, 
                                    float _simx1,
                                    float _simx2,
                                    float _ddh, 
                                    float _globalweight );

        static float getIntraPSWeight() { return intrapsweight; }
        Clique() { 
            type = UNKNOWN; 
            energie = 0.0; 
            blobs = std::vector<Site *>(); 
            labelscount = std::map<int,uint>();  
            subjectscount = std::map<int, std::set<std::string> >();
        }
        

};

double getOverlap(Point3df bbmin1, Point3df bbmax1, Point3df bbmin2, Point3df bbmax2, uint *no_overlap);

void BuildMaximalOrderCliques ( std::vector<Site *> &sites, std::vector<std::vector<int> > &cliquesDuSite, std::vector<Clique> &cliques);
void BuildDataDrivenCliques ( std::vector<Site *> &sites, std::vector<std::vector<int> > &cliquesDuSite, std::vector<Clique> &cliques);
void BuildGlobalClique ( std::vector<Site *> &sites, std::vector<std::vector<int> > &cliquesDuSite, std::vector<Clique> &cliques );
void BuildLowerScaleCliques ( std::vector<Site *> &sites, std::vector<std::vector<int> > &cliquesDuSite, std::vector<Clique> &cliques );






#endif

