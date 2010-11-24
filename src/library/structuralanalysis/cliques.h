#ifndef AIMS_CLIQUES_H
#define AIMS_CLIQUES_H

#include <math.h>
#include <float.h>
#include <cortical_surface/structuralanalysis/sites.h>

enum typesCliques {
    DATADRIVEN, BESTLOWERSCALE, INTRAPRIMALSKETCH, SIMILARITY, UNKNOWN, GLOBAL
};

class Clique{
    public:
        static float ddweight, intrapsweight, simweight, lsweight, ddx1, ddx2, simx1, simx2, lsx1, lsx2, ddh;

        double energie;
        short type;
        std::vector< Site * > blobs;
        double similarity, distance;
        std::map<int, uint> labelscount;
        std::map<int, std::set<std::string> > subjectscount;
        double pot;

        double computeEnergy ( bool save, uint CLIQUESNBSUJETS ) {
            double energy = 0.0;
            switch ( type ) {
                case DATADRIVEN:
                    ASSERT( blobs.size() == 1 );
                    if ( blobs[0]->label != 0 ) {
                        float measure = blobs[0]->t;
                        if ( measure > ddx1 )
                            energy = ddh;
                        else if ( measure < ddx2 )
                            energy = 1.0;
                        else
                            energy = (1.0 - ddh) / (ddx2 - ddx1) * ( measure - ddx2 ) + 1.0;
                    }
                    else 
                        energy = 0.0;
                    energy *= ddweight;
                    energy *= CLIQUESNBSUJETS;

                break;
                case BESTLOWERSCALE:
                    ASSERT( blobs.size() == 1 );
                    if ( blobs[0]->label != 0 ) {
                        double mean_scale = (blobs[0]->tmax + blobs[0]->tmin) / 2.0;
                        if ( mean_scale < lsx1 )
                            energy = 0.0;
                        else if ( mean_scale > lsx2 )
                            energy = 1.0;
                        else
                            energy = 1.0 / (lsx2 - lsx1) * ( mean_scale - lsx2 ) + 1.0;
                    }
                    else
                        energy = 0.0;
                    energy *= lsweight;
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
                          
                        assert ( distance < 0.0 || similarity < 0.0 );
//                        std::cout << distance << "/" << similarity << " " << std::flush;
                        
                        // MODE OVERLAP
                        if ( similarity > 0.0 ) {
                            if ( similarity > simx2 )
                                energy = -1.0;
                            else if ( similarity < simx1 )
                                energy = 0.0;
                            else
                                energy = (- 1.0) / (simx2 - simx1) * ( similarity - simx2 ) - 1.0;
                        }
                        // MODE DISTANCE
                        else if ( distance > 0.0 ) {
                            if ( distance < simx2 )
                                energy = -1.0;
                            else if ( distance > simx1 )
                                energy = 0.0;
                            else
                                energy = (- 1.0) / (simx2 - simx1) * ( distance - simx2 ) - 1.0;
                        }
                    }
                    else {
                        energy = 0.0;
                    }
                    energy *= simweight;
                break;
                case GLOBAL:
                    energy = 0.0;
                break;
            }
            if (save) energie = energy;
            return energy;
        }

        double updateEnergy ( uint node, int old, bool save, uint CLIQUESNBSUJETS ) {

            double energy = 0.0;
            float _intrapsweight;

            switch ( type ) {

                case DATADRIVEN:
                    if ( old == 0 && blobs[0]->label != 0 )
                        energy = computeEnergy ( false, CLIQUESNBSUJETS );
                    else if ( old != 0 && blobs[0]->label == 0 )
                        energy = - energie;
                break;
                case BESTLOWERSCALE:
                    if ( old == 0 && blobs[0]->label != 0 )
                        energy = computeEnergy ( false, CLIQUESNBSUJETS );
                    else if ( old != 0 && blobs[0]->label == 0 )
                        energy = -energie;
                    break;
                case SIMILARITY:
                    ASSERT( (uint)blobs.size()==2 );
                        // ASSERT(((uint)blobs[0]->index == (uint)node || (uint)blobs[1]->index == (uint)node));
                        // if ((uint)blobs[0]->index == (uint)node) index = 0;
                        // else if ((uint)blobs[1]->index == (uint)node) index = 1;
                    if ( fabs(energie) < DBL_EPSILON ) {
                        if ( (uint) blobs[0]->label == (uint) blobs[1]->label && (uint) blobs[0]->label != 0 )
                            energy = computeEnergy ( false, CLIQUESNBSUJETS );
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
                    if ( save ) {
                        labelscount[blobs[i]->label]++;
                        labelscount[old]--;
                    }
                break;
                case GLOBAL:
                    uint j, k;
                    for ( j = 0 ; j < blobs.size() && (uint) blobs[j]->index != (uint) node ; j++ )
                    {}
                    ASSERT( j < blobs.size() );
                    k = 0;
                    while ( k < blobs.size() && !(blobs[k]->subject == blobs[j]->subject && blobs[k]->label == old ) ) {
                        k++;
                    }
                    if ( save ) {
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
                                    float _lsx1,
                                    float _lsx2,
                                    float _ddh );

        static float getIntraPSWeight() { return intrapsweight; }
        
        Clique() {
            type = UNKNOWN;
            energie = 0.0;
            pot = 100.0;
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

