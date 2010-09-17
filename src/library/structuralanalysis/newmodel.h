#ifndef AIMS_NEWMODEL_H
#define AIMS_NEWMODEL_H
#include <cortical_surface/structuralanalysis/minimization.h>
#include <cortical_surface/structuralanalysis/cliques.h>


class NewModel: public SurfaceBased_StructuralAnalysis{
    private:

        public:

            void Run( int verbose = 1 );
            void Step( std::vector<int> &random, long double temp, uint &mod );
};




#endif

