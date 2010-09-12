#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include <cortical_surface/structuralanalysis/icm.h>


void ICM::Step(uint &mod){
    uint min; int old;
    mod=0;
    for (uint i=0; i<sites.size();i++){
        min = 0;
        std::vector<long double> globalenergieslabels(labels.size());
        old = sites[i]->label;
        for (uint k=0;k<labels.size();k++){
            sites[i]->label = labels[k];
            globalenergieslabels[k]=energy;
            for (uint m=0;m<cliquesDuSite[i].size();m++){
                float cen = cliques[cliquesDuSite[i][m]].updateEnergy(i,old,false,nbsujets);
                globalenergieslabels[k] += cen;
        }
        if (globalenergieslabels[k]<globalenergieslabels[min])
            min=k;
        }
        if (old != labels[min]){
            sites[i]->label = labels[min];
            mod++;
            for (uint m=0;m<cliquesDuSite[i].size();m++){
                energy += cliques[cliquesDuSite[i][m]].updateEnergy(i,old,true,nbsujets);
            }
        }
        else {
            sites[i]->label=old;
        }
    }
}

void ICM::Run(){
    float previous_energy;
    int ite = 0;
    uint mod;
    do {
        std::cout << "ite:" << ite++ << " " << std::flush;
        previous_energy = energy;
        Step(mod);
        std::cout << " - E = " << energy << " prev:" << previous_energy << std::endl;
    }
    while(fabs(previous_energy-energy) > 0.01);
}

