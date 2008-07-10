#ifndef AIMS_CLIQUES_H
#define AIMS_CLIQUES_H

#include <math.h>
#include "sites.h"
#include "meshdistance.h"

using namespace std;

enum typesCliques {
  DATADRIVEN, BESTLOWERSCALE, INTRAPRIMALSKETCH, SIMILARITY, UNKNOWN
};

class Clique{
  private:
  static float ddweight, intrapsweight, simweight, lsweight, ddx2, ddx1, ddh;

  public :
    int type;
    vector<Site *> blobs;
    double energie,sigma,rec;
    map<int,uint> labelscount;
    float computeEnergy(bool save, uint CLIQUESNBSUJETS) {
      float energy;
      switch (type){
        case DATADRIVEN:
          ASSERT(blobs.size()==1);
          if (blobs[0]->label != 0){
            if (blobs[0]->tValue > ddx2) energy = ddh;
            else if (blobs[0]->tValue < ddx1) energy = ddweight;
            else // energy = ddweight - (sites[blobs[0]].t- ddx1) * (ddweight - ddh)/(double(ddx2 - ddx1));
            { sigma = 4.0; //((3.0*ddx1+5.0*ddx2)/7.0)/sqrt(2*log(1/((5.0*ddh+1.0*ddweight)/6.0)));          // paramètre de la gaussienne
            energy = (ddweight-ddh)/(ddx1-ddx2) * blobs[0]->tValue + (ddh - (ddweight-ddh)/(ddx1-ddx2)*ddx2);
            }
          }
          else {
            energy = 0.0;
          }
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
          energy=0;
          for (uint i=1;i<labelscount.size();i++){
            if (labelscount[i]<=1)
              energy += 0;
            else
              energy += intrapsweight * (labelscount[i]-1);
          }
          energy *= CLIQUESNBSUJETS;
          break;
        case SIMILARITY:
          ASSERT(blobs.size()==2);
          if (blobs[0]->label == blobs[1]->label && blobs[0]->label != 0){
            sigma = 10.0/sqrt(2*log(10.0));          // paramètre de la gaussienne : le premier 10.0 c'est la distance-seuil à laquelle on veut un potentiel égal à 0.1
//             energy = -simweight*exp(-pow(rec,2)/(2*pow(sigma,2)));
            energy=simweight/20.0-simweight;
          }
          else {
            energy = 0.0;
          }
          break;
      }
      if (save) energie = energy;
      return energy;
    }
    float updateEnergy(uint node, int old, bool save, uint CLIQUESNBSUJETS)  {
      float energy=0.0;
      uint index;
      float _intrapsweight;
      switch(type){
        case DATADRIVEN:
          if (old == 0 && blobs[0]->label != 0)
            energy = computeEnergy(false, CLIQUESNBSUJETS);
          else if (old != 0 && blobs[0]->label == 0)
            energy = -energie;
          break;
        case BESTLOWERSCALE:
          if (old == 0 && blobs[0]->label != 0)
            energy = computeEnergy(false, CLIQUESNBSUJETS);
          else if (old != 0 && blobs[0]->label == 0)
            energy = -energie;
          break;
        case SIMILARITY:
          ASSERT((uint)blobs.size()==2 && ((uint)blobs[0]->index == (uint)node || (uint)blobs[1]->index == (uint)node));
          if ((uint)blobs[0]->index == (uint)node) index = 0;
          else if ((uint)blobs[1]->index == (uint)node) index = 1;
          if (energie == 0 && (uint)blobs[0]->label == (uint)blobs[1]->label && (uint)blobs[0]->label != 0)
            energy = computeEnergy(false, CLIQUESNBSUJETS);
          else if (energie != 0 && ((blobs[0]->label == blobs[1]->label && blobs[0]->label == 0) || (blobs[0]->label != blobs[1]->label)))
            energy = -energie;
          break;
        case INTRAPRIMALSKETCH:
          _intrapsweight = intrapsweight;
          uint i;
          for (i=0;i<blobs.size() && (uint)blobs[i]->index != (uint)node;i++)
          {}
          if (old == blobs[i]->label)
            energy = 0.0;
          else {
            if (old == 0) energy = 0.0;
            else if (labelscount[old] > 1) energy += -_intrapsweight;
            if (blobs[i]->label == 0) energy += 0.0;
            else if (labelscount[blobs[i]->label] > 0) energy += _intrapsweight;
          }
          energy *= CLIQUESNBSUJETS;
          if (save){
            labelscount[blobs[i]->label]++;
            labelscount[old]--;
          }
          break;
      }
      if (save) energie += energy;
     
      return energy;
    }
    void updateLabelsCount();
    static void setParameters(float _ddweight, float _intrapsweight, float _simweight, float _lsweight, float _ddx2, float _ddx1, float _ddh);
    static float getIntraPSWeight(){return intrapsweight;}
    Clique(){ type = UNKNOWN; energie = 0.0; blobs = vector<Site *>(); labelscount = map<int,uint>();  }
};

vector<Clique> ConstruireCliques(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, AimsSurfaceTriangle &mesh);

#endif

