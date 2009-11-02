#ifndef AIMS_CLIQUES_H
#define AIMS_CLIQUES_H

#include <math.h>
#include <cortical_surface/structuralanalysis/sites.h>
#include <cortical_surface/structuralanalysis/meshdistance.h>
// #include "cathier/aims_wrap.h"
// #include "cathier/triangle_mesh_geodesic_map.h"
// #include "cathier/math_functions.h"

using namespace std;

enum typesCliques {
  DATADRIVEN, BESTLOWERSCALE, INTRAPRIMALSKETCH, SIMILARITY, UNKNOWN, DATADRIVEN2
};

class Clique{
  public:
  static float ddweight, intrapsweight, simweight, lsweight, ddx2, ddx1, ddh;

//   public :
    int type;
    vector<Site *> blobs;
    double energie,sigma,rec;
    map<int,uint> labelscount;
    
    float computeEnergy(bool save, uint CLIQUESNBSUJETS) {
//       cout << ddweight << " " << simweight << " " << ddh << " " << ddx1 << " " << ddx2 << ";" ;
      float energy=-1.0;
      switch (type){
        case DATADRIVEN:
          ASSERT(blobs.size()==1);
//           cout << "DD" << blobs[0]->label << " " ;
          if (blobs[0]->label != 0){
//             cout << blobs[0]->label << " " << blobs[0]->t ;
//            if (blobs[0]->t > ddx2) energy = ddh;
            if (blobs[0]->t < ddx1) energy = 1.0; 
            else { 
            energy =  pow(ddx2/2.0,2)/(pow(0.5*ddx2,2)+pow(blobs[0]->t -ddx1,2));
            }
//             cout << " " << energy << " " << ddweight << " " << CLIQUESNBSUJETS << "/";

          }
          else {
            energy = 0.0;
          }
          energy *= ddweight;
          energy *= CLIQUESNBSUJETS;

        break;
//         case BESTLOWERSCALE:
//           ASSERT(blobs.size()==1);
//           if (blobs[0]->label != 0)
//             energy = lsweight * blobs[0]->trep;
//           else
//             energy = 0.0;
//           energy *= CLIQUESNBSUJETS;
//           break;
        case INTRAPRIMALSKETCH:
          energy=0;
          for (uint i=1;i<labelscount.size();i++){
            if (labelscount[i]<=1)
              energy += 0;
            else
              energy += intrapsweight * (labelscount[i]-1);
          }
          energy *= CLIQUESNBSUJETS;
//           energy = 0.0;
//           ASSERT (energy*energy < 0.00001);
          break;
        case SIMILARITY:
          ASSERT(blobs.size()==2);
//           cout << "SIM" << blobs[0]->label << " " << blobs[1]->label << " " << simweight << " ";
          if (blobs[0]->label == blobs[1]->label && blobs[0]->label != 0){
//             sigma = 19.0/sqrt(2*log(10.0));          // paramètre de la gaussienne : le premier 10.0 c'est la distance-seuil à laquelle on veut un potentiel égal à 0.1
//             energy = -simweight*exp(-pow(rec,2)/(2*pow(sigma,2)));
//             cout << endl << rec << " : " << energy << endl;
//             energy=simweight/20.0-simweight;            
//             if (rec > 20.0)
//               energy = rec/5.0 - 4.0;
//             else if (rec < 20.0)
              energy = -rec;
//             cout << "((" << rec << "=>" << energy << ")) " ;
            energy *= simweight;
//             cout << energy << "/";
//             energy = 0.0;
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
//       uint index;
      float _intrapsweight;
      switch(type){
        case DATADRIVEN:
          if (old == 0 && blobs[0]->label != 0)
            energy = computeEnergy(false, CLIQUESNBSUJETS);
          else if (old != 0 && blobs[0]->label == 0)
            energy = -energie;
          break;
//         case BESTLOWERSCALE:
//           if (old == 0 && blobs[0]->label != 0)
//             energy = computeEnergy(false, CLIQUESNBSUJETS);
//           else if (old != 0 && blobs[0]->label == 0)
//             energy = -energie;
//           break;
        case SIMILARITY:
          ASSERT((uint)blobs.size()==2);
//           ASSERT(((uint)blobs[0]->index == (uint)node || (uint)blobs[1]->index == (uint)node));
//           if ((uint)blobs[0]->index == (uint)node) index = 0;
//           else if ((uint)blobs[1]->index == (uint)node) index = 1;
          if (energie*energie < 0.000000000001){
            if ((uint)blobs[0]->label == (uint)blobs[1]->label && (uint)blobs[0]->label != 0)
              energy = computeEnergy(false, CLIQUESNBSUJETS);
          }
          else if (energie*energie > 0.000000000001){
            if ((blobs[0]->label != blobs[1]->label) || (blobs[1]->label == 0 || blobs[0]->label == 0))
              energy = -energie;
          } 
          else printf("eerr %lf\n", (double)energie);
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
//           energy = 0.0;
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

double getOverlap(Point3df bbmin1, Point3df bbmax1, Point3df bbmin2, Point3df bbmax2, uint *no_overlap);

vector<Clique> ConstruireCliques(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, map<string, AimsSurfaceTriangle> &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons);

// vector<Clique> ConstruireCliquesSimple(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite, map<string, AimsSurfaceTriangle> &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons);

vector<Clique> ConstruireCliquesLastChance(vector<Site *> &sites, vector<vector<int> > &cliquesDuSite);

#endif

