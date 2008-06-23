#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "cluster.h"

using namespace aims;
using namespace carto;
using namespace std;

SWC::SWC(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon){
  MinimizationSetup(primal,mesh,lat,lon);
}

map<int ,int> SWC::getCompConn(vector<uint> &indicesCliques){
  map<int,int> comp;
  uint lcomp=1;
  uint nbcomp=0;
  for (uint i=0;i<indicesCliques.size();i++)
    for (uint j=0;j<cliques[indicesCliques[i]].blobs.size();j++)
      comp[cliques[indicesCliques[i]].blobs[j]->index]=0;
  
  for (uint i=0;i<indicesCliques.size();i++){
    if (cliques[indicesCliques[i]].type == SIMILARITY){
      if (comp[cliques[indicesCliques[i]].blobs[0]->index] == 0){
        if (comp[cliques[indicesCliques[i]].blobs[1]->index] == 0){
          comp[cliques[indicesCliques[i]].blobs[0]->index] = lcomp;
          comp[cliques[indicesCliques[i]].blobs[1]->index] = lcomp;
          lcomp++;
          nbcomp++;
        }
        else {
          comp[cliques[indicesCliques[i]].blobs[0]->index] = comp[cliques[indicesCliques[i]].blobs[1]->index];
        }
      }
      else if (comp[cliques[indicesCliques[i]].blobs[1]->index] == 0){
        comp[cliques[indicesCliques[i]].blobs[1]->index] = comp[cliques[indicesCliques[i]].blobs[0]->index];
      }
      else if (comp[cliques[indicesCliques[i]].blobs[0]->index] == comp[cliques[indicesCliques[i]].blobs[1]->index]){
        // clique cyclique ne pas s'inquiéter mais répéter 100 fois rapidement
      }
      else { // fusion
        nbcomp--;
        for (uint j=0;j<comp.size();j++)
          if (comp[j] == comp[cliques[indicesCliques[i]].blobs[1]->index])
            comp[j] = comp[cliques[indicesCliques[i]].blobs[0]->index];
      }
    }
  }
  return comp;
}

vector<uint> SWC::getCliquesTurnedOn(float temp){
  vector<uint> turnedOn;
  for (uint i=0;i<cliques.size();i++){
    if (cliques[i].type == SIMILARITY){
      if (cliques[i].blobs[0]->label == cliques[i].blobs[1]->label){
        long double tirage = ((long double)UniformRandom() * 1.0);
        long double qe= exp(-(cliques[i].rec/(temp *cliques[i].blobs[0]->t*cliques[i].blobs[1]->t)));
        if (cliques[i].blobs[0]->subject == cliques[i].blobs[1]->subject) qe /= 2.0;
        if (tirage<qe) turnedOn.push_back(i);
      }
    }
  }
  return turnedOn;
}


void SWC::Run(){

  FILE * f1;   f1 = fopen ("/home/grg/recuit.txt","w");
  float  temp= 10.0;
  uint ite=0;
  while (temp>0.0000001){
    cout << " T=" << temp << " it="<< ite++ << " " << flush ;

    vector<uint> turnedOn(getCliquesTurnedOn(temp));
    map<int,int> comp(getCompConn(turnedOn));
    map<int,int>::iterator it;
    for (it=comp.begin();it!=comp.end();it++){
      sites[it->first]->label = it->second;
    }
    for (uint i0=0;i0<sites.size();i0++){
      fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
    }
    fprintf(f1, "\n");

    temp *= 0.99;
    cout << endl;
  }
  fclose(f1);
}

void SWC::Step(vector<int> &random, long double temp, uint &mod){
  
  long double somme=0.0;
  int old;
  mod=0;
  for (uint i=0; i<random.size();i++){
    somme=0.0;
    vector<long double> globalenergieslabels(labels.size());
    old = sites[random[i]]->label;
    for (uint k=0;k<labels.size();k++){
      sites[random[i]]->label = labels[k];
      globalenergieslabels[k]=energy;
      for (uint m=0;m<cliquesDuSite[random[i]].size();m++){
        long double cen = cliques[cliquesDuSite[random[i]][m]].updateEnergy(random[i],old,false,nbsujets);
        globalenergieslabels[k] += cen;
      }
      somme += exp(-globalenergieslabels[k]/temp);
      if (isnan(somme)) cout << "#####################################" << endl;
    }

    long double somme2=0.0;

    for (uint k=0;k<labels.size();k++){
      somme2 += exp(-globalenergieslabels[k]/temp)/somme;
      globalenergieslabels[k] = somme2;

    }
      
    long double tirage = ((long double)UniformRandom() * 1.0);
    uint acc;
    for (acc=0;acc<globalenergieslabels.size() && globalenergieslabels[acc]<tirage;acc++) ; //cout << globalenergieslabels[acc] << " " ;

    if (old != labels[acc]) {
      sites[random[i]]->label=labels[acc];
      for (uint m=0;m<cliquesDuSite[random[i]].size();m++){
        energy += cliques[cliquesDuSite[random[i]][m]].updateEnergy(random[i],old,true,nbsujets);
      }
      mod++;
    }
    else {
      sites[random[i]]->label = old;
    }
  }
}


