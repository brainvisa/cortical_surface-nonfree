#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "anneal.h"

using namespace aims;
using namespace carto;
using namespace std;

Anneal::Anneal(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon){
  MinimizationSetup(primal,mesh,lat,lon);
}

void Anneal::Step(vector<int> &random, long double temp, uint &mod){
  
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

void Anneal::Run(){
  energy = getTotalEnergy();
  cout << "energie initiale : " << energy << endl;
  
  vector<int> indices_start;
  for(uint i=0;i<sites.size();i++)
    indices_start.push_back(i);
  
  long double temp=1000.0;

  uint mod=1, ite=0, nb_under_threshold=0;
  cout << "\a" << flush;   cout.precision(2);   cout.setf(ios_base::fixed, ios_base::floatfield);
  FILE * f1;   f1 = fopen ("/home/grg/recuit.txt","w");
  
  while (nb_under_threshold<5 || mod!=0){
    if (mod!=0) nb_under_threshold=0;
    else nb_under_threshold++;
    cout << " T=" << temp << " it="<< ite++ << " " ;
    for (uint i0=0;i0<sites.size();i0++){
      fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
    }
    fprintf(f1, "\n");
    vector<int> indices(indices_start);
    vector<int> random;
    for (uint i=0;i<sites.size();i++){
      int index = (int)(UniformRandom() * indices.size());
      random.push_back(indices[index]);
      indices.erase(indices.begin()+index);
    }
    
    Step(random, temp, mod);

    float Eintra, Edd, Esim, Els, energy2;
    Eintra = getTypeEnergy(INTRAPRIMALSKETCH);
    Edd = getTypeEnergy(DATADRIVEN);
    Esim = getTypeEnergy(SIMILARITY);
    Els = getTypeEnergy(BESTLOWERSCALE);
    energy2 = Eintra + Edd + Esim + Els;
    cout << "E:" << energy2 << " (" << Eintra << ";" << Edd << ";" << Esim << ";" << Els << ") chg:" << mod << " " << flush;
    ShortSummaryLabels();
    cout << " E=" << energy << endl;
    temp = temp*0.98;
  }

  fclose(f1);
}
