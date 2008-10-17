#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "anneal.h"

using namespace aims;
using namespace carto;
using namespace std;

Anneal::Anneal(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons){
  MinimizationSetup(primal,meshes,lats,lons);
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
      
      uint nclsim1=0,nclsim2=0,nbips1=0,nbips2=0;

      for (uint n=0;n<cliquesDuSite[random[i]].size();n++){
        uint aux = cliquesDuSite[random[i]][n];
        if (cliques[aux].type == SIMILARITY && cliques[aux].blobs[0]->label == cliques[aux].blobs[1]->label){
          if (cliques[aux].blobs[0]->label == sites[random[i]]->label)
            nclsim1++;
          else if (cliques[aux].blobs[0]->label== old)
            nclsim2++;
        }
        long double cen = cliques[aux].updateEnergy(random[i],old,false,nbsujets);
        globalenergieslabels[k] += cen;
      }
      for (uint n=0;n<ipscliques.size();n++){
        if (cliques[ipscliques[n]].type == INTRAPRIMALSKETCH && cliques[ipscliques[n]].blobs[0]->subject != sites[random[i]]->subject){
          nbips1 += cliques[ipscliques[n]].labelscount[labels[k]];
          nbips2 += cliques[ipscliques[n]].labelscount[old];
        }
      }
      globalenergieslabels[k] += 1.0 * (nbips1-nclsim1 - (nbips2-nclsim2));
      cout << nbips1 << ";" << nclsim1 << ";" << nbips2 << ";" << nclsim2 << "=" <<globalenergieslabels[k] << " " ;
      
//       for (uint m=0;m<cliquesDuSite[random[i]].size();m++){

//       }
      somme += exp(-globalenergieslabels[k]/temp);
      if (isnan(somme)) cout << "#####################################" << endl;
    }

    long double somme2=0.0;

//     if (temp<0.9) cout << "dist:[";
    for (uint k=0;k<labels.size();k++){
      somme2 += exp(-globalenergieslabels[k]/temp)/somme;
      globalenergieslabels[k] = somme2;
//       if (temp<0.9){
//         cout << somme2 << " " ;
//       }
    }
//     if (temp<0.9) cout << "]";
      
    long double tirage = ((long double)UniformRandom() * 1.0);
    uint acc;
    for (acc=0;acc<globalenergieslabels.size() && globalenergieslabels[acc]<tirage;acc++){} //cout << globalenergieslabels[acc] << " " ;

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
  //   Initialization();
  for (uint i=0;i<sites.size();i++)
    sites[i]->label = 0; //(int)sites[i]->tValue;
  for (uint k=0;k<cliques.size();k++){
    cliques[k].updateLabelsCount();
    cliques[k].computeEnergy(true,nbsujets);
  }
  energy = getTotalEnergy();

  cout << "energie initiale : " << energy << endl;
  ShortSummaryLabels();


  for (uint i=0;i<cliques.size();i++)
    if (cliques[i].type == INTRAPRIMALSKETCH)
      ipscliques.push_back(i);
  
  cout << ipscliques.size() << " cliques intraps" << endl;
  
  
  vector<int> indices_start;
  for(uint i=0;i<sites.size();i++)
    indices_start.push_back(i);
  
  long double temp=300.0;

  uint mod=1, ite=0, nb_under_threshold=0,test=1;

  cout.precision(2);       

//   cout.setf(ios_base::fixed, ios_base::floatfield);

//     blobsnodes[sites[i]->subject]
  
  FILE * f1;   f1 = fopen ("/home/grg/recuit.txt","w");

  while ((nb_under_threshold<5 || mod!=0) && test!=0){
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
    ASSERT(random.size() == sites.size());
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
    if (temp< 0.9){
      cin >> test;
    }
  }

  for (uint k=0;k<cliques.size();k++){
    cliques[k].updateLabelsCount();
    cliques[k].computeEnergy(true,nbsujets);
  }
  energy = getTotalEnergy();

  cout << "energie initiale : " << energy << endl;
  ShortSummaryLabels();
  
  fclose(f1);
}
