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
    vector<long double> globalenergieslabels(labels.size()), expenergies(labels.size()),total(labels.size());
    old = sites[random[i]]->label;
    
//     cout << "e:" << energy;

    for (uint k=0;k<labels.size();k++){
//       cout << " old:" << old << " l:" <<labels[k];
      sites[random[i]]->label = labels[k];
      globalenergieslabels[k]=energy;
//       cout << "GLOBAL AVANT: " << globalenergieslabels[k] << " ";
      
      int nclsim1=0,nclsim2=0,nbips1=0,nbips2=0;
      
      for (uint n=0;n<cliquesDuSite[random[i]].size();n++){
        uint aux = cliquesDuSite[random[i]][n];
        if (cliques[aux].type == DATADRIVEN){
          globalenergieslabels[k]+=cliques[aux].updateEnergy(random[i],old,false,nbsujets);
//           cout << " dd:"<<cliques[aux].updateEnergy(random[i],old,false,nbsujets);
        }
        else if (cliques[aux].type == SIMILARITY){
          globalenergieslabels[k]+=cliques[aux].updateEnergy(random[i],old,false,nbsujets);
//           cout << " sim:"<< cliques[aux].updateEnergy(random[i],old,false,nbsujets);
          uint index=0;
          if (cliques[aux].blobs[0]->index==random[i]) index = 1;
          else if (cliques[aux].blobs[1]->index==random[i]) index = 0;
          else ASSERT(false);
          if (cliques[aux].blobs[index]->label == labels[k] && labels[k] != 0) nclsim1++;
          if (cliques[aux].blobs[index]->label == old && old != 0) nclsim2++;
        }
      }
      for (uint n=0;n<ipscliques.size();n++){
        uint aux = ipscliques[n];
        if (cliques[aux].blobs[0]->subject != sites[random[i]]->subject){
          if (labels[k] !=0)
            nbips1+=cliques[aux].labelscount[labels[k]];
          if (old != 0)
            nbips2+=cliques[aux].labelscount[old];
        }
      }
      
      
          
//       ASSERT(nbips1>=nclsim1 || (cout << "["<<nbips1 << ">=" << nclsim1 << "]"<<endl && false));
//       ASSERT(nbips2>=nclsim2 || (cout << "["<< nbips2 << ">=" << nclsim2 << "]"<< endl && false));
      total[k] = (nbips1-nclsim1 - (nbips2-nclsim2));
//       cout << " nips1:" << nbips1 << " nclsim1:" << nclsim1 << " nips2:" << nbips2 << " nclsim2:"<<nclsim2 << " total[k]:" << total[k];
      
      globalenergieslabels[k] += 4.0 * total[k];
//       cout << " gel:" << globalenergieslabels[k];
      

//       }
//       if (globalenergieslabels[k]/temp > 700.0) globalenergieslabels[k] = 700.0*temp;
//      /* else */if (globalenergieslabels[k]/temp < -700.0) { cout << "g:"<< globalenergieslabels[k] << ";gt" << globalenergieslabels[k]/temp << ";e" << exp(-globalenergieslabels[k]/temp) << "=>" ; globalenergieslabels[k] = -700.0*temp; cout << exp(-globalenergieslabels[k]/temp) << " ";}
//       if( E >= 0 && E < 0.1 ) E = 0.1;
      somme += exp(-globalenergieslabels[k]/temp);
      if (isnan(somme)) cout << "#####################################" << endl;
    }
//     cout << "|||||||";

    long double somme2=0.0;
    uint acc;
    if (somme > exp(700.0)) {
      acc=0;
      for (uint k=0;k<labels.size();k++)
        if (globalenergieslabels[k]<globalenergieslabels[acc]) acc = k;
      
    } 
    else {
      if (somme > exp(700.0)) cout << "dist:[";
      for (uint k=0;k<labels.size();k++){
        somme2 += exp(-globalenergieslabels[k]/temp)/somme;
        expenergies[k] = somme2;
        if (somme > exp(700.0)){
          cout << somme2 << " " ;
        }
      }
      if (somme > exp(700.0)) cout << "]";
        
      long double tirage = ((long double)UniformRandom() * 1.0);

      for (acc=0;acc<expenergies.size() && expenergies[acc]<tirage;acc++){} //cout << globalenergieslabels[acc] << " " ;
    }
    if (old != labels[acc]) {
//       cout << " acc:" << acc;
      sites[random[i]]->label=labels[acc];
      for (uint m=0;m<cliquesDuSite[random[i]].size();m++){
         energy +=cliques[cliquesDuSite[random[i]][m]].updateEnergy(random[i],old,true,nbsujets);
      }
//       cout << "GLOBAL : " << globalenergieslabels[acc] << " ";
//       energy = globalenergieslabels[acc];
      energy += 4.0*total[acc];
//       cout << " E!:" <<energy;
//       cout << "gel:"<<energy << " " ;
      mod++;
//       cin >> mod;
    }
    else {
      sites[random[i]]->label = old;
    }
//     energy = getTotalEnergy();

  }
  FILE * f;   f = fopen ("/home/grg/energy.txt","a");      fprintf(f, "%3lf\n", (float)energy); fclose(f);

}

void Anneal::Run(){
  //   Initialization();
//   for (uint i=0;i<sites.size();i++)
//     sites[i]->label = 0; //(int)sites[i]->tValue;
  for (uint k=0;k<cliques.size();k++){
    cliques[k].updateLabelsCount();
    cliques[k].computeEnergy(true,nbsujets);
  }

  for (uint i=0;i<cliques.size();i++)
    if (cliques[i].type == INTRAPRIMALSKETCH)
      ipscliques.push_back(i);
  
  cout << ipscliques.size() << " cliques intraps" << endl;
  
  energy = getTotalEnergy();

  cout << "energie initiale : " << energy << endl;
  ShortSummaryLabels();



  
  
  vector<int> indices_start;
  for(uint i=0;i<sites.size();i++)
    indices_start.push_back(i);
  
  long double temp=300.0;

  uint mod=1, ite=0, nb_under_threshold=0,test=1;

  cout.precision(2);       

//   cout.setf(ios_base::fixed, ios_base::floatfield);

//     blobsnodes[sites[i]->subject]
  
  FILE * f1;   f1 = fopen ("/home/grg/recuit.txt","w");      
  cin >> mod;

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

//     float Eintra, Edd, Esim, Els, energy2;
//     Eintra = getTypeEnergy(INTRAPRIMALSKETCH);
//     Edd = getTypeEnergy(DATADRIVEN);
//     Esim = getTypeEnergy(SIMILARITY);
//     Els = getTypeEnergy(BESTLOWERSCALE);
//     energy2 = getTotalEnergy();
//     cout << "E:" << energy2 << " (" << Eintra << ";" << Edd << ";" << Esim << ";" << Els << ") chg:" << mod << " " << flush;
//     cout << "E:" << Edd+Esim+Esimil << " (" << Edd << ";" << Esim<< ";" << Esimil << ") 
    cout << " chg:" << mod << " " << flush;


    ShortSummaryLabels();
    cout << " E=" << energy << endl; //"/" << getTotalEnergy() << endl;
    temp = temp*0.98;
//     if (energy<-7800.0){
//       cin >> test;
//     }
  }

  for (uint k=0;k<cliques.size();k++){
    cliques[k].updateLabelsCount();
    cliques[k].computeEnergy(true,nbsujets);
  }
  energy = getTotalEnergy();

  cout << "energie finale : " << energy << endl;
  ShortSummaryLabels();
  
  fclose(f1);
//   for (uint i=0;i<sites.size();i++)
//     sites[i]->label = (uint)sites[i]->t;
}
