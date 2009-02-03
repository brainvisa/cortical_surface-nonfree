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
  set<uint>::iterator it;
  for (uint i=0; i<random.size();i++){
    somme=0.0;
    
    old = sites[random[i]]->label;
    
    vector<int> zoneLab;
//     zoneLab.push_back(0);
    uint no_overlap=0;
//     cin >> no_overlap ;
//     cout << sites[random[i]]->boundingbox_min[0] << " " << sites[random[i]]->boundingbox_min[1] << " " << sites[random[i]]->boundingbox_max[0] << " " << sites[random[i]]->boundingbox_max[1] << " "  << endl;
//     cout << sites[random[i]]->nodes_list.size() << endl;
//     for (int k=0;k<labelsZones.size();k++){
//       Point3df bbmin1 = sites[random[i]]->boundingbox_min, bbmax1 = sites[random[i]]->boundingbox_max;
//       no_overlap=0;
//       double rec = getOverlap(bbmin1, bbmax1, Point3df(labelsZones[k].first[0],labelsZones[k].first[1],0.0), Point3df(labelsZones[k].second[0],labelsZones[k].second[1] ,0.0), &no_overlap);
//       if (no_overlap == 0){
//         zoneLab.push_back(k+1);
//       }
//     }
    
    for (it=listeZones[random[i]].begin();it!=listeZones[random[i]].end();it++)
      zoneLab.push_back(*it);
    vector<long double> globalenergieslabels(zoneLab.size()), expenergies(zoneLab.size()),total(zoneLab.size());
// cout << endl;
//     cout << zoneLab.size() << endl;
    for (uint k=0;k<zoneLab.size();k++){
      sites[random[i]]->label = zoneLab[k];
      globalenergieslabels[k]=energy;
      
      int nclsim1=0,nclsim2=0,nbips1=0,nbips2=0;
      
      for (uint n=0;n<cliquesDuSite[random[i]].size();n++){
        uint aux = cliquesDuSite[random[i]][n];
        if (cliques[aux].type == DATADRIVEN){
          globalenergieslabels[k]+=cliques[aux].updateEnergy(random[i],old,false,nbsujets);
        }
        else if (cliques[aux].type == SIMILARITY){
          globalenergieslabels[k]+=cliques[aux].updateEnergy(random[i],old,false,nbsujets);
          uint index=0;
          if (cliques[aux].blobs[0]->index==(uint)random[i]) index = 1;
          else if (cliques[aux].blobs[1]->index==(uint)random[i]) index = 0;
          else ASSERT(false);
          if (cliques[aux].blobs[index]->label == zoneLab[k] && zoneLab[k] != 0) nclsim1++;
          if (cliques[aux].blobs[index]->label == old && old != 0) nclsim2++;
        }
      }
      for (uint n=0;n<ipscliques.size();n++){
        uint aux = ipscliques[n];
        if (cliques[aux].blobs[0]->subject != sites[random[i]]->subject){
          if (zoneLab[k] !=0)
            nbips1+=cliques[aux].labelscount[zoneLab[k]];
          if (old != 0)
            nbips2+=cliques[aux].labelscount[old];
        }
      }
      
      
          
      total[k] = (nbips1-nclsim1 - (nbips2-nclsim2));

      
      globalenergieslabels[k] += Clique::intrapsweight * total[k];
      
      somme += exp(-globalenergieslabels[k]/temp);
      if (isnan(somme)) cout << "#####################################" << endl;
    }

    long double somme2=0.0;
    uint acc;
    if (somme > exp(700.0)) {
      acc=0;
      for (uint k=0;k<zoneLab.size();k++)
        if (globalenergieslabels[k]<globalenergieslabels[acc]) acc = k;
      
    } 
    else {
      if (somme > exp(700.0)) cout << "dist:[";
      for (uint k=0;k<zoneLab.size();k++){
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
    if (old != zoneLab[acc]) {
      sites[random[i]]->label=zoneLab[acc];
      for (uint m=0;m<cliquesDuSite[random[i]].size();m++){

         energy +=cliques[cliquesDuSite[random[i]][m]].updateEnergy(random[i],old,true,nbsujets);
      }
      

      energy += Clique::intrapsweight*total[acc];
      mod++;
    }
    else {
      sites[random[i]]->label = old;
    }

  }


}

void Anneal::Run(int verbose){
  //   Initialization();
//   set<uint>::iterator it;
//   map<uint,uint> values;
//   //ATTENTION TOUT CE BORDEL CA VA POUR LES SIMULATIONS SEULEMENT
//   values[10]=1;
//   values[500]=2;
//   values[1000]=3;
//   values[2000]=4;
//   values[3000]=5;
//   for (uint i=0;i<sites.size();i++){
//     int found=-1;
//     sites[i]->label=0;
//     set<uint> nl(sites[i]->nodes_list);
//     for (it = nl.begin();it!=nl.end() && found == -1;it++){
//       if (*it==10 || *it==500 || *it==1000 || *it==2000 || *it==3000)
//         found = *it;
//     }
//     if (found != -1) 
//       if (sites[i]->t>5.875) sites[i]->label = values[found];
//   }
  for (uint i=0;i<sites.size();i++)
//     cout << sites[i]->rank << flush;
    sites[i]->label = 0; //(int)sites[i]->tValue;

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
  
  FILE * f1;   f1 = fopen (recuitpath.data(),"w");      
  FILE * f;   f = fopen (energypath.data(),"a"); 
  fprintf(f, "== DEBUT NOUVEAU RECUIT ==\n"); 
  
//   cin >> test;
//   test=0;

  while (nb_under_threshold<5 || mod!=0){
//    while (temp>200.0){
    if (mod!=0) nb_under_threshold=0;
    else nb_under_threshold++;
    cout << " T=" << temp << " it="<< ite++ << " " ;
 
    for (uint i0=0;i0<sites.size();i0++){
      fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
    }
    fprintf(f, "%3lf\n", (float)energy); 
    
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

    cout << " chg:" << mod << " " << flush;


   if (verbose == 1) ShortSummaryLabels();
// double everif= getTotalEnergy();
    cout << " E=" << energy << endl; //<< " Everif=" << everif << endl;
// if (pow(everif-energy,2)<0.01) cout << "ok"<<endl; else printf("no %.3lf\n",(double)(everif-energy));
    
    temp = temp*0.99;

  }

  for (uint k=0;k<cliques.size();k++){
    cliques[k].updateLabelsCount();
    cliques[k].computeEnergy(true,nbsujets);
  }
  energy = getTotalEnergy();

  cout << "energie finale : " << energy << endl;
  ShortSummaryLabels();
  fprintf(f, "%3lf\nFIN RECUIT\n", (float)energy); 
  
  fclose(f1);fclose(f);
//   for (uint i=0;i<sites.size();i++)
//     sites[i]->label = (uint)sites[i]->t;
}
