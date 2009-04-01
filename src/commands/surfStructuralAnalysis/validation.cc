#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "validation.h"

using namespace aims;
using namespace carto;
using namespace std;


vector<int> StructuralAnalysis_Validation::creerHisto(vector<double> &samples, uint histosize, float *mini, float *step){
  
  vector<int> histo(histosize);
  for (uint j=0;j<histosize;j++)
    histo[j]=0;
  float maxi=-10000000.0;
  *step=0.0;
  *mini=10000000.0;
  for (uint j=0; j<samples.size();j++){
    if (samples[j]>maxi) maxi = samples[j];
    if (samples[j]<*mini) *mini = samples[j];
  }
  *step = (maxi-*mini)/(float)histosize;
  cout << "mini:"<<*mini<< " maxi:" << maxi << " step:" << *step << " nombre d'échantillons:" << samples.size() << endl;
  for (uint j=0;j<samples.size();j++){
    if(samples[j]==maxi) histo[histosize-1]++;
    else histo[(samples[j]-*mini)/(maxi-*mini)*histosize]++;
  }
  cout << endl;


  return histo;
}

vector<int> StructuralAnalysis_Validation::creerHisto2(vector<double> &samples, double step, float *mini){
  

  float maxi=-10000000.0;
  uint histosize;
  *mini=10000000.0;
  for (uint j=0; j<samples.size();j++){
    if (samples[j]>maxi) maxi = samples[j];
    if (samples[j]<*mini) *mini = samples[j];
  }
  histosize = (maxi-*mini)/(float)step;
  vector<int> histo(histosize);
  for (uint j=0;j<histosize;j++)
    histo[j]=0;
  cout << "mini:"<<*mini<< " maxi:" << maxi << " step:" << step << " histosize:" << histosize << " nombre d'échantillons:" << samples.size() << endl;
  for (uint j=0;j<samples.size();j++){
    if(samples[j]==maxi) histo[histosize-1]++;
    else histo[(samples[j]-*mini)/(maxi-*mini)*histosize]++;
  }
  cout << endl;


  return histo;
}

void StructuralAnalysis_Validation::printHisto(vector<int> &histo, float mini, float step, int type, FILE *f){
  if (type==HORIZONTAL){
    cout.precision(3);
    for (uint j=0;j<histo.size();j++)
      cout << (float)(j)*(float)step+(float)mini << "  ";
    cout << endl;
    for (uint j=0;j<histo.size();j++)
      cout << histo[j] << "  ";
    cout << endl;
    cout.precision(5);
  }
  else if (type==VERTICAL){
    for (uint j=0;j<histo.size();j++){
//             cout << (float)(j)*(float)step+(float)mini << " " << histo[j] << endl;
      fprintf(f, "%lf %d\n", (float)(j)*(float)step+(float)mini , histo[j]);
    }
    cout << endl;
  }
}

void StructuralAnalysis_Validation::printFile(vector<double> &samples, FILE *f){

  for (uint j=0;j<samples.size();j++){
//             cout << (float)(j)*(float)step+(float)mini << " " << histo[j] << endl;
//             if (samples[j] > 10000.0) cout << "attention outlier " ;
    fprintf(f, "%lf\n", samples[j]);
  }
  cout << endl;

}
// vector<float> StructuralAnalysis_Validation::getPseudoSamplesPermut(vector<uint> &listeBlobs) {
// vector<float> samples;
// float compac=0.0, Ttest=0.0, sum=0.0, energy=0.0;
// vector<uint> permut(listeBlobs.size());
// for (uint k=0;k<permut.size();k++)
//   permut[k]=0;
// bool stop=false;
// if (listeBlobs.size()==0) stop=true;
// 
// 
// 
// while(!stop){
//   uint site_courant=permut.size()-1;
//   
// 
//   permut[site_courant]++;
//   while (permut[site_courant]==2){
//     permut[site_courant]=0;
//       site_courant--;
//       permut[site_courant]++;
//   }
//   compac=0.0; Ttest=0.0; sum=0.0; energy=0.0;
//   for (uint j=0;j<permut.size();j++)
//     if (permut[j]==0) compac += sites[listeBlobs[j]]->t;
//     else compac -= ssb->sites[listeBlobs[j]]->t;
//   
//   compac /= listeBlobs.size();
//   for (uint j=0;j<permut.size();j++){
//     if (permut[j]!=0) 
//       ssb->sites[listeBlobs[j]]->t = -ssb->sites[listeBlobs[j]]->t;
//     sum += pow(ssb->sites[listeBlobs[j]]->t-compac,2);
//   }
// 
//   energy = ssb->getLabelEnergy(ssb->sites[listeBlobs[0]]->label);
//   for (uint j=0;j<permut.size();j++)
//     if (permut[j]!=0) 
//       ssb->sites[listeBlobs[j]]->t = -ssb->sites[listeBlobs[j]]->t;
// 
// 
//   Ttest = sqrt((float)(listeBlobs.size()*(listeBlobs.size()-1)))*compac/sqrt(sum);
// //   cout << compac << " " << Ttest << " " << energy <<  endl;
//   samples.push_back(Ttest);
// 
//   stop=true;
//   for (uint j=0;j<permut.size()&&stop==true;j++)
//     if (permut[j]!=1) 
//       stop = false;
//   
// }
// return samples;
// }


vector<double> StructuralAnalysis_Validation::getPseudoSamplesBootstrap(vector<uint> &listeBlobs, uint type) {
  vector<double> samples;
  double compac=0.0, Ttest=0.0, sum=0.0, energy=0.0;

// vector<uint> permut(listeBlobs.size());
  vector<float> old_t(listeBlobs.size());
  for (uint j=0;j<listeBlobs.size();j++)
    old_t[j] = ssb->sites[listeBlobs[j]]->tValue;

  for (uint i=0;i<5000;i++){

    for (uint j=0;j<listeBlobs.size();j++){
      ssb->sites[listeBlobs[j]]->tValue = old_t[(float)UniformRandom() * listeBlobs.size()];
    }
    compac=0.0; Ttest=0.0; sum=0.0; energy=0.0;
    for (uint j=0;j<listeBlobs.size();j++)
      compac += ssb->sites[listeBlobs[j]]->tValue;


    compac /= listeBlobs.size();
    for (uint j=0;j<listeBlobs.size();j++)
      sum += pow(ssb->sites[listeBlobs[j]]->tValue-compac,2);


//   energy = ssb->getLabelEnergy(ssb->sites[listeBlobs[0]]->label);
    for (uint j=0;j<listeBlobs.size();j++)
      ssb->sites[listeBlobs[j]]->tValue = old_t[j];


    Ttest = sqrt((float)(listeBlobs.size()*(listeBlobs.size()-1)))*compac/sqrt(sum);
    uint nbblobsrec= 0;
    double rec=0.0;
    for (uint k=0;k<listeBlobs.size()-1;k++)
      for (uint m=k+1;m<listeBlobs.size();m++){
      Point3df bbmax1=ssb->sites[listeBlobs[k]]->boundingbox_max, bbmax2=ssb->sites[listeBlobs[m]]->boundingbox_max;
      Point3df bbmin1=ssb->sites[listeBlobs[k]]->boundingbox_min, bbmin2=ssb->sites[listeBlobs[m]]->boundingbox_min;

      uint no_overlap=1;
      float reco = getOverlap(bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap);
      if (no_overlap==0) {
        rec += reco;
      }
      nbblobsrec++;

      }
//   cout << compac << " " << Ttest << " " << energy <<  endl;
      if (type == 0)
        samples.push_back(compac);
      else if (type ==1)
        samples.push_back(Ttest);
      else if (type== 2){
        samples.push_back(compac);
        samples.push_back(rec/nbblobsrec);
      }



  }


  return samples;

}
vector<double> StructuralAnalysis_Validation::getPseudoSamplesFullBootstrap2(vector<uint> &listeBlobs) {
  vector<double> samples;
  double compac=0.0, rec=0.0;
  vector<double> t,sim;
  set<double> tset,simset;
  set<double>::iterator it;
  set<uint> activblob;
  cout << " T=" << flush;
  for (uint j=0;j<listeBlobs.size();j++)
    activblob.insert(listeBlobs[j]);
  for (uint j=0;j<listeBlobs.size();j++){
    cout << ssb->sites[listeBlobs[j]]->tValue << " " ;
    tset.insert(ssb->sites[listeBlobs[j]]->tValue);
    for (uint k=0;k<ssb->cliquesDuSite[listeBlobs[j]].size();k++)
      if (ssb->cliques[ssb->cliquesDuSite[listeBlobs[j]][k]].type == SIMILARITY) {
      uint blob;
      if (ssb->cliques[ssb->cliquesDuSite[listeBlobs[j]][k]].blobs[0]->index == listeBlobs[j]) blob=1;
      else if (ssb->cliques[ssb->cliquesDuSite[listeBlobs[j]][k]].blobs[1]->index == listeBlobs[j]) blob=0;
      else assert(false);
      if (activblob.find(ssb->cliques[ssb->cliquesDuSite[listeBlobs[j]][k]].blobs[blob]->index) != activblob.end())
        simset.insert(ssb->cliques[ssb->cliquesDuSite[listeBlobs[j]][k]].rec);
      }
  }
  for (it=tset.begin();it!=tset.end();it++)
    t.push_back(*it);
  cout << " sim=";
  for (it=simset.begin();it!=simset.end();it++){
    cout << *it << " ";
    sim.push_back(*it);
  }
  cout << endl;
  for (uint i=0;i<5000;i++){
    compac = 0.0; rec=0.0;
    vector<double> sampt(t.size()), samprec(sim.size());
  
    for (uint j=0;j<sampt.size();j++){
      sampt[j] = t[(float)UniformRandom() * listeBlobs.size()];
      compac += sampt[j];
    }
    for (uint j=0;j<samprec.size();j++){
      samprec[j] = sim[(float)UniformRandom() * sim.size()];
      rec += samprec[j];
    }
    samples.push_back(compac/(double)t.size());
    samples.push_back(rec/(double)sim.size());


  }

  return samples;
}



vector<double> StructuralAnalysis_Validation::getPseudoSamplesFullBootstrap(vector<uint> &listeBlobs, uint type) {
  vector<double> samples;
  double compac=0.0, Ttest=0.0, sum=0.0, energy=0.0;

// vector<uint> permut(listeBlobs.size());
  vector<float> old_t(listeBlobs.size());

  map<uint,float> old_sim;
  map<uint,float>::iterator it,jt;
  for (uint j=0;j<listeBlobs.size();j++){
    old_t[j] = ssb->sites[listeBlobs[j]]->tValue;
    for (uint k=0;k<ssb->cliquesDuSite[listeBlobs[j]].size();k++)
      if (ssb->cliques[ssb->cliquesDuSite[listeBlobs[j]][k]].type == SIMILARITY) {
      pair<uint,float> clsim;
      clsim.first = ssb->cliquesDuSite[listeBlobs[j]][k];
      clsim.second = ssb->cliques[ssb->cliquesDuSite[listeBlobs[j]][k]].rec;
      old_sim.insert(clsim);
      }
  }

  for (uint i=0;i<500;i++){

    for (uint j=0;j<listeBlobs.size();j++){
      ssb->sites[listeBlobs[j]]->tValue = old_t[(float)UniformRandom() * listeBlobs.size()];
    }
    double rec=0.0;
    uint nbblobsrec=0;
    for (it=old_sim.begin();it!=old_sim.end();it++){
      uint random = (float)UniformRandom() * old_sim.size(),j=0;
      for (jt=old_sim.begin();j<random && jt!=old_sim.end();jt++,j++){}
      ssb->cliques[(*it).first].rec=(*jt).second;
      rec += (*jt).second;
      nbblobsrec++;
    }
    compac=0.0; Ttest=0.0; sum=0.0; energy=0.0;
    for (uint j=0;j<listeBlobs.size();j++)
      compac += ssb->sites[listeBlobs[j]]->tValue;


    compac /= listeBlobs.size();
    for (uint j=0;j<listeBlobs.size();j++)
      sum += pow(ssb->sites[listeBlobs[j]]->tValue-compac,2);


    energy = ssb->getLabelEnergy(ssb->sites[listeBlobs[0]]->label);
    for (uint j=0;j<listeBlobs.size();j++)
      ssb->sites[listeBlobs[j]]->tValue = old_t[j];
    for (it=old_sim.begin();it!=old_sim.end();it++){
      ssb->cliques[(*it).first].rec=(*it).second;
    }


    Ttest = sqrt((float)(listeBlobs.size()*(listeBlobs.size()-1)))*compac/sqrt(sum);
//   cout << compac << " " << Ttest << " " << energy <<  endl;
    if (type==0)
      samples.push_back(energy);
    else if (type ==2){
      samples.push_back(compac);
      samples.push_back(rec/nbblobsrec);
    }



  }
  return samples;

}

// void StructuralAnalysis_Validation::getRandomLabelsEnergies(long double nb, FILE *f){
//   vector<uint> old;
// //   vector<double> samples(nb);
//   double sample;
//   
//   
//   
//   
//   for (uint i=0;i<ssb->sites.size();i++)
//     old.push_back(ssb->sites[i]->label);
//   for (uint i=0;i<nb;i++){
//     if (i%1000000==0) cout << i << "/" << nb << endl;
//     for (uint j=0;j<ssb->sites.size();j++){
//       ssb->sites[j]->label = ssb->labels[(float)UniformRandom() * ssb->labels.size()];
//     }
//     sample = ssb->getTotalEnergy();
//     fprintf(f, "%3lf\n", (float)sample); 
//   }
//   
//   for (uint i=0;i<ssb->sites.size();i++)
//     ssb->sites[i]->label=old[i];
// }

// long double fonctionnelle(long double x){
//   if (x<0) return 0;
//   if (x==0) return 0;
//   else {long double res=1, aux=x;
//     while (aux>1){
//       res=res*aux;
// //       cout << aux<<"-"<<res << " " << flush;
//       aux--;
//     }
//     return res;
//   }
//   
// }

vector<int> StructuralAnalysis_Validation::getCompConn(vector<uint> &indicesCliques, set<uint> &listeSites){
  vector<int> comp(ssb->sites.size());
  uint blob0,blob1;
  Site *s0, *s1;
  int lcomp,nbcomp,aux;
  int label0,label1;
  set<uint>::iterator it;
  for (uint i=0;i<ssb->sites.size();i++)
    comp[i]=-1;
  lcomp=0;
  nbcomp=0;
  for (it = listeSites.begin();it != listeSites.end();it++)
    comp[*it]=0;
  for (uint i=0;i<indicesCliques.size();i++)
    if (ssb->cliques[indicesCliques[i]].type == SIMILARITY){
    s0 = ssb->cliques[indicesCliques[i]].blobs[0];
    s1 = ssb->cliques[indicesCliques[i]].blobs[1];
    blob0=s0->index;
    blob1=s1->index;
    comp[blob0]=0;
    comp[blob1]=0;
    }
    for (uint i=0;i<indicesCliques.size();i++){
      if (ssb->cliques[indicesCliques[i]].type == SIMILARITY){

        s0 = ssb->cliques[indicesCliques[i]].blobs[0];
        s1 = ssb->cliques[indicesCliques[i]].blobs[1];
        blob0=s0->index;
        blob1=s1->index;

        if (comp[blob0] == 0){
          if (comp[blob1] == 0){
            lcomp++;
            nbcomp++;
            comp[blob0] = lcomp;
            comp[blob1] = lcomp;
  //             cout << "a"<< nbcomp << "-" ;
          }
          else if (comp[blob1] > 0){
            comp[blob0] = comp[blob1];
          }
        }
        else if (comp[blob1] == 0){
          comp[blob1] = comp[blob0];
        }
        else if (comp[blob0] == comp[blob1]){
    // clique cyclique ne pas s'inquiéter mais répéter 100 fois rapidement
        }
        else if (comp[blob0] > 0 && comp[blob1] > 0){ // fusion
          label0 = comp[blob0];
          label1 = comp[blob1];

          if (label0>label1) {
            aux=label0;
            label0=label1;
            label1=aux;
          }

  //           cout << nbcomp << ";" << lcomp << "(" << label1 << "-" << label0 << ")-";
          for (uint i0=0;i0<comp.size();i0++){
            if (comp[i0] == label1)
              comp[i0] = label0;
          }
          for (uint i0=0;i0<comp.size();i0++){
            if (comp[i0] == lcomp){
              comp[i0] = label1;
  //               cout << "o" ;
            }
          }
          lcomp--;
          nbcomp--;
        }
      }
    }
    vector<uint> zeros;
    for (uint i=0;i<comp.size();i++)
      if (comp[i]==0) zeros.push_back(i);
    for (uint i=0;i<zeros.size();i++){
      lcomp++;
      nbcomp++;
      comp[zeros[i]]=lcomp;
    }
//   cout << "nbcomp:" << nbcomp << endl;
    return comp;
}

vector<set<uint> > StructuralAnalysis_Validation::getCompConnVector(vector<int> &comp){
  vector<set<uint> > cc;
  uint cpt=0,cpt2,nbsites=0;
  for (uint i=0;i<comp.size();i++)
    if (comp[i]>=0) nbsites++;
//   cout << nbsites << endl;
  for (uint i0=0;cpt<nbsites;i0++){
    cpt2=0;
//     cout << i0 << " " << flush;
    cc.push_back(set<uint>());
    for (uint i=0;i<comp.size();i++)
      if (comp[i]==(int)i0) {
      cpt2++;cpt++;
      cc[i0].insert(i);
      }
  }
//   cout << endl;
  return cc;
}


uint StructuralAnalysis_Validation::nbcombinaisons(set<uint> &graphe, uint card){
  uint nbfinal=0;
  set<uint>::iterator it2;
  set<unsigned short int>::iterator it;
  set<set<unsigned short int> >::iterator itt;
  set<set<unsigned short int> > composantes;
  set<set<unsigned short int> > composantes_np1;
  uint taille =2;

  for (it2=graphe.begin();it2!=graphe.end();it2++){
      // cout << "\b\b\b\b\b\b\b\b\b\b\b\b\bit="<< *it << flush;
    set<unsigned short int> composante;
    composante.insert(*it2);
    composantes_np1.insert(composante);
  }
  cout << "nb taille 1 :" << composantes_np1.size()<< endl;

  while (taille <=card && composantes_np1.size() < 1000000){
    composantes.clear();
    composantes = set<set<unsigned short int> >(composantes_np1);
    composantes_np1.clear();
//     cout << "composantes:"<<composantes.size() << endl;
    uint aux=0;
    for (itt = composantes.begin();itt!=composantes.end() && composantes_np1.size() < 1000000;itt++){
      if (aux++%1000==0) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << aux << "/" << composantes.size()<<"("<< composantes_np1.size()<<")" << flush;
      set<unsigned short int> voisins;
//       set<unsigned short int> currcomp(*itt);
      for (it = (*itt).begin();it!= (*itt).end(); it++){
        for (uint i=0;i<ssb->cliquesDuSite[*it].size() ;i++){
          uint currclique = ssb->cliquesDuSite[*it][i];
          if (ssb->cliques[currclique].type==SIMILARITY){
//         cout << "last="<<last << " i="<<i << " "<< flush; ;
            if ((*itt).find(ssb->cliques[currclique].blobs[0]->index) == (*itt).end() && graphe.find(ssb->cliques[currclique].blobs[0]->index)!=graphe.end())
              voisins.insert(ssb->cliques[currclique].blobs[0]->index);
            if ((*itt).find(ssb->cliques[currclique].blobs[1]->index) == (*itt).end() && graphe.find(ssb->cliques[currclique].blobs[1]->index)!=graphe.end())
              voisins.insert(ssb->cliques[currclique].blobs[1]->index);
          }
        }
      }
//       cout << "voisins.size() = " << voisins.size() << endl;
      for (it = voisins.begin();it !=voisins.end();it++){
        set<unsigned short int> comp_np1((*itt));
        comp_np1.insert(*it);
        composantes_np1.insert(comp_np1);
      }

    }
    cout << "nb taille="<<taille << " :" << composantes_np1.size()<< endl;
    if (card <5){
    for (itt =composantes_np1.begin();itt!=composantes_np1.end();itt++){
      for (it = (*itt).begin();it!= (*itt).end(); it++)
        cout << *it << "-" << flush;
      cout << " " << flush;
    }
    cout << endl;
    }
    taille++;
  }
//   cout << endl;


  return composantes_np1.size();
}

double StructuralAnalysis_Validation::WalshTest(vector<double> &samplesdist, int r){
  
  std::sort(samplesdist.begin(), samplesdist.end());
  double c =(double) ceil(sqrt(2*samplesdist.size())); uint k=r+c; double b2 = 1.0/0.05;
  double a = (1.0 + sqrt(b2) * sqrt((c-b2)/(c-1)))/(c-b2-1.0);
//   Point2df res(samplesdist[0]-(1+a)*samplesdist[1]+a*samplesdist[k-1],samplesdist[samplesdist.size()-1]-(1+a)*samplesdist[samplesdist.size()-2]+a*samplesdist[samplesdist.size()-k]);
// Xr - (1+a)Xr+1 + aXk < 0
  return samplesdist[r-1] - (1+a)*samplesdist[r] + a*samplesdist[k-1];
  return samplesdist[samplesdist.size()-r]-(1+a)*samplesdist[samplesdist.size()-r-1]+a*samplesdist[samplesdist.size()-k];
//   return res;
}

vector<double> StructuralAnalysis_Validation::getCaracSample(vector<uint> &composante){
  double tmoy=0.0, rec=0.0, sum=0.0, compac=0.0, Ttest=0.0, compaccent;
  uint nbblobsrec=0;
      

  for (uint k=0;k<composante.size();k++){
    tmoy += ssb->sites[composante[k]]->t;
    compac += ssb->sites[composante[k]]->tValue;
    compaccent += ssb->sites[composante[k]]->t2;
  }
  compac /= composante.size();
  for (uint k=0;k<composante.size();k++)
    sum += pow(ssb->sites[composante[k]]->tValue-compac,2);
  Ttest = sqrt((float)(composante.size()*(composante.size()-1)))*compac/sqrt(sum);
  for (uint k=0;k<composante.size()-1;k++)
    for (uint m=k+1;m<composante.size();m++){
    Point3df bbmax1=ssb->sites[composante[k]]->boundingbox_max, bbmax2=ssb->sites[composante[m]]->boundingbox_max;
    Point3df bbmin1=ssb->sites[composante[k]]->boundingbox_min, bbmin2=ssb->sites[composante[m]]->boundingbox_min;

    uint no_overlap=1;
    float reco = getOverlap(bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap);
    if (no_overlap==0) {
      rec += reco;
    }
    nbblobsrec++;
    }
  vector<double> sample;
  sample.push_back(compac);
  sample.push_back(Ttest);
  sample.push_back(rec/(double)nbblobsrec);
  sample.push_back(ssb->getClusterEnergy(composante));
  return sample;

}

vector<double> StructuralAnalysis_Validation::getBackup(vector<uint> &composante){
  double tmoy=0.0, rec=0.0, sum=0.0, compac=0.0, Ttest=0.0;
  uint nbblobsrec=0;
  vector<double> sample;

  for (uint k=0;k<composante.size();k++){
    tmoy += ssb->sites[composante[k]]->t;
    sample.push_back(ssb->sites[composante[k]]->tValue);
    compac += ssb->sites[composante[k]]->tValue;
  }
  compac /= composante.size();
  for (uint k=0;k<composante.size();k++)
    sum += pow(ssb->sites[composante[k]]->tValue-compac,2);
  Ttest = sqrt((float)(composante.size()*(composante.size()-1)))*compac/sqrt(sum);
//         Ttest = compac/(sqrt(sum)/sqrt((float)composante.size()));
  for (uint k=0;k<composante.size()-1;k++)
    for (uint m=k+1;m<composante.size();m++){
    Point3df bbmax1=ssb->sites[composante[k]]->boundingbox_max, bbmax2=ssb->sites[composante[m]]->boundingbox_max;
    Point3df bbmin1=ssb->sites[composante[k]]->boundingbox_min, bbmin2=ssb->sites[composante[m]]->boundingbox_min;

    uint no_overlap=1;
    float reco = getOverlap(bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap);
    if (no_overlap==0) {
      rec += reco;
      sample.push_back(reco);
      nbblobsrec++;
    }
    }
    sample.push_back(compac);
    sample.push_back(sqrt(sum));
    return sample;

}


void StructuralAnalysis_Validation::ValidTestLastChance(){
  
  set<uint> existinglabels;
  set<uint>::iterator it;
  for (uint i=0;i<ssb->sites.size();i++)
    if (ssb->sites[i]->label!=0)
      existinglabels.insert(ssb->sites[i]->label);
  uint cpt=0;
  vector<set<uint> > activblobs(existinglabels.size());
  for (it = existinglabels.begin();it!=existinglabels.end();it++){
    for (uint j=0;j<ssb->sites.size();j++)
      if (ssb->sites[j]->label==(int)*it)
        activblobs[cpt].insert(j);
    cpt++;
  }

  vector<uint> cliquesV;
  set<uint> sitesV;
  for (uint j=0;j<ssb->cliques.size();j++)
    cliquesV.push_back(j);
  for (uint j=0;j<ssb->sites.size();j++)
    sitesV.insert(j);
  vector<int> ccc= getCompConn(cliquesV,sitesV);
  vector<set<uint> > cc = getCompConnVector(ccc);
  uint startsite;
  set<uint> activblobsglobal;
  for (uint i=0;i<activblobs.size();i++)
    for (it=activblobs[i].begin();it!=activblobs[i].end();it++)
      activblobsglobal.insert(*it);
  cout << "ABG=" << activblobsglobal.size() << endl;

  for (uint i=0;i<activblobs.size();i++){
    cout << activblobs[i].size() << endl;
    // on va tirer au sort des clusters de taille activblobs[i].size()
    vector<int> select_cc,tirage(ssb->sites.size());
    for (uint k=0;k<cc.size();k++)
      if (cc[k].size()>=i) select_cc.push_back(k);
    for (uint k=0;k<ssb->sites.size();k++){
      tirage[k]=-1;
      if (cc[ccc[k]].size()>=i) tirage[k] = ccc[k];
    }
    vector<double> samplesT, samplesTtest,samplesRec;
    vector<vector<double> > samplesCarac;

    for (uint j=0;j<5000;j++){
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << j << flush;
                //création du cluster
      uint startsite;
      do{
        startsite = (float)UniformRandom() * ssb->sites.size();
      }
      while(tirage[startsite]<1 && cc[tirage[startsite]].size()<activblobs[i].size() && activblobsglobal.find(*it)!=activblobsglobal.end());
      vector<uint> composante;
      composante.push_back(startsite);
                
      set<uint> dejapris = set<uint>(activblobsglobal);
                
//                 set<uint> dejapris;
//                 set<uint> autorized(cc[select_cc[cci]]);
//                 autorized.erase(*it);
//                 for (it = activblobs[i].begin();it!=activblobs[i].end();it++)
//                   autorized.erase(*it);
                
      dejapris.insert(*it);
                

//                 cout << "o" << flush;
      set<uint> voisins;
      cout << "start:"<<*it << " " << flush;
      while (composante.size()<activblobs[i].size()){
                      
        uint test = composante[composante.size()-1];
        for (uint m=0;m<ssb->cliquesDuSite[test].size();m++){
          uint currclique = ssb->cliquesDuSite[test][m];
          if (ssb->cliques[currclique].type == SIMILARITY){
            uint blob0 = ssb->cliques[currclique].blobs[0]->index, blob1 = ssb->cliques[currclique].blobs[1]->index;
//                             if (activblobs[i].find(blob0) != activblobs[i].end()) cout << "ACTIV0 " << flush;
//                             if (activblobs[i].find(blob1) != activblobs[i].end()) cout << "ACTIV1 " << flush;
            if (dejapris.find(blob0) == dejapris.end()) // && activblobs[i].find(blob0) == activblobs[i].end() )
              voisins.insert(blob0);
            if (dejapris.find(blob1) == dejapris.end()) //&& activblobs[i].find(blob1) == activblobs[i].end())
              voisins.insert(blob1);
//                             if (dejapris.find(blob0) == dejapris.end() && dejapris.find(blob1) == dejapris.end())
//                               cout << "a" << endl;
          }
        }
        if (voisins.size()==0) {
//                             cout << "RESTART " << flush;
          do{
            startsite = (float)UniformRandom() * ssb->sites.size();
          }
          while(tirage[startsite]<1 && cc[tirage[startsite]].size()<activblobs[i].size()&&activblobsglobal.find(*it)!=activblobsglobal.end());
          composante.clear();
          composante.push_back(startsite);
//                             cout << " start:" << *it << " " << flush;
          dejapris.clear();
          dejapris = set<uint>(activblobsglobal);

          dejapris.insert(*it);
          voisins.clear();
        }
        else {
          uint random = (float)UniformRandom() * voisins.size();
    //                         if (voisins.size()==0) cout << "r=" <<random << endl;
          it=voisins.begin();
          for (uint n=0;n<random;n++,it++){}
    //                         if (voisins.size()==0) cout << "it=" <<*it << endl;
          if (dejapris.find(*it)==dejapris.end()){
            composante.push_back(*it);
            dejapris.insert(*it);
            voisins.erase(it);
    //                             ssb->sites[*it]->label=i;
                            //                       }
          }
    //                       else cout << "ALARM" << endl;
        }
      }
      vector<double> sample(getCaracSample(composante));
      samplesCarac.push_back(sample);
      samplesT.push_back(sample[0]);
      samplesRec.push_back(sample[2]);

      samplesTtest.push_back(sample[1]);


    }
    vector<uint> activblob;
    for (it = activblobs[i].begin();it!=activblobs[i].end();it++)
      activblob.push_back(*it);

    vector<double> sample(getCaracSample(activblob));
    vector<double> samplesTBootstrap(getPseudoSamplesBootstrap(activblob,0));
    vector<double> samplesTtestBootstrap(getPseudoSamplesBootstrap(activblob,1));
    vector<double> samples2DBootstrap(getPseudoSamplesFullBootstrap2(activblob));

          //CREATION DES HISTOGRAMMES
    float mini,step;
    cout << "histo" << endl;
    FILE *ft,*frec,*fttest,*factiv,*fttactiv,*f2D, *fD, *f2Dactiv;
    ft=fopen("/home/grg/histo_t.txt","w");
    fttest=fopen("/home/grg/histo_ttest.txt","w");
    factiv=fopen("/home/grg/histo_tactiv.txt","w");
    fttactiv=fopen("/home/grg/histo_ttestactiv.txt","w");
    f2Dactiv=fopen("/home/grg/histo_2Dactiv.txt","w");
    frec=fopen("/home/grg/histo_rec.txt","w");
    uint histosize=40;
    f2D=fopen("/home/grg/histo_2D.txt","w");
    fD=fopen("/home/grg/histo_dist.txt","w");
    Point2df barycentre(0.0,0.0);
    for (uint j=0;j<samplesCarac.size();j++){
      barycentre[0] += samplesCarac[j][0];
      barycentre[1] += samplesCarac[j][2];
      fprintf(f2D,"%f %f\n", samplesCarac[j][0], samplesCarac[j][2]);
    }
    for (uint j=0;j<samples2DBootstrap.size()/2;j=j+2){
      fprintf(f2Dactiv,"%f %f\n", samples2DBootstrap[j], samples2DBootstrap[j+1]);
    }
          
    barycentre /= samplesCarac.size();
    vector<double> samplesdist;
    for (uint j=0;j<samplesCarac.size();j++){
      samplesdist.push_back(sqrt(pow(barycentre[0]-samplesCarac[j][0],2) + pow(barycentre[1]- samplesCarac[j][2],2)));
    }
    cout << "histoDist" << endl;
    vector<int> histodist(creerHisto(samplesdist,histosize,&mini,&step));
    printHisto(histodist,mini,step,VERTICAL,fD);

    cout << "histoRec" << endl;
    vector<int> histoRec(creerHisto(samplesRec,histosize,&mini,&step));
    printHisto(histoRec,mini,step,VERTICAL,frec);

    cout << "histoT" << endl;
    vector<int> histoT(creerHisto(samplesT,histosize,&mini,&step));
    printHisto(histoT,mini,step,VERTICAL,ft);

    cout << "histoTtest" << endl;
    vector<int> histoTtest(creerHisto(samplesTtest,histosize,&mini,&step));
    printHisto(histoTtest,mini,step,VERTICAL,fttest);

    cout << "histoTActiv" << endl;
    vector<int> histoTActiv(creerHisto(samplesTBootstrap,histosize,&mini,&step));
    printHisto(histoTActiv,mini,step,VERTICAL,factiv);

    cout << "histoTtestActiv" << endl;
    vector<int> histoTtestActiv(creerHisto(samplesTtestBootstrap,histosize,&mini,&step));
    printHisto(histoTtestActiv,mini,step,VERTICAL,fttactiv);

    fclose(ft);fclose(fttest);fclose(factiv);fclose(fttactiv);fclose(fD); fclose(f2D);fclose(f2Dactiv);fclose(frec);

    cout << "label " << ssb->sites[activblob[0]]->label << " tmoy=" << sample[0] << " ttest =" << sample[1] <<" rec=" << sample[2] << endl;

    vector<double> samplestout(samplesT);
    for (uint i=0;i<samplesTBootstrap.size();i++)
      samplestout.push_back(samplesTBootstrap[i]);
    std::sort(samplestout.begin(),samplestout.end());
          
    cout << " wn=" << WalshTest(samplestout,250) << endl;
    set<double> setBoot;
    set<double>::iterator dit;
    for (uint j=0;j<samplesTBootstrap.size();j++)
      setBoot.insert(samplesTBootstrap[j]);
    if (setBoot.size()>10){
      dit=setBoot.begin();
      for (uint j=0;j<setBoot.size();j++,dit++){
        cout << samplesT[j] << " " << *dit << endl;
      }
      for (uint j=setBoot.size();j<samplesT.size();j++){
        cout << samplesT[j] << endl;
      }
    }
    cin >> mini;
    

  }
  


}




void StructuralAnalysis_Validation::ValidTest(){
  for (uint i=0;i<ssb->cliques.size();i++){
    if (ssb->cliques[i].type == SIMILARITY){
      assert(ssb->cliques[i].blobs[0]->index <ssb->sites.size());
      assert(ssb->cliques[i].blobs[1]->index <ssb->sites.size());
    }
  }
  cout << endl;
    
  vector<uint> cliquesV;
  set<uint> sitesV;
  for (uint j=0;j<ssb->cliques.size();j++)
    cliquesV.push_back(j);
  for (uint j=0;j<ssb->sites.size();j++)
   sitesV.insert(j);
  vector<int> ccc= getCompConn(cliquesV,sitesV);
  vector<set<uint> > cc = getCompConnVector(ccc);
    
  for (uint k=0;k<cc.size();k++){
    cout << cc[k].size() << " ";
  }

  vector<int> old;
  for (uint k=0;k<ssb->sites.size();k++){
    old.push_back(ssb->sites[k]->label);
    ssb->sites[k]->label = 0;
  }

    

  for (uint i=2;i<=ssb->nbsujets;i++) {
    cout << i << "blobs dans une composante "<< endl;
    set<uint>::iterator it;
    vector<int> select_cc,tirage(ssb->sites.size());
    for (uint k=0;k<cc.size();k++)
      if (cc[k].size()>=i) select_cc.push_back(k);
    for (uint k=0;k<ssb->sites.size();k++){
      tirage[k]=-1;
      if (cc[ccc[k]].size()>=i) tirage[k] = ccc[k];
    }
    vector<double> samples, samplesT, samplesRec, samplesTtest, samplesNeg; vector<vector<double> > backupTtest;
    vector<vector<double> > samplesCarac;
    cout << "neg:" << flush;
    for (uint j=0;j<10000;j++){
//                 cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << j << flush;
                //création du cluster
      uint startsite;
      do{
        startsite = (float)UniformRandom() * ssb->sites.size();
      }
      while(tirage[startsite]<1);
//                 it = cc[select_cc[cci]].begin();
//                 for (uint n=0;n<startsite;n++,it++){}
      vector<uint> composante;
      composante.push_back(startsite);
      set<uint> dejapris;
      dejapris.insert(startsite);

//                 cout << "o" << flush;
      set<uint> voisins;
      while (composante.size()<i){
                      
        uint test = composante[composante.size()-1];
        for (uint m=0;m<ssb->cliquesDuSite[test].size();m++){
          uint currclique = ssb->cliquesDuSite[test][m];
          if (ssb->cliques[currclique].type == SIMILARITY){
            if (dejapris.find(ssb->cliques[currclique].blobs[0]->index) == dejapris.end())
              voisins.insert(ssb->cliques[currclique].blobs[0]->index);
            if (dejapris.find(ssb->cliques[currclique].blobs[1]->index) == dejapris.end())
              voisins.insert(ssb->cliques[currclique].blobs[1]->index);
          }
        }
        if (voisins.size()==0) {
//                             cout << "RESTART " << flush;
          do{
            startsite = (float)UniformRandom() * ssb->sites.size();
          }
          while(tirage[startsite]<1);
//                             for (uint k=0;k<composante.size();k++)
//                               ssb->sites[composante[k]]->label=0;
          composante.clear();
          composante.push_back(startsite);
//                             cout << " start:" << startsite << " " << flush;
          dejapris.clear();

          dejapris.insert(startsite);
          voisins.clear();
        }
        else {
          uint random = (float)UniformRandom() * voisins.size();
          it=voisins.begin();
          for (uint n=0;n<random;n++,it++){}
          if (dejapris.find(*it)==dejapris.end()){
            composante.push_back(*it);
            dejapris.insert(*it);
//                             ssb->sites[*it]->label=i;
            voisins.erase(it);
          }
        }
//                       else cout << "ALARM" << endl;
      }

//                 cout << "o" << flush;
      vector<double> sample(getCaracSample(composante));
      samplesCarac.push_back(sample);
//                 cout << samplesCarac[samplesCarac.size()-1][0] << " " << samplesCarac[samplesCarac.size()-1][2] << "("<< samplesCarac[samplesCarac.size()-1][3] <<")"<< endl;
      samplesT.push_back(sample[0]);
      samplesRec.push_back(sample[2]);
                
      samplesTtest.push_back(sample[1]);
      vector<double> backup(getBackup(composante));
      backup.push_back(sample[1]);
      backupTtest.push_back(backup);
      long double E =ssb->getClusterEnergy(composante);
      samples.push_back(E);
      if (E<0.0){ cout << sample[0] << "-" << sample[2] << " " << flush;
        samplesNeg.push_back(sample[0]);}
                
//                 for (uint k=0;k<composante.size();k++){
//                   ssb->sites[composante[k]]->label=0;
//                 }
    }


    cout << endl;
          //CREATION DES HISTOGRAMMES
    float mini,step;
    cout << "histo" << endl;
    FILE *f,*ft,*fttest,*frec,*f2D,*fD, *fttestRaw,*fNeg;
    f=fopen("/home/grg/histo.txt","w");
    ft=fopen("/home/grg/histo_t.txt","w");
    fttest=fopen("/home/grg/histo_ttest.txt","w");
    frec=fopen("/home/grg/histo_rec.txt","w");
    f2D=fopen("/home/grg/histo_2D.txt","w");
    fD=fopen("/home/grg/histo_dist.txt","w");
    fttestRaw=fopen("/home/grg/plot_ttest.txt","w");
    fNeg=fopen("/home/grg/histo_neg.txt","w");

    uint histosize=40;
    Point2df barycentre(0.0,0.0);
    for (uint j=0;j<samplesCarac.size();j++){
      barycentre[0] += samplesCarac[j][0];
      barycentre[1] += samplesCarac[j][2];
      fprintf(f2D,"%f %f\n", samplesCarac[j][0], samplesCarac[j][2]);
    }
    barycentre /= samplesCarac.size();
    vector<double> samplesdist;
    for (uint j=0;j<samplesCarac.size();j++){
      samplesdist.push_back(sqrt(pow(barycentre[0]-samplesCarac[j][0],2) + pow(barycentre[1]- samplesCarac[j][2],2)));
    }
    cout << "histoDist" << endl;
    vector<int> histodist(creerHisto(samplesdist,histosize,&mini,&step));
    printHisto(histodist,mini,step,VERTICAL,fD);
          
    vector<int> histo(creerHisto(samples,histosize,&mini,&step));
    printHisto(histo,mini,step,VERTICAL,f);
    histosize=40;
    cout << "histoT" << endl;
    vector<int> histoT(creerHisto(samplesT,histosize,&mini,&step));
    printHisto(histoT,mini,step,VERTICAL,ft);

    cout << "histoRec" << endl;
    vector<int> histoRec(creerHisto(samplesRec,histosize,&mini,&step));
    printHisto(histoRec,mini,step,VERTICAL,frec);


    cout << "histoTtest" << endl;
    vector<int> histoTtest(creerHisto(samplesTtest,histosize,&mini,&step));
    printHisto(histoTtest,mini,step,VERTICAL,fttest);
    printFile(samplesTtest,fttestRaw);

    cout << "histoTneg" << endl;
    vector<int> histoNeg(creerHisto(samplesNeg,histosize,&mini,&step));
    printHisto(histoNeg,mini,step,VERTICAL,fNeg);

//           std::sort(histodist.begin(), histodist.end());

    vector<double> integ(histosize);
    double sum2=0.0;
    for (uint j=0;j<integ.size();j++){
      integ[j] = (float)1.0 * (float)histoTtest[j] + (float)sum2;
            // cout << integ[j] << " " << flush;
      sum2 += histoTtest[j];
    }
    float alpha = 0.0;
    float percent = 99.0;
    cout << percent*sum2/100.0 << endl;
    for (uint j=0;j<integ.size() && integ[j]<percent*sum2/100.0;j++) alpha=j;
    double diff = samplesTtest.size()-sum2;
    cout << sum2 << " diff (=0):" << diff << " alpha("<<percent <<"%) : " << alpha*step+mini << endl;

    for (int j=1;j<(int)ssb->labels.size();j++){
      vector<uint> composante;
      for (uint k=0;k<ssb->sites.size();k++){
        if (old[k]==j){
          composante.push_back(k);
        }
      }
//             cout << composante.size() << endl;
      if (composante.size()==i){
        vector<double> sample(getCaracSample(composante));
        vector<double> samplesdist2(samples);
        long double energysample = ssb->getClusterEnergy(composante);
        double dist = sqrt(pow(barycentre[0]-samples[0],2) + pow(barycentre[1]- sample[2],2));
        std::sort(samplesdist2.begin(),samplesdist2.end());
//               cout << "label " << j << " tmoy=" << sample[0] << " ttest =" << sample[1] <<" rec=" << sample[2] << " dist=" << dist << " wn="<< WalshTest(samplesdist2,1) << endl;
        samplesdist2.push_back(energysample); //dist);
        std::sort(samplesdist2.begin(),samplesdist2.end());
        int r;
//               for (r=1;r>=0 && samplesdist2[samplesdist2.size()-r] != energysample;r++)
//                 {}
//               assert (r!=-1);
        for (r=0;r<(int)samplesdist2.size() && pow(samplesdist2[r] - energysample,2) >0.001;r++)
        {}
        cout << r << endl;
        assert (r<(int)samplesdist2.size());
        r++;
            
              
              
        cout << "label " << j << " tmoy=" << sample[0] << " ttest =" << sample[1] <<" rec=" << sample[2] << " dist=" << dist << "("<< samplesdist2[samplesdist2.size()-1] <<") wn="<< WalshTest(samplesdist2, r) << " r=" << r << endl;
        cout << "label:" << j << " ttest=" << sample[1];
                
        vector<double> labelbackup(getBackup(composante));
        for (uint n=0;n<i;n++)
          cout << " t=" << labelbackup[n];
        for (uint n=i;n<labelbackup.size()-2;n++)
          cout << " rec=" <<labelbackup[n];
        cout << " moy="<< labelbackup[labelbackup.size()-2] << " std=" << labelbackup[labelbackup.size()-1];
        cout << endl<<endl;
//               for (uint k=0;k<3;k++){
//                 while (samplesdist2[samplesdist2.size()-1-k] == sample[1]) k++;
//                 cout << "rank:" << k+1 << " ttest=" << samplesdist2[samplesdist2.size()-1-k] ;
//                 uint m;
//                 for (m=0;m<backupTtest.size()&& backupTtest[m][backupTtest[m].size()-1] != samplesdist2[samplesdist2.size()-1-k];m++){}
//                 for (uint n=0;n<i;n++)
//                   cout << " t=" << backupTtest[m][n];
//                 for (uint n=i;n<backupTtest[m].size()-3;n++)
//                   cout << " rec=" <<backupTtest[m][n];
//                 cout << " moy="<< backupTtest[m][backupTtest[m].size()-3] << " std=" << backupTtest[m][backupTtest[m].size()-2];
//                 cout << endl<<endl;
//               }
      }
    }



    fclose(f);fclose(ft);fclose(fttest);fclose(frec);fclose(f2D);fclose(fD);fclose(fttestRaw);fclose(fNeg);
    cin >> mini;
  }
  for (uint k=0;k<ssb->sites.size();k++)
    ssb->sites[k]->label=old[k];

}

void StructuralAnalysis_Validation::ValidAround(){
  set<uint> existinglabels;
  set<uint>::iterator it;
  for (uint i=0;i<ssb->sites.size();i++)
    if (ssb->sites[i]->label!=0)
      existinglabels.insert(ssb->sites[i]->label);
  uint cpt=0;
  vector<set<uint> > activblobs(existinglabels.size());
  for (it = existinglabels.begin();it!=existinglabels.end();it++){
    for (uint j=0;j<ssb->sites.size();j++)
      if (ssb->sites[j]->label==(int)*it)
        activblobs[cpt].insert(j);
    cpt++;
  }

  vector<uint> cliquesV;
  set<uint> sitesV;
  for (uint j=0;j<ssb->cliques.size();j++)
    cliquesV.push_back(j);
  for (uint j=0;j<ssb->sites.size();j++)
    sitesV.insert(j);
  vector<int> ccc= getCompConn(cliquesV,sitesV);
  vector<set<uint> > cc = getCompConnVector(ccc);
  uint startsite;
  set<uint> activblobsglobal;
  for (uint i=0;i<activblobs.size();i++)
    for (it=activblobs[i].begin();it!=activblobs[i].end();it++)
      activblobsglobal.insert(*it);
  cout << "ABG=" << activblobsglobal.size() << endl;

  
  for (uint i=0;i<activblobs.size();i++){
          set<uint> forbidden,autorized(activblobs[i]);
          
            
          cout << endl << activblobs[i].size() << endl;
          for (it=activblobs[i].begin();it!=activblobs[i].end();it++)
            cout << *it << " ";
          cout << endl;
          
          // on va tirer au sort des clusters de taille activblobs[i].size()
          vector<int> select_cc,tirage(ssb->sites.size());
          for (uint k=0;k<cc.size();k++)
            if (cc[k].size()>=activblobs[i].size()) select_cc.push_back(k);
          for (uint k=0;k<ssb->sites.size();k++){
            tirage[k]=-1;
            if (cc[ccc[k]].size()>=activblobs[i].size()) tirage[k] = ccc[k];
          }
          vector<double> samples,samplesT, samplesTtest,samplesRec,samplesNeg;
          vector<vector<double> > samplesCarac;
          set<uint>::iterator it;
          uint j=0;
          for (it = activblobs[i].begin();it!=activblobs[i].end();it++,j++){
            uint test = *it;
            for (uint m=0; m<ssb->cliquesDuSite[test].size();m++){
              uint currclique = ssb->cliquesDuSite[test][m];
              if (ssb->cliques[currclique].type == SIMILARITY){
                autorized.insert(ssb->cliques[currclique].blobs[0]->index);
                autorized.insert(ssb->cliques[currclique].blobs[1]->index);
              }
            }
          }
//           autorized.insert(0);

//           for (uint j=0;j<ssb->sites.size();j++){
//             if (autorized.find(j)==autorized.end()) {
//               forbidden.insert(j);
//             }
//           }
          cout << "restent :"<< ssb->sites.size()-forbidden.size() << endl;
//           cout << nbcombinaisons(autorized,activblobs[i].size()) << "possibilités"<< endl;
//           for (uint j=0;j<ssb->sites.size();j++)
//             if (forbidden.find(j)==forbidden.end())
//               cout << j << "(" << ccc[j]<< ") ";
//           cout << endl;

          for (uint j=0;j<10000;j++){
//                 cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << j << flush;
                //création du cluster
                uint startsite;
                
                set<uint> dejapris = set<uint>(forbidden);
                do{ startsite = (float)UniformRandom() * ssb->sites.size(); }
                while(tirage[startsite]<1 || dejapris.find(startsite)!=dejapris.end());
                vector<uint> composante;
                composante.push_back(startsite);
                          
//                 if (startsite==*(activblobs[i].begin())) cout << startsite << " " << flush;

                          
                dejapris.insert(startsite);
                          

                set<uint> voisins;
//                 cout << "start:"<<*it << " " << flush;
                while (composante.size()<activblobs[i].size()){
                      
                  uint test = composante[composante.size()-1];
                  for (uint m=0;m<ssb->cliquesDuSite[test].size();m++){
                    uint currclique = ssb->cliquesDuSite[test][m];
                    if (ssb->cliques[currclique].type == SIMILARITY){
                      
                      if (dejapris.find(ssb->cliques[currclique].blobs[0]->index) == dejapris.end())
                        voisins.insert(ssb->cliques[currclique].blobs[0]->index);
                      if (dejapris.find(ssb->cliques[currclique].blobs[1]->index) == dejapris.end())
                        voisins.insert(ssb->cliques[currclique].blobs[1]->index);
//                       cout << ssb->cliques[currclique].blobs[0]->index << "-" << ssb->cliques[currclique].blobs[1]->index << " " << flush;
                    }
                  }
                  assert(voisins.size()!=0);
//                   if (voisins.size()==0) {
//                     dejapris.clear();
//                     dejapris = set<uint>(forbidden);
//                     do{
//                       startsite = (float)UniformRandom() * ssb->sites.size();
//                     }
//                     while(tirage[startsite]<1 || dejapris.find(startsite) != dejapris.end());
//                     composante.clear();
//                     composante.push_back(startsite);   
//                     dejapris.insert(startsite);
//                     voisins.clear();
//                   }
//                   else {
                    uint random = (float)UniformRandom() * voisins.size();
                    it=voisins.begin();
                    for (uint n=0;n<random;n++,it++){}
                    if (dejapris.find(*it)==dejapris.end()){
                      composante.push_back(*it);
                      dejapris.insert(*it);
                      voisins.erase(it);
                    }
//                   }
                }
                
                vector<double> sample(getCaracSample(composante));
                samplesCarac.push_back(sample);
                samples.push_back(sample[3]);
                samplesT.push_back(sample[0]);
                samplesRec.push_back(sample[2]);
                if (sample[3]<0.0) { samplesNeg.push_back(sample[0]);
                }
                
          }
      
          float mini,step;
          FILE *f,*ft,*frec,*fNeg,*fAct,*fActRec;
          f=fopen("/home/grg/histo.txt","w");
          ft=fopen("/home/grg/histo_t.txt","w");
          frec=fopen("/home/grg/histo_rec.txt","w");
          fNeg=fopen("/home/grg/histo_neg.txt","w");
          fAct=fopen("/home/grg/barres_activ.txt","w");
          fActRec=fopen("/home/grg/barres_activrec.txt","w");
          cout << endl<< "nombre de neg:"<< samplesNeg.size()<< endl;
          uint histosize=40;

          
          cout << "histoT" << endl;
          vector<int> histoT(creerHisto(samplesT,histosize,&mini,&step));
          printHisto(histoT,mini,step,VERTICAL,ft);

          cout << "histoRec" << endl;
          vector<int> histoRec(creerHisto(samplesRec,histosize,&mini,&step));
          printHisto(histoRec,mini,step,VERTICAL,frec);

          cout << "histoTneg" << endl;
          vector<int> histoNeg(creerHisto(samplesNeg,histosize,&mini,&step));
          printHisto(histoNeg,mini,step,VERTICAL,fNeg);
          
          cout << "histo" << endl;
          vector<int> histo(creerHisto(samples,histosize,&mini,&step));
          printHisto(histo,mini,step,VERTICAL,f);

          int maxiheight=-10000, maxiheight2=-10000;
          for (uint j=0;j<histosize;j++){
            if (histoT[j]>maxiheight) maxiheight=histoT[j];
            if (histoRec[j]>maxiheight2) maxiheight2=histoRec[j];
          }
          maxiheight+=10;
          maxiheight2+=10;
          cout << "maxiheight =" << maxiheight << endl;
          
          for (uint j=0;j<activblobs.size();j++){
            if (activblobs[j].size() == activblobs[i].size()){
            vector<uint> composante;
            vector<double> samplesTri(samplesT);
            
            for (it=activblobs[j].begin();it!=activblobs[j].end();it++)
              composante.push_back(*it);
            vector<double> sample(getCaracSample(composante));
            samplesTri.push_back(sample[0]);
            std::sort(samplesTri.begin(),samplesTri.end());
            double r;
            for (r=0;r<(int)samplesTri.size() && pow(samplesTri[r] - sample[0],2) >0.001;r++){}
            assert (r!=samplesTri.size());
            cout << "perc. " << j <<":" << r/(double)samplesTri.size()*100.0 << " " << sample[0] << endl;
            fprintf(fAct,"%f %f\n", sample[0], (double)maxiheight);
            }
            if (activblobs[j].size() == activblobs[i].size()){
              vector<uint> composante;
              vector<double> samplesTri(samplesRec);
            
              for (it=activblobs[j].begin();it!=activblobs[j].end();it++)
                composante.push_back(*it);
              vector<double> sample(getCaracSample(composante));
              samplesTri.push_back(sample[2]);
              std::sort(samplesTri.begin(),samplesTri.end());
              double r;
              for (r=0;r<(int)samplesTri.size() && pow(samplesTri[r] - sample[2],2) >0.001;r++){}
              assert (r!=samplesTri.size());
              cout << "perc. " << ssb->sites[*(activblobs[j].begin())]->label <<":" << r/(double)samplesTri.size()*100.0 << " " << sample[2] << endl;
              fprintf(fActRec,"%f %f\n", sample[2], (double)maxiheight2);
            }
          }
            
//           vector<double> sample(getCaracSample(composante));
//           vector<double> samplesdist2(samples);
//           long double energysample = ssb->getClusterEnergy(composante);
//           cout << "label " << j << " tmoy=" << sample[0] << " ttest =" << sample[1] <<" rec=" << sample[2] <<  endl;
//           cout << "label:" << j << " ttest=" << sample[1];
          fclose(f);fclose(ft);fclose(frec);fclose(fNeg);fclose(fAct);fclose(fActRec);
          cin>>mini;

  }




}
