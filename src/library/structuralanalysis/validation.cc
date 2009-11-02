#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include <cortical_surface/structuralanalysis/validation.h>

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
  double c =(double) ceil(sqrt(2.F*samplesdist.size())); uint k=r+c; double b2 = 1.0/0.05;
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
    tmoy += ssb->sites[composante[k]]->tValue;
    compac += ssb->sites[composante[k]]->t;
//     compaccent += ssb->sites[composante[k]]->t2;
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
  vector<vector<float> > results(activblobs.size());
  uint activindex=0;
  set<uint> processedsizes;
  cout << "ABG=" << activblobsglobal.size() << endl;
  
  
  for (uint i=0;i<activblobs.size();i++) {
    if (processedsizes.find(activblobs[i].size())==processedsizes.end()){
          processedsizes.insert(activblobs[i].size());
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

          cout << "restent :"<< ssb->sites.size()-forbidden.size() << endl;

          for (uint j=0;j<10000;j++){
//                 cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << j << flush;
                //création du cluster
                uint startsite;
                
                set<uint> dejapris = set<uint>(forbidden);
                do{ startsite = (float)UniformRandom() * ssb->sites.size(); }
                while(tirage[startsite]<1 || dejapris.find(startsite)!=dejapris.end());
                vector<uint> composante;
                composante.push_back(startsite);
                          
                          
                dejapris.insert(startsite);
                          

                set<uint> voisins;
                while (composante.size()<activblobs[i].size()){
                      
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
                  assert(voisins.size()!=0);

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
//           FILE *f,*ft,*frec,*fNeg,*fAct,*fActRec;
//           f=fopen("/volatile/operto/histo.txt","w");
//           ft=fopen("/volatile/operto/histo_t.txt","w");
//           frec=fopen("/volatile/operto/histo_rec.txt","w");
//           fNeg=fopen("/volatile/operto/histo_neg.txt","w");
//           fAct=fopen("/volatile/operto/barres_activ.txt","w");
//           fActRec=fopen("/volatile/operto/barres_activrec.txt","w");
          cout << endl<< "nombre de neg:"<< samplesNeg.size()<< endl;
          uint histosize=40;

          
          cout << "histoT" << endl;
          vector<int> histoT(creerHisto(samplesT,histosize,&mini,&step));
//           printHisto(histoT,mini,step,VERTICAL,ft);

          cout << "histoRec" << endl;
          vector<int> histoRec(creerHisto(samplesRec,histosize,&mini,&step));
//           printHisto(histoRec,mini,step,VERTICAL,frec);

          cout << "histoTneg" << endl;
          vector<int> histoNeg(creerHisto(samplesNeg,histosize,&mini,&step));
//           printHisto(histoNeg,mini,step,VERTICAL,fNeg);
          
          cout << "histo" << endl;
          vector<int> histo(creerHisto(samples,histosize,&mini,&step));
//           printHisto(histo,mini,step,VERTICAL,f);

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
                results[activindex].push_back(ssb->sites[*(activblobs[j].begin())]->label);
                results[activindex].push_back(r/(double)samplesTri.size()*100.0);
//                 fprintf(fAct,"%f %f\n", sample[0], (double)maxiheight);
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
                results[activindex].push_back(r/(double)samplesTri.size()*100.0);
                activindex++;
//                 fprintf(fActRec,"%f %f\n", sample[2], (double)maxiheight2);
            }
          }
            
//           vector<double> sample(getCaracSample(composante));
//           vector<double> samplesdist2(samples);
//           long double energysample = ssb->getClusterEnergy(composante);
//           cout << "label " << j << " tmoy=" << sample[0] << " ttest =" << sample[1] <<" rec=" << sample[2] <<  endl;
//           cout << "label:" << j << " ttest=" << sample[1];
//           fclose(f);fclose(ft);fclose(frec);fclose(fNeg);fclose(fAct);fclose(fActRec);
//           cin>>mini;
    }

  }
  cout << "fin " << endl;
  for (uint i=0;i<ssb->sites.size();i++){
    ssb->sites[*it]->t_rankperc = 0.0;
    ssb->sites[*it]->sim_rankperc = 0.0;
    ssb->sites[*it]->significance = 0.0;
  }
  for (uint i=0;i<results.size();i++){
      cerr << i << ":" << flush;
      cerr << "label " << results[i][0] << ": t-perc=" << results[i][1] << " sim-perc=" << results[i][2] << endl;
      uint j;
      for (j=0;ssb->sites[*(activblobs[j].begin())]->label != results[i][0] && j<activblobs.size();j++){}
      for (it=activblobs[j].begin();it!=activblobs[j].end();it++){
        ssb->sites[*it]->t_rankperc = results[i][1];
        ssb->sites[*it]->sim_rankperc = results[i][2];
        ssb->sites[*it]->significance = (results[i][1]+results[i][2])/2.0;
      }
        
      
  }



}
