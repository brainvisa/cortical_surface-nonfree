#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "cluster.h"
#include <aims/mesh/texture.h>
#include <aims/mesh/surfacegen.h>
#include <aims/mesh/surfaceOperation.h>

using namespace aims;
using namespace carto;
using namespace std;

SWC::SWC(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons){
  MinimizationSetup(primal,meshes,lats,lons);
}

vector<int> SWC::getCompConn(vector<uint> &indicesCliques, set<uint> &listeSites){
  vector<int> comp(sites.size());
  uint blob0,blob1;
  Site *s0, *s1;
  int lcomp,nbcomp,aux;
  int label0,label1;
  set<uint>::iterator it;
  for (uint i=0;i<sites.size();i++)
    comp[i]=-1;
  lcomp=0;
  nbcomp=0;
  for (it = listeSites.begin();it != listeSites.end();it++)
    comp[*it]=0;
  for (uint i=0;i<indicesCliques.size();i++)
    if (cliques[indicesCliques[i]].type == SIMILARITY){
      s0 = cliques[indicesCliques[i]].blobs[0];
      s1 = cliques[indicesCliques[i]].blobs[1];
      blob0=s0->index;
      blob1=s1->index;
      comp[blob0]=0;
      comp[blob1]=0;
    }
  for (uint i=0;i<indicesCliques.size();i++){
    if (cliques[indicesCliques[i]].type == SIMILARITY){

      s0 = cliques[indicesCliques[i]].blobs[0];
      s1 = cliques[indicesCliques[i]].blobs[1];
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

vector<set<uint> > SWC::getCompConnVector(vector<int> &comp){
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

long double frec(long double rec){
  if (rec>20.0) return 0.0;
  else if (rec<0.0) return 1.0;
  else return -1.0/20.0*rec+1.0;
}
long double ft(long double t){
  if (t>9.0) return 1.0;
  else if (t<2.0) return 0.0;
  else return 1.0/7.0*t-2.0/7.0;
}

vector<uint> SWC::getCliquesTurnedOn(float temp, vector<uint> &indicesCliques){
  vector<uint> turnedOn;

  for (uint i=0;i<indicesCliques.size();i++){
    if (cliques[indicesCliques[i]].type == SIMILARITY){
//       if(cliques[indicesCliques[i]].blobs[0]->label==cliques[indicesCliques[i]].blobs[1]->label){
      long double tirage = ((long double)UniformRandom() * 1.0);
//         long double qe= exp(-((cliques[indicesCliques[i]].rec)*(temp+1)));
      long double frect = frec(cliques[indicesCliques[i]].rec)*ft(cliques[indicesCliques[i]].blobs[0]->t)*ft(cliques[indicesCliques[i]].blobs[1]->t);
      long double qe = (1-exp(-frect*temp)); //(1-exp(-1.0));
//         if ((cliques[indicesCliques[i]].blobs[0]->t+cliques[indicesCliques[i]].blobs[1]->t)/2.0 < 5.0) qe=0.0;
//         qe += 0.1 * (-1.0/20.0*cliques[indicesCliques[i]].rec+1.0);
//         long double qe = exp(-(cliques[indicesCliques[i]].rec+1.0));
//         cout << qe << " " ;
      if (frect<0.01) {
      }
      else if (frect>0.99) {
        turnedOn.push_back(indicesCliques[i]);
      }
      else if (tirage<qe){
//           cout << cliques[indicesCliques[i]].rec << "/"<< temp << " " << flush;

          turnedOn.push_back(indicesCliques[i]);
        }
//       }
    }
  }
//   cout << " ts:" << turnedOn.size() << "/" <<indicesCliques.size() <<  endl;
  return turnedOn;
}

long double SWC::getCompacite(set<uint> &comp, bool verb){
  set<string> subj;
  uint penal = 0;
  long double compac=0.0,rec=0.0;
  set<uint> auxcliques;
  set<uint>::iterator it;
  float t=0.0;
  for (it=comp.begin();it!=comp.end();it++){
    for (uint i=0;i<cliquesDuSite[*it].size();i++){
      uint aux=cliquesDuSite[*it][i],index0,index1;
      index0=cliques[aux].blobs[0]->index;
      index1=cliques[aux].blobs[1]->index;
      if ((comp.find(index0) != comp.end() && index0 != *it) || (comp.find(index1) != comp.end() && index1 != *it))
        auxcliques.insert(cliquesDuSite[*it][i]);
    }
        
    compac += 1.0-ft(sites[*it]->tValue);
    if (subj.find(sites[*it]->subject)!=subj.end()) penal=penal+1;
    subj.insert(sites[*it]->subject);
  }
  t = (float) compac;
//   cout << "compsize:" << comp.size() << endl ;
//   cout << "compac:" << compac << endl;
  for (it=auxcliques.begin();it!=auxcliques.end();it++)
    rec += frec(cliques[*it].rec); // /10.0+0.9;
  
  compac += -rec;
//   cout << "rec:" << rec<< endl;
  compac += penal*100.0;
  
  if (verb) cout << "[" << t << ";" << -rec << ";" << subj.size() << ";" << penal <<"]";
//   compac *=10.0;
//   compac /= (float)comp.size();
//   compac *= ((float)subj.size()/(float)nbsujets);
//   compac /= (float)penal;
//   cout << "= "<< compac << " - " ;
  
  return compac;

}



void SWC::Run2(){
  float temp=10.0;

  vector<uint> indicesCliques, turnedOn;
  set<uint> indicesSet, listeSites;
  set<uint>::iterator it,it1,it2,it3;
  FILE * f1;   f1 = fopen (recuitpath.data(),"w");
  int ite=0,acc;
  for (uint i=0;i<cliques.size();i++)
    indicesCliques.push_back(i);
  for (uint i=0;i<sites.size();i++)
    listeSites.insert(i);
  

  vector<int> auxcomp(getCompConn(indicesCliques, listeSites));
  vector<set<uint> > ssgraphes;
  vector<set<uint> > allCompConn(getCompConnVector(auxcomp));
  for (uint i= 0;i<allCompConn.size();i++)
    if (allCompConn[i].size()>0)
      ssgraphes.push_back(allCompConn[i]);
  for (uint i= 0;i<ssgraphes.size();i++)
    cout << getCompacite(ssgraphes[i]) << " " ;
  cout << endl;

  for (uint i=0;i<sites.size();i++){
    if (sites[i]->t < 0) sites[i]->label = 0;
    else sites[i]->label=sites[i]->t;
  }
  for (uint i0=0;i0<sites.size();i0++){
    fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
  }
  fprintf(f1, "\n");
  

  
  while (temp>0.001){

    cout << " T=" << temp << " it="<< ite++ << " ssg:" << ssgraphes.size() << " "  <<endl;
    allCompConn.clear();
    vector<set<uint> > allCC;
    for (uint i=0;i<ssgraphes.size();i++){
      cout << ssgraphes[i].size() << "(" << getCompacite(ssgraphes[i],false) << ") ";
      indicesCliques.clear();
      indicesSet.clear();
      for (it1=ssgraphes[i].begin();it1!=ssgraphes[i].end();it1++)
        for (it2=ssgraphes[i].begin();it2!=ssgraphes[i].end();it2++)
          for (uint m=0;m<cliquesDuSite[*it1].size();m++)
            if (cliques[cliquesDuSite[*it1][m]].type==SIMILARITY)
              if ((cliques[cliquesDuSite[*it1][m]].blobs[0]->index==*it1 && cliques[cliquesDuSite[*it1][m]].blobs[1]->index==*it2) || (cliques[cliquesDuSite[*it1][m]].blobs[1]->index==*it1 && cliques[cliquesDuSite[*it1][m]].blobs[0]->index==*it2))
              indicesSet.insert(cliquesDuSite[*it1][m]);
      for (it=indicesSet.begin();it!=indicesSet.end();it++)
        indicesCliques.push_back(*it);

      turnedOn.clear();
      turnedOn = getCliquesTurnedOn(temp,indicesCliques);
      auxcomp.clear();


      auxcomp=getCompConn(turnedOn, ssgraphes[i]);

      vector<set<uint> > compconn(getCompConnVector(auxcomp));
      for (uint j=0;j<compconn.size();j++)
        allCC.push_back(compconn[j]);
    }
    cout << endl;
    cout << allCC.size() << " composantes connexes - " ;
    
    for (uint j=0;j<allCC.size();j++)
      if (allCC[j].size()>1)
        allCompConn.push_back(allCC[j]);
    cout << allCompConn.size() << " composantes de plus de un élement" << endl;
    if (allCompConn.size() == 0){
      temp=temp*0.99;
            continue;
    }
    vector<long double> compacList,compacDist;
    long double som = 0.0;
    for (uint j=0;j<allCompConn.size();j++){
      cout << "compac:";
      long double aux = getCompacite(allCompConn[j]);
      cout << aux;
      compacList.push_back(exp(-aux/(10.0*temp)));
      som += exp(-aux/(10.0*temp));
      cout << "("<< allCompConn[j].size() <<")" << endl;
    }
    long double  som2=0.0;
    for (uint j=0;j<allCompConn.size();j++){
      som2+=compacList[j]/som;
      cout << som2 << ";" ;
      compacDist.push_back(som2);
    }
    cout << endl;
    
    
//     ssgraphes.clear();
//     for (uint i= 0;i<allCompConn.size();i++)
//       if (allCompConn[i].size()>0)
//         ssgraphes.push_back(allCompConn[i]);
//
      
//     for (uint i0=0;i0<sites.size();i0++){
//       int cc=-1;
//       for (uint j=0;cc==-1 && j<allCompConn.size();j++)
//         if (allCompConn[j].find(i0)!=allCompConn[j].end())
//           cc = j;
//       if (cc==-1) cc=0; else cc++;
//       fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, cc+1);
//     }

    
    
//     for (uint j=0;j<allCompConn.size();j++)
//       if (allCompConn[j].size()>0)
//         cout << allCompConn[j].size() << " ";
   
    uint tir = (uint)((double)UniformRandom() * allCompConn.size());
    cout << "tirage:" << tir <<" @ " ;
    indicesSet.clear();
    int activeSG=-1;
    for (uint i=0;i<ssgraphes.size()&&activeSG==-1;i++){
      if (ssgraphes[i].find(*(allCompConn[tir].begin()))!=ssgraphes[i].end())
        activeSG=i;
    }
    cout << "activeSG:" << activeSG << endl;

    ASSERT(activeSG!=-1);
    for (it = allCompConn[tir].begin() ; it != allCompConn[tir].end() ; it ++){
      ASSERT(ssgraphes[activeSG].find(*it) != ssgraphes[activeSG].end());
    }
    set<uint> graphesAdjacents;
    set<uint> sitesAdjacents;
    uint auxind;
    for (it=allCompConn[tir].begin();it!=allCompConn[tir].end();it++)
      for (uint m=0;m<cliquesDuSite[*it].size();m++){
        if (cliques[cliquesDuSite[*it][m]].type == SIMILARITY){
          if (cliques[cliquesDuSite[*it][m]].blobs[0]->index == *it)
            auxind=1;
          else auxind=0;
          if (allCompConn[tir].find(cliques[cliquesDuSite[*it][m]].blobs[auxind]->index) == allCompConn[tir].end())
            sitesAdjacents.insert(cliques[cliquesDuSite[*it][m]].blobs[auxind]->index);
          for (uint j=0;j<ssgraphes.size();j++)
            if (j!=(uint)activeSG)
              if (ssgraphes[j].find(cliques[cliquesDuSite[*it][m]].blobs[auxind]->index)!=ssgraphes[j].end())
                graphesAdjacents.insert(j);
        }
      }
//     cout << graphesAdjacents.size() << endl;
    
    
    long double compacSSGraph=0.0, compacCompConn=0.0;
    vector<long double> compacDistrib;
    vector<long double> compacListe;
    compacSSGraph = getCompacite(ssgraphes[activeSG]);
    cout << compacSSGraph <<" (" << ssgraphes[activeSG].size() << ") / ";
    compacCompConn = getCompacite(allCompConn[tir]);
    cout << compacCompConn << endl;
    long double somme=0.0;
    compacDistrib.push_back(exp(-compacSSGraph/temp));
    compacListe.push_back(compacSSGraph);
    somme += exp(-compacSSGraph/temp);
//     cout << "sum:" <<somme << endl;
    compacDistrib.push_back(exp(-compacCompConn/temp));
    compacListe.push_back(compacCompConn);
    somme += exp(-compacCompConn/temp);
//     cout << "sum:" << somme << endl;

    cout << "graphes adj:" << graphesAdjacents.size() << " " << "sites adj:" << sitesAdjacents.size() << endl;
    for (it=graphesAdjacents.begin();it!=graphesAdjacents.end();it++){
      set<uint> tempgraph(ssgraphes[*it]);
      for (it1=allCompConn[tir].begin();it1!=allCompConn[tir].end();it1++)
        tempgraph.insert(*it1);
      float compac = getCompacite(tempgraph);
      cout << compac << "(" << tempgraph.size() << ") " ;
      compacDistrib.push_back(exp(-compac/temp));
      compacListe.push_back(compac);
      somme += exp(-compac/temp);
//       cout << somme << endl;
    }
    cout << " - " ;
    for (it=sitesAdjacents.begin();it!=sitesAdjacents.end();it++){
      set<uint> tempgraph(allCompConn[tir]);
      tempgraph.insert(*it);
      float compac = getCompacite(tempgraph);
      cout << compac << "(" << tempgraph.size() << ") " ;
      compacDistrib.push_back(exp(-compac/temp));
      compacListe.push_back(compac);
      somme += exp(-compac/temp);
    }
    cout << endl;
// //     cout << "somme : "<< somme << " " ;
//     long double somme2=0.0;
// 
// //     cout << somme << endl;
// 
//     for (uint i=0;i<compacDistrib.size();i++){
//       somme2 += compacDistrib[i]/somme;
//       compacDistrib[i] = somme2;
//       cout << somme2 << ";" ;
//     }
//
    cout << "=" ;


    
    float tirage = (double)UniformRandom();
// 
//     
    int graphAdjFus=0, indiceSite=0;
    if (false){
    for (acc=0;(uint)acc<compacDistrib.size() && compacDistrib[acc]<tirage;acc++) {}
    }
    else{
      cout << "mode ICM:" << endl;
//       temp=0.0;
      uint mini=0;
      for (uint i=0;i<compacListe.size();i++){
//         cout << compacListe[i] << " ";
        if (compacListe[i]<compacListe[mini]) mini=i;
      }
      acc=mini;
    }
//     cout << "acc:" << acc << endl;
    if (acc == 0) cout << "no fusion"<< endl;
    else if (acc == 1) cout << "new subgraph" << endl;
    else if (acc > 1 && (uint)acc < graphesAdjacents.size()+2){
      it = graphesAdjacents.begin();
      for (uint i=0;i<(uint)(acc-2);i++,it++) {}
      graphAdjFus = *it;
      cout << "fusion avec sgraph d'indice " << graphAdjFus << endl;
    }
    else if ((uint)acc>= graphesAdjacents.size()+2){
      it=sitesAdjacents.begin();
      for (uint i=0;i<acc-graphesAdjacents.size()-2;i++)
        it++;
      indiceSite = *it;
      cout << "new subgraph avec rajout du site d'indice " << indiceSite << endl;

    }
    if (acc == 1) {
      //CREATION NEW SUBGRAPH
      set<uint> old;
      for (it = ssgraphes[activeSG].begin();it != ssgraphes[activeSG].end();it++)
        if (allCompConn[tir].find(*it) == allCompConn[tir].end()) old.insert(*it);
      ASSERT(allCompConn[tir].size() + old.size() == ssgraphes[activeSG].size());
//       cout << "old.size:"<< old.size() << " allComp[tir].size:" << allCompConn[tir].size() <<  endl;
      if (old.size()!=0){
      ssgraphes[activeSG].clear();
      for (it = old.begin();it!=old.end();it++)
        ssgraphes[activeSG].insert(*it);
      ssgraphes.push_back(allCompConn[tir]);
      }
//       cout << "fin new sub" << endl;
    }
    else if (acc >1 && (uint)acc < graphesAdjacents.size()+2){
      // FUSION AVEC LE SSGRAPH adjacents d'indice acc-2
      for (it = allCompConn[tir].begin();it!=allCompConn[tir].end();it++)
        ssgraphes[graphAdjFus].insert(*it);
      set<uint> old;
      for (it = ssgraphes[activeSG].begin();it != ssgraphes[activeSG].end();it++)
        if (allCompConn[tir].find(*it) == allCompConn[tir].end()) old.insert(*it);
      ASSERT(allCompConn[tir].size() + old.size() == ssgraphes[activeSG].size());

//       cout << "old.size:"<< old.size() << " allComp[tir].size:" << allCompConn[tir].size() <<  endl;

      if (old.size()!=0){
      ssgraphes[activeSG].clear();
      for (it = old.begin();it!=old.end();it++)
        ssgraphes[activeSG].insert(*it);
      }
      else{
//         cout << "erase" << ssgraphes.size() << endl;
        ssgraphes.erase(ssgraphes.begin()+activeSG);
        cout << ssgraphes.size() << endl;
      }
    }
    else if ((uint)acc >= graphesAdjacents.size() + 2){
      // FUSION AVEC UN SITE ADJACENT D'INDICE acc - graphesAdjacents.size() -2
      graphAdjFus=-1;
      for (uint i=0;i<ssgraphes.size()&&graphAdjFus==-1;i++)
        if (ssgraphes[i].find(indiceSite)!=ssgraphes[i].end())
          graphAdjFus=i;
      ASSERT(graphAdjFus != -1);
      set<uint> old;
      for (it = ssgraphes[activeSG].begin();it != ssgraphes[activeSG].end();it++)
        if (allCompConn[tir].find(*it) == allCompConn[tir].end()) old.insert(*it);
      ASSERT(allCompConn[tir].size() + old.size() == ssgraphes[activeSG].size());
//       cout << "old.size:"<< old.size() << " allComp[tir].size:" << allCompConn[tir].size() << " " << activeSG << " " << graphAdjFus <<  endl;
      if (old.size()!=0){
        ssgraphes[activeSG].clear();
        for (it = old.begin();it!=old.end();it++)
          ssgraphes[activeSG].insert(*it);
        ssgraphes.push_back(allCompConn[tir]);
        ssgraphes[ssgraphes.size()-1].insert(indiceSite);
      }
      else {
        ssgraphes[activeSG].insert(indiceSite);

      }
      // FIN CREATION NOUVEAU SOUS GRAPHE AVEC LA COMPOSANTE CONNEXE TIREE AU SORT
      // RAJOUT SITE ADJACENT ET SUPPRESSION DU SITE DU SOUS-GRAPHE AUQUEL IL APPARTENAIT

      old.clear();
      for (it = ssgraphes[graphAdjFus].begin();it != ssgraphes[graphAdjFus].end();it++)
        if (*it!=(uint)indiceSite) old.insert(*it);
      ASSERT(1 + old.size() == ssgraphes[graphAdjFus].size());

//       cout << "old.size:"<< old.size() << endl ; //" allComp[tir].size:" << allCompConn[tir].size() <<  endl;


      if (old.size()!=0){
        ssgraphes[graphAdjFus].clear();
        for (it = old.begin();it!=old.end();it++)
          ssgraphes[graphAdjFus].insert(*it);
      }
      else{
//         cout << "erase" << ssgraphes.size() << endl;
        ssgraphes.erase(ssgraphes.begin()+graphAdjFus);
        cout << ssgraphes.size() << endl;
      }
      uint chksum3=0;
      for (uint i=0;i<ssgraphes.size();i++)
        chksum3+=ssgraphes[i].size();
        ASSERT(chksum3==sites.size());
      cout << chksum3 << "/" << sites.size() << endl;
    }
// 
    uint chksum=0;
    for (uint i=0;i<sites.size();i++)
      sites[i]->label=0;
    uint lab=1;
    for (uint i=0;i<ssgraphes.size();i++){
      chksum+=ssgraphes[i].size();
//       cout << "ssg.size:" << ssgraphes[i].size() << endl;
      float compac = getCompacite(ssgraphes[i],false);
//       cout << "compacite :" << compac << endl;
      if (compac > 480.0) compac=480.0;
      if (ssgraphes[i].size()>1){
        for (it=ssgraphes[i].begin();it!=ssgraphes[i].end();it++){
          sites[*it]->label= lab; //(int) (compac+20.0);
        }
        lab++;
      }
        
//       }

    }
    ASSERT(chksum==sites.size() || (cout << chksum << "/" << sites.size() << endl && false));


    for (uint i0=0;i0<sites.size();i0++){
      fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
    }
    fprintf(f1, "\n");
    temp=temp*0.99;
    cout << "====================================================================================================================" << endl;
    cout << "====================================================================================================================" << endl;
    cin >> ite;
  }

  for (uint i=0;i<sites.size();i++)
    sites[i]->label=0;
  for (uint i=0;i<ssgraphes.size();i++){
    float compac = getCompacite(ssgraphes[i],false);
    if (compac > 100.0) compac=100.0;
    if (ssgraphes[i].size()>1){
      for (it=ssgraphes[i].begin();it!=ssgraphes[i].end();it++){
        sites[*it]->label= (int) (compac+20.0);
      }
    }
    else{
      for (it=ssgraphes[i].begin();it!=ssgraphes[i].end();it++){
        sites[*it]->label= 120.0;
      }
    }
      
  }
  for (uint i0=0;i0<sites.size();i0++){
    fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
  }
  fprintf(f1, "\n");
  
  fclose(f1);

}

void SWC::Run(){
  float temp=10.0;
  vector<uint> indicesCliques,labelscpt;
//   long double tirage, qe;
  set<uint> nodes;
  uint cpt=0,ite=0;

  FILE * f1;   f1 = fopen (recuitpath.data(),"w");

  while (temp>0.1){

    cout << " T=" << temp << " it="<< ite++ << " " ;

    vector<uint> indicesCliques;
    // attention IL FAUT INITIALISER INDICESCLIQUES LA PREMIERE FOIS EN DEHORS DE LA BOUCLE
    indicesCliques = getCliquesTurnedOn(temp,indicesCliques);
    // attention à initialiser listeSites;
    set<uint> listeSites;
    vector<int> comp(getCompConn(indicesCliques, listeSites));

      
    for (uint i=0;i<sites.size();i++)
      if (comp[i]<0) sites[i]->label=0;
      else  sites[i]->label = comp[i];
      
    cpt=0;
    cout << sites.size() << " ";
    for (uint i0=0;cpt<sites.size();i0++){
      vector<uint> conn;
      set<string> subj;
      float moy=0;
      for (uint i=0;i<comp.size();i++)
        if (comp[i]==labels[i0] || (i0==0 && comp[i]==-1)){
        conn.push_back(i);
        subj.insert(sites[i]->subject);
        moy+=sites[i]->t;
        }
        moy/=conn.size();
        cpt+=conn.size();
        moy*=subj.size();
      long double tirage = ((long double)UniformRandom() * 1.0);
      long double qe= exp(-(1/moy));
      if (tirage>qe) //remise à zéro
        for (uint i=0;i<conn.size();i++)
          comp[conn[i]] = 0;
        cout << "L" << labels[i0] << ":" << conn.size() << ";" << moy << ";" << tirage << ">" << qe << "-";
    }

          
    for (uint i=0;i<sites.size();i++)
      if (comp[i]<0) sites[i]->label=0;
    else  sites[i]->label = comp[i];

    for (uint i0=0;i0<sites.size();i0++){
      fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
    }
    fprintf(f1, "\n");
    temp=temp*0.99; 
    cout << endl;
  }

  
  for (uint i=0;i<cliques.size();i++)
    cliques[i].updateLabelsCount();
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


