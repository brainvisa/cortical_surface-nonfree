#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "cluster.h"
#include <aims/mesh/texture.h>

using namespace aims;
using namespace carto;
using namespace std;

SWC::SWC(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon){
  MinimizationSetup(primal,mesh,lat,lon);
}

vector<int> SWC::getCompConn(vector<uint> &indicesCliques){
  vector<int> comp(sites.size());
  for (uint i=0;i<sites.size();i++)
    comp[i]=-1;
  int lcomp=0;
  int nbcomp=0;
  int label0,label1;
  for (uint i=0;i<indicesCliques.size();i++)
    if (cliques[indicesCliques[i]].type == SIMILARITY){
      Site *b0,*b1;
      b0 = cliques[indicesCliques[i]].blobs[0];
      b1 = cliques[indicesCliques[i]].blobs[1];
      uint blob0=b0->index, blob1=b1->index;
      comp[blob0]=0;
      comp[blob1]=0;
    }
  
  for (uint i=0;i<indicesCliques.size();i++){
    if (cliques[indicesCliques[i]].type == SIMILARITY){
      Site *b0,*b1;
      b0 = cliques[indicesCliques[i]].blobs[0];
      b1 = cliques[indicesCliques[i]].blobs[1];
      uint blob0=b0->index, blob1=b1->index;
//       if (!(comp[blob0]>0 && comp[blob1]>0)) cout << indicesCliques[i] << "/" << indicesCliques.size() << ":" << comp[blob0] << "-" << comp[blob1] << ":";

      if (comp[blob0] == 0){
        if (comp[blob1] == 0){
          lcomp++;
          nbcomp++;
          comp[blob0] = lcomp;
          comp[blob1] = lcomp;

          cout << "a"<< nbcomp << "-" ;
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
        int aux;
        if (label0>label1) {
          aux=label0;
          label0=label1;
          label1=aux;
        }

        cout << nbcomp << ";" << lcomp << "(" << label1 << "-" << label0 << ")-";
        for (uint i0=0;i0<comp.size();i0++){
//           if (cliques[indicesCliques[i0]].type == SIMILARITY){
          if (comp[i0] == label1)
            comp[i0] = label0;
//           if (comp[i0] == label1)
//             comp[i0] = label0;
          }
        for (uint i0=0;i0<comp.size();i0++){
//           if (cliques[indicesCliques[i0]].type == SIMILARITY){
          if (comp[i0] == lcomp){
            comp[i0] = label1;
            cout << "o" ;
          }
//           if (comp[i0] == lcomp)
//             comp[i0] = label1;
          }
        lcomp--;
        nbcomp--;
        
      }
    
    if (!(comp[blob0]==comp[blob1])){
      cout << "(" ;
    for (int i0=0;i0<100;i0++){
//       cout << i0 << ":";
      int cnt=0;
      for (uint j=0;j<comp.size();j++)
        if (comp[j]==i0) cnt++;
      if (cnt>0) cout << i0 << ":" << cnt << "-";
    }
    cout << ")" ;
    }
    }
  }
  set<uint> nodes;
  cout << nbcomp << "!" << endl;
  for (uint i0=0;i0<indicesCliques.size();i0++)
    if (cliques[indicesCliques[i0]].type == SIMILARITY){
      if (comp[cliques[indicesCliques[i0]].blobs[0]->index] > nbcomp || comp[cliques[indicesCliques[i0]].blobs[1]->index] > nbcomp) cout << "!!!!!!!!!!!!" << endl;
      nodes.insert(cliques[indicesCliques[i0]].blobs[0]->index);
      nodes.insert(cliques[indicesCliques[i0]].blobs[1]->index);
    }
  uint cpt=0,cpt2;
  vector<uint> labelscpt;
  for (uint i0=0;i0<100;i0++)
    labelscpt.push_back(0);
  for (uint i0=0;cpt<nodes.size();i0++){
    cout << i0 << ":";
    cpt2=0;
    for (uint i=0;i<comp.size();i++)
      if ((uint)comp[i]==i0) {
      cpt++;cpt2++;
      cout << sites[i]->index << "-" <<sites[i]->tValue << "-" << sites[i]->node  << "\t\t;\t";
      }
      labelscpt[cpt2]++;
    cout << endl;
  }
  cout << endl;
  for (uint i0=0;i0<100;i0++)
    cout << labelscpt[i0] << ";";
  cout << endl;
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

void SWC::getListeTriangles(){
  vector<AimsVector<uint,3> > triangles;
  vector<set<uint> > poly(sites.size());
  for (uint i=0;i<cliques.size();i++){
    if (cliques[i].type == SIMILARITY){
      poly[cliques[i].blobs[0]->index].insert(cliques[i].blobs[1]->index);
      poly[cliques[i].blobs[1]->index].insert(cliques[i].blobs[0]->index);
    }
  }
  
  set<uint>::iterator it,jt,kt;
  uint k=0;
  set<uint> vertices;
  vector<uint> corres(sites.size());
  map<int,int> corres2;
  for (uint i=0;i<poly.size();i++){
    for (it=poly[i].begin();it!=poly[i].end();it++)
      for (jt=poly[i].begin();jt!=poly[i].end();jt++){
        if (poly[*it].find(*jt)!=poly[*it].end()){
          if (i<*it && *it<*jt){
            AimsVector<uint,3> v;
            vertices.insert(i);
            vertices.insert(*it);
            vertices.insert(*jt);
            for (k=0,kt=vertices.begin();kt!=vertices.end() && *kt!=i;kt++,k++){}
            v[0]=k;
            corres[i] = k;
            corres2[k]= i;
            for (k=0,kt=vertices.begin();kt!=vertices.end() && *kt!=*it;kt++,k++){}
            v[1]=k;
            corres[*it]=k;
            corres2[k]=*it;
            for (k=0,kt=vertices.begin();kt!=vertices.end() && *kt!=*jt;kt++,k++){}
            v[2]=k;
            corres[*jt]=k;
            corres2[k]=*jt;
            triangles.push_back(v);
          }
        }
      }
  }
  
//   for (uint i=0;i<triangles.size();i++)
//     cout << triangles[i][0] <<" " << triangles[i][1] << " " << triangles[i][2] << "-";
  cout << "nombre de vertex :" << vertices.size() << endl;

  cout << "nombre de triangles :" << triangles.size() << endl;
  Reader<AimsSurfaceTriangle> r("/home/grg/V10_sphere.mesh");
  AimsSurfaceTriangle sphere;
  r.read(sphere);
  AimsSurfaceTriangle mesh;
  for (kt=vertices.begin();kt!=vertices.end();kt++){
    AimsVector<float,3> p(20.0*cos(sites[*kt]->gravitycenter[1])*cos(sites[*kt]->gravitycenter[0]),20.0*sin(sites[*kt]->gravitycenter[1])*cos(sites[*kt]->gravitycenter[0]),20.0*sin(sites[*kt]->gravitycenter[0]));
    p = sphere[0].vertex()[sites[*kt]->gravitycenter[2]];
    AimsVector<float,3> q(p);
    q /= p.norm();
    q *= sites[*kt]->tValue * 10.0;
    p += q;
    mesh[0].vertex().push_back(p);
  }
  for (uint i=0;i<triangles.size();i++)
    mesh[0].polygon().push_back(AimsVector<uint,3>(triangles[i][0],triangles[i][1],triangles[i][2]));
  mesh[0].updateNormals();
  Writer<AimsSurfaceTriangle> w("/home/grg/testNewMesh.mesh");
  w.write(mesh);
  vector<uint> clik;
  for (uint i=0;i<cliques.size();i++)
    if (cliques[i].type == SIMILARITY)
      clik.push_back(i);
  map<int,int>::iterator itt;
  TimeTexture<short> tex(1,0);
  for (uint i=0;i<mesh[0].vertex().size();i++)
    tex[0].push_back(0);

  vector<int> cc = getCompConn(clik);
  for (uint i=0;i<sites.size();i++)
    sites[i]->label = 0;
  for (uint i=0;i<sites.size();i++){
    if (cc[i]==-1) cc[i] = 0;
    if (cc[i]<20) sites[i]->label = cc[i];
    tex[0].item(corres[i])=cc[i];
  }
  for (uint i=0;i<cliques.size();i++)
    cliques[i].updateLabelsCount();
  
  Writer<TimeTexture<short> > w2("/home/grg/testNewTex.tex");
  w2.write(tex);

  FILE * f1;   f1 = fopen ("/home/grg/recuit.txt","w");
  for (uint i0=0;i0<sites.size();i0++){
    fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
  }
  fprintf(f1, "\n");
  fclose(f1);
  
}

void SWC::Run(){
  getListeTriangles();

  
//   FILE * f1;   f1 = fopen ("/home/grg/recuit.txt","w");
//   float  temp= 10.0;
//   uint ite=0;
//   while (temp>0.0000001){
//     cout << " T=" << temp << " it="<< ite++ << " " << flush ;
// 
//     vector<uint> turnedOn(getCliquesTurnedOn(temp));
//     map<int,int> comp(getCompConn(turnedOn));
//     uint nbcc=0;
//     
//     map<int,int>::iterator it;
//     for (it=comp.begin();it!=comp.end();it++)
//       if (nbcc<it->second) nbcc=it->second;
//     uint tirage = (uint)((long double)UniformRandom() * nbcc);
//     
//     for (it=comp.begin();it!=comp.end();it++){
//       sites[it->first]->label = it->second;
//     }
//     for (uint i0=0;i0<sites.size();i0++){
//       fprintf(f1, "%s %d %d %d-", sites[i0]->subject.data(), sites[i0]->index, sites[i0]->graph_index, sites[i0]->label);
//     }
//     fprintf(f1, "\n");
// 
//     temp *= 0.99;
//     cout << endl;
//   }
//   fclose(f1);
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


