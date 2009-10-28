#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "minimization.h"

using namespace aims;
using namespace carto;
using namespace std;

float max(float a, float b){
  if (a>b) return a; else return b;
  return a;
}
float min(float a, float b){
  if (a<b) return a; else return b;
  return a;
}

struct ltstr_vec
{
  bool operator()(const vector<uint> p1, const vector<uint> p2)
  {
    assert(p1.size()==p2.size());
    for (uint i=0;i<p1.size();i++){
      if (p1[p1.size()-1-i]<p2[p1.size()-1-i])
        return true;
      else if (p1[p1.size()-1-i]>p2[p1.size()-1-i])
        return false;
    }
  }
};


void SurfaceBased_StructuralAnalysis::MinimizationSetup(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons){}

void SurfaceBased_StructuralAnalysis::MinimizationSetup(Graph &primal){
  cout << "Construction du vecteur de sites ..." << flush;

  sites = ConstruireSites(primal); 
  cout << "done (" << sites.size() << " sites)" << endl;

  set<string> subjects;

  cout << endl << "  done" << endl;
  for (uint i=0;i<sites.size();i++)
    subjects.insert(sites[i]->subject);
  nbsujets = subjects.size();
  cout << "Construction des cliques ... " << flush;
  cliques = ConstruireCliquesLastChance(sites,cliquesDuSite); 
  
  
  uint nb_cl_sim=0, nb_cl_dd=0, nb_cl_intraps=0, nb_cl_lower=0;
  for (uint i=0;i<cliques.size();i++){
    if (cliques[i].type == SIMILARITY) nb_cl_sim++;
    else if (cliques[i].type == DATADRIVEN) nb_cl_dd++;
    else if (cliques[i].type == BESTLOWERSCALE) nb_cl_lower++;
    else if (cliques[i].type == INTRAPRIMALSKETCH) nb_cl_intraps++;
  }
  cout << " done (" << nb_cl_sim << " cliques de similarité ; " << nb_cl_dd << " cliques datadriven ; " << nb_cl_lower << " cliques lower ; " << nb_cl_intraps << " cliques intraps ; " << cliques.size() << " cliques en tout)" << endl;
 
  uint i=0;
 
  for (uint it=0;it<3;it++)
    for (float zonelat = -20;zonelat<180.0;zonelat += 36.0)
      for ( float zonelon = -20;zonelon<360.0;zonelon += 72.0){
  
        pair<Point2df, Point2df> zone;
        zone.first = Point2df(max(zonelat,0.0), max(zonelon, 0.0));
        zone.second = Point2df(min(zonelat + 76.0, 180.0), min(zonelon + 112.0,360.0));
//         cout << i << " " << zone.first[0] << ";" << zone.first[1] << " " << zone.second[0] << ";" << zone.second[1] << endl;
        labelsZones.push_back(zone);
        i++;
      }
  cout << labelsZones.size() << " zones" << endl;
  for (i=0;i<labelsZones.size()+1;i++)
    labels.push_back(i);
  vector<int> zonescount, labelscount;
  for (i=0;i<labelsZones.size()+1;i++){
    zonescount.push_back(0);
    labelscount.push_back(0);
  }
  zonesListesBlobs = vector<set<uint> >(labelsZones.size()+1);
  listeZones = vector<set<uint> > (sites.size());
  for (uint j=0;j<sites.size();j++){
    uint count=0;
    listeZones[j].insert(0);
    zonesListesBlobs[0].insert(j);
    for (uint k=0;k<labelsZones.size();k++){
      Point3df bbmin1 = sites[j]->boundingbox_min, bbmax1 = sites[j]->boundingbox_max;
      uint no_overlap=0;
      getOverlap(bbmin1, bbmax1, Point3df(labelsZones[k].first[0],labelsZones[k].first[1],0.0), Point3df(labelsZones[k].second[0],labelsZones[k].second[1] ,0.0), &no_overlap);
      if (no_overlap == 0){
        zonescount[k+1]++;
        zonesListesBlobs[k+1].insert(j);
        listeZones[j].insert(k+1);
        count++;
      }
    }
    labelscount[count]++;
    assert(listeZones[j].size()!=0);
  }
  for (i=0;i<labelsZones.size()+1;i++)
    cout << zonescount[i] << " ";
  cout << endl;
  for (i=0;i<labelsZones.size()+1;i++)
    cout << labelscount[i] << " ";
  cout << endl;




    
}

SurfaceBased_StructuralAnalysis::SurfaceBased_StructuralAnalysis(Graph &primal){
  MinimizationSetup(primal);
}

SurfaceBased_StructuralAnalysis::SurfaceBased_StructuralAnalysis(Graph &primal, map<string, AimsSurfaceTriangle> &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons){
  MinimizationSetup(primal,meshes,lats,lons);

//   Initialization();
//   for (uint k=0;k<cliques.size();k++){
//     cliques[k].updateLabelsCount();
//     cliques[k].computeEnergy(true,nbsujets);
//   }
// 
//   for (uint i=0;i<cliques.size();i++)
//     if (cliques[i].type == INTRAPRIMALSKETCH)
//       ipscliques.push_back(i);
// 
//   cout << ipscliques.size() << " cliques intraps" << endl;
// 
//   energy = getTotalEnergy();
// 
//   cout << "energie initiale : " << energy << endl;
//   SummaryLabels();
}


// void SurfaceBased_StructuralAnalysis::Initialization(){
//   vector<uint> indices;
//   for (uint i=0;i<cliques.size();i++)
//     indices.push_back(i);
//   map<int,int> comp(getCompConn(indices));
// 
//   for (uint i=0;i<sites.size();i++)
//     sites[i]->label = comp[i];
// }

void SurfaceBased_StructuralAnalysis::Initialization(){
  energy=0.0;
  cout << "Initialisation :" << endl;
//   for (uint i=0;i<sites.size();i++)
//     sites[i]->label = labels[0];
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
  SummaryLabels();
  
}

long double SurfaceBased_StructuralAnalysis::getClusterEnergy(vector<uint> &composante){
  vector<int> old;
  for (uint i=0;i<sites.size();i++){
    old.push_back(sites[i]->label);
    sites[i]->label=0;
  }
  for (uint i=0;i<composante.size();i++)
    sites[composante[i]]->label = 1;
  
  long double energy = getTotalEnergy();
  for (uint i=0;i<sites.size();i++)
    sites[i]->label = old[i];
  return energy;

}

long double SurfaceBased_StructuralAnalysis::getLabelEnergy(int label, int type){
  bool test = true;
  long double energydd=0.0, energysim=0.0, energy=0.0;
  uint nclsim=0,nbips=0;
  vector<int> bysub(nbsujets);
  
  for (uint i=0;i<cliques.size();i++){
    cliques[i].updateLabelsCount();
    cliques[i].computeEnergy(true,nbsujets);
    if (type==UNKNOWN || cliques[i].type == type){
      test = true;
      for (uint j=0;j<cliques[i].blobs.size() && test == true;j++)
        if (cliques[i].blobs[j]->label != label ) test = false;
      if (test) { 
          if (cliques[i].type==DATADRIVEN) {energydd += cliques[i].energie; }
          else if (cliques[i].type==SIMILARITY)  { nclsim++; energysim += cliques[i].energie;
          }

//         if (false && cliques[i].type==SIMILARITY) cout << "(("<<cliques[i].rec << "(" << cliques[i].blobs[0]->label << "(" <<  cliques[i].blobs[0]->node << ")" << "-" << cliques[i].blobs[1]->label <<"(" <<  cliques[i].blobs[1]->node << ")" << ")=>" << cliques[i].energie << ")) ";
      }
//       if (cliques[i].type == INTRAPRIMALSKETCH && cliques[i].labelscount[label] > 1 && label != 0) {energy += (cliques[i].labelscount[label]-1)*Clique::getIntraPSWeight();}
    }
  }
    uint i,j;
    for (uint n=0;n<ipscliques.size();n++){
      bysub[n] = cliques[ipscliques[n]].labelscount[label];
//       cout << bysub[n] << " " ;
    }
//     cout << endl;
    uint nb=0;
    for (i=0;i<nbsujets-1;i++){
      for (j=i+1;j<nbsujets;j++){
        nb += bysub[i]*bysub[j];
      }
    }
    uint nb2=0;
    for (i=0;i<nbsujets;i++)
      if (bysub[i]>1) nb2 += bysub[i]-1;
    
//     cout << nclsim << "|"<< energydd << " " << energysim << "|";
    nbips += nb;
// cout << Clique::intrapsweight*(nbips-nclsim) << " ";
    ASSERT(nbips>=nclsim || (cout << nbips << ">=" << nclsim << endl && false));
//   energy += Clique::intrapsweight*(nbips-nclsim);
  energy += Clique::intrapsweight*nb2*nbsujets;
  energy += energydd;
  energy += energysim;

  return energy;
}

long double SurfaceBased_StructuralAnalysis::getTypeEnergy(int type){ // RETOURNE L'ENERGIE PAR TYPE DE CLIQUE
  long double energy=0.0;
  for (uint i=0;i<cliques.size();i++)
    if (cliques[i].type == type)
      energy += cliques[i].energie;
    return energy;
}

long double SurfaceBased_StructuralAnalysis::getTotalEnergy(){
  long double energy=0.0;
  int nclsim=0,nbips=0,nb2=0;
  vector<int> bysub(nbsujets);

  for (uint i=0;i<cliques.size();i++){
    cliques[i].updateLabelsCount();
    if (cliques[i].type == DATADRIVEN){
      energy += cliques[i].computeEnergy(true, nbsujets);
    }
    else if (cliques[i].type == SIMILARITY){
      energy += cliques[i].computeEnergy(true,nbsujets);
      if (cliques[i].blobs[0]->label == cliques[i].blobs[1]->label && cliques[i].blobs[0]->label != 0)
        nclsim++;
        
    }
  }
    
  for (uint k=1;k<labels.size();k++){

    uint i,j;
    for (uint n=0;n<ipscliques.size();n++){
      bysub[n] = cliques[ipscliques[n]].labelscount[labels[k]];
    // cout << bysub[n] << " " ;
    }
    // cout << endl;
    uint nb=0;
    for (i=0;i<nbsujets-1;i++){
      for (j=i+1;j<nbsujets;j++){
        nb += bysub[i]*bysub[j];
      }
    }
    //  cout << "nb"<< nb << " ";
    nbips += nb;
    for (i=0;i<nbsujets;i++)
      if (bysub[i]>1) nb2 += bysub[i]-1;
  }

    
    
  //   Esimil = 4.0*(nbips-nclsim);
  ASSERT(nbips>=nclsim || (cout << nbips << ">=" << nclsim << endl && false));
//   energy += Clique::intrapsweight*(nbips-nclsim);
  energy += Clique::intrapsweight*nb2*nbsujets;
  
  return energy;
}



void SurfaceBased_StructuralAnalysis::SummaryLabels(){
  long double Edd,Esim,energy,Esub,Eintra;
  cout << labels[0] << ":";
  uint nblabel=0;
  for (uint il=0;il<cliques.size();il++){
    nblabel += cliques[il].labelscount[0];
  }
  cout << nblabel << " - ";
  vector<uint> nblab;

  long double Etotal=0.0;
  cout << endl << endl;
  FILE * f;
  if (energypath!=""){
   f = fopen (energypath.data(),"a");
   fprintf(f, "== SUMMARYLABELS ==\n");
  }
  for (uint lab=1;lab<labels.size();lab++){
    Edd=0.0; Esim=0.0; Esub=0.0; Eintra=0.0;
    
    energy=0.0;
    int nclsim=0,nbips=0;
    vector<int> bysub(nbsujets);
  
    for (uint i=0;i<cliques.size();i++){
      cliques[i].updateLabelsCount();
      if (cliques[i].type == DATADRIVEN && cliques[i].blobs[0]->label==labels[lab]){
        energy += cliques[i].computeEnergy(true, nbsujets);
        Edd += cliques[i].computeEnergy(true, nbsujets);
      }
      else if (cliques[i].type == SIMILARITY && cliques[i].blobs[0]->label == cliques[i].blobs[1]->label && cliques[i].blobs[0]->label== labels[lab]){
        energy += cliques[i].computeEnergy(true,nbsujets);
        Esim += cliques[i].computeEnergy(true, nbsujets);
        if (cliques[i].blobs[0]->label == cliques[i].blobs[1]->label && cliques[i].blobs[0]->label != 0)
        nclsim++;
          
      }
    }

    uint k,j;
    for (uint n=0;n<ipscliques.size();n++){
      bysub[n] = cliques[ipscliques[n]].labelscount[labels[lab]];
    }
    uint nb=0;
    for (k=0;k<nbsujets-1;k++){
      for (j=k+1;j<nbsujets;j++){
        nb += bysub[k]*bysub[j];
      }
    }
    nbips += nb;

    uint nb2=0.0;
    for (uint i=0;i<nbsujets;i++)
      if (bysub[i]>1) nb2 += bysub[i]-1;


    ASSERT(nbips>=nclsim || (cout << nbips << ">=" << nclsim << endl && false));
//     energy += Clique::intrapsweight*(nbips-nclsim);
//     Esub += Clique::intrapsweight*(nbips-nclsim);
    energy += Clique::intrapsweight*nb2*nbsujets;
    Eintra += Clique::intrapsweight*nb2*nbsujets;
    Etotal +=energy;

    ASSERT(pow(Edd+Esim+Esub+Eintra-energy,2)<0.01);
    ASSERT(Esub<0.0001);

    nblabel=0;
    for (uint il=0;il<cliques.size();il++){
      nblabel += cliques[il].labelscount[labels[lab]];
    }
    if (nblabel != 0){cout << "label " << labels[lab] << " : " << " (" << Edd << ";" << Esim << ";" << Esub<<";" <<Eintra<<") " << energy << " " << flush;

    cout << nblabel << " - " << endl;
    // fprintf(f,"label %d : (%3lf;%3lf) %i;\n", labels[lab],(double)Edd,(double)Esim,nblabel);
    }
    nblab.push_back(nblabel);

  }
  cout << " ";
  for (uint il=0;il<nblab.size();il++)
    if (nblab[il] != 0) cout << "<<" << nblab[il] << ">>-";
    else cout << nblab[il] << "-" ;
  for (uint i=0;i<sites.size();i++){
    for (uint j=0;j<nblab.size();j++){
      if (sites[i]->label != 0)
        sites[i]->label_occur_number = nblab[sites[i]->label-1];
      else 
        sites[i]->label_occur_number = 0;
    }
  }
  cout <<"\b ";
  cout << "Etot=" << Etotal << endl;

  fclose(f);
}

void SurfaceBased_StructuralAnalysis::ShortSummaryLabels(){
  long double Edd,Esim,energy,Esub,Eintra;
  cout << labels[0] << ":";
  uint nblabel=0;
  for (uint il=0;il<cliques.size();il++){
    nblabel += cliques[il].labelscount[0];
  }
  cout << nblabel << " - ";
  vector<uint> nblab;

  long double Etotal=0.0;

  for (uint lab=1;lab<labels.size();lab++){
    Edd=0.0; Esim=0.0; Esub=0.0; Eintra=0.0;
    
    energy=0.0;
    int nclsim=0,nbips=0;
    vector<int> bysub(nbsujets);
  
    for (uint i=0;i<cliques.size();i++){
      cliques[i].updateLabelsCount();
      if (cliques[i].type == DATADRIVEN && cliques[i].blobs[0]->label==labels[lab]){
        energy += cliques[i].computeEnergy(true, nbsujets);
        Edd += cliques[i].computeEnergy(true, nbsujets);
      }
      else if (cliques[i].type == SIMILARITY && cliques[i].blobs[0]->label == cliques[i].blobs[1]->label && cliques[i].blobs[0]->label== labels[lab]){
        energy += cliques[i].computeEnergy(true,nbsujets);
        Esim += cliques[i].computeEnergy(true, nbsujets);
        if (cliques[i].blobs[0]->label == cliques[i].blobs[1]->label && cliques[i].blobs[0]->label != 0)
        nclsim++;
          
      }
    }
      
    
  
    uint k,j;
    for (uint n=0;n<ipscliques.size();n++){
      bysub[n] = cliques[ipscliques[n]].labelscount[labels[lab]];
    }
    uint nb=0;
    for (k=0;k<nbsujets-1;k++){
      for (j=k+1;j<nbsujets;j++){
        nb += bysub[k]*bysub[j];
      }
    }
    nbips += nb;

    uint nb2=0;
    for (uint k=0;k<nbsujets;k++)
      if (bysub[k]>1) nb2 +=bysub[k]-1;
    

    ASSERT(nbips>=nclsim || (cout << nbips << ">=" << nclsim << endl && false));
//     energy += Clique::intrapsweight*(nbips-nclsim);
//     Esub += Clique::intrapsweight*(nbips-nclsim);
    energy += Clique::intrapsweight*nb2*nbsujets;
    Eintra += Clique::intrapsweight*nb2*nbsujets;
    Etotal +=energy;

    ASSERT(pow(Edd+Esim+Esub+Eintra-energy,2)<0.01);
    ASSERT(Esub<0.0001);
    
    nblabel=0;
    for (uint il=0;il<cliques.size();il++){
      nblabel += cliques[il].labelscount[labels[lab]];
    }
    if (nblabel != 0){
    cout << "L"<<labels[lab] << "(" << energy << "=" << Edd << "+" << Esim << "+" << Esub << "):" << flush;
    cout << nblabel << " - ";
    }
    nblab.push_back(nblabel);
  }
    cout << " ";
    for (uint il=0;il<nblab.size();il++)
      cout << nblab[il] << "-" ;
    cout <<"\b ";
    cout << "Etot=" << Etotal << " "  ;
  
}



void SurfaceBased_StructuralAnalysis::setModelParameters(float _ddweight, float _intrapsweight, float _simweight, float _lsweight, float _ddx2, float _ddx1, float _ddh, float _ddweight2, float _dd2x1, float _dd2x2){
  Clique::setParameters(_ddweight, _intrapsweight, _simweight, _lsweight, _ddx2, _ddx1, _ddh, _ddweight2, _dd2x1, _dd2x2);
}

void SurfaceBased_StructuralAnalysis::StoreToGraph(Graph &primal){
  std::set<Vertex *>::iterator iv, jv;
  vector<float> bc1, bc2;
  float tmin_1, tmax_1, trep, tvalue1;
  int index1;
  int node, label_occur_number;
  string subject1, subject2;
  map<float, vector<pair<float, uint > > >::iterator meshIt;
  vector<pair<float, uint> >::iterator yIt;
  for (iv=primal.vertices().begin() ; iv!=primal.vertices().end(); ++iv){
    string test;
    (*iv)->getProperty("index", index1);
    (*iv)->getProperty( "subject", subject1);
    (*iv)->getProperty("label",test);
    (*iv)->getProperty( "subject", subject1 );
    (*iv)->getProperty( "gravity_center", bc1);
    (*iv)->getProperty( "tmin", tmin_1);
    (*iv)->getProperty( "tmax", tmax_1);
    (*iv)->getProperty( "node", node);
    (*iv)->getProperty( "trep", trep);
    (*iv)->getProperty( "tValue", tvalue1);
    for (uint i=0;i<sites.size();i++)
      if (sites[i]->graph_index == index1 && sites[i]->subject == subject1){
        std::ostringstream s;
        s << sites[i]->label ;
        (*iv)->setProperty("label", s.str());
        (*iv)->setProperty("name", s.str());
        node = sites[i]->node;
        (*iv)->setProperty( "node", node);
        label_occur_number = sites[i]->label_occur_number;
        (*iv)->setProperty( "label_occur_number", label_occur_number);
        
        
      }
  }

}



long double frec0(long double rec){
  if (rec>20.0) return 0.0;
  else if (rec<0.0) return 1.0;
  else return -1.0/20.0*rec+1.0;
}
long double ft0(long double t){
  if (t>16.0) return 1.0;
  else if (t<8.0) return 0.0;
  else return 1.0/8.0*t-8.0/8.0;
}

// long double SurfaceBased_StructuralAnalysis::getCompacite(set<uint> &comp, bool verb){
//   set<string> subj;
//   uint penal = 0;
//   long double compac=0.0,rec=0.0;
//   set<uint> auxcliques;
//   set<uint>::iterator it;
//   float t=0.0;
//   for (it=comp.begin();it!=comp.end();it++){
//     for (uint i=0;i<cliquesDuSite[*it].size();i++){
//       uint aux=cliquesDuSite[*it][i],index0,index1;
//       index0=cliques[aux].blobs[0]->index;
//       index1=cliques[aux].blobs[1]->index;
//       if ((comp.find(index0) != comp.end() && index0 != *it) || (comp.find(index1) != comp.end() && index1 != *it))
//         auxcliques.insert(cliquesDuSite[*it][i]);
//     }
//         
//     compac += sites[*it]->t;
//     if (subj.find(sites[*it]->subject)!=subj.end()) penal=penal+1;
//     subj.insert(sites[*it]->subject);
//   }
//   compac /= comp.size();
//   t = (float) compac;
// 
//   if (verb) cout << "[" << t << ";" << -rec << ";" << subj.size() << ";" << penal <<"]";
// 
//   return compac;
// 
// }




   
  










