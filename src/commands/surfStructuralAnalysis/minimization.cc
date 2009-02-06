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


void SurfaceBased_StructuralAnalysis::MinimizationSetup(Graph &primal, map<string, AimsSurfaceTriangle > &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons){
//   cout << "Building alternate representation of input mesh..." << flush;
//   map<float, vector<pair<float, uint> > > altmesh = getAlternateMesh(mesh, lat, lon);
//   cout << " done " << endl;
  cout << "Construction du vecteur de sites ..." << flush;
  string dir = "/home/grg/data/nmr_surface/";

//   vector<vector<Clique> > allcliques;


  sites = ConstruireSites(primal); //altmesh);
  cout << "done (" << sites.size() << " sites)" << endl;

//   vector<uint> histo,histo_t;
//   for (uint i=0;i<200;i++){
//     histo.push_back(0);
//     histo_t.push_back(0);
//   }
//   for (uint i=0;i<sites.size();i++){
//     if (sites[i]->tValue<-5.0)
//       histo[0]++;
//     else if (sites[i]->tValue>194.0)
//       histo[199]++;
//     else
//     histo[(uint)(8*(sites[i]->tValue+5.0))]++;
//     
//     if (sites[i]->t<-5.0)
//       histo_t[0]++;
//     else if (sites[i]->t>194.0)
//       histo_t[199]++;
//     else
//       histo_t[(uint)(8*(sites[i]->t+5.0))]++;
//   }
//   float cpt = -5.0;
//   for (uint i=0;i<200;i++){
//     cout <<  cpt << " " << histo[i] << " " << histo_t[i] << endl;
//     cpt += 0.125;
//   }
//   cin >> dir;

  set<string> subjects;

  cout << endl << "  done" << endl;
  for (uint i=0;i<sites.size();i++)
    subjects.insert(sites[i]->subject);
  nbsujets = subjects.size();
  cout << "Construction des cliques ... " << flush;
  cliques = ConstruireCliquesLastChance(sites,cliquesDuSite,meshes, lats,lons);
  
  
// === HISTOGRAMME ==========================================
  FILE * f;   f = fopen ("/home/grg/histo.txt","w"); 
  
  vector<uint> histo,histo_t;
  for (uint i=0;i<1600;i++){
    histo.push_back(0);
    histo_t.push_back(0);
  }
  for (uint i=0;i<cliques.size();i++){
    if (cliques[i].rec>40.0)
      histo[1599]++;
    else
      histo[(uint)(cliques[i].rec*cliques[i].rec)]++;
  }
  float cp = 0.0;
  for (uint i=0;i<1600;i++){
    fprintf(f, "%3lf %d\n", sqrt(cp), histo[i]); 
//     cout <<  sqrt(cp) << " " << histo[i] << endl;
    cp += 1;
  }
  fclose(f);

  for (uint i=0;i<sites.size();i++){
    sites[i]->label = 1;
  }

  for (uint i=0;i<cliques.size();i++)
    cliques[i].computeEnergy(false,nbsujets);


  
  cin >> dir;
// ===  FIN HISTOGRAMME ==========================================


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
        cout << i << " " << zone.first[0] << ";" << zone.first[1] << " " << zone.second[0] << ";" << zone.second[1] << endl;
        labelsZones.push_back(zone);
        i++;
      }
  cout << labelsZones.size() << " zones" << endl;
  for (i=0;i<labelsZones.size()+1;i++)
//   for (i=0;i<10;i++)
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
    for (int k=0;k<labelsZones.size();k++){
//       cout << "miaou"<<k << endl;
      Point3df bbmin1 = sites[j]->boundingbox_min, bbmax1 = sites[j]->boundingbox_max;
      uint no_overlap=0;
      getOverlap(bbmin1, bbmax1, Point3df(labelsZones[k].first[0],labelsZones[k].first[1],0.0), Point3df(labelsZones[k].second[0],labelsZones[k].second[1] ,0.0), &no_overlap);
//       cout << "rec="<< rec << endl ;
      if (no_overlap == 0){
//         cout << k << endl;
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

SurfaceBased_StructuralAnalysis::SurfaceBased_StructuralAnalysis(Graph &primal, map<string, AimsSurfaceTriangle> &meshes, map<string, TimeTexture<float> > &lats, map<string, TimeTexture<float> > &lons){
  MinimizationSetup(primal,meshes,lats,lons);
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
  cout << "Initalisation :" << endl;
  for (uint i=0;i<sites.size();i++)
    sites[i]->label = labels[0];
  for (uint k=0;k<cliques.size();k++){
    cliques[k].updateLabelsCount();
    cliques[k].computeEnergy(true,nbsujets);
  }
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
//     cout << nclsim << "|"<< energydd << " " << energysim << "|";
    nbips += nb;
// cout << Clique::intrapsweight*(nbips-nclsim) << " ";
    ASSERT(nbips>=nclsim || (cout << nbips << ">=" << nclsim << endl && false));
  energy += Clique::intrapsweight*(nbips-nclsim);
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
  int nclsim=0,nbips=0;
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
//       cout << bysub[n] << " " ;
    }
//     cout << endl;
    uint nb=0;
    for (i=0;i<nbsujets-1;i++){
      for (j=i+1;j<nbsujets;j++){
        nb += bysub[i]*bysub[j];
      }
    }
//     cout << "nb"<< nb << " ";
    nbips += nb;
  }
    
//   Esimil = 4.0*(nbips-nclsim);
  ASSERT(nbips>=nclsim || (cout << nbips << ">=" << nclsim << endl && false));
  energy += Clique::intrapsweight*(nbips-nclsim);
  
  return energy;
}

// double SurfaceBased_StructuralAnalysis::getTotalEnergyLastChance(uint site, uint newlabel){
//   double energy=0.0;
//   uint old = sites[site]->label;
//   sites[site]->label = newlabel;
// //   cout << "totalenergycalcul" << endl;
//   for (uint i=0;i<cliques.size();i++){
//     long double cen=0.0;
//     if (cliques[i].type == DATADRIVEN){
//       uint ms=cliques[i].blobs[0]->index;
//       if (sites[ms]->label!=0){
//         if (sites[ms]->t>15.0)  cen = 10.0;
//         else if (sites[ms]->t<5.0) cen = 0.001;
//         else cen = sites[ms]->t * (-5.0*0.999/10.0) + (10.0+5.0*0.999)/10.0;
//       }
//     }
//     else if (cliques[i].type == SIMILARITY){
//       uint ms1 = cliques[i].blobs[0]->index;
//       uint ms2 = cliques[i].blobs[1]->index;
//       if (sites[ms1]->label == sites[ms2]->label && sites[ms1]->label != 0){
//         cen = cliques[i].rec / 2.0 - 1.0;
//       }
//     }
//     energy+=cen;
//     
//   }
//   sites[site]->label = old;
//   return energy;
// }



void SurfaceBased_StructuralAnalysis::SummaryLabels(){
  long double Edd,Esim,energy,Esub;
  cout << labels[0] << ":";
  uint nblabel=0;
  for (uint il=0;il<cliques.size();il++){
    nblabel += cliques[il].labelscount[0];
  }
  cout << nblabel << " - ";
  vector<uint> nblab;

  long double Etotal=0.0;
  cout << endl << endl;
  FILE * f;   f = fopen (energypath.data(),"a"); 
  fprintf(f, "== SUMMARYLABELS ==\n"); 

  for (uint lab=1;lab<labels.size();lab++){
    Edd=0.0; Esim=0.0; Esub=0.0;
    
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

    ASSERT(nbips>=nclsim || (cout << nbips << ">=" << nclsim << endl && false));
    energy += Clique::intrapsweight*(nbips-nclsim);
    Esub += Clique::intrapsweight*(nbips-nclsim);
    Etotal +=energy;

    ASSERT(pow(Edd+Esim+Esub-energy,2)<0.01);
//     cout << "L"<<labels[lab] << "(" << energy << "=" << Edd << "+" << Esim << "+" << Esub << "):" << flush;
    
  

    
  
  nblabel=0;
    for (uint il=0;il<cliques.size();il++){
      nblabel += cliques[il].labelscount[labels[lab]];
    }
    if (nblabel != 0){cout << "label " << labels[lab] << " : " << " (" << Edd << ";" << Esim << ";" << Esub<<") " << energy << " " << flush;

    cout << nblabel << " - " << endl;
fprintf(f,"label %d : (%3lf;%3lf) %i;\n", labels[lab],(double)Edd,(double)Esim,nblabel);
    }
    nblab.push_back(nblabel);
    


    

  }
    cout << " ";
    for (uint il=0;il<nblab.size();il++)
      if (nblab[il] != 0) cout << "<<" << nblab[il] << ">>-";
      else cout << nblab[il] << "-" ;
    cout <<"\b ";
    cout << "Etot=" << Etotal << endl;
    
  


  
  fclose(f);
}

void SurfaceBased_StructuralAnalysis::ShortSummaryLabels(){
  long double Edd,Esim,energy,Esub;
  cout << labels[0] << ":";
  uint nblabel=0;
  for (uint il=0;il<cliques.size();il++){
    nblabel += cliques[il].labelscount[0];
  }
  cout << nblabel << " - ";
  vector<uint> nblab;

  long double Etotal=0.0;

  for (uint lab=1;lab<labels.size();lab++){
    Edd=0.0; Esim=0.0; Esub=0.0;
    
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

    ASSERT(nbips>=nclsim || (cout << nbips << ">=" << nclsim << endl && false));
    energy += Clique::intrapsweight*(nbips-nclsim);
    Esub += Clique::intrapsweight*(nbips-nclsim);
    Etotal +=energy;

    ASSERT(pow(Edd+Esim+Esub-energy,2)<0.01);
    
    
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
  int node;
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
//         cout << node << " ";
        (*iv)->getProperty( "node", node);
//         cout << node << " " ;
        
      }
  }

}



// void SurfaceBased_StructuralAnalysis::Validation() {
//   
//   for (uint i=1;i<labels.size();i++){
//     
//     set<uint> listeBlobs;
//     set<uint>::iterator it,it1;
// 
// //     for (uint j=0;j<sites.size();j++){
// //       if (sites[j]->label == labels[i]){
// //         for (it=listeZones[j].begin();it!=listeZones[j].end();it++){
// //           if (*it==0) continue;
// //           for (it1=zonesListesBlobs[*it].begin();it1!=zonesListesBlobs[*it].end();it1++){
// //             listeBlobs.insert(*it1);
// //           }
// //         }
// //       }
// //     }
//     uint zone=i;
//       for (it1=zonesListesBlobs[zone].begin();it1!=zonesListesBlobs[zone].end();it1++)
//         listeBlobs.insert(*it1);
// 
// 
//     cout << listeBlobs.size() << " blobs pour ce label" << endl;
//     long double testcount=1.0;
//     for (it=listeBlobs.begin();it!=listeBlobs.end();it++){
// //       cout << testcount << " ";
//       testcount *= 2.0;// listeZones[*it].size();
//     }
//     cout << "label " << labels[i] << ":" << testcount << " possibilités" << endl;
//     
//     for (uint k=0;k<cliques.size();k++){
//       cliques[k].updateLabelsCount();
//       cliques[k].computeEnergy(true,nbsujets);
//     }
//     ipscliques.clear();
//     for (uint k=0;k<cliques.size();k++)
//       if (cliques[k].type == INTRAPRIMALSKETCH)
//         ipscliques.push_back(k);
//     
//     cout << ipscliques.size() << " cliques intraps" << endl;
//     
// 
//     uint old;
//     vector<uint> permut(listeBlobs.size());
//     vector<uint> listeBlobsV;
//     vector<vector<uint> > listeLabV(listeZones.size());
//     for (uint k=0;k<listeZones.size();k++){
// //       for (it=listeZones[k].begin();it!=listeZones[k].end();it++)
//         listeLabV[k].push_back(0);
//         listeLabV[k].push_back(zone);
//     }
//     for (it=listeBlobs.begin();it!=listeBlobs.end();it++)
//       listeBlobsV.push_back(*it);
//     for (uint k=0;k<permut.size();k++)
//       permut[k]=0;
//     for (uint j=0;j<listeBlobsV.size();j++)
//       sites[j]->label=listeLabV[j][permut[j]];
//     double energy = getTotalEnergy();
//     cout << "energie initiale : " << energy << endl;
//     long double count =1;
//     bool stop=false;
//     if (listeBlobs.size()==0) stop=true;
//     while(!stop){
//           uint site_courant=permut.size()-1;
//           
//           // permutation suivante
//           old = sites[site_courant]->label;
//           permut[site_courant]++;
// //           for (int j=permut.size()-1;j>=0;j--)
//           while (permut[site_courant]==listeZones[site_courant].size()){
//             permut[site_courant]=0;
// 
//             // mise à jour intermédiaire
//             old = sites[site_courant]->label;
//             sites[site_courant]->label = listeLabV[site_courant][permut[site_courant]];
//             int nclsim1=0,nclsim2=0,nbips1=0,nbips2=0;
// 
//             for (uint n=0;n<cliquesDuSite[site_courant].size();n++){
//               uint aux = cliquesDuSite[site_courant][n];
//               if (cliques[aux].type == DATADRIVEN){
//                 energy +=cliques[aux].updateEnergy(site_courant,old,false,nbsujets);
//               }
//               else if (cliques[aux].type == SIMILARITY){
//                 energy +=cliques[aux].updateEnergy(site_courant,old,false,nbsujets);
//                 uint index=0;
//                 if (cliques[aux].blobs[0]->index==(uint)site_courant) index = 1;
//                 else if (cliques[aux].blobs[1]->index==(uint)site_courant) index = 0;
//                 else ASSERT(false);
//                 if (cliques[aux].blobs[index]->label == listeLabV[site_courant][permut[site_courant]] && listeLabV[site_courant][permut[site_courant]] != 0) nclsim1++;
//                 if (cliques[aux].blobs[index]->label == old && old != 0) nclsim2++;
//               }
//             }
//             for (uint n=0;n<ipscliques.size();n++){
//               uint aux = ipscliques[n];
//               if (cliques[aux].blobs[0]->subject != sites[site_courant]->subject){
//                 if (listeLabV[site_courant][permut[site_courant]] !=0)
//                   nbips1+=cliques[aux].labelscount[listeLabV[site_courant][permut[site_courant]]];
//                 if (old != 0)
//                   nbips2+=cliques[aux].labelscount[old];
//               }
//             }
//     
//             energy += Clique::intrapsweight * (nbips1-nclsim1 - (nbips2-nclsim2));
// 
//             site_courant--;
//             old = sites[site_courant]->label;
//             permut[site_courant]++;
//             
//           }
//           
//           assert(sites[site_courant]->label!= listeLabV[site_courant][permut[site_courant]]);
//           sites[site_courant]->label= listeLabV[site_courant][permut[site_courant]];
//           
//           stop=true;
//           for (uint j=0;j<permut.size()&&stop==true;j++)
//             if (permut[j]!=listeZones[j].size()-1) 
//               stop = false;
// 
// //           for (uint j=0;j<permut.size();j++)
// //             cout << listeZones[j][permut[j]] << " ";
// //           cout << endl;
//           
//           
//           
//           // mise à jour énergie
//           int nclsim1=0,nclsim2=0,nbips1=0,nbips2=0;
// 
//           for (uint n=0;n<cliquesDuSite[site_courant].size();n++){
//             uint aux = cliquesDuSite[site_courant][n];
//             if (cliques[aux].type == DATADRIVEN){
//               energy +=cliques[aux].updateEnergy(site_courant,old,false,nbsujets);
//             }
//             else if (cliques[aux].type == SIMILARITY){
//               energy +=cliques[aux].updateEnergy(site_courant,old,false,nbsujets);
//               uint index=0;
//               if (cliques[aux].blobs[0]->index==(uint)site_courant) index = 1;
//               else if (cliques[aux].blobs[1]->index==(uint)site_courant) index = 0;
//               else ASSERT(false);
//               if (cliques[aux].blobs[index]->label == listeLabV[site_courant][permut[site_courant]] && listeLabV[site_courant][permut[site_courant]] != 0) nclsim1++;
//               if (cliques[aux].blobs[index]->label == old && old != 0) nclsim2++;
//             }
//           }
//           for (uint n=0;n<ipscliques.size();n++){
//             uint aux = ipscliques[n];
//             if (cliques[aux].blobs[0]->subject != sites[site_courant]->subject){
//               if (listeLabV[site_courant][permut[site_courant]] !=0)
//                 nbips1+=cliques[aux].labelscount[listeLabV[site_courant][permut[site_courant]]];
//               if (old != 0)
//                 nbips2+=cliques[aux].labelscount[old];
//             }
//           }
//   
//           energy += Clique::intrapsweight * (nbips1-nclsim1 - (nbips2-nclsim2));
// //           cout << energy << " " << flush;
// //           cin >> nclsim1;
// //           vector<uint> histo;
// //           for (uint i=0;i<100;i++){
// //             histo.push_back(0);
// //           }
// //           for (uint i=0;i<cliques.size();i++){
// //             if (cliques[i].rec>40.0)
// //               histo[1599]++;
// //             else
// //               histo[(uint)(cliques[i].rec*cliques[i].rec)]++;
// //           }
//         
//           count++;
//           cout << count << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << flush;
//         }
// 
// //       sites[*it]->label = old;
//       }
//     }

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

long double SurfaceBased_StructuralAnalysis::getCompacite(set<uint> &comp, bool verb){
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
        
    compac += sites[*it]->t;
    if (subj.find(sites[*it]->subject)!=subj.end()) penal=penal+1;
    subj.insert(sites[*it]->subject);
  }
  compac /= comp.size();
  t = (float) compac;
//   cout << "compsize:" << comp.size() << endl ;
//   cout << "compac:" << compac << endl;
//   for (it=auxcliques.begin();it!=auxcliques.end();it++)
//     rec += frec0(cliques[*it].rec); // /10.0+0.9;
  
//   compac += -rec;
//   cout << "rec:" << rec<< endl;
//   compac += penal*100.0;
  
  if (verb) cout << "[" << t << ";" << -rec << ";" << subj.size() << ";" << penal <<"]";
//   compac *=10.0;
//   compac /= (float)comp.size();
//   compac *= ((float)subj.size()/(float)nbsujets);
//   compac /= (float)penal;
//   cout << "= "<< compac << " - " ;
  
  return compac;

}


vector<float> SurfaceBased_StructuralAnalysis::getPseudoSamplesPermut(vector<uint> &listeBlobs) {
vector<float> samples;
float compac=0.0, Ttest=0.0, sum=0.0, energy=0.0;
vector<uint> permut(listeBlobs.size());
for (uint k=0;k<permut.size();k++)
  permut[k]=0;
bool stop=false;
if (listeBlobs.size()==0) stop=true;



while(!stop){
  uint site_courant=permut.size()-1;
  

  permut[site_courant]++;
  while (permut[site_courant]==2){
    permut[site_courant]=0;
      site_courant--;
      permut[site_courant]++;
  }
  compac=0.0; Ttest=0.0; sum=0.0; energy=0.0;
  for (uint j=0;j<permut.size();j++)
    if (permut[j]==0) compac += sites[listeBlobs[j]]->t;
    else compac -= sites[listeBlobs[j]]->t;
  
  compac /= listeBlobs.size();
  for (uint j=0;j<permut.size();j++){
    if (permut[j]!=0) 
      sites[listeBlobs[j]]->t = -sites[listeBlobs[j]]->t;
    sum += pow(sites[listeBlobs[j]]->t-compac,2);
  }

  energy = getLabelEnergy(sites[listeBlobs[0]]->label);
  for (uint j=0;j<permut.size();j++)
    if (permut[j]!=0) 
      sites[listeBlobs[j]]->t = -sites[listeBlobs[j]]->t;


  Ttest = sqrt(listeBlobs.size()*(listeBlobs.size()-1))*compac/sqrt(sum);
//   cout << compac << " " << Ttest << " " << energy <<  endl;
  samples.push_back(Ttest);

  stop=true;
  for (uint j=0;j<permut.size()&&stop==true;j++)
    if (permut[j]!=1) 
      stop = false;
  
}
return samples;
}


vector<float> SurfaceBased_StructuralAnalysis::getPseudoSamplesBootstrap(vector<uint> &listeBlobs) {
vector<float> samples;
float compac=0.0, Ttest=0.0, sum=0.0, energy=0.0;

// vector<uint> permut(listeBlobs.size());
vector<float> old_t(listeBlobs.size());
for (uint j=0;j<listeBlobs.size();j++)
  old_t[j] = sites[listeBlobs[j]]->t;

for (uint i=0;i<1000;i++){
  
  for (uint j=0;j<listeBlobs.size();j++){
    sites[listeBlobs[j]]->t = old_t[(float)UniformRandom() * listeBlobs.size()];
  }
  compac=0.0; Ttest=0.0; sum=0.0; energy=0.0;
  for (uint j=0;j<listeBlobs.size();j++)
    compac += sites[listeBlobs[j]]->t;
   
  
  compac /= listeBlobs.size();
  for (uint j=0;j<listeBlobs.size();j++)
    sum += pow(sites[listeBlobs[j]]->t-compac,2);
  

  energy = getLabelEnergy(sites[listeBlobs[0]]->label);
  for (uint j=0;j<listeBlobs.size();j++)
    sites[listeBlobs[j]]->t = old_t[j];


  Ttest = sqrt(listeBlobs.size()*(listeBlobs.size()-1))*compac/sqrt(sum);
//   cout << compac << " " << Ttest << " " << energy <<  endl;
  samples.push_back(energy);

  
  
}


return samples;

}

vector<float> SurfaceBased_StructuralAnalysis::getPseudoSamplesFullBootstrap(vector<uint> &listeBlobs) {
vector<float> samples;
float compac=0.0, Ttest=0.0, sum=0.0, energy=0.0;

// vector<uint> permut(listeBlobs.size());
vector<float> old_t(listeBlobs.size());

map<uint,float> old_sim;
map<uint,float>::iterator it,jt;
for (uint j=0;j<listeBlobs.size();j++){
  old_t[j] = sites[listeBlobs[j]]->t;
  for (uint k=0;k<cliquesDuSite[listeBlobs[j]].size();k++)
    if (cliques[cliquesDuSite[listeBlobs[j]][k]].type == SIMILARITY) {
      pair<uint,float> clsim;
      clsim.first = cliquesDuSite[listeBlobs[j]][k];
      clsim.second = cliques[cliquesDuSite[listeBlobs[j]][k]].rec;
      old_sim.insert(clsim);
    }
}

for (uint i=0;i<1000;i++){
  
  for (uint j=0;j<listeBlobs.size();j++){
    sites[listeBlobs[j]]->t = old_t[(float)UniformRandom() * listeBlobs.size()];
  }
  for (it=old_sim.begin();it!=old_sim.end();it++){
    uint random = (float)UniformRandom() * old_sim.size(),j=0;
    for (jt=old_sim.begin();j<random && jt!=old_sim.end();jt++,j++){}
    cliques[(*it).first].rec=(*jt).second;
  }
  compac=0.0; Ttest=0.0; sum=0.0; energy=0.0;
  for (uint j=0;j<listeBlobs.size();j++)
    compac += sites[listeBlobs[j]]->t;
   
  
  compac /= listeBlobs.size();
  for (uint j=0;j<listeBlobs.size();j++)
    sum += pow(sites[listeBlobs[j]]->t-compac,2);
  

  energy = getLabelEnergy(sites[listeBlobs[0]]->label);
  for (uint j=0;j<listeBlobs.size();j++)
    sites[listeBlobs[j]]->t = old_t[j];
  for (it=old_sim.begin();it!=old_sim.end();it++){
    cliques[(*it).first].rec=(*it).second;
  }


  Ttest = sqrt(listeBlobs.size()*(listeBlobs.size()-1))*compac/sqrt(sum);
//   cout << compac << " " << Ttest << " " << energy <<  endl;
  samples.push_back(energy);

  
  
}


return samples;

}




void SurfaceBased_StructuralAnalysis::Validation(int type) {
  uint test=0,test1;
  for (int i=1;i<labels.size();i++){
    vector<uint> listeBlobs;
    for (uint j=0;j<sites.size();j++){
      if (sites[j]->label == i){
        listeBlobs.push_back(j);
      }
    }
//     cout << "label "<< i << " : " << getCompacite(listeBlobs,true) << endl;
    if  (listeBlobs.size()==0) continue;
    
    cout << "Label " << i << ":" ;
    cout << listeBlobs.size() << " ";
    float compac=0.0, Ttest=0.0, sum=0.0, energy=0.0;
    for (uint j=0;j<listeBlobs.size();j++)
      compac += sites[listeBlobs[j]]->t;
    compac /= listeBlobs.size();
    for (uint j=0;j<listeBlobs.size();j++)
      sum += pow(sites[listeBlobs[j]]->t-compac,2);
    Ttest = sqrt(listeBlobs.size()*(listeBlobs.size()-1))*compac/sqrt(sum);
    cout.precision(5);
    cout << compac << " T=" << Ttest << " E=" << getLabelEnergy(i) << endl;
    int histosize;
    vector<float> samples;
    if (type == PERMUT){
      cout << "PERMUT" << endl;
      samples=getPseudoSamplesPermut(listeBlobs);
      histosize =10;
    }
    else if (type == BOOTSTRAP){
      cout << "BOOTSTRAP" << endl;
      samples=getPseudoSamplesFullBootstrap(listeBlobs);
      histosize =100;
    }

    while(histosize>9){
      cin >> histosize;
      vector<uint> histo(histosize);
      vector<float> integ(histosize);
      float mini=10000000.0, maxi=-10000000.0, step=0.0;
      sum=0.0;
      for (uint j=0; j<samples.size();j++){
        if (samples[j]>maxi) maxi = samples[j];
        if (samples[j]<mini) mini = samples[j];
      }
      step = (maxi-mini)/(float)histosize;
      cout << "mini:"<<mini<< " maxi:" << maxi << " step:" << step << " nombre d'échantillons:" << samples.size() << endl;
      for (uint j=0;j<samples.size();j++){
        histo[(samples[j]-mini)/(maxi-mini)*(histosize-1)]++;
      }
        
      for (uint j=0;j<histosize;j++){
        integ[j]= (float)step * (float)histo[j] + (float)sum;
        sum += histo[j];
        cout << j*step+mini << " ";
      }
      cout << endl;
      for (uint j=0;j<histosize;j++)
        cout << j*step+mini << " " << histo[j] << endl;;
      cout << endl;
  //     for (uint j=0;j<histosize;j++)
  //       cout << integ[j] << " ";
  //     cout << endl;
      uint alpha = 0;
      uint percent = 99;
      for (uint j=0;integ[j]<percent*sum/100.0;j++) alpha=j;
  
      cout << "diff (=0):" << samples.size() - sum << " alpha("<<percent <<"%) : " << alpha*step+mini;
  
      cout << endl<<endl;
    
  

    }
  }
 

}


   
  










