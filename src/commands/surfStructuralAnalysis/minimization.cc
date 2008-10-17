#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "minimization.h"

using namespace aims;
using namespace carto;
using namespace std;

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
//   allcliques.push_back(cliques);

  


//   for (uint i=0;i<cliques.size();i++){
//     if (cliques[i].type==SIMILARITY){
//       if (cliques[i].blobs[0]->rank <5 && cliques[i].blobs[1]->rank <5){
// //         cout << cliques[i].blobs[0]->rank << " " << cliques[i].blobs[0]->subject << " - " << cliques[i].blobs[1]->rank << " " << cliques[i].blobs[1]->subject << " : " << cliques[i].rec << endl;
//         float distance = 0.0, ecart = 0.0;
//         vector<float> vectdist;
//         for (uint j=0;j<allcliques.size();j++){
//           for (uint k=0;k<allcliques[j].size();k++){
//             if ((allcliques[j][k].blobs[0]->index == cliques[i].blobs[0]-> index && allcliques[j][k].blobs[1]->index == cliques[i].blobs[1]-> index) || (allcliques[j][k].blobs[0]->index == cliques[i].blobs[1]-> index && allcliques[j][k].blobs[1]->index == cliques[i].blobs[0]-> index)){
// //               cout << /*allcliques[j][k].blobs[0]->rank << " " << allcliques[j][k].blobs[0]->subject << " - " << allcliques[j][k].blobs[1]->rank << " " << allcliques[j][k].blobs[1]->subject << " : " <<*/ allcliques[j][k].rec << " " ;
//               distance += allcliques[j][k].rec;
//               vectdist.push_back(allcliques[j][k].rec);
//             }
//           }
//         }
//         distance /= allcliques.size();
//         cout << distance << " ";
//         
//         for (uint j=0;j<vectdist.size();j++){
//           ecart += (vectdist[j] - distance) *(vectdist[j] - distance) ;
//         }
//         ecart /= vectdist.size();
//         cout << sqrt(ecart) << endl;
//       }
//         
//     }
//   }


    
  
  
  
  
  uint nb_cl_sim=0, nb_cl_dd=0, nb_cl_intraps=0, nb_cl_lower=0;
  for (uint i=0;i<cliques.size();i++){
    if (cliques[i].type == SIMILARITY) nb_cl_sim++;
    else if (cliques[i].type == DATADRIVEN) nb_cl_dd++;
    else if (cliques[i].type == BESTLOWERSCALE) nb_cl_lower++;
    else if (cliques[i].type == INTRAPRIMALSKETCH) nb_cl_intraps++;
  }
  cout << " done (" << nb_cl_sim << " cliques de similarité ; " << nb_cl_dd << " cliques datadriven ; " << nb_cl_lower << " cliques lower ; " << nb_cl_intraps << " cliques intraps ; " << cliques.size() << " cliques en tout)" << endl;
  for (int i=0;i<20;i++)
    labels.push_back(i);
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

double SurfaceBased_StructuralAnalysis::getLabelEnergy(int label, int type){
  bool test = true;
  double energy=0.0;
  for (uint i=0;i<cliques.size();i++){
    if (type==UNKNOWN || cliques[i].type == type){
      test = true;
      for (uint j=0;j<cliques[i].blobs.size() && test == true;j++)
        if (cliques[i].blobs[j]->label != label ) test = false;
      if (test) energy += cliques[i].energie;
      if (cliques[i].type == INTRAPRIMALSKETCH && cliques[i].labelscount[label] > 1 && label != 0) {energy += (cliques[i].labelscount[label]-1)*Clique::getIntraPSWeight();}
    }
  }
  return energy;
}

double SurfaceBased_StructuralAnalysis::getTypeEnergy(int type){ // RETOURNE L'ENERGIE PAR TYPE DE CLIQUE
  double energy=0.0;
  for (uint i=0;i<cliques.size();i++)
    if (cliques[i].type == type)
      energy += cliques[i].energie;
    return energy;
}

double SurfaceBased_StructuralAnalysis::getTotalEnergy(){
  double energy=0.0;
  for (uint i=0;i<cliques.size();i++)
    energy += cliques[i].energie;
  return energy;
}

double SurfaceBased_StructuralAnalysis::getTotalEnergyLastChance(uint site, uint newlabel){
  double energy=0.0;
  uint old = sites[site]->label;
  sites[site]->label = newlabel;
//   cout << "totalenergycalcul" << endl;
  for (uint i=0;i<cliques.size();i++){
    long double cen=0.0;
    if (cliques[i].type == DATADRIVEN){
      uint ms=cliques[i].blobs[0]->index;
      if (sites[ms]->label!=0){
        if (sites[ms]->t>15.0)  cen = 10.0;
        else if (sites[ms]->t<5.0) cen = 0.001;
        else cen = sites[ms]->t * (-5.0*0.999/10.0) + (10.0+5.0*0.999)/10.0;
      }
    }
    else if (cliques[i].type == SIMILARITY){
      uint ms1 = cliques[i].blobs[0]->index;
      uint ms2 = cliques[i].blobs[1]->index;
      if (sites[ms1]->label == sites[ms2]->label && sites[ms1]->label != 0){
        cen = cliques[i].rec / 2.0 - 1.0;
      }
    }
    energy+=cen;
    
  }
  sites[site]->label = old;
  return energy;
}



void SurfaceBased_StructuralAnalysis::SummaryLabels(){
  float Eintra,Edd,Els,Esim,Etot;
  for (uint i=1;i<labels.size();i++){
    Eintra = getLabelEnergy(labels[i], INTRAPRIMALSKETCH);
    Edd = getLabelEnergy(labels[i], DATADRIVEN);
    Els = getLabelEnergy(labels[i], BESTLOWERSCALE);
    Esim = getLabelEnergy(labels[i], SIMILARITY);
    Etot = Eintra + Edd + Els + Esim;
    cout << "label " << labels[i] << " : " << Etot << " (" << Eintra << ";" << Edd << ";" << Esim << ";" << Els << ") " << flush;

    
    uint nblabel=0;
    for (uint il=0;il<cliques.size();il++)
      nblabel += cliques[il].labelscount[i];
    cout << nblabel << ";";
    for (uint j=0;j<sites.size();j++)
      if (sites[j]->label==labels[i]){
        cout << sites[j]->index << "(" << sites[j]->subject << ")-";
      }
    cout<<  endl;
  }
}

void SurfaceBased_StructuralAnalysis::ShortSummaryLabels(){
  float Eintra,Edd,Els,Esim,Etot;
  cout << labels[0] << ":";
  uint nblabel=0;
  for (uint il=0;il<cliques.size();il++)
    nblabel += cliques[il].labelscount[0];
  cout << nblabel << "-";
  for (uint i=1;i<labels.size();i++){
    Eintra = getLabelEnergy(labels[i], INTRAPRIMALSKETCH);
    Edd = getLabelEnergy(labels[i], DATADRIVEN);
    Els = getLabelEnergy(labels[i], BESTLOWERSCALE);
    Esim = getLabelEnergy(labels[i], SIMILARITY);
    Etot = Eintra + Edd + Els + Esim;
    cout << "L"<<labels[i] << "(" << Etot << "=" << Eintra << "+" << Edd << "+" << Esim << "+" << Els << "):" << flush;

    
    nblabel=0;
    vector<uint> nblab;
    for (uint il=0;il<cliques.size();il++)
      nblabel += cliques[il].labelscount[i];
    cout << nblabel << " - ";
    nblab.push_back(nblabel);
    if (nblabel==10 && Eintra<0.001) {
      cout << "DETAIL";

      for (uint jl=0;jl<cliques.size();jl++){
        if (cliques[jl].type == SIMILARITY){
          if (cliques[jl].blobs[0]->label==labels[i] && cliques[jl].blobs[1]->label==labels[i])
            cout << cliques[jl].blobs[0]->index << "-" << cliques[jl].blobs[0]->subject << ";"<< cliques[jl].blobs[1]->index << "-" << cliques[jl].blobs[1]->subject << "=" << cliques[jl].energie <<"(" << cliques[jl].rec << ");" << flush;
        }
      }
    }
    cout << " ";
    for (uint il=0;il<nblab.size();il++)
      cout << nblab[il] << "-" ;
    cout <<"\b ";
  }
}



void SurfaceBased_StructuralAnalysis::setModelParameters(float _ddweight, float _intrapsweight, float _simweight, float _lsweight, float _ddx2, float _ddx1, float _ddh){
  Clique::setParameters(_ddweight, _intrapsweight, _simweight, _lsweight, _ddx2, _ddx1, _ddh);
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




