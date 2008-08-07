#include <aims/getopt/getopt2.h>
#include <aims/math/random.h>
#include "minimization.h"

using namespace aims;
using namespace carto;
using namespace std;

void SurfaceBased_StructuralAnalysis::MinimizationSetup(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon){
  cout << "Building alternate representation of input mesh..." << flush;
  map<float, vector<pair<float, uint> > > altmesh = getAlternateMesh(mesh, lat, lon);
  cout << " done " << endl;
  cout << "Construction du vecteur de sites ..." << flush;
  sites = ConstruireSites(primal, mesh, lat, lon); //altmesh);
  cout << "done (" << sites.size() << " sites)" << endl;
  set<string> subjects;
  
  cout << endl << "  done" << endl;
  for (uint i=0;i<sites.size();i++)
    subjects.insert(sites[i]->subject);
  nbsujets = subjects.size();
  cout << "Construction des cliques ... " << flush;
  cliques = ConstruireCliques(sites,cliquesDuSite,mesh);
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

SurfaceBased_StructuralAnalysis::SurfaceBased_StructuralAnalysis(Graph &primal, AimsSurfaceTriangle &mesh, TimeTexture<float> &lat, TimeTexture<float> &lon){
  MinimizationSetup(primal,mesh,lat,lon);
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

void SurfaceBased_StructuralAnalysis::SummaryLabels(){
  float Eintra,Edd,Els,Esim,Etot;
  for (uint i=0;i<labels.size();i++){
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
        cout << sites[j]->index << "(" << sites[j]->node << ")-";
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



