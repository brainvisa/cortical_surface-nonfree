/* Copyright (c) 1995-2007 CEA
 *
 *  This software and supporting documentation were developed by
 *      CEA/DSV/SHFJ
 *      4 place du General Leclerc
 *      91401 Orsay cedex
 *      France
 *
 * This software is governed by the CeCILL license version 2 under 
 * French law and abiding by the rules of distribution of free software.
 * You can  use, modify and/or redistribute the software under the 
 * terms of the CeCILL license version 2 as circulated by CEA, CNRS
 * and INRIA at the following URL "http://www.cecill.info". 
 *  
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 * 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license version 2 and that you accept its terms.
 */
#include <iostream>
#include <cstdlib>
#include <aims/data/data_g.h>
#include <aims/io/io_g.h>
#include <iomanip>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/texture.h>
#include <aims/getopt/getopt2.h>
#include <aims/graph/graphmanip.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/primalsketch/scalespace.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>
#include <aims/primalsketch/primalSketch.h>
#include <cortical_surface/structuralanalysis/representation.h>
#include <cortical_surface/structuralanalysis/cliques.h>
#include <cortical_surface/structuralanalysis/iograph.h>
#include <cortical_surface/structuralanalysis/anneal.h>
#include <cortical_surface/structuralanalysis/blobs.h>

using namespace aims;
using namespace carto;
using namespace std;


class SubjectData{
  public :
    string subject;
    AimsSurfaceTriangle mesh;
    TimeTexture<float> tex;
    TimeTexture<float> lat;
    TimeTexture<float> lon;    
};

//##############################################################################

vector<int> set2vector(set<int> &s){
  
  vector<int> v;
  set<int>::iterator it;
  for (it=s.begin();it!=s.end();it++)
    v.push_back(*it);     
  return v;  
  
}

//##############################################################################

set<int> vector2set(vector<int> &v){
  
  set<int> s;
  for (uint i=0;i<v.size();i++)
    s.insert(v[i]);     
  return s;  
  
}

//##############################################################################


vector<string> getVectorStringFromGraph(Graph &graph, string graph_property){
  vector<string> v;
  if( graph.hasProperty( graph_property ) )  {
    Object slist = graph.getProperty( graph_property ); // note the different getProperty() method
    cout << "node with 'sujets' property:\n";
    Object oit = slist->objectIterator();  // iterator on the list
    while( oit->isValid() ) {
      Object s = oit->currentValue(); // the list element, type Object
      string ss = s->getString(); // extract as std::string or convert to string      
      cout << ss << ", ";
      v.push_back(ss);
      oit->next();
    }
    cout << endl;    
  }     
  
  return v;
  
}

//##############################################################################

void recoverGroupData ( Graph &graph,
                        map<string, SubjectData> &data){
  vector<string> listSujets, listGraphPaths, listMeshPaths, listTexPaths, listLatPaths, listLonPaths;;
  listSujets = getVectorStringFromGraph(graph, "sujets");
  listGraphPaths = getVectorStringFromGraph(graph, "indiv_graphs");
  listTexPaths = getVectorStringFromGraph(graph, "textures");
  listLatPaths = getVectorStringFromGraph(graph, "latitudes");
  listLonPaths = getVectorStringFromGraph(graph, "longitudes");
  listMeshPaths = getVectorStringFromGraph(graph, "meshes");
    
  
  for (uint i = 0 ; i < listSujets.size() ; i++){
    pair<string, SubjectData> pSubjData;
    pSubjData.first = listSujets[i];
    Reader<TimeTexture<float> > texRdr ( listTexPaths[i] ) ;
    texRdr.read(pSubjData.second.tex);
    
    Reader<AimsSurfaceTriangle> meshRdr (listMeshPaths[i] ) ;
    meshRdr.read(pSubjData.second.mesh);
    
    Reader<TimeTexture<float> > latRdr ( listLatPaths[i] );
    latRdr.read(pSubjData.second.lat);
  
    Reader<TimeTexture<float> > lonRdr ( listLonPaths[i] );
    lonRdr.read(pSubjData.second.lon);
    
    // Checking the data
    cout << " subject : " << pSubjData.second.subject << endl;
    cout << "  texture : " << pSubjData.second.tex[0].nItem() << " values" << endl;
    cout << "  mesh : " << pSubjData.second.mesh[0].vertex().size() << " nodes" << endl;
    cout << "  lat : " << pSubjData.second.lat[0].nItem() << " values" << endl;
    cout << "  lon : " << pSubjData.second.lon[0].nItem() << " values" << endl;
    
    data.insert(pSubjData); 
  }                          
                          
}

//##############################################################################

void readGroupGraph ( Graph &graph,
                      map<string, SubjectData> &data,
                      vector<surf::ScaleSpaceBlob *> &ssblobs, 
                      vector<surf::SSBClique> &cliques,
                      vector<Vertex *> &listVertex){

  cout << "Recovering the data..." << endl;
  recoverGroupData(graph, data);
    
  set<Vertex *>::iterator iv;
  // Recovering the scale-space blobs
  string sujet;
  int index;
  int newindex=0;
  cout << " Recovering the scale-space blobs..." << endl;
  for (iv = graph.vertices().begin() ; iv != graph.vertices().end(); ++iv){
    string label;
    
    float t, tmin, tmax, tvalue;
    vector<int> representation;

    if ((*iv)->getSyntax() == "ssb"){
      (*iv)->getProperty("index", index);
      (*iv)->getProperty( "subject", sujet);
      (*iv)->getProperty("label", label);      
      (*iv)->getProperty( "tmin", tmin);
      (*iv)->getProperty( "tmax", tmax);
      (*iv)->getProperty( "t", t);
//       (*iv)->getProperty( "tValue", tvalue);
//       (*iv)->getProperty( "rank", rank);
//       (*iv)->setProperty("label", label);
      (*iv)->getProperty( "nodes", representation);      

      ssblobs.push_back(new surf::ScaleSpaceBlob());
      surf::ScaleSpaceBlob *s=ssblobs[ssblobs.size()-1];
      listVertex.push_back(*iv);


//       (*iv)->getProperty("label", label);

      s->index = newindex++;
      s->label = atoi(label.data());    

      s->graph_index = index;
      s->subject = sujet;
      s->tmin = tmin;
      s->tmax = tmax;
      s->t = t;
      s->nodes = vector2set(representation);      
      (*iv)->setProperty( "sites_index", (int)(ssblobs.size()-1));
      index =0;
      (*iv)->getProperty( "sites_index", index);
//       cout << "index:" << index << " " << flush;
    }

  }


  cout << "ssblobs.size() :"<< ssblobs.size() << endl;

  
  // Recovering the links...
  Edge *e;
  Vertex::iterator jv;
  Edge::iterator kv;
  
  cout << " Recovering the similarity cliques..." << endl;
  for (iv = graph.vertices().begin() ; iv!=graph.vertices().end(); ++iv){
    if ((*iv)->getSyntax() == "ssb"){
//         (*iv)->getProperty( "sites_index", index);
//         cout << "idx:" << index << " " << flush;
        (*iv)->getProperty("index", index);
        (*iv)->getProperty( "subject", sujet);
        for (jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++){
        e = *jv;
        if (e->getSyntax() == "b2b"){
          float similarity;
          e->getProperty("similarity", similarity);
          for (kv = e->begin() ; kv != e->end() ; kv++){
            if ((*kv)->getSyntax() == "ssb"){
              int indexB2;
              string sujetB2;
              (*kv)->getProperty("index", indexB2);
              (*kv)->getProperty("subject", sujetB2);
              if (!(indexB2 == index && sujetB2 == sujet)){
                int blobs_index1, blobs_index2;
                (*iv)->getProperty("sites_index", blobs_index1);
                (*kv)->getProperty("sites_index", blobs_index2);
//                 cout << cliques.size() << " " << ssblobs.size() << " " << blobs_index1 << " " << blobs_index2 << flush;
                if (blobs_index1 < blobs_index2){
                  cliques.push_back(surf::SSBClique(ssblobs[blobs_index1], ssblobs[blobs_index2], similarity));
                }
                
              }

            }
          }
        }
      }
    }
  }

  cout << "ssbcliques.size() :"<< cliques.size() << endl << endl;


}

//##############################################################################

void convertSSBlobsToSites(vector<surf::ScaleSpaceBlob *> &ssblobs, vector<Site *> &sites){
//    int newindex= 0;
  for (uint i = 0 ; i < ssblobs.size() ; i++){
//     if (ssblobs[i]->t>2.0){
     sites.push_back(new Site());
     Site *s = sites[sites.size() - 1];
     s->index = ssblobs[i]->index;
     s->graph_index = ssblobs[i]->index;
     s->subject = ssblobs[i]->subject;
     s->label = ssblobs[i]->label;
     s->tmin = ssblobs[i]->tmin;
     s->tmax = ssblobs[i]->tmax;
     s->t = ssblobs[i]->t;
     s->nodes_list = ssblobs[i]->nodes;
//     }
  }    
    
}

//##############################################################################

void getCliquesFromSSBCliques ( vector<surf::SSBClique> &ssbcliques, 
                                vector<Site *> &sites,
                                vector<Clique> &cliques,
                                vector<vector<int> > &cliquesDuSite){
  
 cliquesDuSite = vector<vector<int> >(sites.size()); 
 
 for (uint i = 0 ; i < ssbcliques.size() ; i++){
   Clique simc;
   simc.type = SIMILARITY;
   simc.rec = ssbcliques[i].similarity;
   surf::ScaleSpaceBlob *ssb1, *ssb2; 
   int iSSB1, iSSB2;   
   ssb1 = ssbcliques[i].ssb1;
   ssb2 = ssbcliques[i].ssb2;
   iSSB1 = ssb1->index;
   iSSB2 = ssb2->index;

   cliquesDuSite[sites[iSSB1]->index].push_back(i);
   cliquesDuSite[sites[iSSB2]->index].push_back(i);
   
   simc.blobs.push_back(sites[iSSB1]);
   simc.blobs.push_back(sites[iSSB2]);
   cliques.push_back(simc);
   
 }
  for (uint i=0; i<sites.size();i++){
  
    for (uint n=0;n<cliquesDuSite[i].size();n++){
        uint aux = cliquesDuSite[i][n];
        if (cliques[aux].type == SIMILARITY){
          if ( cliques[aux].blobs[0]->index == (uint) i ) {}
            else if (cliques[aux].blobs[1]->index ==(uint) i ) {}
            else {
                cout << i << " " << aux << " " << cliques[aux].type << " " << cliques[aux].blobs.size() << " " << cliques[aux].blobs[0]->index << " " << cliques[aux].blobs[1]->index << endl;
                ASSERT(false);
              }          
        }
    }
  
  }
  cout << "ok" << endl;
  
}

//##############################################################################

void computeSitesBoundingBoxes(vector<Site *> &sites, map<string, SubjectData> &data){

  for (uint i = 0 ; i < sites.size() ; i++) {
    set<int> sTemp(sites[i]->nodes_list);
    pair<Point2df, Point2df> bb = getBoundingBox(sTemp, data[sites[i]->subject].lat, data[sites[i]->subject].lon);
    sites[i]->boundingbox_min = Point3df(bb.first[0], bb.first[1], 0);
    sites[i]->boundingbox_max = Point3df(bb.second[0], bb.second[1], 0);
  }

}


//##############################################################################

int main( int argc, const char **argv ){
  try {

    int mode=0;
    string indivGraphPaths = "",
           groupGraphPath = "",
           meshPaths = "",
           texPaths = "",
           latPaths = "",
           lonPaths = "",
           flatPaths = "",
           sujets = "";

    AimsApplication app( argc, argv, "surfLabelsTex2Graph" );
    app.addOption( meshPaths, "-m", "mesh");
    app.addOption( texPaths, "-t", "texture");
    app.addOption( indivGraphPaths, "-g", "indiv graphs");
    app.addOption( groupGraphPath, "-G", "group graph");
    app.addOption( mode, "-M", "mode : 0 (compute the primal sketches) - 1 (take previously computed primal sketches)",1);
    app.addOption( sujets, "-s", "sujet");
    app.addOption( latPaths, "--lat", "latitude");
    app.addOption( lonPaths, "--lon", "longitude");
    app.addOption( flatPaths, "--flat", "flat",1);
    app.initialize();
    if (mode == 2){
        cout << "Reading the group graph..." << endl;
        Graph graph; 
        Reader<Graph> rdrGroupGraph(groupGraphPath);
        rdrGroupGraph.read(graph);
        cout << "Recovering the sites and cliques..." << endl;
        // Getting the sites and cliques
        vector<surf::ScaleSpaceBlob *> ssblobs;
        vector<surf::SSBClique> ssbcliques;
        vector<Vertex *> listVertex;
        map<string, SubjectData> data;
        
        readGroupGraph(graph, data, ssblobs, ssbcliques, listVertex);
        
        // faire l'analyse = étiquetter
        Anneal swc;
//         vector<Site *> sites;
//         vector<Clique> cliques;
//         vector<vector<int> > cliquesDuSite;
        
        // Sauvegarder les labels dans les graphes individuels
        //  sur le modèle de comment on faisait avant
        
        convertSSBlobsToSites(ssblobs, swc.sites);
        cout << swc.sites.size() << " sites generated" << endl;
        
        getCliquesFromSSBCliques(ssbcliques, swc.sites, swc.cliques, swc.cliquesDuSite);
        
        ConstruireCliquesDataDriven(swc.sites, swc.cliquesDuSite, swc.cliques);        
        ConstruireCliquesIntraPS(swc.sites, swc.cliquesDuSite, swc.cliques);
        cout << swc.cliques.size() << " cliques created" << endl;

        
        set<string> subjects;

        cout << endl << "  done" << endl;
        for (uint i=0;i< swc.sites.size();i++)
          subjects.insert(swc.sites[i]->subject);
        swc.nbsujets = subjects.size();
        
        uint nb_cl_sim=0, nb_cl_dd=0, nb_cl_intraps=0, nb_cl_lower=0;
        for (uint i=0;i<swc.cliques.size();i++){
          if (swc.cliques[i].type == SIMILARITY) nb_cl_sim++;
          else if (swc.cliques[i].type == DATADRIVEN) nb_cl_dd++;
          else if (swc.cliques[i].type == BESTLOWERSCALE) nb_cl_lower++;
          else if (swc.cliques[i].type == INTRAPRIMALSKETCH) nb_cl_intraps++;
        }
        cout << " done (" << nb_cl_sim << " cliques de similarité ; " << nb_cl_dd << " cliques datadriven ; " << nb_cl_lower << " cliques lower ; " << nb_cl_intraps << " cliques intraps ; " << swc.cliques.size() << " cliques en tout)" << endl;
        computeSitesBoundingBoxes(swc.sites, data);
        swc.prepareLabelsZones();
        
        
  
        float _ddweight=0.8, _intrapsweight = 4.0, _simweight=1.0, _lsweight=1.0, _ddx1 = 8.0, _ddx2 = 4.0, _ddh=0.0001;
        cout << _ddweight << "-" << _intrapsweight << "-" << _simweight << "-" << _lsweight << "-" << _ddx2 << "-" << _ddx1 << "-" << _ddh << endl;
        swc.setModelParameters(_ddweight, _intrapsweight, _simweight, _lsweight, _ddx2, _ddx1, _ddh);
        swc.run = 1;
        swc.Initialization();
      
        swc.Run(1);
        swc.SummaryLabels();
 
    }
      
      
  

      // On vient d'écrire le graphe pour l'analyse (à créer : une fonction qui lit
      //   ce graphe et qui récupère les sites et cliques)
          
      // Eventuellement une fonction pour séparer le graphe en sous-graphes correspondants
      //   aux différents sujets
    
     
    return EXIT_SUCCESS;
  }
  catch( carto::user_interruption & )
  {
  }
  catch( exception & e )
  {
    cerr << e.what() << endl;
  }
}
