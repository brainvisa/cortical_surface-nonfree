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

// Function that builds a collection of surf::GreyLevelBlob and surf::ScaleSpaceBlob 
// objects from a previously computed Primal Sketch
void construireBlobs(PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch, 
                     vector<surf::GreyLevelBlob *> &blobs, vector<surf::ScaleSpaceBlob *> &ssblobs, bool initNull = true){
  
//     // Inititalization of the results vectors "blobs" and "ssblobs"
    if (initNull){
      blobs.clear();
      ssblobs.clear();
    }
    
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*> listBlobs 
         = sketch.BlobSet();
    uint iBlob=blobs.size(), iSSblob=ssblobs.size();
    
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*>::iterator itSSB;
    list<GreyLevelBlob<SiteType<AimsSurface<3, Void> >::type > *>::iterator itGLB;
    
    ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *ssb;
    set< SiteType<AimsSurface<3, Void> >::type, 
      ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> >::iterator itPoints;
    
    
    for (itSSB = listBlobs.begin() ; itSSB != listBlobs.end() ; itSSB++){
      
      // For each scale-space blob, we create a surf::ScaleSpaceBlob in "ssblobs" containing 
      //    various surf::GreyLevelBlob objects (being themselves contained in a general resulting 
      // "blobs" vector).
      ssb = *itSSB;
      ssblobs.push_back(new surf::ScaleSpaceBlob());
      surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size() - 1];
//       ssblob->index = iSSblob;
      ssblob->subject = sketch.Subject();
      ssblob->tmin = 999.0;
      ssblob->tmax = -999.0;
      
      for (itGLB = ssb->glBlobs.begin(); itGLB != ssb->glBlobs.end(); itGLB++){
        
        // For each grey-level blob, we create a Blob
        blobs.push_back(new surf::GreyLevelBlob());
        surf::GreyLevelBlob *blob = blobs[blobs.size()-1];
        
        // Each surf::GreyLevelBlob has a specific index iBlob, and a surf::ScaleSpaceBlob has an isurf::ScaleSpaceBlob
//         blob->index = iBlob++;
        blob->ssb_parent = ssblob;
//         blob->subject = sketch.Subject();
        blob->t = (*itGLB)->measurements.t;
        blob->scale = (*itGLB)->GetScale();
        
        // The surf::GreyLevelBlob's nodeslist contains its corresponding nodes indices on the 
        //    mesh it was extracted from.
        set<SiteType<AimsSurface<3, Void> >::type, 
             ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> > listePoints 
                 = (*itGLB)->GetListePoints();
        for (itPoints = listePoints.begin() ; itPoints != listePoints.end() ; itPoints++)
          (blob->nodes).insert((*itPoints).second);
        
        ssblob->blobs.insert(blob);
        
        if (blob->scale < ssblob->tmin)
          ssblob->tmin=blob->scale;
        if (blob->scale > ssblob->tmax)
          ssblob->tmax=blob->scale;
      }
      
      ssblob->t = ssb->GetMeasurements().t;
      
      iSSblob++;
    }
    
    cout << " iBlob : " << iBlob << 
            " iSSblob : " << iSSblob <<
            " blobs.size : " << blobs.size() << 
            " ssblobs.size : " << ssblobs.size() << endl;
}

//##############################################################################

// Function that creates a Scale Space from a texture, a mesh and two coordinates textures
ScaleSpace<AimsSurface<3, Void>, Texture<float> > getScaleSpace(
                                  TimeTexture<float> &laTexture, 
                                  AimsSurfaceTriangle &laMesh, 
                                  TimeTexture<float> &lat, 
                                  TimeTexture<float> &longit){
                                    
      float moy=0.0;
      for (uint i=0;i<laTexture.nItem();i++)
        if (laTexture.item(i) == laTexture.item(i))
          moy += laTexture.item(i);
         
      moy /= laTexture.nItem();
      for (uint i=0;i<laTexture.nItem();i++)
        if (laTexture.item(i) != laTexture.item(i))
          laTexture.item(i) = moy;

      cout << "Size=" << laMesh[0].vertex().size() << endl;
      cout << "Ssmoother creation" << endl;

      FiniteElementSmoother<3, float> *smooth;
      smooth=new FiniteElementSmoother<3, float>(0.01, &(laMesh[0]));
      
      ScaleSpace< AimsSurface<3, Void>, Texture<float> > scale_space(
            &(laMesh[0]), &(laTexture[0]), smooth);
            
      cout << "Scale-space creation" << endl;
      
      vector<Point3df> *coordinates;
      coordinates=new vector<Point3df>();
      
      for (uint i=0;i<lat[0].nItem();i++)
        (*coordinates).push_back(Point3df(lat[0].item(i), longit[0].item(i),i));
      
      scale_space.PutCoordinates(coordinates);
      scale_space.GenerateDefaultScaleSpace(8.0);
      return scale_space;
}


//##############################################################################

void setupData( map<string, SubjectData> &data,
                vector<string> &listTexPaths,
                vector<string> &listMeshPaths,
                vector<string> &listLatPaths,
                vector<string> &listLonPaths,
                vector<string> &listSujets){
  
  assert( listLatPaths.size() == listSujets.size() );
  assert( listLonPaths.size() == listSujets.size() );
  assert( listMeshPaths.size() == listSujets.size() );
  assert( listTexPaths.size() == listSujets.size() );
  
  for ( uint i = 0 ; i < listSujets.size() ; i++ ) {
    
    pair<string, SubjectData > pSubjData;
//     SubjectData subjData;
    pSubjData.first = listSujets[i];
//     pSubjData.second = subjData;
    pSubjData.second.subject = listSujets[i];
    
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

// Creates an Aims Graph for only ONE subject with scale-space blobs, grey-level 
//  blobs and links between both types

void ConstruireIndividualGraph( Graph *graph,
                                vector<surf::GreyLevelBlob *> &blobs,
                                vector<surf::ScaleSpaceBlob *> &ssblobs,
                                string meshPath,
                                string texPath,
                                string latPath,
                                string lonPath,
                                string sujet,
                                string repMeshPath ){
  // From the two vectors, we build a graph containing the ssb, the glb and
  //  the links between ssb and glb                                 

  cerr << "Graph construction..." << endl;
  vector<float> resolution,bbmin2D,bbmax2D;
  vector<int> bbmin, bbmax;
  resolution.push_back(1.0); resolution.push_back(1.0); resolution.push_back(1.0);
  //   bbmin.push_back(mesh[0].minimum()[0]-1); bbmin.push_back(mesh[0].minimum()[1]-1); bbmin.push_back(mesh[0].minimum()[2]-1);
  //   bbmax.push_back(mesh[0].maximum()[0]+1); bbmax.push_back(mesh[0].maximum()[1]+1); bbmax.push_back(mesh[0].maximum()[2]+1);
  bbmin.push_back(-10); bbmin.push_back(-10); bbmin.push_back(-10);
  bbmax.push_back(10); bbmax.push_back(10); bbmax.push_back(10); 
  graph->setProperty( "filename_base", "*");

  graph->setProperty("voxel_size", resolution);
  graph->setProperty("boundingbox_min", bbmin);
  graph->setProperty("boundingbox_max", bbmax);
  graph->setProperty("mesh", meshPath);
  graph->setProperty("sujet", sujet);
  graph->setProperty("texture", texPath);
  graph->setProperty("latitude", latPath);
  graph->setProperty("longitude", lonPath);
  
  AimsSurfaceTriangle repMesh;
  Reader<AimsSurfaceTriangle> rdrMesh ( repMeshPath );
  rdrMesh.read(repMesh);
  AimsSurfaceTriangle *objects = new AimsSurfaceTriangle();
  vector<set<int> > nodes_lists;
  TimeTexture<float> lat,lon;
  Reader<TimeTexture<float> > rdrLat(latPath), rdrLon(lonPath);
  rdrLat.read(lat);
  rdrLon.read(lon);
  
  
  // Extracting mesh patches for the graph
  cout << "Extracting mesh patches for the graph... (from " << repMeshPath << ")" << endl;
  cout << " vertex : " << repMesh[0].vertex().size() << endl;
  cout << " polygon : " << repMesh[0].polygon().size() << endl;
  cout << " size : " << repMesh.size() << endl;
  *objects = getBlobsSphericalMeshes( blobs, repMesh[repMesh.size()-1], lat[0], lon[0], nodes_lists);
  cout << " done" << endl;
  
  
  Vertex *vert;
  carto::rc_ptr<AimsSurfaceTriangle> ptr;
  aims::GraphManip manip;
  vector<Vertex *> listVertSSB( ssblobs.size() ), listVertGLB( blobs.size() );
    
      
  // Let's add the scale-space blobs
  
  cout << "Adding scale-space blobs..." << endl;
  
  for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
        
    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("ssb");
    
//     vert->setProperty("index", i );
    vert->setProperty( "label", "0");
    vert->setProperty( "t", ssblobs[i]->t);
    vert->setProperty( "subject", ssblobs[i]->subject);
    vert->setProperty( "tmin", ssblobs[i]->tmin);
    vert->setProperty( "tmax", ssblobs[i]->tmax);
    ssblobs[i]->index = i;
//     vert->setProperty( "tValue", 100.0);
     

        
    listVertSSB[  i  ] = vert;
  }
  cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added... done" << endl; 
    
  // Let's add the grey-level blobs
  
  cout << "Adding grey-level blobs..." << endl;
  
  int iNbGLB = 0;
  for (int i = 0 ; i < (int) blobs.size() ; i++) {
        
    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index
    
    cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
    vert = graph->addVertex("glb");
    
//     vert->setProperty("index", i );
    vert->setProperty("t", blobs[i]->t);
    vert->setProperty( "scale", blobs[i]->scale);
    vert->setProperty( "nodes", blobs[i]->nodes);
    blobs[i]->index = i;
    listVertGLB[ i ] = vert;
        
    // We associate the proper mesh patch from "objects" to the vertex
    ptr=carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
    (*ptr)[0]=(*objects)[i];
    manip.storeAims(*graph, vert, "glb", ptr);
    vert->setProperty("glb_label", i);
    
    
  }
  cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbGLB << " blobs added... done" << endl; 
   
  // Let's add the links between scale-space and grey-level blobs
  
  cout << "Adding links between both types..." << endl;
  
  uint iNbLinks=0;
  for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
    
    set<surf::GreyLevelBlob *>::iterator itB1;
    set<surf::GreyLevelBlob *> &listGLB = ssblobs[i]->blobs;
    
    for (itB1 = listGLB.begin(); itB1 != listGLB.end() ; itB1++) {
      
      Vertex *v1, *v2;
      
      v1 = listVertSSB[ i ];
      v2 = listVertGLB[(*itB1)->index];
      graph->addEdge(v1,v2,"s2g");
      iNbLinks++;
//       edge->setProperty("ssb_index", ssblobs[i]->index);
//       edge->setProperty("glb_index", (*itB1)->index);
      
    }
  }
    cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " links added... done" << endl; 

}

//##############################################################################

void FromRawTexturesToIndividualGraphsViaPrimalSketches ( string sujets,
                                                          string indivGraphPaths,
                                                          string meshPaths,
                                                          string texPaths,
                                                          string latPaths,
                                                          string lonPaths,
                                                          string repMeshPaths){
  
  map<string, SubjectData> data;
  vector<string> listSujets = splitGraphFile(sujets);
  vector<string> listGraphPaths = splitGraphFile(indivGraphPaths);
  vector<string> listMeshPaths = splitGraphFile(meshPaths);
  vector<string> listTexPaths = splitGraphFile(texPaths);
  vector<string> listLatPaths = splitGraphFile(latPaths);
  vector<string> listLonPaths = splitGraphFile(lonPaths);
  vector<string> listRepMeshPaths = splitGraphFile(repMeshPaths);
      
  cerr << "  split string sujets -> " << listSujets.size() << " sujets" << endl << endl;

  cout << "Reading the data..." << endl;
  setupData(data, listTexPaths, listMeshPaths, listLatPaths, listLonPaths, listSujets);
  cout << "done" << endl<< endl;
  
  
  vector<surf::GreyLevelBlob *> blobs;
  vector<surf::ScaleSpaceBlob *> ssblobs;

  
  // Processing every subject...
  for (uint i = 0 ; i < listSujets.size() ; i++) {
    
    string sujet = listSujets[i];
    
    // Checking the data
    cout << " subject : " << sujet << endl;
    cout << "  texture : " << data[sujet].tex[0].nItem() << " values" << endl;
    cout << "  mesh : " << data[sujet].mesh[0].vertex().size() << " nodes" << endl;
    cout << "  lat : " << data[sujet].lat[0].nItem() << " values" << endl;
    cout << "  lon : " << data[sujet].lon[0].nItem() << " values" << endl;        
    
    // Generating a scale-space
    ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss(
        getScaleSpace( data[sujet].tex, data[sujet].mesh, data[sujet].lat, data[sujet].lon) );

    // Constructing a primal-sketch        
    PrimalSketch<AimsSurface<3, Void>, Texture<float> > sketch(sujet, &ss, SURFACE);

    // Launching the computation of the PS (tmin, tmax, statfile, intersection_criterium)
    sketch.ComputePrimalSketch(1.0, 8.0, "", 10);

    
    // Getting the blobs from PS structure
    cout << "Blobs vectors construction..." << endl;
    construireBlobs(sketch, blobs, ssblobs);      
    cout << "blobs.size() = " << blobs.size() << endl;
    cout << "ssblobs.size() = " << ssblobs.size() << endl;
    
    // Converting the blobs into an Aims Graph
    Graph *tmpGraph = new Graph("BlobsArg");
    ConstruireIndividualGraph(tmpGraph, blobs, ssblobs, listMeshPaths[i], listTexPaths[i], listLatPaths[i], listLonPaths[i], sujet, listRepMeshPaths[i] );
    
    // Storing the graph on the hard disk
    cout << "Writing graph .. " << listGraphPaths[i] << endl;
    Writer<Graph> wtrGraph(listGraphPaths[i]);
    wtrGraph.write(*tmpGraph); 
    cout << "done" << endl << endl;
    
  // Next subject...
  }
  
  // Deleting the blobs vectors
  for (uint i = 0 ; i < blobs.size() ; i++)
    delete(blobs[i]);
  for (uint i = 0 ; i < ssblobs.size() ; i++)
    delete(ssblobs[i]);
    
}


//##############################################################################



int main( int argc, const char **argv ){
  try {

    string indivGraphPaths = "",
           groupGraphPath = "",
           meshPaths = "",
           texPaths = "",
           latPaths = "",
           lonPaths = "",
           repMeshPaths = "",
           sujets = "";

    AimsApplication app( argc, argv, "surfMesh2Graph" );
    app.addOption( meshPaths, "-m", "mesh");
    app.addOption( texPaths, "-t", "texture");
    app.addOption( indivGraphPaths, "-g", "indiv graphs");
    app.addOption( groupGraphPath, "-G", "group graph");
    app.addOption( sujets, "-s", "sujet");
    app.addOption( latPaths, "--lat", "latitude");
    app.addOption( lonPaths, "--lon", "longitude");
    app.addOption( repMeshPaths, "--repM", "repMesh",1);
    app.initialize();


    FromRawTexturesToIndividualGraphsViaPrimalSketches
        ( sujets, indivGraphPaths, meshPaths, texPaths, latPaths, lonPaths, repMeshPaths);
      
    
    
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
