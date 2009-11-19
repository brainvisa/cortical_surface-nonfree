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
#include <cortical_surface/structuralanalysis/representation.h>
#include <cortical_surface/structuralanalysis/iograph.h>




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
  graph->setProperty("representation_mesh", repMeshPath);

  AimsSurfaceTriangle repMesh;
  Reader<AimsSurfaceTriangle> rdrMesh ( repMeshPath );
  rdrMesh.read(repMesh);
  AimsSurfaceTriangle *objects = new AimsSurfaceTriangle();
  vector<set<int> > nodes_lists;
  TimeTexture<float> lat,lon,tex;
  Reader<TimeTexture<float> > rdrLat(latPath), rdrLon(lonPath), rdrTex(texPath);
  rdrLat.read(lat);
  rdrLon.read(lon);
  rdrTex.read(tex);
  TimeTexture<short> texshort(1,tex[0].nItem());
  for (uint i = 0 ; i < tex[0].nItem() ; i ++){
    texshort[0].item(i) = tex[0].item(i);
  }


  // Extracting mesh patches for the graph
  cout << "Extracting mesh patches for the graph... (from " << repMeshPath << ")" << endl;
  cout << " vertex : " << repMesh[0].vertex().size() << endl;
  cout << " polygon : " << repMesh[0].polygon().size() << endl;
  cout << " size : " << repMesh.size() << endl;
  *objects = getLabelObjectsOnASphere( texshort, repMesh[repMesh.size()-1], lat[0], lon[0], nodes_lists);
  cout << " done" << endl;


  Vertex *vert;
  carto::rc_ptr<AimsSurfaceTriangle> ptr;
  aims::GraphManip manip;
//   vector<Vertex *> listVertSSB( ssblobs.size() ), listVertGLB( blobs.size() );


  // Let's add the scale-space blobs

//   cout << "Adding scale-space blobs..." << endl;
//
//   for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
//
//     // For every scale-space blob, we create a vertex in the Aims graph : we define
//     //   its properties and store a link between the created vertex and the blob index
//
//     cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
//     vert = graph->addVertex("ssb");
//
// //     vert->setProperty("index", i );
//     vert->setProperty( "label", "0");
//     vert->setProperty( "t", ssblobs[i]->t);
//     vert->setProperty( "subject", ssblobs[i]->subject);
//     vert->setProperty( "tmin", ssblobs[i]->tmin);
//     vert->setProperty( "tmax", ssblobs[i]->tmax);
//     ssblobs[i]->index = i;
// //     vert->setProperty( "tValue", 100.0);
//
//     listVertSSB[  i  ] = vert;
//   }
//   cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added... done" << endl;

  // Let's add the grey-level blobs

  cout << "Adding grey-level blobs..." << endl;

  int iNbGLB = 0;
//   for (int i = 0 ; i < (int) blobs.size() ; i++) {

    // For every scale-space blob, we create a vertex in the Aims graph : we define
    //   its properties and store a link between the created vertex and the blob index

//     cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
//     vert = graph->addVertex("glb");

//     vert->setProperty("index", i );
//     vert->setProperty("t", blobs[i]->t);
//     vert->setProperty( "scale", blobs[i]->scale);
//     vert->setProperty( "nodes", blobs[i]->nodes);
//     blobs[i]->index = i;
//     listVertGLB[ i ] = vert;


  for (int i = 0 ; i < (int) (*objects).size() ; i++){
    vert = graph->addVertex("glb");
    // We associate the proper mesh patch from "objects" to the vertex
    ptr=carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
    (*ptr)[0]=(*objects)[i];
    manip.storeAims(*graph, vert, "glb", ptr);
    vert->setProperty("glb_label", i);


  }
  cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added... done" << endl;

//   // Let's add the links between scale-space and grey-level blobs
//
//   cout << "Adding links between both types..." << endl;
//
//   uint iNbLinks=0;
//   for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
//
//     set<surf::GreyLevelBlob *>::iterator itB1;
//     set<surf::GreyLevelBlob *> &listGLB = ssblobs[i]->blobs;
//
//     for (itB1 = listGLB.begin(); itB1 != listGLB.end() ; itB1++) {
//
//       Vertex *v1, *v2;
//
//       v1 = listVertSSB[ i ];
//       v2 = listVertGLB[(*itB1)->index];
//       graph->addEdge(v1,v2,"s2g");
//       iNbLinks++;
// //       edge->setProperty("ssb_index", ssblobs[i]->index);
// //       edge->setProperty("glb_index", (*itB1)->index);
//
//     }
//   }
//     cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " links added... done" << endl;

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
//     ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss(
//         getScaleSpace( data[sujet].tex, data[sujet].mesh, data[sujet].lat, data[sujet].lon) );

    // Constructing a primal-sketch
//     PrimalSketch<AimsSurface<3, Void>, Texture<float> > sketch(sujet, &ss, SURFACE);

    // Launching the computation of the PS (tmin, tmax, statfile, intersection_criterium)
//     sketch.ComputePrimalSketch(1.0, 8.0, "", 10);


    // Getting the blobs from PS structure
//     cout << "Blobs vectors construction..." << endl;
//     construireBlobs(sketch, blobs, ssblobs);
//     cout << "blobs.size() = " << blobs.size() << endl;
//     cout << "ssblobs.size() = " << ssblobs.size() << endl;

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

    int mode=0;
    string indivGraphPaths = "",
           groupGraphPath = "",
           meshPaths = "",
           texPaths = "",
           latPaths = "",
           lonPaths = "",
           repMeshPaths = "",
           sujets = "";

    AimsApplication app( argc, argv, "surfLabelsTex2Graph" );
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
