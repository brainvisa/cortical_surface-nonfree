

#include <cortical_surface/structuralanalysis/region.h>
#include <cortical_surface/structuralanalysis/texturetoblobs.h>

#include <aims/getopt/getopt2.h>
#include <aims/primalsketch/primalSketch.h>

#include <aims/graph/graphmanip.h>
#include <aims/mesh/surfaceOperation.h>
#include <cortical_surface/structuralanalysis/representation.h>
#include <cortical_surface/structuralanalysis/blobs.h>
#include <cortical_surface/surfacereferential/gyri/mesh_operations.h>
#include <aims/primalsketch/primalSketchUtil.h>
#include <aims/distancemap/meshdistance.h>


using namespace aims;
using namespace carto;
using namespace std;

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

void storeCoordinatesInScaleSpace ( SubjectData &regionData, ScaleSpace<AimsSurface<3, Void>, Texture<float> >  &ss ) {
    vector<Point3df> *coordinates;
    coordinates = new vector<Point3df>();

    if ( regionData.coordinates == LATLON_2D ) {

        for ( uint i = 0 ; i < regionData.lat->nItem() ; i++ )
            (*coordinates).push_back(
                Point3df( regionData.lat->item(i), regionData.lon->item(i), i ) );

    }
    else {

        for ( uint i = 0 ; i < regionData.lat->nItem() ; i++ )
            (*coordinates).push_back(
                Point3df( regionData.lat->item(i), -1.0, i ) );

    }
    ss.PutCoordinates(coordinates);
}

void TextureToBlobs::PrimalSketchRegionMode (   vector<surf::GreyLevelBlob *> &blobs,
                                                vector<surf::ScaleSpaceBlob *> &ssblobs,
                                                surf::Region &region,
                                                SubjectData &regionData,
                                                string scaleSpacePath, string blobsPath,
                                                bool recover,
                                                float scale_max ) {

    // Creating the smoother...
    cout << endl << "  ══ Smoother creation... " << endl;
    FiniteElementSmoother<3, float> smooth ( 0.01, regionData.mesh, regionData.weightLapl );
    ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss ( regionData.mesh, regionData.tex, &smooth );

    if ( regionData.coordinates == LATLON_2D || regionData.coordinates == LAT_1D ) {
        storeCoordinatesInScaleSpace( regionData, ss );
    }


    if ( recover ) {
        TimeTexture<float> scale_space;
        // Extracting The Gyrus-Specific Data From The Scale-Space Texture
        Reader< TimeTexture<float> > rdrScaleSpace ( scaleSpacePath ) ;
        rdrScaleSpace.read ( scale_space );
        assert( scale_space.size() > 1 );
        TimeTexture<float> regionScaleSpace;
        for ( uint i = 1 ; i < scale_space.size() - 1 ; i++ )
            regionScaleSpace[i-1] = region.getLocalFromGlobalTexture( scale_space[i] );

        // Recovering Previously Computed Scale-Space
        cerr << "   ══ Recovering scale-space " << endl;
        ss.uploadPreviouslyComputedScaleSpace(regionScaleSpace);

    }
    else {

        // Generating A Scale-Space And Writing It Onto Hard Disk
        cerr << "══ Computing scale-space " << endl;

        ss.GenerateDefaultScaleSpace( 128.0 );

        cout << "══ Writing scale-space (after rebuilding the whole texture from the gyrus-specific data)..." << endl;
        TimeTexture<float> regionScaleSpace = ss.getScaleSpaceTexture( );

        TimeTexture<float> scale_space;
        for ( uint i = 0 ; i < regionScaleSpace.size() ; i++ )
            scale_space[i] = region.getGlobalFromLocalTexture( regionScaleSpace[i] );

        Writer<TimeTexture<float> > wtrScaleSpace ( scaleSpacePath );
        wtrScaleSpace.write ( scale_space );
    }

    TimeTexture<float> regionBlobsTex;

    TextureToBlobs::PrimalSketch ( regionData, blobs, ssblobs, &ss, regionBlobsTex, scale_max );

    // Update Coordinates With Conversion Gyrus/Global Indices
    for ( uint i = 0 ; i < blobs.size() ; i ++ ) {
        set<int>::iterator it;
        set<int> nodes_aux (blobs[i]->nodes);
        map<int, vector<float> > coordinates_aux (blobs[i]->coordinates);
        map<int, vector<float> > raw_coordinates_aux (blobs[i]->raw_coordinates);

        blobs[i]->nodes.clear();
        blobs[i]->coordinates.clear();
        blobs[i]->raw_coordinates.clear();
        for ( it =  nodes_aux.begin() ; it != nodes_aux.end() ; it ++ ) {
            blobs[i]->nodes.insert( region.nodes[*it] );
            blobs[i]->coordinates[region.nodes[*it]] = vector<float>(coordinates_aux[*it]);
            blobs[i]->raw_coordinates[region.nodes[*it]] = vector<float>(raw_coordinates_aux[*it]);
        }


    }

    TimeTexture<float> blobs_tex;
    cout << "══ Writing blobs texture (after rebuilding the whole texture from the gyrus-specific data)..." << endl;
    for ( uint i = 0 ; i < regionBlobsTex.size() ; i++ )
        blobs_tex[i] = region.getGlobalFromLocalTexture( regionBlobsTex[i] );
    Writer<TimeTexture<float> > wtrBlobs ( blobsPath );
    wtrBlobs.write ( blobs_tex );


}

void TextureToBlobs::PrimalSketchGlobalMode (   vector<surf::GreyLevelBlob *> &blobs,
                                    vector<surf::ScaleSpaceBlob *> &ssblobs,
                                    SubjectData &subject,
                                    string scaleSpacePath,
                                    string blobsPath,
                                    bool recover,
                                    float scale_max ) {

    FiniteElementSmoother<3, float> smooth ( 0.01, subject.mesh, subject.weightLapl );
    ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss ( subject.mesh, subject.tex, &smooth );

    if ( subject.coordinates == LATLON_2D || subject.coordinates == LAT_1D ) {
        storeCoordinatesInScaleSpace( subject, ss );
    }
    cout << subject.mesh->vertex().size() << " " << subject.tex->nItem() << endl;

    if ( recover ) {
        TimeTexture<float> scale_space;

        Reader< TimeTexture<float> > rdrScaleSpace ( scaleSpacePath ) ;
        rdrScaleSpace.read ( scale_space );

        // Recovering Previously Computed Scale-Space

//        cout << "  Checking that the scale-space has more than just one texture..." << flush;
        assert( scale_space.size() > 1 );
//        cout << "OK (" << scale_space.size() << ") ( tex[0] is supposed to be the original texture (scale : 0) )" << endl;
//        uint nb_scales = scale_space.size() - 1;
        TimeTexture<float> scaleSpaceTex;
//        cout << "Filtering in the first " << nb_scales << " scales (to allow using bigger pre-computed scale-spaces)..." << flush;
        for ( uint i = 1 ; i < scale_space.size() ; i++ )
            scaleSpaceTex[i-1] = scale_space[i];
//
//        cout << "OK" << endl;
        cerr << "   ══ Recovering scale-space " << endl;

        ss.uploadPreviouslyComputedScaleSpace(scaleSpaceTex);

    }
    else {

        // Generating A Scale-Space And Writing It Onto Hard Disk
        cerr << "══ Computing scale-space " << endl;
        ss.GenerateDefaultScaleSpace( 128.0 );
        TimeTexture<float> scale_space;
        cout << "══ Writing scale-space..." << endl;
        scale_space = ss.getScaleSpaceTexture( );
        Writer<TimeTexture<float> > wtrScaleSpace ( scaleSpacePath );
        wtrScaleSpace.write ( scale_space );

    }

    cout << "══ Computing primal sketch..." << endl;
    TimeTexture<float> blobs_tex;
    TextureToBlobs::PrimalSketch ( subject, blobs, ssblobs, &ss, blobs_tex, scale_max );
    cout << "══ Writing blobs texture..." << endl;
    Writer<TimeTexture<float> > wtrBlobs ( blobsPath );
    wtrBlobs.write ( blobs_tex );

}

void TextureToBlobs::PrimalSketch ( SubjectData &subject,
                    std::vector<surf::GreyLevelBlob *> &blobs,
                    std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                    ScaleSpace<AimsSurface<3, Void>, Texture<float> > *ss,
                    TimeTexture<float> &blobs_texture,
                    float scale_max ) {



    // Constructing A Primal-Sketch
    aims::PrimalSketch<AimsSurface<3, Void>, Texture<float> > sketch ( subject.subject_id, SURFACE );
    cout << "  ══ setting scale-space..." << endl;
    sketch.SetScaleSpace(ss);

    // Scale_max == -1.0 Enforces Scale_max To Be Autodetermined (Max Scale From Ss)
    if ( scale_max == -1.0 ) {
        set<float>::iterator it;
        set<float> scales = ss->GetScaleList();
        scale_max = 0.0;
        for (it = scales.begin() ; it != scales.end() ; it ++ )
            if (*it > scale_max)
                scale_max = *it;
    }

    cout << "  ══ Computing primal sketch..." << endl;
    // Launching The Computation Of The PS (tmin, tmax, statfile, intersection_criterium)
    sketch.ComputePrimalSketch( 1.0, scale_max, "", 1 );

    // Writing A Blob Texture
    cerr << "   ══ Getting scale-space blob texture... " << endl;
    blobs_texture = GetSSBlobTexture( & sketch );

    // Getting The Blobs From The Primal Sketch Structure
    cerr << "    ▪ Blobs vectors construction..." << endl;
    getBlobsFromPrimalSketch ( subject, sketch, blobs, ssblobs );


    cout << blobs.size() << " grey-level blobs / " << ssblobs.size() << " scale-space blobs" << endl;



}


//##############################################################################

// Function that builds a collection of surf::GreyLevelBlob and surf::ScaleSpaceBlob
// objects from a previously computed Primal Sketch
void TextureToBlobs::getBlobsFromPrimalSketch ( SubjectData & subject,
                     aims::PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch,
                     vector<surf::GreyLevelBlob *> &blobs,
                     vector<surf::ScaleSpaceBlob *> &ssblobs,
                     bool initNull ) {

    // Initialization of the results vectors "blobs" and "ssblobs"
    if (initNull){
        blobs.clear();
        ssblobs.clear();
    }

    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*> listBlobs
         = sketch.BlobSet();

    std::list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*>::iterator itSSB, itSSBaux, itSSBaux2;
    std::list<GreyLevelBlob<SiteType<AimsSurface<3, Void> >::type > *>::iterator itGLB;

    ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *ssb;
    std::set< SiteType<AimsSurface<3, Void> >::type,
        ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> >::iterator itPoints;

    std::map< ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*, surf::ScaleSpaceBlob * > ssbMap;

    for (itSSB = listBlobs.begin() ; itSSB != listBlobs.end() ; itSSB++){

        // For each scale-space blob, we create a surf::ScaleSpaceBlob in "ssblobs" containing
        //    various surf::GreyLevelBlob objects (being themselves contained in a general resulting
        // "blobs" vector).
        ssb = *itSSB;

        ssblobs.push_back(new surf::ScaleSpaceBlob());
        surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size() - 1];
        ssblob->subject = sketch.Subject();
        ssblob->tmin = 999.0;
        ssblob->tmax = -999.0;

        // We save a link between the pointer and the corresponding index in ssblobs
        std::pair< ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*, surf::ScaleSpaceBlob * > p;
        p.first = ssb;
        p.second = ssblob;
        ssbMap.insert(p);

        for ( itGLB = ssb->glBlobs.begin() ; itGLB != ssb->glBlobs.end() ; itGLB++ ){

            // For each grey-level blob, we create a Blob
            blobs.push_back( new surf::GreyLevelBlob() );
            surf::GreyLevelBlob *blob = blobs[blobs.size()-1];

            blob->ssb_parent = ssblob;
            blob->t = (*itGLB)->measurements.t;
            blob->scale = (*itGLB)->GetScale();

            // The surf::GreyLevelBlob's nodeslist contains its corresponding nodes indices on the
            //    mesh it was extracted from.
            std::set<SiteType<AimsSurface<3, Void> >::type,
                 ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> > listePoints
                     = (*itGLB)->GetListePoints();
            for ( itPoints = listePoints.begin() ; itPoints != listePoints.end() ; itPoints++ ) {
                (blob->nodes).insert( (*itPoints).second );

                if ( subject.coordinates == LATLON_2D ) {
                    (blob->coordinates)[(*itPoints).second] = vector<float>(2);

        			(blob->coordinates)[(*itPoints).second][0] = subject.lat->item((*itPoints).second);
                    (blob->coordinates)[(*itPoints).second][1] = subject.lon->item((*itPoints).second);
                }

                (blob->raw_coordinates)[(*itPoints).second] = vector<float>(3);
                (blob->raw_coordinates)[(*itPoints).second][0] = subject.mesh->vertex()[(*itPoints).second][0];
                (blob->raw_coordinates)[(*itPoints).second][1] = subject.mesh->vertex()[(*itPoints).second][1];
                (blob->raw_coordinates)[(*itPoints).second][2] = subject.mesh->vertex()[(*itPoints).second][2];
            }

            ssblob->blobs.insert(blob);

            if ( blob->scale < ssblob->tmin )
                ssblob->tmin = blob->scale;
            if ( blob->scale > ssblob->tmax )
                ssblob->tmax = blob->scale;
        }

        ssblob->t = ssb->GetMeasurements().t;

    }


    // Now that every SSB has been created, we can create links (bifurcations)
    // between them

    // BIFURCATIONS
    std::list<Bifurcation<SiteType<AimsSurface<3, Void> >::type> *> bifurcations = sketch.BifurcationList();
    std::list<Bifurcation<SiteType<AimsSurface<3, Void> >::type> *>::iterator bifurcIt;
    for ( bifurcIt = bifurcations.begin() ; bifurcIt != bifurcations.end() ; bifurcIt ++ ) {
        list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *> topBlobs, bottomBlobs, topBottomBlobs, bottomTopBlobs;
        topBlobs = (*bifurcIt)->TopBlobs();
        bottomBlobs = (*bifurcIt)->BottomBlobs();

        for (itSSBaux = topBlobs.begin() ; itSSBaux != topBlobs.end() ; itSSBaux++) {
            for (itSSBaux2 = bottomBlobs.begin() ; itSSBaux2 != bottomBlobs.end() ; itSSBaux2++){
                ssbMap[*itSSBaux2]->topBlobs.insert( ssbMap[*itSSBaux] );
                ssbMap[*itSSBaux]->bottomBlobs.insert( ssbMap[*itSSBaux2] );
            }
        }
    }

    std::cout << " blobs.size : " << blobs.size() <<
            " ssblobs.size : " << ssblobs.size() << std::endl;

}

//##############################################################################

void TextureToBlobs::GreyLevelBlobsFromTexture ( SubjectData &subject,
                    vector<surf::GreyLevelBlob *> &blobs,
                    vector<surf::ScaleSpaceBlob *> &ssblobs,
                    string blobsPath ) {

    FiniteElementSmoother<3, float> smooth ( 0.01, subject.mesh, subject.weightLapl );
    ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss ( subject.mesh, subject.tex, &smooth );
    ss.AddScale(1.0, *(subject.tex));

    if ( subject.coordinates == LATLON_2D || subject.coordinates == LAT_1D )
        storeCoordinatesInScaleSpace( subject, ss );

    cout << "══ Computing primal sketch..." << endl;
    TimeTexture<float> blobs_tex;
    PrimalSketch( subject, blobs, ssblobs, &ss, blobs_tex, 1.0 );

    Writer<TimeTexture<float> > wtrBlobs ( blobsPath );
    wtrBlobs.write ( blobs_tex );
}



//##############################################################################

//void TextureToBlobs::RecoverBlobsFromIndivGraph( Graph *graph,
//                            vector<surf::GreyLevelBlob *> &blobs,
//                            vector<surf::ScaleSpaceBlob *> &ssblobs,
//                            bool initNull ){
//    if (initNull){
//      blobs.clear();
//      ssblobs.clear();
//    }
//
//    set<Vertex *>::iterator iv;
//    Edge *e;
//    Vertex::iterator jv;
//    Edge::iterator kv;
//
//    vector<Vertex *> listVertSSB, listVertGLB;
//    map<int, set<int> > listGLBindices;
//    int iNbLinks = 0;
//    int iNbGLB = 0;
//    int iNbSSB = 0;
//
//    // A first pass to get the Grey-level Blobs
//    for (iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv){
//      if ((*iv)->getSyntax() == "glb"){
//        int index;
//        float scale, t;
//        vector<int> nodes_list;
//
//        blobs.push_back(new surf::GreyLevelBlob());
//        surf::GreyLevelBlob *blob = blobs[blobs.size()-1];
//
//        (*iv)->getProperty("scale", scale);
//        (*iv)->getProperty("nodes", nodes_list);
//        (*iv)->getProperty("t", t);
//
//        blob->index = iNbGLB++;
//        index = blob->index;
//        (*iv)->setProperty("index", (int) index );
////         (*iv)->getProperty("index", index );
////         cout << "idx: " << index << " " << flush;
//
//        blob->nodes = vector2set(nodes_list);
//        blob->scale = scale;
//        blob->t = t;
//        blob->ssb_parent = NULL;
//
//      }
//    }
//
//
//    // Another pass to get the Scale-space Blobs
//    for (iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv){
//      if ((*iv)->getSyntax() == "ssb"){
//        int index;
//        float tmax, tmin, t;
//
//        ssblobs.push_back(new surf::ScaleSpaceBlob());
//        surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size()-1];
//
//        (*iv)->getProperty("tmax", tmax);
//        (*iv)->getProperty("tmin", tmin);
//        (*iv)->getProperty("t", t);
//
//        ssblob->index = iNbSSB++;
//        index = ssblob->index;
//        (*iv)->setProperty("index", (int) index );
////         cout << "idx2:" << index << flush;
//
//        ssblob->tmax = tmax;
//        ssblob->tmin = tmin;
//        ssblob->t = t;
//        string sujet;
//        graph->getProperty("sujet", sujet);
//        ssblob->subject = sujet;
//
//        if (listGLBindices.find(index) == listGLBindices.end())
//          listGLBindices[index] = set<int>();
//        for (jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++){
//          e = *jv;
//          if (e->getSyntax() == "s2g"){
//            for (kv = e->begin() ; kv != e->end() ; kv++){
//              if ((*kv)->getSyntax() == "ssb"){
//
//              }
//              else if ((*kv)->getSyntax() == "glb"){
//                int iGLBindex;
//                (*kv)->getProperty("index", iGLBindex);
////                 cout << "igl:" << iGLBindex << " " << flush;
//                listGLBindices[index].insert(iGLBindex);
//                iNbLinks++;
//              }
//            }
//          }
//        }
//
//      }
//    }
//    cout << "Rebuilding links..." << endl;
//    for (uint i = ssblobs.size() - iNbSSB ; i < ssblobs.size() ; i++){
//      int index = ssblobs[i]->index;
//      set<int>::iterator it;
//      for (it = listGLBindices[index].begin() ; it != listGLBindices[index].end() ; it++){
//        ssblobs[i]->blobs.insert(blobs[blobs.size() - iNbGLB + *it]);
//        blobs[blobs.size() - iNbGLB + *it]->ssb_parent = ssblobs[i];
//
//      }
//    }
//
//
//    cout << iNbGLB << " blobs added" << endl;
//    cout << iNbSSB << " ssblobs added " << endl;
//    cout << iNbLinks << " links added" << listGLBindices.size() << endl;
//
//    cout << "Checking links..." << endl;
////    for (uint i = blobs.size() - iNbGLB ; i < blobs.size() ; i++){
////        assert(blobs[i]->ssb_parent != NULL);
////
////    }
//
//}

void TextureToBlobs::getGreyLevelBlobsFromIndividualGraph ( Graph *graph,
                            SubjectData &subject,
                            std::vector <surf::GreyLevelBlob *> &blobs,
                            bool initNull ){
    if ( initNull ){
        blobs.clear();
    }
    std::set <Vertex *>::iterator iv;

    int iNbGLB = 0;

    for ( iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv ) {
        if ( (*iv)->getSyntax() == "glb" ) {
            int index, label;
            float scale, t;
            vector<int> nodes_list;
            vector<float> latitudes, longitudes;

            blobs.push_back( new surf::GreyLevelBlob() );
            surf::GreyLevelBlob *blob = blobs[blobs.size()-1];

            (*iv)->getProperty( "scale", scale );

            (*iv)->getProperty( "t", t );
            (*iv)->getProperty( "label", label );
            (*iv)->getProperty( "nodes", nodes_list );
            (*iv)->getProperty( "x", latitudes );
            (*iv)->getProperty( "y", longitudes );
            assert(nodes_list.size() == latitudes.size());

            for ( uint i = 0 ; i < nodes_list.size() ; i++ ) {
                
                (blob->raw_coordinates)[nodes_list[i]] = vector<float>(3);
                (blob->raw_coordinates)[nodes_list[i]][0] = subject.mesh->vertex()[nodes_list[i]][0];
                (blob->raw_coordinates)[nodes_list[i]][1] = subject.mesh->vertex()[nodes_list[i]][1];
                (blob->raw_coordinates)[nodes_list[i]][2] = subject.mesh->vertex()[nodes_list[i]][2];

                if ( subject.coordinates == LATLON_2D ) {
                    ( blob->coordinates)[nodes_list[i]] = vector<float>(2);
                    ( blob->coordinates)[nodes_list[i]][0] = latitudes[i];
                    ( blob->coordinates)[nodes_list[i]][1] = longitudes[i];
                }

            }

            blob->index = iNbGLB++;
            blob->label = label;
            index = blob->index;
            (*iv)->setProperty( "index", (int) index );

            blob->nodes = vector2set(nodes_list);
            blob->scale = scale;
            blob->t = t;
            blob->ssb_parent = NULL;
        }
    }
}


void TextureToBlobs::getScaleSpaceBlobsFromIndividualGraph ( Graph *graph,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            map<int, set<int> > &listGLBindices,
                            bool initNull ){
    if ( initNull )
        ssblobs.clear();

    std::set< Vertex * >::iterator iv;
    Edge *e;
    Vertex::iterator jv;
    Edge::iterator kv;

    listGLBindices = std::map<int, std::set<int> >();
    int iNbLinks = 0;
    int iNbSSB = 0;

    for ( iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv ) {
        if ( (*iv)->getSyntax() == "ssb" ) {
            int index;
            float tmax,
                  tmin,
                  t;
            std::string subject_id, label;

            ssblobs.push_back( new surf::ScaleSpaceBlob() );
            surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size()-1];

            (*iv)->getProperty( "tmax", tmax);
            (*iv)->getProperty( "tmin", tmin);
            (*iv)->getProperty( "t", t);
            (*iv)->getProperty( "subject", subject_id);
            (*iv)->getProperty( "label", label);

            ssblob->index = iNbSSB++;
            index = ssblob->index;
            (*iv)->setProperty( "index", (int) index );
            (*iv)->setProperty( "sites_index", (int)(ssblobs.size()-1) );

            ssblob->tmax = tmax;
            ssblob->tmin = tmin;
            ssblob->t = t;
            ssblob->subject = subject_id;
            ssblob->label = atoi(label.data());

            if ( listGLBindices.find( index ) == listGLBindices.end() )
                listGLBindices[index] = set<int>();

            for ( jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++ ) {
                e = *jv;

                if ( e->getSyntax() == "s2g" ) {
                    for ( kv = e->begin() ; kv != e->end() ; kv++ ) {
                        if ( (*kv)->getSyntax() == "ssb" ) {

                        }
                        else if ( (*kv)->getSyntax() == "glb" ) {
                            int iGLBindex;
                            (*kv)->getProperty( "index", iGLBindex );
                            listGLBindices[index].insert( iGLBindex );
                            iNbLinks++;
                        }
                    }
                }
            }
        }
    }


}


void TextureToBlobs::RecoverBlobsFromIndividualGraph( Graph *graph,
                            SubjectData &subject,
                            std::vector<surf::GreyLevelBlob *> &blobs,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            bool initNull ){

    std::map <int, std::set<int> > listGLBindices;
    int iNbGLB = blobs.size();
    int iNbSSB = ssblobs.size();

    getGreyLevelBlobsFromIndividualGraph( graph, subject, blobs, initNull );
    iNbGLB = blobs.size() - iNbGLB;

    getScaleSpaceBlobsFromIndividualGraph( graph, ssblobs, listGLBindices, initNull );
    iNbSSB = ssblobs.size() - iNbSSB;

    for ( uint i = ssblobs.size() - iNbSSB ; i < ssblobs.size() ; i++ ) {
        int index = ssblobs[i]->index;
        std::set<int>::iterator it;
        for ( it = listGLBindices[index].begin() ; it != listGLBindices[index].end() ; it++ ) {
            ssblobs[i]->blobs.insert( blobs[blobs.size() - iNbGLB + *it] );
            blobs[blobs.size() - iNbGLB + *it]->ssb_parent = ssblobs[i];

        }
    }

    cout << iNbGLB << " blobs added" << endl;
    cout << iNbSSB << " ssblobs added " << endl;

    for (uint i = blobs.size() - iNbGLB ; i < blobs.size() ; i++)
        assert(blobs[i]->ssb_parent != NULL);
}


void TextureToBlobs::RecoverBlobsFromGLBOnly( Graph *graph,
                            SubjectData &subject,
                            std::vector<surf::GreyLevelBlob *> &blobs,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            bool initNull ){

    int iNbGLB = blobs.size();
    getGreyLevelBlobsFromIndividualGraph ( graph, subject, blobs, initNull );
    iNbGLB = blobs.size() - iNbGLB;

    int iNbSSB = 0;

    // We Add A Fictitious Ssb Per Glb
    for ( uint i = 0 ; i < iNbGLB ; i++ ) {
        ssblobs.push_back( new surf::ScaleSpaceBlob() );
        surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size()-1];

        ssblob->index = iNbSSB++;
        ssblob->t = blobs[blobs.size() - iNbGLB + i]->t;
        ssblob->subject = subject.subject_id;

        ssblob->blobs.insert( blobs[blobs.size() - iNbGLB + i] );
        ssblob->getNodesFromBlob( * (ssblob->blobs.begin()) );
        blobs[blobs.size() - iNbGLB + i]->ssb_parent = ssblob;
    }

    assert( iNbSSB == iNbGLB );

    cout << iNbGLB << " blobs added" << endl;
    cout << iNbSSB << " ssblobs added " << endl;

    for (uint i = blobs.size() - iNbGLB ; i < blobs.size() ; i++)
        assert(blobs[i]->ssb_parent != NULL);
}

////##############################################################################


// Creates an Aims Graph for only ONE subject with scale-space blobs, grey-level
//  blobs and links between both types

void TextureToBlobs::AimsGraph (   Graph *graph,
								   SubjectData & subject,
                                   vector<surf::Blob *> &blobs ) {

    // From the two vectors, we build a graph containing the ssb, the glb and
    //  the links between ssb and glb
    vector<float> resolution, bbmin2D, bbmax2D;
    vector<int> bbmin, bbmax;
    resolution.push_back(1.0); resolution.push_back(1.0); resolution.push_back(1.0);

    bbmin.push_back(-10); bbmin.push_back(-10); bbmin.push_back(-10);
    bbmax.push_back(10); bbmax.push_back(10); bbmax.push_back(10);
    graph->setProperty( "filename_base", "*");

    graph->setProperty("voxel_size", resolution);
    graph->setProperty("boundingbox_min", bbmin);
    graph->setProperty("boundingbox_max", bbmax);
    graph->setProperty("mesh", subject.paths.meshPath);
    graph->setProperty("subject", subject.subject_id);
    graph->setProperty("texture", subject.paths.texPath);
    if ( subject.paths.latPath != "" )
        graph->setProperty("latitude", subject.paths.latPath);
    if ( subject.paths.lonPath != "" )
        graph->setProperty("longitude", subject.paths.lonPath);

    vector<set<int> > nodes_lists;
    Vertex *vert;
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;

    // Initializing blobs indices
//    for (int i = 0 ; i < (int) blobs.size() ; i++)
//        blobs[i]->index = i;

    cout << "════ Extracting meshes for the grey-level blobs..." << endl;

	cout << "      ░░░ mode Blobs Meshes From Mesh ░░░     " << endl;
    for ( uint i = 0 ; i < blobs.size() ; i++ )
        blobs[i]->getAimsSphereAtMaxNode ( *(subject.tex) );

    cout << "════ Adding blobs..." << endl;

    for ( int i = 0 ; i < (int) blobs.size() ; i++ ) {

        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index

        cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
        vert = graph->addVertex("glb");

        vert->setProperty( "index", blobs[i]->index );
        vert->setProperty( "nodes", blobs[i]->nodes );
        vert->setProperty( "label", "0");

        if ( blobs[i]->coordinates.size() != 0 ) {
            vector<float> latitudes, longitudes;
            set<int>::iterator it;
            for ( it = blobs[i]->nodes.begin() ; it != blobs[i]->nodes.end() ; it++ ) {
                latitudes.push_back( blobs[i]->coordinates[*it][0] );
                if ( blobs[i]->coordinates[*it].size() >= 2 )
                    longitudes.push_back( blobs[i]->coordinates[*it][1] );
            }
            vert->setProperty( "x", latitudes );
            if ( latitudes.size() == longitudes.size() ) {
                vert->setProperty( "y", longitudes );
            }
            else {
                assert(longitudes.size() == 0);
            }
        }

        // We associate the proper mesh patch from "objects" to the vertex
        ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
        (*ptr)[0]=blobs[i]->mesh;
        manip.storeAims(*graph, vert, "glb", ptr);
    }
    cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added in total (SSB and GLB)" << endl;



}

//##############################################################################

// Creates an Aims Graph for only ONE subject with scale-space blobs, grey-level
//  blobs and links between both types
//
void TextureToBlobs::AimsGraph (   Graph *graph,
								   SubjectData & subject,
                                   vector<surf::GreyLevelBlob *> &blobs,
                                   vector<surf::ScaleSpaceBlob *> &ssblobs ) {
    // From the two vectors, we build a graph containing the ssb, the glb and
    //  the links between ssb and glb
    vector<float> resolution, bbmin2D, bbmax2D;
    vector<int> bbmin, bbmax;
    resolution.push_back(1.0); resolution.push_back(1.0); resolution.push_back(1.0);

    bbmin.push_back(-10); bbmin.push_back(-10); bbmin.push_back(-10);
    bbmax.push_back(10); bbmax.push_back(10); bbmax.push_back(10);
    graph->setProperty( "filename_base", "*");

    graph->setProperty("voxel_size", resolution);
    graph->setProperty("boundingbox_min", bbmin);
    graph->setProperty("boundingbox_max", bbmax);
    graph->setProperty("mesh", subject.paths.meshPath);
    graph->setProperty("subject", subject.subject_id);
    graph->setProperty("texture", subject.paths.texPath);

    if ( subject.paths.latPath != "" )
        graph->setProperty("latitude", subject.paths.latPath);
    if ( subject.paths.lonPath != "" )
        graph->setProperty("longitude", subject.paths.lonPath);

    vector<set<int> > nodes_lists;
    vector< pair<surf::GreyLevelBlob *, surf::GreyLevelBlob *> > blobsIndices;
//    TimeTexture<float> lat, lon;
    Vertex *vert;
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    vector<Vertex *> listVertSSB( ssblobs.size() ), listVertGLB( blobs.size() );

    // Initializing blobs indices
    for (int i = 0 ; i < (int) ssblobs.size() ; i++)
        ssblobs[i]->index = i;
    for (int i = 0 ; i < (int) blobs.size() ; i++)
        blobs[i]->index = i;


    // Let's add the scale-space blobs
    cout << "════ Adding scale-space blobs..." << endl;

    for ( int i = 0 ; i < (int) ssblobs.size() ; i++ ) {

        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index

        cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
        vert = graph->addVertex("ssb");

        vert->setProperty( "label", "0");
        vert->setProperty( "t", ssblobs[i]->t);
        vert->setProperty( "subject", ssblobs[i]->subject);
        vert->setProperty( "tmin", ssblobs[i]->tmin);
        vert->setProperty( "tmax", ssblobs[i]->tmax);
        listVertSSB[  i  ] = vert;
    }
    cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added... done" << endl;



    // Let's add the grey-level blobs
    cout << "════ Extracting meshes for the grey-level blobs..." << endl;

    std::cout << "      ░░░ mode AimsSphereAtMaxNode ░░░     " << std::endl;
    for ( uint i = 0 ; i < blobs.size() ; i++ )
//        blobs[i]->getAimsEllipsoidAtMaxNode ( *(subject.tex) );
        blobs[i]->getAimsSphereAtMaxNode( *(subject.tex), 0.3);


    cout << "════ Adding grey-level blobs..." << endl;

    for ( int i = 0 ; i < (int) blobs.size() ; i++ ) {

        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index

        std::cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << std::flush ;
        vert = graph->addVertex("glb");

        vert->setProperty( "t", blobs[i]->t);
        vert->setProperty( "scale", blobs[i]->scale );
        vert->setProperty( "subject", blobs[i]->ssb_parent->subject );
        vert->setProperty( "nodes", blobs[i]->nodes );

        if ( subject.coordinates == LATLON_2D ) {
            std::vector<float> latitudes, longitudes;
            std::set<int>::iterator it;
            for ( it = blobs[i]->nodes.begin() ; it != blobs[i]->nodes.end() ; it++ ) {
                latitudes.push_back( (float) blobs[i]->coordinates[*it][0] );
                if ( blobs[i]->coordinates[*it].size() >= 2 )
                    longitudes.push_back( (float) blobs[i]->coordinates[*it][1] );
            }
            vert->setProperty( "x", latitudes );
            
            if ( latitudes.size() == longitudes.size() ) {
                vert->setProperty( "y", longitudes );
            }
            else {
                if ( longitudes.size() != 0 ) {
                    cout << latitudes.size() << " " << longitudes.size() << endl;
                }
                assert(longitudes.size() == 0);
            }
        }
        else if ( subject.coordinates == LAT_1D ) {
            //TODO
            cout << "To Be Implemented" << endl;

        }

        listVertGLB[ i ] = vert;

        // We associate the proper mesh patch from "objects" to the vertex
        ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
        (*ptr)[0] = blobs[i]->mesh;
        manip.storeAims(*graph, vert, "glb", ptr);

    }
    cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added in total (SSB and GLB)" << endl;


    // Let's add the links between scale-space and grey-level blobs
    cout << "════ Adding links between SSB and GLB (GLB forming SSB)..." << endl;

    uint iNbLinks=0;
    for (int i = 0 ; i < (int) ssblobs.size() ; i++) {

        set<surf::GreyLevelBlob *>::iterator itB1;
        set<surf::GreyLevelBlob *> &listGLB = ssblobs[i]->blobs;

        for (itB1 = listGLB.begin(); itB1 != listGLB.end() ; itB1++) {

            Vertex *v1, *v2;

            v1 = listVertSSB[ i ];
            v2 = listVertGLB[(*itB1)->index];

            graph->addEdge(v1, v2, "s2g");

            iNbLinks++;
        }
    }
    cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " links added... done" << endl;

    // Let's add the links between grey-level blobs
    cout << "════ Collecting the list of GLB relations..." << endl;
    vector< pair<surf::GreyLevelBlob *, surf::GreyLevelBlob *> > blobsPairs;

    for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
        set<surf::GreyLevelBlob *>::iterator itB;
        set<surf::GreyLevelBlob *> &unsortedListGLB = ssblobs[i]->blobs;
        set<surf::GreyLevelBlob *, ltBlobs> listGLB;
        set<surf::GreyLevelBlob *, ltBlobs>::iterator itB1, itB2;
        for ( itB = unsortedListGLB.begin() ; itB != unsortedListGLB.end() ; itB ++ )
            listGLB.insert(*itB);
        ASSERT( unsortedListGLB.size() == listGLB.size() );
        itB1 = listGLB.begin();
        itB2 = itB1;
        if ( itB2 != listGLB.end() )
            itB2++;
        else
            ASSERT(false);

        while ( itB2 != listGLB.end() ) {
            surf::GreyLevelBlob *glb1, *glb2;
            glb1 = (*itB1);
            glb2 = (*itB2);
            pair<surf::GreyLevelBlob *, surf::GreyLevelBlob *> p;
            p.first = glb1;
            p.second = glb2;
            blobsPairs.push_back(p);
            itB1++, itB2++;
        }
    }
    cout << "  " << blobsPairs.size() << " paires de blobs" << endl;

    cout << "════ Extracting meshes for the GLB relations..." << endl;
    AimsSurfaceTriangle *relations = new AimsSurfaceTriangle();
    *relations = getG2GRelationsMeshes( blobsPairs, NODES_BARYCENTERS );

    cout << "════ Adding links between GLB from same SSB (with corresponding meshes for visualization)..." << endl;

    iNbLinks=0;
    for ( uint i = 0 ; i < blobsPairs.size() ; i++ ) {
        cout << "\b\b\b\b\b\b\b\b\b\b\b" << i << flush;
        // We associate the proper mesh patch from "objects" to the vertex

        Vertex *v1, *v2;
        surf::GreyLevelBlob *glb1, *glb2;
        glb1 = blobsPairs[i].first;
        glb2 = blobsPairs[i].second;
        v1 = listVertGLB[glb1->index];
        v2 = listVertGLB[glb2->index];
        Edge *e = graph->addEdge(v1, v2, "g2g");

        ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
        (*ptr)[0]=(*relations)[i];
        manip.storeAims(*graph, e, "g2g", ptr);

        iNbLinks++;
    }
    cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " links added... done" << endl;


    // Let's add the links between scale-space blobs (bifurcations)
    cout << "════ Collecting the list of existing bifurcations..." << endl;
    vector< pair< surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> > bifurcPairs;

    for (int i = 0 ; i < (int) ssblobs.size() ; i++) {

        set<surf::ScaleSpaceBlob *>::iterator itB1;
        set<surf::ScaleSpaceBlob *> &listSSB = ssblobs[i]->topBlobs;
        for (itB1 = listSSB.begin(); itB1 != listSSB.end() ; itB1++) {

            if ( ssblobs[i]->index < (*itB1)->index ) {
                pair<surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> p;
                p.first = ssblobs[i];
                p.second = (*itB1);
                bifurcPairs.push_back(p);
            }
        }
        listSSB = ssblobs[i]->bottomBlobs;
        for (itB1 = listSSB.begin() ; itB1 != listSSB.end() ; itB1++) {

            if ( ssblobs[i]->index < (*itB1)->index ) {
                pair<surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> p;
                p.first = ssblobs[i];
                p.second = (*itB1);
                bifurcPairs.push_back(p);
            }
        }
    }

    cout << "════ Extracting meshes for the bifurcations relations..." << endl;
    AimsSurfaceTriangle *bifurcations = new AimsSurfaceTriangle();
    *bifurcations = getBifurcationRelationsMeshes( bifurcPairs, NODES_BARYCENTERS );
    cout << "  " << bifurcPairs.size() << " paires de bifurcations" << endl;
    cout << "  " << bifurcations->size() << " bifurcations meshes" << endl;

    cout << "════ Adding links between SSB (bifurcations)..." << endl;
    iNbLinks=0;

    for ( uint i = 0 ; i < bifurcPairs.size() ; i++ ) {

        Vertex *v1, *v2;
        ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
        (*ptr)[0] = (*bifurcations)[i];
        surf::ScaleSpaceBlob *ssb1, *ssb2;
        ssb1 = bifurcPairs[i].first;
        ssb2 = bifurcPairs[i].second;
        v1 = listVertSSB[ ssb1->index ];
        v2 = listVertSSB[ ssb2->index ];
        Edge *e = graph->addEdge(v1, v2, "bifurcation");
        manip.storeAims(*graph, e, "bifurcation", ptr);
        e->setProperty("type", "test");
        iNbLinks++;

    }
    cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " bifurcations added... done" << endl;


}

//##############################################################################

void TextureToBlobs::AimsGroupGraph ( Graph *graph,
                        map<string, SubjectData *> data,
//                       vector<Blob *> &blobs,
                            vector<surf::ScaleSpaceBlob *> &ssblobs,
                            vector<surf::SSBClique> &cliques ) {

    std::cerr << "Construction du graphe de groupe..." << std::endl;
    std::vector<float> resolution,bbmin2D,bbmax2D;
    std::vector<int> bbmin, bbmax;
    resolution.push_back(1.0); resolution.push_back(1.0); resolution.push_back(1.0);
    bbmin.push_back(-10); bbmin.push_back(-10); bbmin.push_back(-10);
    bbmax.push_back(10); bbmax.push_back(10); bbmax.push_back(10);
    graph->setProperty( "filename_base", "*");
    std::map <std::string, SubjectData *>::iterator it;


    std::vector<std::string> listTexPaths, listMeshPaths, listLatPaths, listLonPaths, listIndivGraphPaths, listSubjects;
    std::map <std::string,int> subj_names;
    uint i = 0;
    for ( it = data.begin() ; it != data.end() ; it ++ ) {
        SubjectData *subject = (*it).second;
        subj_names[(*it).first] = i++;
        listTexPaths.push_back(subject->paths.texPath);
        listMeshPaths.push_back(subject->paths.meshPath);
        listLatPaths.push_back(subject->paths.latPath);
        listLonPaths.push_back(subject->paths.lonPath);
        listIndivGraphPaths.push_back(subject->paths.graphPath);
        listSubjects.push_back(subject->subject_id);
    }

    graph->setProperty( "voxel_size", resolution );
    graph->setProperty( "boundingbox_min", bbmin );
    graph->setProperty( "boundingbox_max", bbmax );
    graph->setProperty( "meshes", listMeshPaths );
    graph->setProperty( "subjects", listSubjects );
    graph->setProperty( "textures", listTexPaths );
    graph->setProperty( "latitudes", listLatPaths );
    graph->setProperty( "longitudes", listLonPaths );
    graph->setProperty( "indiv_graphs", listIndivGraphPaths );

    Vertex *vert;
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    std::vector<Vertex *> listVertSSB( ssblobs.size() ), listVertGLB( ssblobs.size() );

    for ( int i = 0 ; i < (int) ssblobs.size() ; i++ ) {

        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index

        cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
        vert = graph->addVertex("ssb");
        ssblobs[i]->index = i;
        vert->setProperty( "index", i );
        vert->setProperty( "label", "0" );

        vert->setProperty( "t", ssblobs[i]->t ); // subj_names[ssblobs[i]->subject]);
        vert->setProperty( "subject", ssblobs[i]->subject );
        vert->setProperty( "tmin", ssblobs[i]->tmin );
        vert->setProperty( "tmax", ssblobs[i]->tmax );
        vert->setProperty( "tValue", 100.0 );
        vert->setProperty( "nodes", set2vector((*(ssblobs[i]->blobs.begin()))->nodes) );

         // We associate the proper mesh patch from "objects" to the vertex
         ptr = carto::rc_ptr<AimsSurfaceTriangle>( new AimsSurfaceTriangle );
         (*ptr)[0] = ssblobs[i]->mesh;
         manip.storeAims( *graph, vert, "ssb", ptr );
         vert->setProperty( "ssb_label", i );
         listVertSSB[ i ] = vert;
    }

    std::cout << " nodes " << std::endl << "════ Extracting meshes for the interblobs relations..." << std::endl;
    AimsSurfaceTriangle *relations = new AimsSurfaceTriangle();
    *relations = getB2BRelationsMeshes( cliques, NODES_BARYCENTERS );

    std::cout << "Building cliques in the Aims group graph..." << endl;
    for ( uint i = 0 ; i < cliques.size() ; i++ ) {

        // For every clique, we get the two corresponding vertices from listVertices
        Vertex *v1, *v2;
        int aux1, aux2;
        for ( aux1 = 0 ;
              aux1 < (int) ssblobs.size() &&
                 !(cliques[i].ssb1->index == ssblobs[aux1]->index &&
                     cliques[i].ssb1->subject == ssblobs[aux1]->subject );
             aux1++ )
        {  }
        for ( aux2 = 0 ;
              aux2 < (int) ssblobs.size() &&
                 !(cliques[i].ssb2->index == ssblobs[aux2]->index &&
                 cliques[i].ssb2->subject == ssblobs[aux2]->subject ) ;
             aux2++ )
        {  }
        assert( aux1 != ssblobs.size() && aux2 != ssblobs.size() );
        v1 = listVertSSB[aux1];
        v2 = listVertSSB[aux2];
        Edge *edge = graph->addEdge( v1, v2, "b2b" );

        edge->setProperty( "similarity", cliques[i].similarity );
        edge->setProperty( "distance", cliques[i].distance );

        ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
        (*ptr)[0]=(*relations)[i];
        manip.storeAims( *graph, edge, "b2b", ptr );
        edge->setProperty( "b2b_label", i );
    }
}

void TextureToBlobs::ReadAimsGroupGraph (   Graph &graph,
                                            GroupData &data,
                                            vector<surf::ScaleSpaceBlob *> &ssblobs,
                                            vector<surf::SSBClique> &cliques,
                                            vector<Vertex *> &listVertex   ) {

    std::cout << "Recovering the data..." << std::endl;

    std::set<Vertex *>::iterator iv;
    std::string sujet;
    int index;
    int newindex=0;

    // Recovering The Scale-Space Blobs

    std::cout << " Recovering the scale-space blobs..." << std::endl;
    std::map<int, std::set<int> > listGLBindices;

    TextureToBlobs::getScaleSpaceBlobsFromIndividualGraph ( &graph, ssblobs, listGLBindices, false );

    std::cout << "ssblobs.size() :"<< ssblobs.size() << std::endl;

    // Recovering The Clique Links...
    Edge *e;
    Vertex::iterator jv;
    Edge::iterator kv;

    std::cout << " Recovering the similarity cliques..." << std::endl;

    for ( iv = graph.vertices().begin() ; iv != graph.vertices().end() ; ++iv ) {

        if ( (*iv)->getSyntax() == "ssb" ){

            (*iv)->getProperty( "sites_index", index );
            (*iv)->getProperty( "subject", sujet );

            // Iterating On The Edges Of The Current Vertex
            for ( jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++ ) {
                e = *jv;

                if ( e->getSyntax() == "b2b" ){

                    float similarity;
                    e->getProperty( "similarity", similarity );

                    for ( kv = e->begin() ; kv != e->end() ; kv++ ){

                        if ( (*kv)->getSyntax() == "ssb" ){

                            int indexB2;
                            string sujetB2;
                            (*kv)->getProperty( "sites_index", indexB2 );
                            (*kv)->getProperty( "subject", sujetB2 );

                            if ( !(indexB2 == index && sujetB2 == sujet) ){

                                int blobs_index1, blobs_index2;
                                (*iv)->getProperty( "sites_index", blobs_index1 );
                                (*kv)->getProperty( "sites_index", blobs_index2 );

                                assert ( (index==blobs_index1 && indexB2 == blobs_index2) ||
                                  (index == blobs_index2 && indexB2 == blobs_index1));

                                if ( blobs_index1 < blobs_index2 ){
                                  cliques.push_back( surf::SSBClique( ssblobs[ blobs_index1 ] , ssblobs[ blobs_index2 ], similarity) );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "ssbcliques.size() :"<< cliques.size() << std::endl << std::endl;
}

//set<int> TextureToBlobs::getFilteringNodes( SubjectData & subject ) {
//
//    set<int> nodes;
//
//    if ( subject.coordinates == LAT_1D ) {
//
//        cout << "Filtering : using only a latitude texture" << endl;
//        // TODO
//        cout << "TO BE REIMPLEMENTED..." << endl;
//        for ( uint i = 0 ; i < subject.mesh->vertex().size() ; i++ )
//            nodes.insert(i);
//
//    }
//    else if ( subject.coordinates == LATLON_2D ) {
//
//        cerr << "Filtering : using a latitude and a longitude textures" << endl;
////        set<uint>::iterator it;
////        for ( it = region.nodes )
////            nodes.insert( subjData.gyrusVertices[i] );
//
//    }
//
//    return nodes;
//}
//
//
//void TextureToBlobs::filterBlobs (  SubjectData & subject,
//									vector<surf::ScaleSpaceBlob *> &ssblobs,
//                                    vector<surf::GreyLevelBlob *> & filteredBlobs,
//                                    vector<surf::ScaleSpaceBlob *> & filteredSsblobs ) {
//
////    if ( subjData.getFilterMode() == NO_FILTER ) {
//        cerr << "No coordinates textures provided" << endl;
//        // if no coordinates textures, that means no bounding box specified,
//        //   then no filtering
//        set<int> filteringNodes;
//        for ( uint i = 0 ; i < subject.mesh->vertex().size() ; i++ )
//            filteringNodes.insert(i);
//        filteringBlobs ( ssblobs, filteredBlobs, filteredSsblobs, filteringNodes );
//
////    }
////    else if ( subjData.getFilterMode() == GYRUS ){
////
////        cout << " Filtering by gyrus..." << endl;
////        set<int> filteringNodes = getFilteringNodes( );
////        cout<< filteringNodes.size() << " filtering nodes" << endl;
////
////        filteringBlobs ( ssblobs, filteredBlobs, filteredSsblobs, filteringNodes );
////
////    }
//
//}
