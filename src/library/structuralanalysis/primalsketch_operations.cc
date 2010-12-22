

#include <cortical_surface/structuralanalysis/primalsketch_operations.h>


void storeCoordinatesInScaleSpace ( SubjectData &regionData, ScaleSpace<AimsSurface<3, Void>, Texture<float> >  &ss ) 
{
    std::vector<Point3df> *coordinates;
    coordinates = new vector<Point3df>();

    if ( regionData.coordinates == LATLON_2D ) 
    {
        for ( uint i = 0 ; i < regionData.lat->nItem() ; i++ )
            (*coordinates).push_back(
                Point3df( regionData.lat->item(i), regionData.lon->item(i), i ) );

    }
    else 
    {
        for ( uint i = 0 ; i < regionData.lat->nItem() ; i++ )
            (*coordinates).push_back(
                Point3df( regionData.lat->item(i), -1.0, i ) );

    }
    ss.PutCoordinates(coordinates);
}

void TextureToBlobs::PrimalSketchRegionMode (   //std::vector<surf::GreyLevelBlob *> &blobs,
                                                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                                surf::Region &region,
                                                SubjectData &regionData,
                                                std::string scaleSpacePath, string blobsPath,
                                                bool recover,
                                                float scale_max ) 
{

    // Creating the smoother...
    std::cout << std::endl << "  ══ Smoother creation... " << std::endl;
    FiniteElementSmoother<3, float> smooth ( 0.01, regionData.mesh, regionData.weightLapl );
    ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss ( regionData.mesh, regionData.tex, &smooth );

    if ( regionData.coordinates == LATLON_2D || regionData.coordinates == LAT_1D ) 
    {
        storeCoordinatesInScaleSpace( regionData, ss );
    }


    if ( recover ) 
    {
        TimeTexture<float> scale_space;
        // Extracting The Gyrus-Specific Data From The Scale-Space Texture
        Reader< TimeTexture<float> > rdrScaleSpace ( scaleSpacePath ) ;
        rdrScaleSpace.read ( scale_space );
        assert( scale_space.size() > 1 );
        TimeTexture<float> regionScaleSpace;
        for ( uint i = 1 ; i < scale_space.size() - 1 ; i++ )
            regionScaleSpace[i-1] = region.getLocalFromGlobalTexture( scale_space[i] );

        // Recovering Previously Computed Scale-Space
        std::cerr << "      Recovering scale-space " << std::flush;
        ss.uploadPreviouslyComputedScaleSpace(regionScaleSpace);

    }
    else 
    {

        // Generating A Scale-Space And Writing It Onto Hard Disk
        std::cerr << "   Computing scale-space " << std::endl;

        ss.GenerateDefaultScaleSpace( 128.0 );

        std::cout << "══ Writing scale-space (after rebuilding the whole texture from the gyrus-specific data)..." << std::endl;
        TimeTexture<float> regionScaleSpace = ss.getScaleSpaceTexture( );

        TimeTexture<float> scale_space;
        for ( uint i = 0 ; i < regionScaleSpace.size() ; i++ )
            scale_space[i] = region.getGlobalFromLocalTexture( regionScaleSpace[i] );

        Writer<TimeTexture<float> > wtrScaleSpace ( scaleSpacePath );
        wtrScaleSpace.write ( scale_space );
    }

    TimeTexture<float> regionBlobsTex;

    TextureToBlobs::PrimalSketch ( regionData, ssblobs, &ss, regionBlobsTex, scale_max );
    std::vector<surf::GreyLevelBlob *> blobs = TextureToBlobs::recoverGreyLevelBlobs ( ssblobs );
    // Update Coordinates With Conversion Gyrus/Global Indices
    for ( uint i = 0 ; i < blobs.size() ; i ++ ) {
        std::set<int>::iterator it;
        std::set<int> nodes_aux (blobs[i]->nodes);
        std::map<int, vector<float> > coordinates_aux (blobs[i]->coordinates);
        std::map<int, vector<float> > raw_coordinates_aux (blobs[i]->raw_coordinates);

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
    std::cout << "══ Writing blobs texture (after rebuilding the whole texture from the gyrus-specific data)..." << std::endl;
    for ( uint i = 0 ; i < regionBlobsTex.size() ; i++ )
        blobs_tex[i] = region.getGlobalFromLocalTexture( regionBlobsTex[i] );
    Writer<TimeTexture<float> > wtrBlobs ( blobsPath );
    wtrBlobs.write ( blobs_tex );


}

void TextureToBlobs::PrimalSketchGlobalMode (   //std::vector<surf::GreyLevelBlob *> &blobs,
                                    std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                    SubjectData &subject,
                                    std::string scaleSpacePath,
                                    std::string blobsPath,
                                    bool recover,
                                    float scale_max ) 
{
    FiniteElementSmoother<3, float> smooth ( 0.01, subject.mesh, subject.weightLapl );
    ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss ( subject.mesh, subject.tex, &smooth );

    if ( subject.coordinates == LATLON_2D || subject.coordinates == LAT_1D ) 
    {
        storeCoordinatesInScaleSpace( subject, ss );
    }
    std::cout << subject.mesh->vertex().size() << " " << subject.tex->nItem() << std::endl;

    if ( recover ) 
    {
        TimeTexture<float> scale_space;

        aims::Reader< TimeTexture<float> > rdrScaleSpace ( scaleSpacePath ) ;
        rdrScaleSpace.read ( scale_space );

        // Recovering Previously Computed Scale-Space
        // cout << "  Checking that the scale-space has more than just one texture..." << flush;
        assert( scale_space.size() > 1 );
        //cout << "OK (" << scale_space.size() << ") ( tex[0] is supposed to be the original texture (scale : 0) )" << std::endl;
        //uint nb_scales = scale_space.size() - 1;
        TimeTexture<float> scaleSpaceTex;
        //cout << "Filtering in the first " << nb_scales << " scales (to allow using bigger pre-computed scale-spaces)..." << flush;
        for ( uint i = 1 ; i < scale_space.size() ; i++ )
            scaleSpaceTex[i-1] = scale_space[i];

        //cout << "OK" << std::endl;
        std::cerr << "      Recovering scale-space " << std::endl;
        ss.uploadPreviouslyComputedScaleSpace(scaleSpaceTex);
    }
    else 
    {
        // Generating A Scale-Space And Writing It Onto Hard Disk
        std::cerr << "══ Computing scale-space " << std::endl;
        ss.GenerateDefaultScaleSpace( 128.0 );
        TimeTexture<float> scale_space;
        std::cout << "══ Writing scale-space..." << std::endl;
        scale_space = ss.getScaleSpaceTexture( );
        aims::Writer<TimeTexture<float> > wtrScaleSpace ( scaleSpacePath );
        wtrScaleSpace.write ( scale_space );
    }

    std::cout << "══ Computing primal sketch..." << std::endl;
    TimeTexture<float> blobs_tex;
    TextureToBlobs::PrimalSketch ( subject, ssblobs, &ss, blobs_tex, scale_max );
    std::cout << "══ Writing blobs texture..." << std::endl;
    aims::Writer<TimeTexture<float> > wtrBlobs ( blobsPath );
    wtrBlobs.write ( blobs_tex );
}

void TextureToBlobs::PrimalSketch ( SubjectData &subject,
                    //std::vector<surf::GreyLevelBlob *> &blobs,
                    std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                    ScaleSpace<AimsSurface<3, Void>, Texture<float> > *ss,
                    TimeTexture<float> &blobs_texture,
                    float scale_max ) 
{
    // Constructing A Primal-Sketch
    aims::PrimalSketch<AimsSurface<3, Void>, Texture<float> > sketch ( subject.subject_id, SURFACE );
    std::cout << "  ══ Setting scale-space..." << std::endl;
    sketch.SetScaleSpace(ss);

    // Scale_max == -1.0 Enforces Scale_max To Be Autodetermined (Max Scale From Ss)
    if ( scale_max == -1.0 ) 
    {
        set<float>::iterator it;
        set<float> scales = ss->GetScaleList();
        scale_max = 0.0;
        for (it = scales.begin() ; it != scales.end() ; it ++ )
            if (*it > scale_max)
                scale_max = *it;
    }

    std::cout << "  ══ Computing primal sketch..." << std::endl;
    // Launching The Computation Of The PS (tmin, tmax, statfile, intersection_criterium)
    sketch.ComputePrimalSketch ( 1.0, scale_max, "", 1 );

    // Writing A Blob Texture
    std::cerr << "      Getting scale-space blob texture... " << std::endl;
    blobs_texture = GetSSBlobTexture ( & sketch );

    // Getting The Blobs From The Primal Sketch Structure
    std::cerr << "      Blobs vectors construction..." << std::endl;
    getBlobsFromPrimalSketch ( subject, sketch, ssblobs );

    //std::cout << blobs.size() << " grey-level blobs / " << ssblobs.size() << " scale-space blobs" << std::endl;
}

//##############################################################################

void TextureToBlobs::GreyLevelBlobsFromTexture ( SubjectData &subject,
                    //vector<surf::GreyLevelBlob *> &blobs,
                    vector<surf::ScaleSpaceBlob *> &ssblobs,
                    string blobsPath ) 
{
    FiniteElementSmoother<3, float> smooth ( 0.01, subject.mesh, subject.weightLapl );
    ScaleSpace<AimsSurface<3, Void>, Texture<float> > ss ( subject.mesh, subject.tex, &smooth );
    ss.AddScale(1.0, *(subject.tex));

    if ( subject.coordinates == LATLON_2D || subject.coordinates == LAT_1D )
        storeCoordinatesInScaleSpace( subject, ss );

    std::cout << "   Computing primal sketch..." << std::endl;
    TimeTexture<float> blobs_tex;
    TextureToBlobs::PrimalSketch( subject, ssblobs, &ss, blobs_tex, 1.0 );

    Writer<TimeTexture<float> > wtrBlobs ( blobsPath );
    wtrBlobs.write ( blobs_tex );
}

//##############################################################################

// Function that builds a collection of surf::GreyLevelBlob and surf::ScaleSpaceBlob
// objects from a previously computed Primal Sketch

void TextureToBlobs::getBlobsFromPrimalSketch ( SubjectData & subject,
                     aims::PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch,
                     //std::vector<surf::GreyLevelBlob *> &blobs,
                     std::vector<surf::ScaleSpaceBlob *> &ssblobs ) 
{
    std::vector<surf::GreyLevelBlob *> blobs;
    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*> listBlobs
         = sketch.BlobSet();

    std::list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*>::iterator itSSB, itSSBaux, itSSBaux2;
    std::list<GreyLevelBlob<SiteType<AimsSurface<3, Void> >::type > *>::iterator itGLB;

    ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *ssb;
    std::set< SiteType<AimsSurface<3, Void> >::type,
        ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> >::iterator itPoints;

    std::map< ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*, surf::ScaleSpaceBlob * > ssbMap;

    for ( itSSB = listBlobs.begin() ; itSSB != listBlobs.end() ; itSSB++ )
    {
        // For each scale-space blob, we create a surf::ScaleSpaceBlob in "ssblobs" containing
        //    various surf::GreyLevelBlob objects (being themselves contained in a general resulting
        // "blobs" vector).
        ssb = *itSSB;

        ssblobs.push_back(new surf::ScaleSpaceBlob());
        surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size() - 1];
        ssblob->subject = sketch.Subject();
        ssblob->tmin = 999.0;
        ssblob->tmax = -999.0;
        ssblob->scales = std::set<float>();        

        // We save a link between the pointer and the corresponding index in ssblobs
        std::pair< ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*, surf::ScaleSpaceBlob * > p;
        p.first = ssb;
        p.second = ssblob;
        ssbMap.insert(p);

        for ( itGLB = ssb->glBlobs.begin() ; itGLB != ssb->glBlobs.end() ; itGLB++ ) 
        {
            // For each grey-level blob, we create a Blob
            blobs.push_back( new surf::GreyLevelBlob() );
            surf::GreyLevelBlob *blob = blobs[blobs.size()-1];

            blob->ssb_parent = ssblob;
            blob->t = (*itGLB)->measurements.t;
            blob->scale = (*itGLB)->GetScale();
            blob->subject = ssblob->subject;
            ssblob->scales.insert( blob->scale );
            // The surf::GreyLevelBlob's nodeslist contains its corresponding nodes indices on the
            //    mesh it was extracted from.
            std::set<SiteType<AimsSurface<3, Void> >::type,
                 ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> > listePoints
                     = (*itGLB)->GetListePoints();
            for ( itPoints = listePoints.begin() ; itPoints != listePoints.end() ; itPoints++ ) 
            {
                blob->nodes.insert( (*itPoints).second );

                if ( subject.coordinates == LATLON_2D ) 
                {
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
//
//void TextureToBlobs::BlobsFromPrimalSketch ( SubjectData & subject,
//                     aims::PrimalSketch<AimsSurface<3, Void>, Texture<float> > &sketch,
//                     std::vector<surf::GreyLevelBlob *> &blobs,
//                     std::vector<surf::ScaleSpaceBlob *> &ssblobs,
//                     bool initNull ) {
//
//    // Initialization of the results vectors "blobs" and "ssblobs"
//    if (initNull){
//        blobs.clear();
//        ssblobs.clear();
//    }
//
//    list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*> listBlobs
//         = sketch.BlobSet();
//
//    std::list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*>::iterator itSSB, itSSBaux, itSSBaux2;
//    std::list<GreyLevelBlob<SiteType<AimsSurface<3, Void> >::type > *>::iterator itGLB;
//
//    ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *ssb;
//    std::set< SiteType<AimsSurface<3, Void> >::type,
//        ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> >::iterator itPoints;
//
//    std::map< ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*, surf::ScaleSpaceBlob * > ssbMap;
//
//    for (itSSB = listBlobs.begin() ; itSSB != listBlobs.end() ; itSSB++){
//
//        // For each scale-space blob, we create a surf::ScaleSpaceBlob in "ssblobs" containing
//        //    various surf::GreyLevelBlob objects (being themselves contained in a general resulting
//        // "blobs" vector).
//        ssb = *itSSB;
//
//        ssblobs.push_back(new surf::ScaleSpaceBlob());
//        surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size() - 1];
//        ssblob->subject = sketch.Subject();
//        ssblob->tmin = 999.0;
//        ssblob->tmax = -999.0;
//
//        // We save a link between the pointer and the corresponding index in ssblobs
//        std::pair< ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type >*, surf::ScaleSpaceBlob * > p;
//        p.first = ssb;
//        p.second = ssblob;
//        ssbMap.insert(p);
//
//        for ( itGLB = ssb->glBlobs.begin() ; itGLB != ssb->glBlobs.end() ; itGLB++ ){
//
//            // For each grey-level blob, we create a Blob
//            blobs.push_back( new surf::GreyLevelBlob() );
//            surf::GreyLevelBlob *blob = blobs[blobs.size()-1];
//
//            blob->ssb_parent = ssblob;
//            blob->t = (*itGLB)->measurements.t;
//            blob->scale = (*itGLB)->GetScale();
//            blob->subject = ssblob->subject;
//
//            // The surf::GreyLevelBlob's nodeslist contains its corresponding nodes indices on the
//            //    mesh it was extracted from.
//            std::set<SiteType<AimsSurface<3, Void> >::type,
//                 ltstr_p3d<SiteType<AimsSurface<3, Void> >::type> > listePoints
//                     = (*itGLB)->GetListePoints();
//            for ( itPoints = listePoints.begin() ; itPoints != listePoints.end() ; itPoints++ ) {
//                (blob->nodes).insert( (*itPoints).second );
//
//                if ( subject.coordinates == LATLON_2D ) {
//                    (blob->coordinates)[(*itPoints).second] = vector<float>(2);
//
//                    (blob->coordinates)[(*itPoints).second][0] = subject.lat->item((*itPoints).second);
//                    (blob->coordinates)[(*itPoints).second][1] = subject.lon->item((*itPoints).second);
//                }
//
//                (blob->raw_coordinates)[(*itPoints).second] = vector<float>(3);
//                (blob->raw_coordinates)[(*itPoints).second][0] = subject.mesh->vertex()[(*itPoints).second][0];
//                (blob->raw_coordinates)[(*itPoints).second][1] = subject.mesh->vertex()[(*itPoints).second][1];
//                (blob->raw_coordinates)[(*itPoints).second][2] = subject.mesh->vertex()[(*itPoints).second][2];
//            }
//
//            ssblob->blobs.insert(blob);
//
//            if ( blob->scale < ssblob->tmin )
//                ssblob->tmin = blob->scale;
//            if ( blob->scale > ssblob->tmax )
//                ssblob->tmax = blob->scale;
//        }
//
//        ssblob->t = ssb->GetMeasurements().t;
//
//    }
//
//
//    // Now that every SSB has been created, we can create links (bifurcations)
//    // between them
//
//    // BIFURCATIONS
//    std::list<Bifurcation<SiteType<AimsSurface<3, Void> >::type> *> bifurcations = sketch.BifurcationList();
//    std::list<Bifurcation<SiteType<AimsSurface<3, Void> >::type> *>::iterator bifurcIt;
//    for ( bifurcIt = bifurcations.begin() ; bifurcIt != bifurcations.end() ; bifurcIt ++ ) {
//        list<ScaleSpaceBlob<SiteType<AimsSurface<3, Void> >::type > *> topBlobs, bottomBlobs, topBottomBlobs, bottomTopBlobs;
//        topBlobs = (*bifurcIt)->TopBlobs();
//        bottomBlobs = (*bifurcIt)->BottomBlobs();
//
//        for (itSSBaux = topBlobs.begin() ; itSSBaux != topBlobs.end() ; itSSBaux++) {
//            for (itSSBaux2 = bottomBlobs.begin() ; itSSBaux2 != bottomBlobs.end() ; itSSBaux2++){
//                ssbMap[*itSSBaux2]->topBlobs.insert( ssbMap[*itSSBaux] );
//                ssbMap[*itSSBaux]->bottomBlobs.insert( ssbMap[*itSSBaux2] );
//            }
//        }
//    }
//
//    std::cout << " blobs.size : " << blobs.size() <<
//            " ssblobs.size : " << ssblobs.size() << std::endl;
//
//}