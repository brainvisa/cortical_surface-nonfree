


#include <aims/graph/graphmanip.h>
#include <cortical_surface/structuralanalysis/representation.h>
#include <cortical_surface/surfacereferential/gyri/mesh_operations.h>
#include <cortical_surface/structuralanalysis/texturetoblobs.h>
#include <cortical_surface/structuralanalysis/meshdistance.h>

using namespace aims;
using namespace std;


std::vector<int> set2vector ( std::set<int> &s ) {
    std::vector<int> v;
    std::set<int>::iterator it;
    for ( it = s.begin() ; it != s.end() ; it++ )
        v.push_back(*it);
    return v;
}
std::vector<float> set2vector ( std::set<float> &s ) {
    std::vector<float> v;
    std::set<float>::iterator it;
    for ( it = s.begin() ; it != s.end() ; it++ )
        v.push_back(*it);
    return v;
}

//##############################################################################

std::set<int> vector2set(std::vector<int> &v){
    std::set<int> s;
    for ( uint i = 0 ; i < v.size() ; i++ )
        s.insert(v[i]);
    return s;
}
std::set<float> vector2set(std::vector<float> &v){
    std::set<float> s;
    for ( uint i = 0 ; i < v.size() ; i++ )
        s.insert(v[i]);
    return s;
}

//##############################################################################


void TextureToBlobs::DestroyBlobs ( std::vector<surf::ScaleSpaceBlob *> &ssblobs ) {
    std::vector<surf::GreyLevelBlob *> blobs;
    std::set<surf::GreyLevelBlob *>::iterator it;

    for ( uint i = 0 ; i < ssblobs.size() ; i++ ) {
        for ( it = ssblobs[i]->blobs.begin() ; it != ssblobs[i]->blobs.end() ; it ++ ) {
            blobs.push_back( *it );
        }
    }
    for ( uint i = 0 ; i < blobs.size() ; i++ )
        delete ( blobs[i] );
    for ( uint i = 0 ; i < ssblobs.size() ; i++ )
        delete ( ssblobs[i] );
}

void TextureToBlobs::DestroyBlobs ( std::vector<surf::Blob *> &blobs ) {
    for ( uint i = 0 ; i < blobs.size() ; i++ )
        delete ( blobs[i] );
}

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
            std::string subject_id;
            float scale, t;
            std::vector<int> nodes_list;
            std::vector<float> latitudes, longitudes;

            blobs.push_back( new surf::GreyLevelBlob() );
            surf::GreyLevelBlob *blob = blobs[blobs.size()-1];

            (*iv)->getProperty( "subject", subject_id );
            (*iv)->getProperty( "scale", scale );
            (*iv)->getProperty( "t", t );
            (*iv)->getProperty( "label", label );
            (*iv)->getProperty( "nodes", nodes_list );
            (*iv)->getProperty( "x", latitudes );
            (*iv)->getProperty( "y", longitudes );

//            assert(nodes_list.size() == latitudes.size());

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

            blob->subject = subject_id;
            blob->nodes = vector2set(nodes_list);
            blob->scale = scale;
            blob->t = t;
            blob->ssb_parent = NULL;
        }
    }
}


void TextureToBlobs::getScaleSpaceBlobsFromIndividualGraph ( Graph *graph,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            std::map<int, std::set<int> > &listGLBindices,
                            bool initNull ){
    if ( initNull )
        ssblobs.clear();

    std::set< Vertex * >::iterator iv;
    Edge *e;
    Vertex::iterator jv;
    Edge::iterator kv;

    listGLBindices = std::map<int, std::set<int> >();
    int iNbLinks = 0;
    int iNbSSB = ssblobs.size();

    for ( iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv ) {
        if ( (*iv)->getSyntax() == "ssb" ) {
            int index;
            float tmax,
                  tmin,
                  t;
            std::vector<float> scales;
            std::string subject_id, label;

            ssblobs.push_back( new surf::ScaleSpaceBlob() );
            surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size()-1];

            (*iv)->getProperty( "tmax", tmax);
            (*iv)->getProperty( "tmin", tmin);
            (*iv)->getProperty( "scales", scales);
            (*iv)->getProperty( "t", t);
            (*iv)->getProperty( "subject", subject_id);
            (*iv)->getProperty( "label", label);

            ssblob->index = iNbSSB++;
            index = ssblob->index;
            (*iv)->setProperty( "index", (int) index );
            (*iv)->setProperty( "sites_index", (int)(ssblobs.size()-1) );

            ssblob->tmax = tmax;
            ssblob->tmin = tmin;
            ssblob->scales = vector2set(scales);
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
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            bool initNull ){

    std::vector<surf::GreyLevelBlob *> blobs;
    std::set<surf::GreyLevelBlob *>::iterator it;
    for ( uint i = 0 ; i < ssblobs.size() ; i++ ) {
        for ( it = ssblobs[i]->blobs.begin() ; it != ssblobs[i]->blobs.end() ; it ++ )
            blobs.push_back( *it );
    }

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
//    iNbGLB = blobs.size() - iNbGLB;

    int iNbSSB = ssblobs.size();
    // We Add A Fictitious Ssb Per Glb
    for ( uint i = 0 ; i < blobs.size() - iNbGLB ; i++ ) {
        ssblobs.push_back( new surf::ScaleSpaceBlob() );
        surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size()-1];

        ssblob->index = iNbSSB + i;
        ssblob->t = blobs[iNbGLB + i]->t;
        ssblob->subject = subject.subject_id;
        ssblob->tmin = blobs[iNbGLB + i]->scale;
        ssblob->tmax = blobs[iNbGLB + i]->scale;
        ssblob->scales.insert(blobs[iNbGLB + i]->scale);
        ssblob->blobs.insert( blobs[iNbGLB + i] );
        ssblob->getNodesFromBlob( * (ssblob->blobs.begin()) );
        assert(ssblob->nodes.size() != 0 );
        blobs[iNbGLB + i]->ssb_parent = ssblob;
    }


    assert( blobs.size() - iNbGLB == ssblobs.size() - iNbSSB );

    cout << blobs.size() - iNbGLB << " blobs added" << endl;
    cout << ssblobs.size() - iNbSSB << " ssblobs added " << endl;

    for (uint i = blobs.size() - iNbGLB ; i < blobs.size() ; i++)
        assert(blobs[i]->ssb_parent != NULL);
}

//##############################################################################

void TextureToBlobs::BlobsFromLabelTexture ( std::vector<surf::Blob *> &blobs,
                             SubjectData &subject ) {

    std::set<short> labels;
    std::set<short>::iterator it;
    for ( uint i = 0 ; i < subject.tex->nItem() ; i++ )
        labels.insert( subject.tex->item(i) );
    std::cout << labels.size() << " labels detected" << std::endl;

    for ( it = labels.begin() ; it != labels.end() ; it ++ ) {
        if ( *it > 0 ) {
            std::cout << *it << std::endl;
            surf::Blob *blob = new surf::Blob();
            blobs.push_back(blob);

            blob->t = *it;
            blob->subject = subject.subject_id;
            blob->raw_coordinates = std::map<int, std::vector<float> > ();
            for ( uint i = 0 ; i < subject.tex->nItem() ; i++ )
                if ( subject.tex->item(i) == *it ) {
                    blob->nodes.insert(i);
                    blob->raw_coordinates[i] = std::vector<float>();
                    for ( uint j = 0 ; j < 3 ; j ++ )
                        blob->raw_coordinates[i].push_back ( subject.mesh->vertex()[i][j] );
                }
        }
    }

    std::cout << blobs.size() << " blobs extracted" << std::endl;


}

//##############################################################################

std::vector<surf::GreyLevelBlob *> recoverGreyLevelBlobs ( const std::vector<surf::ScaleSpaceBlob *> &ssblobs ) {
    std::vector<surf::GreyLevelBlob *> blobs;
    std::set< surf::GreyLevelBlob *>::iterator it;

    for ( uint i = 0 ; i < ssblobs.size() ; i++ )
        for ( it = ssblobs[i]->blobs.begin() ; it != ssblobs[i]->blobs.end() ; it ++ )
            blobs.push_back( *it );

    return blobs;
}

void defineGraphGlobalProperties ( Graph *graph,
                                   SubjectData & subject ) {
    std::vector<float> resolution, bbmin2D, bbmax2D;
    std::vector<int> bbmin, bbmax;
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

}

void addScaleSpaceBlobsToGraph ( Graph *graph,
                                 const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                 std::vector<Vertex *> &listVertSSB,
                                 SubjectData &subject,
                                 bool buildMeshes = true,
                                 bool storeMeshes = true ) {

    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    Vertex *vert;

    if ( buildMeshes ) {
        std::cout << "════ Extracting meshes for the scale-space blobs..." << std::endl;

        std::cout << "      ░░░ mode AimsSphereAtMaxNode ░░░     " << std::endl;
        for ( uint i = 0 ; i < ssblobs.size() ; i++ )
            ssblobs[i]->getAimsSphereAtMaxNode( *(subject.tex), 0.4);
    }

    std::cout << "════ Adding scale-space blobs..." << std::endl;
    listVertSSB = std::vector<Vertex *> ( ssblobs.size() );
    for ( int i = 0 ; i < (int) ssblobs.size() ; i++ ) {

        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index

        cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
        vert = graph->addVertex("ssb");

        vert->setProperty( "label", "0");
        //assert(ssblobs[i]->nodes.size() != 0);
        vert->setProperty( "nodes", set2vector(ssblobs[i]->nodes));
        vert->setProperty( "t", ssblobs[i]->t);
        vert->setProperty( "subject", ssblobs[i]->subject);
        vert->setProperty( "tmin", ssblobs[i]->tmin);
        vert->setProperty( "tmax", ssblobs[i]->tmax);
        vert->setProperty( "scales", set2vector(ssblobs[i]->scales) );
        listVertSSB[  i  ] = vert;

        if ( storeMeshes ) {
            ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
            (*ptr)[0] = ssblobs[i]->mesh;
            manip.storeAims(*graph, vert, "ssb", ptr);
        }
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added... done" << std::endl;

}

void addGreyLevelBlobsToGraph ( Graph *graph,
                                const std::vector<surf::GreyLevelBlob *> &blobs,
                                std::vector<Vertex *> &listVertGLB,
                                SubjectData &subject,
                                bool buildMeshes = true,
                                bool storeMeshes = true ) {

    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    Vertex *vert;

    if ( buildMeshes ) {
        std::cout << "════ Extracting meshes for the grey-level blobs..." << endl;

        std::cout << "      ░░░ mode AimsSphereAtMaxNode ░░░     " << std::endl;
        for ( uint i = 0 ; i < blobs.size() ; i++ )
            // blobs[i]->getAimsEllipsoidAtMaxNode ( *(subject.tex) );
            blobs[i]->getAimsSphereAtMaxNode( *(subject.tex), 0.4);
    }

    std::cout << "════ Adding grey-level blobs..." << endl;
    listVertGLB = std::vector<Vertex *> ( blobs.size() );
    for ( int i = 0 ; i < (int) blobs.size() ; i++ ) {

        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index

        std::cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << std::flush ;
        vert = graph->addVertex("glb");
        vert->setProperty( "index", blobs[i]->index );
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

        if ( storeMeshes ) {
            // We associate the proper mesh patch from "objects" to the vertex
            ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
            (*ptr)[0] = blobs[i]->mesh;
            manip.storeAims(*graph, vert, "glb", ptr);
        }

    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added in total (SSB and GLB)" << std::endl;
}

void addBlobsToGraph ( Graph *graph,
                       const std::vector<surf::Blob *> &blobs,
                       std::vector<Vertex *> &listVertGLB,
                       SubjectData &subject,
                       bool buildMeshes = true,
                       bool storeMeshes = true ) {

    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    Vertex *vert;

    if ( buildMeshes ) {
        std::cout << "════ Extracting meshes for the grey-level blobs..." << endl;

        //std::cout << "      ░░░ mode AimsSphereAtMaxNode ░░░     " << std::endl;
        //for ( uint i = 0 ; i < blobs.size() ; i++ )
        //    blobs[i]->getAimsSphereAtMaxNode( *(subject.tex), 0.4);

        //std::cout << "      ░░░ mode AimsEllipsoidAtMaxNode ░░░     " << std::endl;
        //for ( uint i = 0 ; i < blobs.size() ; i++ )
        //    blobs[i]->getAimsEllipsoidAtMaxNode ( *(subject.tex) );

        std::cout << "      ░░░ mode AimsMesh ░░░     " << std::endl;
        for ( uint i = 0 ; i < blobs.size() ; i++ )
            blobs[i]->getAimsMesh ( *(subject.mesh) );
    }

    std::cout << "════ Adding grey-level blobs..." << endl;
    listVertGLB = std::vector<Vertex *> ( blobs.size() );
    for ( int i = 0 ; i < (int) blobs.size() ; i++ ) {

        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index

        std::cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << std::flush ;
        vert = graph->addVertex("glb");
        vert->setProperty( "subject", blobs[i]->subject.data() );
        vert->setProperty( "t", blobs[i]->t );
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

        if ( storeMeshes ) {
            // We associate the proper mesh patch from "objects" to the vertex
            ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
            (*ptr)[0] = blobs[i]->mesh;
            manip.storeAims(*graph, vert, "glb", ptr);
        }

    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added in total (SSB and GLB)" << std::endl;
}


void addScaleSpaceToGreyLevelBlobsRelations ( Graph *graph,
                                              const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                              const std::vector<Vertex *> &listVertSSB,
                                              const std::vector<Vertex *> &listVertGLB ) {

    std::cout << "════ Adding links between SSB and GLB (GLB forming SSB)..." << std::endl;
    uint iNbLinks = 0;
    Vertex *v1, *v2;
    for (int i = 0 ; i < (int) ssblobs.size() ; i++) {

        std::set<surf::GreyLevelBlob *>::iterator itB1;
        std::set<surf::GreyLevelBlob *> &listGLB = ssblobs[i]->blobs;

        for (itB1 = listGLB.begin(); itB1 != listGLB.end() ; itB1++) {

            v1 = listVertSSB[ i ];
            v2 = listVertGLB[(*itB1)->index];
            graph->addEdge(v1, v2, "s2g");
            iNbLinks++;
        }
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " links added... done" << std::endl;
}

void addBetweenGreyLevelBlobsRelations ( Graph *graph,
                                         const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                         const std::vector<Vertex *> &listVertGLB,
                                         bool buildAndStoreRelationsMeshes ) {

    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    AimsSurfaceTriangle *relations;

    std::cout << "════ Collecting the list of GLB relations..." << std::endl;
    std::vector< std::pair< surf::GreyLevelBlob *, surf::GreyLevelBlob *> > blobsPairs;

    for (int i = 0 ; i < (int) ssblobs.size() ; i++) {
        std::set<surf::GreyLevelBlob *>::iterator itB;
        std::set<surf::GreyLevelBlob *> &unsortedListGLB = ssblobs[i]->blobs;
        std::set<surf::GreyLevelBlob *, ltSurfBlobs> listGLB;
        std::set<surf::GreyLevelBlob *, ltSurfBlobs>::iterator itB1, itB2;
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
    std::cout << "  " << blobsPairs.size() << " blobs pairs" << std::endl;

    if ( buildAndStoreRelationsMeshes ) {
        std::cout << "════ Extracting meshes for the GLB relations..." << std::endl;
        relations = new AimsSurfaceTriangle();
        *relations = getG2GRelationsMeshes( blobsPairs, NODES_BARYCENTERS );
    }

    std::cout << "════ Adding links between GLB from same SSB (with corresponding meshes for visualization)..." << std::endl;

    uint iNbLinks = 0;
    for ( uint i = 0 ; i < blobsPairs.size() ; i++ ) {
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b" << i << std::flush;
        // We associate the proper mesh patch from "objects" to the vertex

        Vertex *v1, *v2;
        surf::GreyLevelBlob *glb1, *glb2;
        glb1 = blobsPairs[i].first;
        glb2 = blobsPairs[i].second;
        v1 = listVertGLB[glb1->index];
        v2 = listVertGLB[glb2->index];
        Edge *e = graph->addEdge(v1, v2, "g2g");

        if ( buildAndStoreRelationsMeshes ) {
            ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
            (*ptr)[0] = (*relations)[i];
            manip.storeAims(*graph, e, "g2g", ptr);
        }

        iNbLinks++;
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " links added... done" << std::endl;
}

void addBifurcationsRelations ( Graph *graph,
                                const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                const std::vector<Vertex *> &listVertSSB,
                                bool buildAndStoreRelationsMeshes ) {

    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    AimsSurfaceTriangle *bifurcations;

    std::cout << "════ Collecting the list of existing bifurcations..." << std::endl;
    std::vector< std::pair< surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> > bifurcPairs;

    for (int i = 0 ; i < (int) ssblobs.size() ; i++) {

        std::set<surf::ScaleSpaceBlob *>::iterator itB1;
        std::set<surf::ScaleSpaceBlob *> &listSSB = ssblobs[i]->topBlobs;
        for (itB1 = listSSB.begin(); itB1 != listSSB.end() ; itB1++) {

            if ( ssblobs[i]->index < (*itB1)->index ) {
                std::pair<surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> p;
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

    if ( buildAndStoreRelationsMeshes ) {
        std::cout << "════ Extracting meshes for the bifurcations relations..." << std::endl;
        bifurcations = new AimsSurfaceTriangle();
        *bifurcations = getBifurcationRelationsMeshes( bifurcPairs, NODES_BARYCENTERS );
        std::cout << "  " << bifurcPairs.size() << " bifurcations pairs" << std::endl;
        std::cout << "  " << bifurcations->size() << " bifurcations meshes" << std::endl;
    }

    std::cout << "════ Adding links between SSB (bifurcations)..." << std::endl;
    uint iNbLinks = 0;

    for ( uint i = 0 ; i < bifurcPairs.size() ; i++ ) {
        Vertex *v1, *v2;
        surf::ScaleSpaceBlob *ssb1, *ssb2;
        ssb1 = bifurcPairs[i].first;
        ssb2 = bifurcPairs[i].second;
        v1 = listVertSSB[ ssb1->index ];
        v2 = listVertSSB[ ssb2->index ];
        Edge *e = graph->addEdge(v1, v2, "bifurcation");
        e->setProperty("type", "test");

        if ( buildAndStoreRelationsMeshes ) {
            ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
            (*ptr)[0] = (*bifurcations)[i];
            manip.storeAims(*graph, e, "bifurcation", ptr);
        }

        iNbLinks++;
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " bifurcations added... done" << std::endl;

}

void getGraphModeOptions ( const int graph_mode,
                           bool &buildScaleSpaceBlobs,
                           bool &buildScaleSpaceBlobsMeshes,
                           bool &storeScaleSpaceBlobsMeshes,
                           bool &buildGreyLevelBlobs,
                           bool &buildGreyLevelBlobsMeshes,
                           bool &storeGreyLevelBlobsMeshes,
                           bool &buildSSBToGLBRelations,
                           bool &buildGLBRelations,
                           bool &buildAndStoreGLBRelationsMeshes,
                           bool &buildBifurcations,
                           bool &buildAndStoreBifurcationsMeshes ) {

    buildScaleSpaceBlobs = false;
    buildScaleSpaceBlobsMeshes = false;
    storeScaleSpaceBlobsMeshes = false;
    buildGreyLevelBlobs = false;
    buildGreyLevelBlobsMeshes = false;
    storeGreyLevelBlobsMeshes = false;
    buildSSBToGLBRelations = false;
    buildGLBRelations = false;
    buildAndStoreGLBRelationsMeshes = false;
    buildBifurcations = false;
    buildAndStoreBifurcationsMeshes = false;

    switch ( graph_mode ) {

    case NO_SCALESPACEBLOBS_MESHES:
        buildScaleSpaceBlobs = true;
        buildGreyLevelBlobs = true;
        buildGreyLevelBlobsMeshes = true;
        storeGreyLevelBlobsMeshes = true;
        buildSSBToGLBRelations = true;
        buildGLBRelations = true;
        buildAndStoreGLBRelationsMeshes = true;
        buildBifurcations = true;
        buildAndStoreGLBRelationsMeshes = true;
    break;
    case NO_SCALESPACEBLOBS_MESHES_AND_NO_RELATIONS_MESHES:
        buildScaleSpaceBlobs = true;
        buildGreyLevelBlobs = true;
        buildGreyLevelBlobsMeshes = true;
        storeGreyLevelBlobsMeshes = true;
        buildSSBToGLBRelations = true;
        //buildGLBRelations = true;
        //buildBifurcations = true;
    break;
    case DEFAULT:
        buildScaleSpaceBlobs = true;
        buildScaleSpaceBlobsMeshes = true;
        storeScaleSpaceBlobsMeshes = true;
        buildGreyLevelBlobs = true;
        buildGreyLevelBlobsMeshes = true;
        storeGreyLevelBlobsMeshes = true;
        buildSSBToGLBRelations = true;
        buildGLBRelations = true;
        buildAndStoreGLBRelationsMeshes = true;
        buildBifurcations = true;
        buildAndStoreGLBRelationsMeshes = true;
    break;
    }

}

// Creates an Aims Graph for only ONE subject with scale-space blobs, grey-level
//  blobs and links between both types

void TextureToBlobs::AimsGraph ( Graph *graph,
                                 SubjectData & subject,
                                 const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                 int graph_mode ) {

     bool buildScaleSpaceBlobs, buildScaleSpaceBlobsMeshes, storeScaleSpaceBlobsMeshes,
          buildGreyLevelBlobs, buildGreyLevelBlobsMeshes, storeGreyLevelBlobsMeshes,
          buildSSBToGLBRelations, buildGLBRelations, buildAndStoreGLBRelationsMeshes,
          buildBifurcations, buildAndStoreBifurcationsMeshes;

     getGraphModeOptions ( graph_mode, buildScaleSpaceBlobs, buildScaleSpaceBlobsMeshes,
             storeScaleSpaceBlobsMeshes, buildGreyLevelBlobs, buildGreyLevelBlobsMeshes,
             storeGreyLevelBlobsMeshes, buildSSBToGLBRelations, buildGLBRelations,
             buildAndStoreGLBRelationsMeshes, buildBifurcations, buildAndStoreBifurcationsMeshes);


    // First We Recover The Grey-Level Blobs From The Vector of Scale-Space Blobs
        std::vector<surf::GreyLevelBlob *> blobs;
        blobs = recoverGreyLevelBlobs ( ssblobs );
        std::cout << blobs.size() << " blobs recovered " << std::endl;

    // We Define The Graph Global Properties
        defineGraphGlobalProperties ( graph, subject );

    // Initializing Blobs Indices
        for (int i = 0 ; i < (int) ssblobs.size() ; i++)
            ssblobs[i]->index = i;
        for (int i = 0 ; i < (int) blobs.size() ; i++)
            blobs[i]->index = i;

    //=========================================
    // Let's Add Or Not The Scale-Space Blobs
        std::vector<Vertex *> listVertSSB, listVertGLB;
        if ( buildScaleSpaceBlobs ) {
            addScaleSpaceBlobsToGraph ( graph, ssblobs, listVertSSB, subject, buildScaleSpaceBlobsMeshes, storeScaleSpaceBlobsMeshes );
        }

    //=========================================
    // Let's Add Or Not The Grey-Level Blobs
        if ( buildGreyLevelBlobs ) {
            addGreyLevelBlobsToGraph ( graph, blobs, listVertGLB, subject, buildGreyLevelBlobsMeshes, storeGreyLevelBlobsMeshes );
        }

    //===========================================================================
    // Let's Add Or Not The Links Between Scale-Space Blobs And Grey-Level Blobs
        if ( buildSSBToGLBRelations && ( buildScaleSpaceBlobs && buildGreyLevelBlobs ) ) {
            addScaleSpaceToGreyLevelBlobsRelations ( graph, ssblobs, listVertSSB, listVertGLB );
        }

    //=======================================================
    // Let's Add Or Not The Links Between Grey-Level Blobs
        if ( buildGLBRelations && ( buildGreyLevelBlobs ) ) {
            addBetweenGreyLevelBlobsRelations ( graph, ssblobs, listVertGLB, buildAndStoreGLBRelationsMeshes );
        }

    //======================================================================
    // Let's Add Or Not The Links Between Scale-Space Blobs (Bifurcations)
        if ( buildBifurcations && ( buildScaleSpaceBlobs ) ) {
            addBifurcationsRelations ( graph, ssblobs, listVertSSB, buildAndStoreGLBRelationsMeshes );
        }
}

void TextureToBlobs::AimsGraph ( Graph *graph,
                                 SubjectData & subject,
                                 const std::vector<surf::Blob *> &blobs ) {

     // We Define The Graph Global Properties
         defineGraphGlobalProperties ( graph, subject );

     // Initializing Blobs Indices
         for (int i = 0 ; i < (int) blobs.size() ; i++)
             blobs[i]->index = i;

     //=========================================
     // Let's Add Or Not The Blobs
         std::vector<Vertex *> listVert;
         addBlobsToGraph ( graph, blobs, listVert, subject, true, true );

}

//##############################################################################

void defineGroupGraphGlobalProperties ( Graph *graph,
                                        std::map<string, SubjectData *> data ) {
    std::vector<float> resolution, bbmin2D, bbmax2D;
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

}

void addInterSubjectRelations ( Graph *graph,
                                std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                std::vector<surf::Clique> &cliques,
                                std::vector<Vertex *> &listVertSSB ) {

    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;

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
        //if ( ( ssblobs[aux1]->subject == "Vil1" && ssblobs[aux2]->subject == "Yel1" ) ||
        //        ( ssblobs[aux2]->subject == "Vil1" && ssblobs[aux1]->subject == "Yel1" ) )
        //    std::cout << ssblobs[aux1]->index << "-" << ssblobs[aux2]->index << std::endl;
        Edge *edge = graph->addEdge( v1, v2, "b2b" );

        edge->setProperty( "similarity", cliques[i].similarity );
        edge->setProperty( "distance", cliques[i].distance );

        ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
        (*ptr)[0]=(*relations)[i];
        manip.storeAims( *graph, edge, "b2b", ptr );
        edge->setProperty( "b2b_label", i );
    }
}


void TextureToBlobs::AimsGroupGraph ( Graph *graph,
                                      std::map<string, SubjectData *> data,
                                      std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                      std::vector<surf::Clique> &cliques ) {

    std::cerr << "Building Group Graph..." << std::endl;
    defineGroupGraphGlobalProperties ( graph, data );

//    Vertex *vert;
//    carto::rc_ptr<AimsSurfaceTriangle> ptr;
//    aims::GraphManip manip;
//    std::vector<Vertex *> listVertSSB( ssblobs.size() ), listVertGLB( ssblobs.size() );
//
//    for ( int i = 0 ; i < (int) ssblobs.size() ; i++ ) {
//
//        // For every scale-space blob, we create a vertex in the Aims graph : we define
//        //   its properties and store a link between the created vertex and the blob index
//
//        cerr << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << flush ;
//        vert = graph->addVertex("ssb");
//        ssblobs[i]->index = i;
//        vert->setProperty( "index", i );
//        vert->setProperty( "label", "0" );
//
//        vert->setProperty( "t", ssblobs[i]->t ); // subj_names[ssblobs[i]->subject]);
//        vert->setProperty( "subject", ssblobs[i]->subject );
//        vert->setProperty( "tmin", ssblobs[i]->tmin );
//        vert->setProperty( "tmax", ssblobs[i]->tmax );
//        vert->setProperty( "scales", set2vector(ssblobs[i]->scales) );
//        vert->setProperty( "nodes", set2vector((*(ssblobs[i]->blobs.begin()))->nodes) );
//
//         // We associate the proper mesh patch from "objects" to the vertex
//         ptr = carto::rc_ptr<AimsSurfaceTriangle>( new AimsSurfaceTriangle );
//         (*ptr)[0] = ssblobs[i]->mesh;
//         manip.storeAims( *graph, vert, "ssb", ptr );
//         vert->setProperty( "ssb_label", i );
//         listVertSSB[ i ] = vert;
//    }

    SubjectData _;
    std::vector<Vertex *> listVertSSB;
    addScaleSpaceBlobsToGraph ( graph, ssblobs, listVertSSB, _, false, true );
    addInterSubjectRelations ( graph, ssblobs, cliques, listVertSSB );

}

void TextureToBlobs::ReadAimsGroupGraph (   Graph &graph,
                                            vector<surf::ScaleSpaceBlob *> &ssblobs,
                                            vector<surf::Clique> &cliques ) {

    std::cout << "Recovering the data..." << std::endl;

    std::set<Vertex *>::iterator iv;
    std::string sujet;
    int index;

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

                    float similarity, distance;
                    assert (!( e->hasProperty("similarity") && e->hasProperty("distance") ));
                    e->getProperty( "similarity", similarity );
                    e->getProperty( "distance", distance );
                    if ( e->hasProperty("similarity") )
                        distance = -1.0;
                    else if ( e->hasProperty("distance") )
                        similarity = -1.0;

                    for ( kv = e->begin() ; kv != e->end() ; kv++ ){

                        if ( (*kv)->getSyntax() == "ssb" ){

                            int indexB2;
                            std::string sujetB2;
                            (*kv)->getProperty( "sites_index", indexB2 );
                            (*kv)->getProperty( "subject", sujetB2 );

                            if ( !(indexB2 == index && sujetB2 == sujet) ){

                                int blobs_index1, blobs_index2;
                                (*iv)->getProperty( "sites_index", blobs_index1 );
                                (*kv)->getProperty( "sites_index", blobs_index2 );

                                assert ( (index==blobs_index1 && indexB2 == blobs_index2) ||
                                  (index == blobs_index2 && indexB2 == blobs_index1));

                                if ( blobs_index1 < blobs_index2 )
                                    cliques.push_back( surf::Clique( ssblobs[ blobs_index1 ] , ssblobs[ blobs_index2 ], distance, similarity ) );
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "ssbcliques.size() :"<< cliques.size() << std::endl << std::endl;
}


std::vector<uint> TextureToBlobs::getClustersListsFromGLB ( std::vector<surf::GreyLevelBlob *> &blobs,
                                                   GroupData &data,
                                                   float clustering_distance_threshold ) {

    std::vector<uint>  clusters ( blobs.size() );
    for ( uint i = 0 ; i < blobs.size() ; i++ )
        clusters[i] = i;

    std::set<uint> dejapris;
    for ( uint i = 0 ; i < blobs.size() ; i++ ) {
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << blobs.size() << std::flush ;

        //assert( blobs[i]->blobs.size() == 1 );

        uint max_node = blobs[i]->getMaximumNode( * (data[blobs[i]->subject]->tex) );

        std::map<uint, float> distanceMap  = LocalMeshDistanceMap(
                    data[blobs[i]->subject]->mesh,
                    data[blobs[i]->subject]->neighbours,
                    max_node,
                    25.0 ); //clustering_distance_threshold + 1.0 );

        std::map<uint, float>::iterator it;

        std::set<uint> voisins_max_node;
        for ( it = distanceMap.begin() ; it != distanceMap.end() ; it ++ )
            if ( it->second < clustering_distance_threshold )
                voisins_max_node.insert( it->first );

        voisins_max_node.insert( max_node );

        std::set<uint> blobs_voisins;
        for ( uint j = 0 ; j < blobs.size() ; j++ ) {
            uint max_node2 = blobs[j]->getMaximumNode( * (data[blobs[j]->subject]->tex) );

            if ( blobs[i]->subject == blobs[j]->subject &&
                        voisins_max_node.find(max_node2) != voisins_max_node.end() ) {
                blobs_voisins.insert( j);
            }
        }

        assert( blobs_voisins.find(i) != blobs_voisins.end() );

        std::set<uint>::iterator ite;
        uint color = clusters[i];
        for ( ite = blobs_voisins.begin() ; ite != blobs_voisins.end() ; ite ++ ) {
            clusters[*ite] = color;
        }

    }
    std::cout << std::endl;
    return clusters;
}


void TextureToBlobs::buildBlobsFromClustersLists ( std::vector< surf::GreyLevelBlob *> &blobs,
                                   GroupData & data,
                                   std::vector<uint> &clusters,
                                   std::vector<surf::ScaleSpaceBlob *> &clusteredSsblobs,
                                   float clustering_distance_threshold,
                                   std::string outputTextFile,
                                   bool uniqueGLB ) {

    FILE *f1;
    std::map<std::string, SubjectData *>::iterator it;

            if ( outputTextFile != "" ) {
                f1 = fopen ( outputTextFile.data(), "a" );
                it = data.begin();
                fprintf(f1, "charac_clusters[\'%s\'][%.3f] = {}\n", it->first.data(), clustering_distance_threshold );
            }

    std::set<uint> colors;
    for ( uint i = 0 ; i <  clusters.size() ; i++ )
        colors.insert ( clusters[i] );
    std::set<uint>::iterator ite;
    std::cout << colors.size() << " clustered ssblobs" << std::endl;

    // Now processing each cluster that has previously been defined
    for ( ite = colors.begin() ; ite != colors.end() ; ite ++ ) {

        uint color = *ite;
        std::vector<uint> cluster_blobs;
        for ( uint i = 0 ; i <  clusters.size() ; i++ )
            if ( clusters[i] == color )
                cluster_blobs.push_back(i);

        // Characterization of each cluster
        float distance_moyenne = 0.0;
        std::map< float, uint > mapCountScales;

        for ( uint i = 0 ; i < cluster_blobs.size() ; i ++ ) {

            uint max_node1 = blobs[cluster_blobs[i]]->getMaximumNode( * (data[blobs[cluster_blobs[i]]->subject]->tex) );
            Point3df p1 = data[blobs[cluster_blobs[i]]->subject]->mesh->vertex()[max_node1];
            if ( mapCountScales.find(blobs[cluster_blobs[i]]->scale) == mapCountScales.end() ) {
                mapCountScales[blobs[cluster_blobs[i]]->scale] = 1;
            }
            else {
                mapCountScales[blobs[cluster_blobs[i]]->scale]++;
            }

            if ( i < cluster_blobs.size() - 1 ) {
                for ( uint j = i + 1 ; j < cluster_blobs.size() ; j ++ ) {
                    uint max_node2 = blobs[cluster_blobs[j]]->getMaximumNode( * (data[blobs[cluster_blobs[j]]->subject]->tex) );
                    Point3df p2 = data[blobs[cluster_blobs[i]]->subject]->mesh->vertex()[max_node2];

                    Point3df d = p1-p2;
                    distance_moyenne += d.norm();
                }
            }

        }

        std::cout << cluster_blobs.size() << std::endl;
        if ( cluster_blobs.size() > 1 )
            distance_moyenne = distance_moyenne / ( cluster_blobs.size() * (cluster_blobs.size()-1) / 2.0 );
        else
            assert(distance_moyenne == 0.0);

            if ( outputTextFile != "" ) {
                fprintf(f1, "charac_clusters[\'%s\'][%.3f][%d] = {\n", it->first.data(), clustering_distance_threshold, color);
                fprintf(f1, "\'dist_moy\' : %lf,\n ", (double)(distance_moyenne) );
                fprintf(f1, "\'nb_total_blob\' : %d,\n ", cluster_blobs.size() );
                fprintf(f1, "\'map_count_scales\' : {");
                std::map< float, uint >::iterator ite3;
                for ( ite3 = mapCountScales.begin() ; ite3 != mapCountScales.end() ; ite3++ )
                    fprintf(f1, " %lf : %d, ", (double)(ite3->first), ite3->second );
                fprintf(f1, "}\n");
                fprintf(f1, "}\n\n");
            }


        // Creation of Blobs
        clusteredSsblobs.push_back ( new surf::ScaleSpaceBlob() );

        surf::ScaleSpaceBlob *ssb = clusteredSsblobs[ clusteredSsblobs.size() -1 ];
        ssb->subject = blobs[ cluster_blobs[0]]->subject;

        std::set<float> scales;
        surf::GreyLevelBlob nodes_dummyglb;
        //bool uniqueGLB = false;
        if ( uniqueGLB ) {
            // Should there be one single GLB created by merging all existing GLB ?

            ssb->blobs.insert( new surf::GreyLevelBlob() );
            surf::GreyLevelBlob *glb = *(ssb->blobs.begin());
            glb->t = blobs[cluster_blobs[0] ]->t;
            glb->scale = blobs[cluster_blobs[0] ]->scale;

            for ( uint i = 0 ; i < cluster_blobs.size() ; i ++ ) {
                std::set<int>::iterator ite2;
                for ( ite2 = blobs[ cluster_blobs[i] ]->nodes.begin() ;
                    ite2 != blobs[ cluster_blobs[i] ]->nodes.end() ;
                    ite2++ ) {
                        glb->nodes.insert( *ite2 );
                        glb->raw_coordinates[ *ite2 ] = std::vector<float>(3);
                        for ( uint k = 0 ; k < 3 ; k++ )
                            glb->raw_coordinates[ *ite2 ][k] = blobs[ cluster_blobs[i]]->raw_coordinates[*ite2][k];
                        if ( blobs[cluster_blobs[i]]->coordinates.size() != 0 ) {
                            glb->coordinates[ *ite2 ] = std::vector<float>(2);
                            for ( uint k = 0 ; k < 2 ; k++ ) {
                                glb->coordinates[ *ite2 ][k] = blobs[ cluster_blobs[i]]->coordinates[*ite2][k];
                            }
                        }
                }
                scales.insert( blobs[cluster_blobs[i] ]->scale );
            }
            glb->ssb_parent = ssb;
            glb->getAimsSphereAtMaxNode ( * (data[ssb->subject]->tex), 0.4 );
            ssb->getNodesFromBlob(glb);
            ssb->getAimsSphereAtMaxNode ( * (data[ssb->subject]->tex), 0.4 );
        }
        else {
            //   Or should we keep track of original GLB ?

            for ( uint i = 0 ; i < cluster_blobs.size() ; i ++ ) {

                surf::GreyLevelBlob *glb = new surf::GreyLevelBlob();
                ssb->blobs.insert( glb );
                glb->t = blobs[cluster_blobs[i] ]->t;
                glb->scale = blobs[cluster_blobs[i] ]->scale;
                glb->ssb_parent = ssb;

                std::set<int>::iterator ite2;
                for ( ite2 = blobs[ cluster_blobs[i] ]->nodes.begin() ;
                    ite2 != blobs[ cluster_blobs[i] ]->nodes.end() ;
                    ite2++ ) {
                        glb->nodes.insert( *ite2 );
                        nodes_dummyglb.nodes.insert(*ite2);
                        glb->raw_coordinates[ *ite2 ] = std::vector<float>(3);
                        nodes_dummyglb.raw_coordinates[ *ite2 ] = std::vector<float>(3);
                        for ( uint k = 0 ; k < 3 ; k++ ) {
                            glb->raw_coordinates[ *ite2 ][k] = blobs[ cluster_blobs[i]]->raw_coordinates[*ite2][k];
                            nodes_dummyglb.raw_coordinates[ *ite2 ][k] = blobs[ cluster_blobs[i]]->raw_coordinates[*ite2][k];
                        }
                        if ( blobs[cluster_blobs[i]]->coordinates.size() != 0 ) {
                            glb->coordinates[ *ite2 ] = std::vector<float>(2);
                            nodes_dummyglb.coordinates[ *ite2 ] = std::vector<float>(2);
                            for ( uint k = 0 ; k < 2 ; k++ ) {
                                glb->coordinates[ *ite2 ][k] = blobs[ cluster_blobs[i]]->coordinates[*ite2][k];
                                nodes_dummyglb.coordinates[ *ite2 ][k] = blobs[ cluster_blobs[i]]->coordinates[*ite2][k];
                            }
                        }
                }
                glb->getAimsSphereAtMaxNode ( * (data[ssb->subject]->tex), 0.4 );


                scales.insert( blobs[cluster_blobs[i] ]->scale );
            }

            ssb->getNodesFromBlob(&nodes_dummyglb);

        }

        ssb->tmax = *(scales.rbegin());
        ssb->tmin = *(scales.begin());
        ssb->scales = scales;
        std::cout << ssb->tmax << "-" << ssb->tmin << " " << std::flush;
        std::cout << "("<<ssb->blobs.size()<<")" << std::endl;
        ssb->t = data[ssb->subject]->tex->item(ssb->getMaximumNode(*(data[ssb->subject]->tex) ) );
        ssb->label = 0;

        //else if ( type_distance == DISTANCE_LATITUDES ) {
        //    int maxim_node = ssb->getMaximumNode(* (data[ssb->subject]->tex));
        //    AimsSurfaceTriangle *sph;
        //    Point3df p( (*(ssb->blobs.begin()))->coordinates[maxim_node][0],
        //            (*(ssb->blobs.begin()))->coordinates[maxim_node][1], 0.0);
        //    if (p[1] > 300.0)
        //        p[1] -= 360.0;
        //    sph = SurfaceGenerator::sphere( p , 0.3, 10);
        //    ssb->mesh = (*sph)[0];
        //}
    }

    if ( outputTextFile != "" )
        fclose(f1);
}

double TextureToBlobs::getOverlapMeasure( Point2df bbmin1,
                          Point2df bbmax1,
                          Point2df bbmin2,
                          Point2df bbmax2,
                          uint *no_overlap ){

    float overlap_x, overlap_y, aux;
    double rec = 0.0;

    if ( sqrt(pow( bbmin1[0] - bbmax1[0], 2) ) < 0.0001 )  bbmax1[0] += 0.5;
    if ( sqrt(pow( bbmin1[1] - bbmax1[1], 2) ) < 0.0001 )  bbmax1[1] += 0.5;
    if ( sqrt(pow( bbmin2[0] - bbmax2[0], 2) ) < 0.0001 )  bbmax2[0] += 0.5;
    if ( sqrt(pow( bbmin2[1] - bbmax2[1], 2) ) < 0.0001 )  bbmax2[1] += 0.5;

    if (sqrt(pow(bbmin1[1] - bbmax1[1], 2)) > 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) < 300) {

        if ( 360 - bbmax2[1] < bbmin2[1] ) {
            aux = bbmax1[1];
            bbmax1[1] = bbmin1[1] + 360.0;
            bbmin1[1] = aux;
        }
        else {
            aux = bbmin1[1];
            bbmin1[1] = bbmax1[1] - 360.0;
            bbmax1[1] = aux;
        }
    }

    else if (sqrt(pow( bbmin1[1] - bbmax1[1], 2)) < 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) > 300) {

        if ( 360 - bbmax1[1] < bbmin1[1] ) {
            aux = bbmax2[1];
            bbmax2[1] = bbmin2[1] + 360.0;
            bbmin2[1] = aux;
        }
        else {
            aux = bbmin2[1];
            bbmin2[1] = bbmax2[1] - 360.0;
            bbmax2[1] = aux;
        }
    }
    else if (sqrt(pow( bbmin1[1] - bbmax1[1], 2)) > 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) > 300) {

        aux = bbmin1[1];
        bbmin1[1] = bbmax1[1] - 360.0;
        bbmax1[1] = aux;
        aux = bbmin2[1];
        bbmin2[1] = bbmax2[1] - 360.0;
        bbmax2[1] = aux;
    }

    // ON S'OCCUPE DE LA LATITUDE
    if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) > 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) < 150) {

        if ( 180 - bbmax2[0] < bbmin2[0] ) {
            aux = bbmax1[0];
            bbmax1[0] = bbmin1[0] + 180.0;
            bbmin1[0] = aux;
        }
        else {
            aux = bbmin1[0];
            bbmin1[0] = bbmax1[0] - 180.0;
            bbmax1[0] = aux;
        }
    }
    else if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) < 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) > 150){

        if ( 180 - bbmax1[0] < bbmin1[0] ) {
            aux = bbmax2[0];
            bbmax2[0] = bbmin2[0] + 180.0;
            bbmin2[0] = aux;
        }
        else {
            aux = bbmin2[0];
            bbmin2[0] = bbmax2[0] - 180.0;
            bbmax2[0] = aux;
        }

    }

    else if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) > 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) > 150){

        aux = bbmin1[0];
        bbmin1[0] = bbmax1[0] - 360.0;
        bbmax1[0] = aux;
        aux = bbmin2[0];
        bbmin2[0] = bbmax2[0] - 360.0;
        bbmax2[0] = aux;
    }

    // PRÉTRAITEMENTS EFFECTUÉS ON CALCULE LE RECOUVREMENT

    *no_overlap = 0;
    if ( bbmin1[0] <= bbmin2[0] )
        if ( bbmax1[0] < bbmin2[0] )
            *no_overlap = 1;
        else
            overlap_x = ( bbmax2[0] < bbmax1[0] ? bbmax2[0] : bbmax1[0] ) - bbmin2[0];
    else
        if ( bbmax2[0] < bbmin1[0] )
            *no_overlap = 1;
        else
            overlap_x = ( bbmax1[0] < bbmax2[0] ? bbmax1[0] : bbmax2[0] ) - bbmin1[0];

    if ( *no_overlap == 0 ) {

        if ( bbmin1[1] <= bbmin2[1] )
            if ( bbmax1[1] < bbmin2[1] )
                *no_overlap = 1;
            else
                overlap_y = ( bbmax2[1] < bbmax1[1] ? bbmax2[1] : bbmax1[1] ) - bbmin2[1];
        else
            if ( bbmax2[1] < bbmin1[1] )
                *no_overlap = 1;
            else
                overlap_y = ( bbmax1[1] < bbmax2[1] ? bbmax1[1] : bbmax2[1] ) - bbmin1[1];
        if ( *no_overlap == 0 ) {
            rec = overlap_x * overlap_y;
            double div=  ( bbmax1[0] - bbmin1[0] ) * ( bbmax1[1] - bbmin1[1] )
                + ( bbmax2[0] - bbmin2[0] ) * ( bbmax2[1] - bbmin2[1] ) ;

            rec = 2 * rec / div;
        }
    }

    return rec;
}





//##############################################################################



bool TextureToBlobs::isInside2DBox( Point2df p1, Point2df bbmin, Point2df bbmax) {
    uint no_overlap = 2;

    Point2df bbmin1 (p1[0] - 0.0001, p1[1] - 0.0001),
        bbmax1 (p1[0] + 0.0001, p1[1] + 0.0001),
        bbmin2 (bbmin[0], bbmin[1]),
        bbmax2 (bbmax[0], bbmax[1]);

    getOverlapMeasure( bbmin1, bbmax1, bbmin2, bbmax2, &no_overlap );
    if (no_overlap == 0)
        return true;
    else if (no_overlap == 1)
        return false;
    assert(false);
    return false;
}


void TextureToBlobs::filteringBlobs (  std::vector<surf::ScaleSpaceBlob *> & ssblobs,
                       std::vector<surf::GreyLevelBlob *> &filteredBlobs,
                       std::vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
                       std::set<int> &nodes ){
    // Not much used : allows to filter in blobs provided their support intersects
    //    with nodes belonging to a given set

    for (uint i = 0 ; i < ssblobs.size() ; i++ )
        ssblobs [i] -> index = i;

    std::set< uint > filteredIndices;

    // Filtering according to positions
    for (uint i = 0 ; i < ssblobs.size() ; i++ ){
        std::string subject;
        subject = ssblobs[i]->subject;
        bool firstGLB = false;

        surf::ScaleSpaceBlob *ssb;
        ssb = new surf::ScaleSpaceBlob( ssblobs[i] );
        ssb->blobs.clear();
        ssb->topBlobs.clear();
        ssb->bottomBlobs.clear();

        std::set<surf::GreyLevelBlob *>::iterator itB1;
        for (itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() && !firstGLB ; itB1++){
            uint no_overlap = 2;

            std::set<int>::iterator it;
            std::set<int> intersection;
            for ( it = (*itB1)->nodes.begin() ; it != (*itB1)->nodes.end() ; it++ )
                if ( nodes.find(*it) != nodes.end() )
                    intersection.insert(*it);

            if ( intersection.size() != 0 )
                no_overlap = 0;
            else
                no_overlap = 1;
//             }
            if (no_overlap == 0)
                firstGLB = true;
        }

        if (firstGLB) {
            for ( itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() ; itB1++ ) {

                surf::GreyLevelBlob *glb;
                glb = new surf::GreyLevelBlob( *itB1 );
                ASSERT(glb->nodes.size() == glb->raw_coordinates.size());
                glb->ssb_parent = ssb;
                ssb->blobs.insert(glb);
                filteredBlobs.push_back(glb);
            }
            filteredSsblobs.push_back(ssb);
            filteredIndices.insert ( ssb->index );
        }
        else
            delete(ssb);
    }
    cerr << filteredBlobs.size() << " filtered blobs - " << filteredSsblobs.size() << " filtered ssblobs" << endl;


    // Now that the blobs are filtered, we add the correct bifurcations

    for ( uint i = 0 ; i < filteredSsblobs.size() ; i ++ ) {

        set<surf::ScaleSpaceBlob *> auxTop = ssblobs[filteredSsblobs[i]->index]->topBlobs;
        set<surf::ScaleSpaceBlob *>::iterator it;
        for ( it = auxTop.begin() ; it != auxTop.end() ; it ++ ) {
            if (filteredIndices.find((*it)->index) != filteredIndices.end()) {
                uint i1 = 0;
                for ( ; i1 < filteredSsblobs.size() && filteredSsblobs[i1]->index != (*it)->index ; i1 ++ ) {}
                filteredSsblobs[i]->topBlobs.insert( filteredSsblobs[i1] );
            }
        }

        set<surf::ScaleSpaceBlob *> auxBot = ssblobs[filteredSsblobs[i]->index]->bottomBlobs;
        for ( it = auxBot.begin() ; it != auxBot.end() ; it ++ ) {
            if (filteredIndices.find((*it)->index) != filteredIndices.end()){
                uint i1 = 0;
                for ( ; i1 < filteredSsblobs.size() && filteredSsblobs[i1]->index != (*it)->index ; i1 ++ ) {}
                filteredSsblobs[i]->bottomBlobs.insert(filteredSsblobs[i1]);
            }
        }

    }

}
