
#include <aims/graph/graphmanip.h>
#include <cortical_surface/structuralanalysis/representation.h>
#include <cortical_surface/surfacereferential/gyri/mesh_operations.h>
#include <cortical_surface/structuralanalysis/texturetoblobs.h>
#include <cortical_surface/structuralanalysis/meshdistance.h>

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
    ssblobs.clear();
}

void TextureToBlobs::DestroyBlobs ( std::vector<surf::Blob *> &blobs ) {
    for ( uint i = 0 ; i < blobs.size() ; i++ )
        delete ( blobs[i] );
    blobs.clear();
}

void TextureToBlobs::getGreyLevelBlobsFromGraph ( Graph *graph,
                            SubjectData &subject,
                            std::vector <surf::GreyLevelBlob *> &blobs,
                            bool initNull ){
    if ( initNull )
        blobs.clear();
    std::set <Vertex *>::iterator iv;
    int iNbGLB = 0;
    int filtered_out;
    
    for ( iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv ) {
        if ( (*iv)->getSyntax() == "glb" ) {    
            filtered_out = 0;
            if ( (*iv)->hasProperty("filtered_out") )
                (*iv)->getProperty("filtered_out", filtered_out);
            if ( filtered_out == 0 ) {
                int index, label;
                std::string subject_id;
                float scale, t;
                std::vector<int> nodes_list;
                std::vector<float> latitudes, longitudes;
    
                blobs.push_back( new surf::GreyLevelBlob() );
                surf::GreyLevelBlob *blob = blobs[blobs.size()-1];
                (*iv)->getProperty( "index", index );
                (*iv)->getProperty( "subject", subject_id );
                (*iv)->getProperty( "scale", scale );
                (*iv)->getProperty( "t", t );
                (*iv)->getProperty( "nodes", nodes_list );
                (*iv)->getProperty( "x", latitudes );
                (*iv)->getProperty( "y", longitudes );
                blob->index = index;
                blob->subject = subject_id;
                blob->nodes = vector2set ( nodes_list );
                assert(blob->nodes.size() != 0);
                blob->scale = scale;
                blob->t = t;
                blob->ssb_parent = NULL;
                
                for ( uint i = 0 ; i < nodes_list.size() ; i++ ) {
                    (blob->raw_coordinates)[nodes_list[i]] = std::vector<float>(3);
                    (blob->raw_coordinates)[nodes_list[i]][0] = subject.mesh->vertex()[nodes_list[i]][0];
                    (blob->raw_coordinates)[nodes_list[i]][1] = subject.mesh->vertex()[nodes_list[i]][1];
                    (blob->raw_coordinates)[nodes_list[i]][2] = subject.mesh->vertex()[nodes_list[i]][2];
    
                    if ( subject.coordinates == LATLON_2D ) {
                        ( blob->coordinates)[nodes_list[i]] = std::vector<float>(2);
                        ( blob->coordinates)[nodes_list[i]][0] = latitudes[i];
                        ( blob->coordinates)[nodes_list[i]][1] = longitudes[i];
                    }
                }
            }
        }
    }
}


void TextureToBlobs::getScaleSpaceBlobsFromGraph ( Graph *graph,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            std::map< std::string, std::map<int, std::set<int> > > &listGLBindices,
                            bool initNull ){
    if ( initNull )
        ssblobs.clear();
    std::set< Vertex * >::iterator iv;
    Edge *e;
    Vertex::iterator jv;
    Edge::iterator kv;
    int iNbSSB = ssblobs.size();
    int index, filtered_out;
    float tmax, tmin, t;
    std::vector<float> scales, bc;
    std::vector<int> nodes;
    std::string subject_id, label;

    for ( iv = graph->vertices().begin() ; iv != graph->vertices().end() ; ++iv ) {
        if ( (*iv)->getSyntax() == "ssb" ) {            
            filtered_out = 0;
            if ( (*iv)->hasProperty ( "filtered_out" ) ) 
                (*iv)->getProperty ( "filtered_out", filtered_out );
            if ( filtered_out == 0 ) {
                ssblobs.push_back( new surf::ScaleSpaceBlob() );
                surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size()-1];
                
                (*iv)->getProperty( "index", index );
                (*iv)->getProperty( "tmax", tmax );
                (*iv)->getProperty( "tmin", tmin );
                (*iv)->getProperty( "scales", scales );
                (*iv)->getProperty( "t", t );
                (*iv)->getProperty( "subject", subject_id );
                (*iv)->getProperty( "label", label );
                (*iv)->getProperty( "nodes", nodes );
                (*iv)->getProperty( "gravity_center", bc );
                (*iv)->setProperty( "sites_index", (int) ( ssblobs.size() - 1 ) );
                ssblob->index = index;
                ssblob->tmax = tmax;
                ssblob->tmin = tmin;
                ssblob->scales = vector2set(scales);
                ssblob->t = t;
                ssblob->subject = subject_id;
                ssblob->label = atoi(label.data());
                ssblob->nodes = vector2set ( nodes );
                //ssblob->nodes.insert(node);
                ssblob->raw_coordinates = std::map<int, std::vector<float> >();
                //ssblob->raw_coordinates[node] = std::vector<float> ( bc.size() );
                //for ( uint i = 0 ; i < bc.size() ; i++ )
                //    ssblob->raw_coordinates[node][i] = bc[i];
    
                if ( listGLBindices.find(subject_id) == listGLBindices.end() )
                    listGLBindices[subject_id] = std::map<int, std::set<int> >();
    
                if ( listGLBindices[subject_id].find( index ) == listGLBindices[subject_id].end() )
                    listGLBindices[subject_id][index] = std::set<int>();
    
                for ( jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++ ) {
                    e = *jv;
    
                    if ( e->getSyntax() == "s2g" ) {
                        for ( kv = e->begin() ; kv != e->end() ; kv++ ) {
                            if ( (*kv)->getSyntax() == "ssb" ) {
                            }
                            else if ( (*kv)->getSyntax() == "glb" ) {
                                int iGLBindex;
                                (*kv)->getProperty( "index", iGLBindex );
                                listGLBindices[subject_id][index].insert( iGLBindex );
                            }
                        }
                    }
                }
            }
        }
    }
}

surf::GreyLevelBlob *TextureToBlobs::findBlob ( const std::vector<surf::GreyLevelBlob *> &blobs,
                                std::string subject_id,
                                int index ) {
    for ( uint i = 0 ; i < blobs.size() ; i ++ )
        if ( blobs[i]->subject == subject_id && blobs[i]->index == index )
            return blobs[i];
    return NULL;
}

surf::ScaleSpaceBlob *TextureToBlobs::findBlob ( const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                std::string subject_id,
                                int index ) {
    for ( uint i = 0 ; i < ssblobs.size() ; i ++ )
        if ( ssblobs[i]->subject == subject_id && ssblobs[i]->index == index )
            return ssblobs[i];
    return NULL;
}

int TextureToBlobs::findBlobIndex ( const std::vector<surf::GreyLevelBlob *> &blobs,
                                std::string subject_id,
                                int index ) 
{
    for ( uint i = 0 ; i < blobs.size() ; i ++ )
        if ( blobs[i]->subject == subject_id && blobs[i]->index == index )
            return i;
    return -1;
}

int TextureToBlobs::findBlobIndex ( const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                std::string subject_id,
                                int index ) {
    for ( uint i = 0 ; i < ssblobs.size() ; i ++ )
        if ( ssblobs[i]->subject == subject_id && ssblobs[i]->index == index )
            return i;
    return -1;
}

void TextureToBlobs::RecoverBlobsFromGraph( Graph *graph,
                            SubjectData &subject,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            bool initNull )
{
    std::set<int>::iterator it;
    std::vector<surf::GreyLevelBlob *> blobs = recoverGreyLevelBlobs(ssblobs);

    std::map < std::string, std::map< int, std::set<int> > > listGLBindices;
    int iNbGLB = blobs.size();
    int iNbSSB = ssblobs.size();

    getGreyLevelBlobsFromGraph ( graph, subject, blobs, initNull ) ;
    iNbGLB = blobs.size() - iNbGLB;
    std::cout << iNbGLB << " blobs added" << std::endl;
    
    getScaleSpaceBlobsFromGraph ( graph, ssblobs, listGLBindices, initNull ) ;
    iNbSSB = ssblobs.size() - iNbSSB;

    for ( uint i = ssblobs.size() - iNbSSB ; i < ssblobs.size() ; i++ ) 
    {
        int index = ssblobs[i]->index;
        std::string subject_id = ssblobs[i]->subject;
        for ( it = listGLBindices[subject_id][index].begin() ; it != listGLBindices[subject_id][index].end() ; it++ ) 
        {
            surf::GreyLevelBlob *glb = findBlob ( blobs, subject_id, *it);
            assert ( glb != NULL);
            ssblobs[i]->blobs.insert ( glb );
            glb->ssb_parent = ssblobs[i];
        }
    }
    std::cout << iNbSSB << " ssblobs added " << std::endl;
    for (uint i = blobs.size() - iNbGLB ; i < blobs.size() ; i++)
        assert(blobs[i]->ssb_parent != NULL);
}

void TextureToBlobs::RecoverBlobsFromGLBOnly( Graph *graph,
                            SubjectData &subject,
                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                            bool initNull,
                            float thresholdOnT )
{
    std::vector<surf::GreyLevelBlob *> blobs;    
    if ( thresholdOnT < -99.0 ) 
        getGreyLevelBlobsFromGraph ( graph, subject, blobs, initNull );
    else 
    {
        std::vector<surf::GreyLevelBlob *> proto_blobs;
        getGreyLevelBlobsFromGraph ( graph, subject, proto_blobs, initNull );
        for ( uint i = 0 ; i < proto_blobs.size() ; i++ ) 
        {
            assert( proto_blobs[i]->nodes.size() != 0 );
            if ( proto_blobs[i]->t > thresholdOnT )
                blobs.push_back ( new surf::GreyLevelBlob(*(proto_blobs[i])) );
        }
    }

    int iNbSSB = ssblobs.size();
    // We Add A Fictitious Ssb Per Glb
    for ( uint i = 0 ; i < blobs.size() ; i++ ) 
    {
        ssblobs.push_back( new surf::ScaleSpaceBlob() );
        surf::ScaleSpaceBlob *ssblob = ssblobs[ssblobs.size()-1];

        ssblob->index = blobs[i]->index;
        ssblob->t = blobs[i]->t;
        ssblob->subject = subject.subject_id;
        ssblob->tmin = blobs[i]->scale;
        ssblob->tmax = blobs[i]->scale;
        ssblob->scales.insert(blobs[i]->scale);
        ssblob->blobs.insert( blobs[i] );
        ssblob->getNodesFromBlob( * (ssblob->blobs.begin()) );
        assert( blobs[i]->nodes.size() != 0 );
        assert(ssblob->nodes.size() != 0 );
        blobs[i]->ssb_parent = ssblob;
    }

    std::cout << blobs.size() << " blobs added" << std::endl;
    std::cout << ssblobs.size() - iNbSSB << " ssblobs added " << std::endl;

    for (uint i = blobs.size() ; i < blobs.size() ; i++)
        assert(blobs[i]->ssb_parent != NULL);
}

//##############################################################################

void TextureToBlobs::BlobsFromLabelTexture ( std::vector<surf::Blob *> &blobs,
                             SubjectData &subject ) 
{
    std::set<short> labels;
    std::set<short>::iterator it;
    for ( uint i = 0 ; i < subject.tex->nItem() ; i++ )
        labels.insert( subject.tex->item(i) );
    std::cout << labels.size() << " labels detected" << std::endl;

    for ( it = labels.begin() ; it != labels.end() ; it ++ ) 
    {
        if ( *it > 0 ) 
        {
            std::cout << *it << std::endl;
            surf::Blob *blob = new surf::Blob();
            blobs.push_back(blob);

            blob->t = *it;
            blob->subject = subject.subject_id;
            blob->raw_coordinates = std::map<int, std::vector<float> > ();
            for ( uint i = 0 ; i < subject.tex->nItem() ; i++ )
                if ( subject.tex->item(i) == *it ) 
                {
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

std::vector<surf::GreyLevelBlob *> TextureToBlobs::recoverGreyLevelBlobs ( const std::vector<surf::ScaleSpaceBlob *> &ssblobs ) {
    std::vector<surf::GreyLevelBlob *> blobs;
    std::set< surf::GreyLevelBlob *>::iterator it;

    for ( uint i = 0 ; i < ssblobs.size() ; i++ )
        for ( it = ssblobs[i]->blobs.begin() ; it != ssblobs[i]->blobs.end() ; it ++ )
            blobs.push_back( *it );

    return blobs;
}

void defineGraphGlobalProperties ( Graph *graph,
                                   SubjectData & subject ) 
{
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
                                 std::map<std::string, std::map<int,Vertex *> > &listVertSSB,
                                 std::map<std::string, SubjectData *> &data,
                                 bool storeMeshes = true,
                                 int representationMode = SPHERES ) {

    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    Vertex *vert;
    uint node;
    Point3df maxnode;

    if ( representationMode == SPHERES ) {
        std::cout << "════ Extracting meshes for the scale-space blobs..." << std::endl;
        std::cout << "      ░░░ mode AimsSphereAtMaxNode ░░░     " << std::endl;
        for ( uint i = 0 ; i < ssblobs.size() ; i++ )
            ssblobs[i]->getAimsSphereAtMaxNode( *(data[ssblobs[i]->subject]->tex), 0.4);
    }
    else if ( representationMode == CORTICAL_PATCHES ) {
        std::cout << "════ Extracting meshes for the scale-space blobs..." << std::endl;
        std::cout << "      ░░░ mode AimsMesh ░░░     " << std::endl;
        for ( uint i = 0 ; i < ssblobs.size() ; i++ )
            ssblobs[i]->getAimsMesh( *(data[ssblobs[i]->subject]->mesh) );
    }

    std::cout << "════ Adding scale-space blobs..." << std::endl;
    listVertSSB = std::map<std::string, std::map<int,Vertex *> >();
    for ( int i = 0 ; i < (int) ssblobs.size() ; i++ ) {

        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index

        std::cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << std::flush ;
        vert = graph->addVertex("ssb");

        vert->setProperty( "label", "0");
        vert->setProperty( "index", ssblobs[i]->index);
        vert->setProperty( "t", ssblobs[i]->t);
        vert->setProperty( "subject", ssblobs[i]->subject);
        if ( listVertSSB.find(ssblobs[i]->subject) == listVertSSB.end() )
            listVertSSB[ssblobs[i]->subject] = std::map< int, Vertex * >();
        vert->setProperty( "tmin", ssblobs[i]->tmin);
        vert->setProperty( "tmax", ssblobs[i]->tmax);
        vert->setProperty( "scales", set2vector(ssblobs[i]->scales) );

        if ( ssblobs[i]->nodes.size() > 0 ) {
            vert->setProperty( "nodes", set2vector(ssblobs[i]->nodes));
            node = ssblobs[i]->getMaximumNode(*(data[ssblobs[i]->subject]->tex));
            maxnode = data[ssblobs[i]->subject]->mesh->vertex()[ node ];
            std::vector<float> bc(3);
            for ( uint j = 0 ; j < bc.size() ; j++ )
                bc[j] = maxnode[j];
            vert->setProperty ( "gravity_center", bc );
        }

        listVertSSB[ssblobs[i]->subject][ssblobs[i]->index] = vert;
        if ( storeMeshes ) {
            ptr = carto::rc_ptr<AimsSurfaceTriangle> ( new AimsSurfaceTriangle );
            (*ptr)[0] = ssblobs[i]->mesh;
            manip.storeAims(*graph, vert, "ssb", ptr);
        }
    }

    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added... done" << std::endl;

}

void addGreyLevelBlobsToGraph ( Graph *graph,
                                const std::vector<surf::GreyLevelBlob *> &blobs,
                                std::map<std::string, std::map<int, Vertex *> > &listVertGLB,
                                SubjectData &subject,
                                bool storeMeshes = true,
                                int representationMode = SPHERES ) 
{
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    Vertex *vert;
    if ( representationMode == SPHERES ) 
    {
        std::cout << "════ Extracting meshes for the grey-level blobs..." << std::endl;
        std::cout << "      ░░░ mode AimsSphereAtMaxNode ░░░     " << std::endl;
        for ( uint i = 0 ; i < blobs.size() ; i++ )
            blobs[i]->getAimsSphereAtMaxNode( *(subject.tex), 0.4);
    }
    else if ( representationMode == CORTICAL_PATCHES ) 
    {
        std::cout << "════ Extracting meshes for the grey-level blobs..." << std::endl;
        std::cout << "      ░░░ mode AimsMesh ░░░     " << std::endl;
        for ( uint i = 0 ; i < blobs.size() ; i++ )
            blobs[i]->getAimsMesh( *(subject.mesh) );
    }
    std::cout << "════ Adding grey-level blobs..." << std::endl;
    listVertGLB = std::map<std::string, std::map<int, Vertex *> > ();
    for ( int i = 0 ; i < (int) blobs.size() ; i++ ) 
    {
        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << std::flush ;
        vert = graph->addVertex("glb");
        vert->setProperty( "index", blobs[i]->index );
        vert->setProperty( "t", blobs[i]->t);
        vert->setProperty( "scale", blobs[i]->scale );
        vert->setProperty( "subject", blobs[i]->ssb_parent->subject );
        if ( listVertGLB.find( blobs[i]->ssb_parent->subject ) == listVertGLB.end() )
            listVertGLB[blobs[i]->ssb_parent->subject] = std::map<int, Vertex *>();
        vert->setProperty( "nodes", set2vector(blobs[i]->nodes) );
        if ( subject.coordinates == LATLON_2D ) 
        {
            std::vector<float> latitudes, longitudes;
            std::set<int>::iterator it;
            for ( it = blobs[i]->nodes.begin() ; it != blobs[i]->nodes.end() ; it++ ) 
            {
                latitudes.push_back( (float) blobs[i]->coordinates[*it][0] );
                if ( blobs[i]->coordinates[*it].size() >= 2 )
                    longitudes.push_back( (float) blobs[i]->coordinates[*it][1] );
            }
            vert->setProperty( "x", latitudes );

            if ( latitudes.size() == longitudes.size() ) 
                vert->setProperty( "y", longitudes );
            else 
            {
                if ( longitudes.size() != 0 ) 
                    std::cout << latitudes.size() << " " << longitudes.size() << std::endl;
                assert(longitudes.size() == 0);
            }
        }

        listVertGLB[blobs[i]->ssb_parent->subject][ blobs[i]->index ] = vert;

        if ( storeMeshes ) 
        {
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
                       bool storeMeshes = true,
                       int representationMode = CORTICAL_PATCHES ) 
{
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    Vertex *vert;
    if ( representationMode == SPHERES ) 
    {
        std::cout << "════ Extracting meshes for the grey-level blobs..." << std::endl;
        std::cout << "      ░░░ mode AimsSphereAtMaxNode ░░░     " << std::endl;
        for ( uint i = 0 ; i < blobs.size() ; i++ )
            blobs[i]->getAimsSphereAtMaxNode( *(subject.tex), 0.4);
    }
    else if ( representationMode == CORTICAL_PATCHES ) 
    {
        std::cout << "════ Extracting meshes for the grey-level blobs..." << std::endl;
        std::cout << "      ░░░ mode AimsMesh ░░░     " << std::endl;
        for ( uint i = 0 ; i < blobs.size() ; i++ )
            blobs[i]->getAimsMesh( *(subject.mesh) );
    }

    std::cout << "════ Adding grey-level blobs..." << std::endl;
    listVertGLB = std::vector<Vertex *> ( blobs.size() );
    for ( int i = 0 ; i < (int) blobs.size() ; i++ ) 
    {
        // For every scale-space blob, we create a vertex in the Aims graph : we define
        //   its properties and store a link between the created vertex and the blob index
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b" << graph->order() << std::flush ;
        vert = graph->addVertex("glb");
        vert->setProperty( "subject", blobs[i]->subject.data() );
        vert->setProperty( "t", blobs[i]->t );
        vert->setProperty( "nodes", blobs[i]->nodes );
        if ( subject.coordinates == LATLON_2D ) 
        {
            std::vector<float> latitudes, longitudes;
            std::set<int>::iterator it;
            for ( it = blobs[i]->nodes.begin() ; it != blobs[i]->nodes.end() ; it++ ) 
            {
                latitudes.push_back( (float) blobs[i]->coordinates[*it][0] );
                if ( blobs[i]->coordinates[*it].size() >= 2 )
                    longitudes.push_back( (float) blobs[i]->coordinates[*it][1] );
            }
            vert->setProperty( "x", latitudes );

            if ( latitudes.size() == longitudes.size() ) 
                vert->setProperty( "y", longitudes );
            else 
            {
                if ( longitudes.size() != 0 ) 
                    std::cout << latitudes.size() << " " << longitudes.size() << std::endl;
                assert(longitudes.size() == 0);
            }
        }
        listVertGLB[ i ] = vert;

        if ( storeMeshes ) 
        {
            // We associate the proper mesh patch from "objects" to the vertex
            ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
            (*ptr)[0] = blobs[i]->mesh;
            manip.storeAims ( *graph, vert, "glb", ptr );
        }
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << graph->order() << " blobs added in total (SSB and GLB)" << std::endl;
}

void addScaleSpaceToGreyLevelBlobsRelations ( Graph *graph,
                                              const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                              std::map<std::string, std::map<int, Vertex *> > &listVertSSB,
                                              std::map<std::string, std::map<int, Vertex *> > &listVertGLB ) 
{
    std::cout << "════ Adding links between SSB and GLB (GLB forming SSB)..." << std::endl;
    uint iNbLinks = 0;
    Vertex *v1, *v2;
    for (int i = 0 ; i < (int) ssblobs.size() ; i++) {

        std::set<surf::GreyLevelBlob *>::iterator itB1;
        std::set<surf::GreyLevelBlob *> &listGLB = ssblobs[i]->blobs;

        for ( itB1 = listGLB.begin(); itB1 != listGLB.end() ; itB1++ ) {
            assert( ssblobs[i]->subject == (*itB1)->subject );
            v1 = listVertSSB[ssblobs[i]->subject][ssblobs[i]->index];
            v2 = listVertGLB[(*itB1)->subject][(*itB1)->index];
            graph->addEdge ( v1, v2, "s2g" );
            iNbLinks++;
        }
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " links added... done" << std::endl;
}

void addBetweenGreyLevelBlobsRelations ( Graph *graph,
                                         const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                         std::map<std::string, std::map<int, Vertex *> > &listVertGLB,
                                         bool buildAndStoreRelationsMeshes ) 
{
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    AimsSurfaceTriangle *relations;

    std::cout << "════ Collecting the list of GLB relations..." << std::endl;
    std::vector< std::pair< surf::GreyLevelBlob *, surf::GreyLevelBlob *> > blobsPairs;

    for ( int i = 0 ; i < (int) ssblobs.size() ; i++ ) 
    {
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

        while ( itB2 != listGLB.end() ) 
        {
            surf::GreyLevelBlob *glb1, *glb2;
            glb1 = (*itB1);
            glb2 = (*itB2);
            std::pair<surf::GreyLevelBlob *, surf::GreyLevelBlob *> p;
            p.first = glb1;
            p.second = glb2;
            blobsPairs.push_back(p);
            itB1++, itB2++;
        }
    }
    std::cout << "  " << blobsPairs.size() << " blobs pairs" << std::endl;

    if ( buildAndStoreRelationsMeshes ) 
    {
        std::cout << "════ Extracting meshes for the GLB relations..." << std::endl;
        relations = new AimsSurfaceTriangle();
        *relations = getG2GRelationsMeshes( blobsPairs, NODES_BARYCENTERS );
    }

    std::cout << "════ Adding links between GLB from same SSB (with corresponding meshes for visualization)..." << std::endl;

    uint iNbLinks = 0;
    for ( uint i = 0 ; i < blobsPairs.size() ; i++ ) 
    {
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b" << i << std::flush;
        // We associate the proper mesh patch from "objects" to the vertex

        Vertex *v1, *v2;
        surf::GreyLevelBlob *glb1, *glb2;
        glb1 = blobsPairs[i].first;
        glb2 = blobsPairs[i].second;
        v1 = listVertGLB[glb1->subject][glb1->index];
        v2 = listVertGLB[glb2->subject][glb2->index];
        Edge *e = graph->addEdge(v1, v2, "g2g");

        if ( buildAndStoreRelationsMeshes ) 
        {
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
                                std::map<std::string, std::map<int, Vertex *> > &listVertSSB,
                                bool buildAndStoreRelationsMeshes ) 
{
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    AimsSurfaceTriangle *bifurcations;

    std::cout << "════ Collecting the list of existing bifurcations..." << std::endl;
    std::vector< std::pair< surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> > bifurcPairs;

    for ( int i = 0 ; i < (int) ssblobs.size() ; i++ ) 
    {
        std::set<surf::ScaleSpaceBlob *>::iterator itB1;
        std::set<surf::ScaleSpaceBlob *> &listSSB = ssblobs[i]->topBlobs;
        for (itB1 = listSSB.begin(); itB1 != listSSB.end() ; itB1++) 
        {
            if ( ssblobs[i]->index < (*itB1)->index ) 
            {
                std::pair<surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> p;
                p.first = ssblobs[i];
                p.second = (*itB1);
                bifurcPairs.push_back(p);
            }
        }
        listSSB = ssblobs[i]->bottomBlobs;
        for (itB1 = listSSB.begin() ; itB1 != listSSB.end() ; itB1++) 
        {
            if ( ssblobs[i]->index < (*itB1)->index ) 
            {
                std::pair<surf::ScaleSpaceBlob *, surf::ScaleSpaceBlob *> p;
                p.first = ssblobs[i];
                p.second = (*itB1);
                bifurcPairs.push_back(p);
            }
        }
    }

    if ( buildAndStoreRelationsMeshes ) 
    {
        std::cout << "════ Extracting meshes for the bifurcations relations..." << std::endl;
        bifurcations = new AimsSurfaceTriangle();
        *bifurcations = getBifurcationRelationsMeshes( bifurcPairs, NODES_BARYCENTERS );
        std::cout << "  " << bifurcPairs.size() << " bifurcations pairs" << std::endl;
        std::cout << "  " << bifurcations->size() << " bifurcations meshes" << std::endl;
    }

    std::cout << "════ Adding links between SSB (bifurcations)..." << std::endl;
    uint iNbLinks = 0;

    for ( uint i = 0 ; i < bifurcPairs.size() ; i++ ) 
    {
        Vertex *v1, *v2;
        surf::ScaleSpaceBlob *ssb1, *ssb2;
        ssb1 = bifurcPairs[i].first;
        ssb2 = bifurcPairs[i].second;
        v1 = listVertSSB[ssb1->subject][ssb1->index];
        v2 = listVertSSB[ssb2->subject][ssb2->index];
        Edge *e = graph->addEdge(v1, v2, "bifurcation");
        e->setProperty("type", "test");

        if ( buildAndStoreRelationsMeshes ) 
        {
            ptr = carto::rc_ptr<AimsSurfaceTriangle>(new AimsSurfaceTriangle);
            (*ptr)[0] = (*bifurcations)[i];
            manip.storeAims(*graph, e, "bifurcation", ptr);
        }

        iNbLinks++;
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b  " << iNbLinks << " bifurcations added... done" << std::endl;
}

void getGraphModeOptions ( const int graph_mode,
                           const int representation_mode,
                           bool &buildScaleSpaceBlobs,
                           int &scaleSpaceBlobsMeshesRepresentationMode,
                           bool &storeScaleSpaceBlobsMeshes,
                           bool &buildGreyLevelBlobs,
                           int &greyLevelBlobsMeshesRepresentationMode,
                           bool &storeGreyLevelBlobsMeshes,
                           bool &buildSSBToGLBRelations,
                           bool &buildGLBRelations,
                           bool &buildAndStoreGLBRelationsMeshes,
                           bool &buildBifurcations,
                           bool &buildAndStoreBifurcationsMeshes ) {

    buildScaleSpaceBlobs = false;
    scaleSpaceBlobsMeshesRepresentationMode = NONE;
    storeScaleSpaceBlobsMeshes = false;
    buildGreyLevelBlobs = false;
    greyLevelBlobsMeshesRepresentationMode = NONE;
    storeGreyLevelBlobsMeshes = false;
    buildSSBToGLBRelations = false;
    buildGLBRelations = false;
    buildAndStoreGLBRelationsMeshes = false;
    buildBifurcations = false;
    buildAndStoreBifurcationsMeshes = false;

    switch ( graph_mode ) 
    {
        case NO_SCALESPACEBLOBS_MESHES:
            buildScaleSpaceBlobs = true;
            buildGreyLevelBlobs = true;
            greyLevelBlobsMeshesRepresentationMode = representation_mode;
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
            greyLevelBlobsMeshesRepresentationMode = representation_mode;
            storeGreyLevelBlobsMeshes = true;
            buildSSBToGLBRelations = true;
        break;
        case DEFAULT:
            buildScaleSpaceBlobs = true;
            scaleSpaceBlobsMeshesRepresentationMode = representation_mode;
            storeScaleSpaceBlobsMeshes = true;
            buildGreyLevelBlobs = true;
            greyLevelBlobsMeshesRepresentationMode = representation_mode;
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
                                 SubjectData &subject,
                                 const std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                 int graph_mode,
                                 int representation_mode ) 
{
     bool buildScaleSpaceBlobs, storeScaleSpaceBlobsMeshes,
          buildGreyLevelBlobs, storeGreyLevelBlobsMeshes,
          buildSSBToGLBRelations, buildGLBRelations, buildAndStoreGLBRelationsMeshes,
          buildBifurcations, buildAndStoreBifurcationsMeshes;
     int scaleSpaceBlobsMeshesRepresentationMode, greyLevelBlobsMeshesRepresentationMode;

     getGraphModeOptions ( graph_mode, representation_mode, buildScaleSpaceBlobs, scaleSpaceBlobsMeshesRepresentationMode,
             storeScaleSpaceBlobsMeshes, buildGreyLevelBlobs, greyLevelBlobsMeshesRepresentationMode,
             storeGreyLevelBlobsMeshes, buildSSBToGLBRelations, buildGLBRelations,
             buildAndStoreGLBRelationsMeshes, buildBifurcations, buildAndStoreBifurcationsMeshes);

    // First We Recover The Grey-Level Blobs From The Vector of Scale-Space Blobs
        std::vector<surf::GreyLevelBlob *> blobs = TextureToBlobs::recoverGreyLevelBlobs ( ssblobs );
        std::cout << blobs.size() << " blobs recovered " << std::endl;

    // We Define The Graph Global Properties
        defineGraphGlobalProperties ( graph, subject );

    // Initializing Blobs Indices
        for ( int i = 0 ; i < (int) ssblobs.size() ; i++ )
            ssblobs[i]->index = i;
        for ( int i = 0 ; i < (int) blobs.size() ; i++ )
            blobs[i]->index = i;

    // Let's Add Or Not The Scale-Space Blobs
        std::map<std::string, std::map<int, Vertex *> > listVertSSB, listVertGLB;
        std::map<std::string, SubjectData *> data;
        std::pair<std::string, SubjectData *> p;
        p.first = subject.subject_id;
        p.second = &subject;
        data.insert(p);
        if ( buildScaleSpaceBlobs ) 
            addScaleSpaceBlobsToGraph ( graph, ssblobs, listVertSSB, data, storeScaleSpaceBlobsMeshes, scaleSpaceBlobsMeshesRepresentationMode );

    // Let's Add Or Not The Grey-Level Blobs
        if ( buildGreyLevelBlobs ) 
            addGreyLevelBlobsToGraph ( graph, blobs, listVertGLB, subject, storeGreyLevelBlobsMeshes, greyLevelBlobsMeshesRepresentationMode );

    // Let's Add Or Not The Links Between Scale-Space Blobs And Grey-Level Blobs
        if ( buildSSBToGLBRelations && ( buildScaleSpaceBlobs && buildGreyLevelBlobs ) ) 
            addScaleSpaceToGreyLevelBlobsRelations ( graph, ssblobs, listVertSSB, listVertGLB );

    // Let's Add Or Not The Links Between Grey-Level Blobs
        if ( buildGLBRelations && ( buildGreyLevelBlobs ) ) 
            addBetweenGreyLevelBlobsRelations ( graph, ssblobs, listVertGLB, buildAndStoreGLBRelationsMeshes );

    // Let's Add Or Not The Links Between Scale-Space Blobs (Bifurcations)
        if ( buildBifurcations && ( buildScaleSpaceBlobs ) ) 
            addBifurcationsRelations ( graph, ssblobs, listVertSSB, buildAndStoreGLBRelationsMeshes );
}

void TextureToBlobs::AimsGraph ( Graph *graph,
                                 SubjectData & subject,
                                 const std::vector<surf::Blob *> &blobs,
                                 int representation_mode ) 
{
     // We Define The Graph Global Properties
         defineGraphGlobalProperties ( graph, subject );
     // Initializing Blobs Indices
         for (int i = 0 ; i < (int) blobs.size() ; i++)
             blobs[i]->index = i;
     // Let's Add Or Not The Blobs
         std::vector<Vertex *> listVert;
         addBlobsToGraph ( graph, blobs, listVert, subject, true, representation_mode );
}

//##############################################################################

void defineGroupGraphGlobalProperties ( Graph *graph,
                                        std::map<std::string, SubjectData *> data ) 
{
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
    for ( it = data.begin() ; it != data.end() ; it ++ ) 
    {
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
                                std::map<std::string, std::map<int, Vertex *> > &listVertSSB,
                                bool buildAndStoreRelationsMeshes = false ) 
{
    carto::rc_ptr<AimsSurfaceTriangle> ptr;
    aims::GraphManip manip;
    AimsSurfaceTriangle *relations;

    if ( buildAndStoreRelationsMeshes ) 
    {
        std::cout << " nodes " << std::endl << "════ Extracting meshes for the interblobs relations..." << std::endl;
        relations = new AimsSurfaceTriangle();
        *relations = getB2BRelationsMeshes( cliques, NODES_BARYCENTERS );
    }


    std::cout << "Building cliques in the Aims group graph..." << std::endl;
    for ( uint i = 0 ; i < cliques.size() ; i++ ) 
    {
        // For every clique, we get the two corresponding vertices from listVertices
        Vertex *v1, *v2;
        v1 = listVertSSB[cliques[i].ssb1->subject][cliques[i].ssb1->index];
        v2 = listVertSSB[cliques[i].ssb2->subject][cliques[i].ssb2->index];
        Edge *edge = graph->addEdge( v1, v2, "b2b" );
        edge->setProperty( "similarity", cliques[i].similarity );
        edge->setProperty( "distance", cliques[i].distance );

        if ( buildAndStoreRelationsMeshes ) 
        {
            ptr = carto::rc_ptr<AimsSurfaceTriangle> ( new AimsSurfaceTriangle );
            (*ptr)[0] = (*relations)[i];
            manip.storeAims ( *graph, edge, "b2b", ptr );
            edge->setProperty( "b2b_label", i );
        }
    }

    if ( buildAndStoreRelationsMeshes )
        delete relations;
}


void TextureToBlobs::AimsGroupGraph ( Graph *graph,
                                      std::map<std::string, SubjectData *> data,
                                      std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                      std::vector<surf::Clique> &cliques,
                                      bool buildAndStoreRelationsMeshes ) 
{
    std::cerr << "Building Group Graph..." << std::endl;
    defineGroupGraphGlobalProperties ( graph, data );
    std::map<std::string, std::map<int, Vertex *> > listVertSSB;
    addScaleSpaceBlobsToGraph ( graph, ssblobs, listVertSSB, data, true, NONE );
    addInterSubjectRelations ( graph, ssblobs, cliques, listVertSSB, buildAndStoreRelationsMeshes );
}

void TextureToBlobs::ReadAimsGroupGraph (   Graph &graph,
                                            std::vector<surf::ScaleSpaceBlob *> &ssblobs,
                                            std::vector<surf::Clique> &cliques ) 
{
    std::set<Vertex *>::iterator iv;
    int index, indexB2, graph_index, graph_indexB2;
    std::string sujet, sujetB2;
    float similarity, distance;
    
    // Recovering The Scale-Space Blobs
    std::cout << " Recovering the scale-space blobs..." << std::endl;
    std::map<std::string, std::map<int, std::set<int> > > listGLBindices;
    TextureToBlobs::getScaleSpaceBlobsFromGraph ( &graph, ssblobs, listGLBindices, false );
    std::cout << "ssblobs.size() :"<< ssblobs.size() << std::endl;

    // Recovering The Clique Links...
    Edge *e;
    Vertex::iterator jv;
    Edge::iterator kv;

    std::cout << " Recovering the similarity cliques..." << std::endl;
    for ( iv = graph.vertices().begin() ; iv != graph.vertices().end() ; ++iv ) 
    {
        if ( (*iv)->getSyntax() == "ssb" )
        {
            (*iv)->getProperty( "index", index );
            (*iv)->getProperty( "sites_index", graph_index );
            (*iv)->getProperty( "subject", sujet );

            // Iterating On The Edges Of The Current Vertex
            for ( jv = (*iv)->begin() ; jv != (*iv)->end() ; jv++ )
            {
                e = *jv;

                if ( e->getSyntax() == "b2b" ){

                    e->getProperty( "similarity", similarity );
                    e->getProperty( "distance", distance );
                    assert (!( e->hasProperty("similarity") && e->hasProperty("distance") ) || (similarity == -1.0 || distance == -1.0));
                    if ( ! e->hasProperty("similarity") )
                        similarity = -1.0;
                    else if ( ! e->hasProperty("distance") )
                        distance = -1.0;

                    for ( kv = e->begin() ; kv != e->end() ; kv++ )
                    {
                        if ( (*kv)->getSyntax() == "ssb" )
                        {
                            (*kv)->getProperty( "index", indexB2 );
                            (*kv)->getProperty( "sites_index", graph_indexB2 );
                            (*kv)->getProperty( "subject", sujetB2 );

                            if ( !(indexB2 == index && sujetB2 == sujet) )
                            {
                                assert( sujet != sujetB2 );
                                if ( graph_index < graph_indexB2 )
                                    cliques.push_back( surf::Clique( ssblobs[ graph_index ] , ssblobs[ graph_indexB2 ], distance, similarity ) );
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "ssbcliques.size() : "<< cliques.size() << std::endl << std::endl;
}


std::vector<uint> TextureToBlobs::getClustersListsFromGLB ( std::vector<surf::GreyLevelBlob *> &blobs,
                                                            GroupData &data,
                                                            float clustering_distance_threshold ) 
{
    std::vector<uint>  clusters ( blobs.size() );
    for ( uint i = 0 ; i < blobs.size() ; i++ )
        clusters[i] = i;

    std::set<uint> dejapris;
    for ( uint i = 0 ; i < blobs.size() ; i++ ) 
    {
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << blobs.size() << std::flush ;
        uint max_node = blobs[i]->getMaximumNode( * (data[blobs[i]->subject]->tex) );
        std::map<uint, float> distanceMap  = LocalMeshDistanceMap(
                    data[blobs[i]->subject]->mesh,
                    data[blobs[i]->subject]->neighbours,
                    max_node,
                    clustering_distance_threshold + 10.0 );
        // 10.0mm being a distance larger than the maximum internode distance

        std::map<uint, float>::iterator it;
        std::set<uint> voisins_max_node;
        for ( it = distanceMap.begin() ; it != distanceMap.end() ; it ++ )
            if ( it->second < clustering_distance_threshold )
                voisins_max_node.insert( it->first );

        voisins_max_node.insert( max_node );

        std::set<uint> blobs_voisins;
        for ( uint j = 0 ; j < blobs.size() ; j++ ) 
        {
            uint max_node2 = blobs[j]->getMaximumNode( * (data[blobs[j]->subject]->tex) );
            if ( blobs[i]->subject == blobs[j]->subject &&
                        voisins_max_node.find(max_node2) != voisins_max_node.end() ) 
                blobs_voisins.insert( j);
        }

        assert( blobs_voisins.find(i) != blobs_voisins.end() );
        std::set<uint>::iterator ite;
        uint color = clusters[i];
        for ( ite = blobs_voisins.begin() ; ite != blobs_voisins.end() ; ite ++ ) 
            clusters[*ite] = color;
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
                                   bool uniqueGLB ) 
{
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

        // Creation of clustered ssblobs
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
            glb->subject = blobs[cluster_blobs[0]]->subject;

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
            ssb->getNodesFromBlob(glb);
        }
        else 
        {
            //   Or should we keep track of original GLB ?
            for ( uint i = 0 ; i < cluster_blobs.size() ; i ++ ) {
                surf::GreyLevelBlob *glb = new surf::GreyLevelBlob();
                ssb->blobs.insert( glb );
                glb->t = blobs[cluster_blobs[i] ]->t;
                glb->scale = blobs[cluster_blobs[i] ]->scale;
                glb->ssb_parent = ssb;
                glb->subject = blobs[cluster_blobs[i] ]->subject;

                std::set<int>::iterator ite2;
                for ( ite2 = blobs[ cluster_blobs[i] ]->nodes.begin() ;
                    ite2 != blobs[ cluster_blobs[i] ]->nodes.end() ;
                    ite2++ ) 
                {
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
    }

    if ( outputTextFile != "" )
        fclose(f1);
}

double TextureToBlobs::getOverlapMeasure( Point2df bbmin1,
                          Point2df bbmax1,
                          Point2df bbmin2,
                          Point2df bbmax2,
                          uint *no_overlap )
{
    float overlap_x, overlap_y, aux;
    double rec = 0.0;

    if ( sqrt(pow( bbmin1[0] - bbmax1[0], 2) ) < 0.0001 )  bbmax1[0] += 0.5;
    if ( sqrt(pow( bbmin1[1] - bbmax1[1], 2) ) < 0.0001 )  bbmax1[1] += 0.5;
    if ( sqrt(pow( bbmin2[0] - bbmax2[0], 2) ) < 0.0001 )  bbmax2[0] += 0.5;
    if ( sqrt(pow( bbmin2[1] - bbmax2[1], 2) ) < 0.0001 )  bbmax2[1] += 0.5;

    if (sqrt(pow(bbmin1[1] - bbmax1[1], 2)) > 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) < 300) 
    {
        if ( 360 - bbmax2[1] < bbmin2[1] ) 
        {
            aux = bbmax1[1];
            bbmax1[1] = bbmin1[1] + 360.0;
            bbmin1[1] = aux;
        }
        else 
        {
            aux = bbmin1[1];
            bbmin1[1] = bbmax1[1] - 360.0;
            bbmax1[1] = aux;
        }
    }
    else if (sqrt(pow( bbmin1[1] - bbmax1[1], 2)) < 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) > 300) 
    {
        if ( 360 - bbmax1[1] < bbmin1[1] ) 
        {
            aux = bbmax2[1];
            bbmax2[1] = bbmin2[1] + 360.0;
            bbmin2[1] = aux;
        }
        else 
        {
            aux = bbmin2[1];
            bbmin2[1] = bbmax2[1] - 360.0;
            bbmax2[1] = aux;
        }
    }
    else if (sqrt(pow( bbmin1[1] - bbmax1[1], 2)) > 300 && sqrt(pow( bbmin2[1] - bbmax2[1], 2)) > 300) 
    {
        aux = bbmin1[1];
        bbmin1[1] = bbmax1[1] - 360.0;
        bbmax1[1] = aux;
        aux = bbmin2[1];
        bbmin2[1] = bbmax2[1] - 360.0;
        bbmax2[1] = aux;
    }

    // ON S'OCCUPE DE LA LATITUDE
    if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) > 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) < 150) 
    {
        if ( 180 - bbmax2[0] < bbmin2[0] ) 
        {
            aux = bbmax1[0];
            bbmax1[0] = bbmin1[0] + 180.0;
            bbmin1[0] = aux;
        }
        else 
        {
            aux = bbmin1[0];
            bbmin1[0] = bbmax1[0] - 180.0;
            bbmax1[0] = aux;
        }
    }
    else if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) < 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) > 150)
    {
        if ( 180 - bbmax1[0] < bbmin1[0] ) 
        {
            aux = bbmax2[0];
            bbmax2[0] = bbmin2[0] + 180.0;
            bbmin2[0] = aux;
        }
        else 
        {
            aux = bbmin2[0];
            bbmin2[0] = bbmax2[0] - 180.0;
            bbmax2[0] = aux;
        }
    }
    else if (sqrt(pow( bbmin1[0] - bbmax1[0], 2)) > 150 && sqrt(pow( bbmin2[0] - bbmax2[0], 2)) > 150)
    {
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

    if ( *no_overlap == 0 ) 
    {
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
            double div =  ( bbmax1[0] - bbmin1[0] ) * ( bbmax1[1] - bbmin1[1] )
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


//void TextureToBlobs::filteringBlobs (  std::vector<surf::ScaleSpaceBlob *> & ssblobs,
//                       std::vector<surf::GreyLevelBlob *> &filteredBlobs,
//                       std::vector<surf::ScaleSpaceBlob *> & filteredSsblobs,
//                       std::set<int> &nodes ){
//    // Not much used : allows to filter in blobs provided their support intersects
//    //    with nodes belonging to a given set
//
//    for (uint i = 0 ; i < ssblobs.size() ; i++ )
//        ssblobs [i] -> index = i;
//
//    std::set< uint > filteredIndices;
//
//    // Filtering according to positions
//    for (uint i = 0 ; i < ssblobs.size() ; i++ ){
//        std::string subject;
//        subject = ssblobs[i]->subject;
//        bool firstGLB = false;
//
//        surf::ScaleSpaceBlob *ssb;
//        ssb = new surf::ScaleSpaceBlob( ssblobs[i] );
//        ssb->blobs.clear();
//        ssb->topBlobs.clear();
//        ssb->bottomBlobs.clear();
//
//        std::set<surf::GreyLevelBlob *>::iterator itB1;
//        for (itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() && !firstGLB ; itB1++){
//            uint no_overlap = 2;
//
//            std::set<int>::iterator it;
//            std::set<int> intersection;
//            for ( it = (*itB1)->nodes.begin() ; it != (*itB1)->nodes.end() ; it++ )
//                if ( nodes.find(*it) != nodes.end() )
//                    intersection.insert(*it);
//
//            if ( intersection.size() != 0 )
//                no_overlap = 0;
//            else
//                no_overlap = 1;
////             }
//            if (no_overlap == 0)
//                firstGLB = true;
//        }
//
//        if (firstGLB) {
//            for ( itB1 = ssblobs[i]->blobs.begin() ; itB1 != ssblobs[i]->blobs.end() ; itB1++ ) {
//
//                surf::GreyLevelBlob *glb;
//                glb = new surf::GreyLevelBlob( *itB1 );
//                ASSERT(glb->nodes.size() == glb->raw_coordinates.size());
//                glb->ssb_parent = ssb;
//                ssb->blobs.insert(glb);
//                filteredBlobs.push_back(glb);
//            }
//            filteredSsblobs.push_back(ssb);
//            filteredIndices.insert ( findBlobIndex( ssblobs, ssb->subject, ssb->index ) );
//        }
//        else
//            delete(ssb);
//    }
//    std::cerr << filteredBlobs.size() << " filtered blobs - " << filteredSsblobs.size() << " filtered ssblobs" << std::endl;
//
//    // Now that the blobs are filtered, we add the correct bifurcations
//
//    for ( uint i = 0 ; i < filteredSsblobs.size() ; i ++ ) {
//
//        std::set<surf::ScaleSpaceBlob *> auxTop = findBlob(ssblobs, filteredSsblobs[i]->subject, filteredSsblobs[i]->index)->topBlobs;
//        std::set<surf::ScaleSpaceBlob *>::iterator it;
//        for ( it = auxTop.begin() ; it != auxTop.end() ; it ++ ) {
//            int bi = findBlobIndex ( filteredSsblobs, (*it)->subject, (*it)->index );
//            if ( bi != -1 )
//                filteredSsblobs[i]->topBlobs.insert( filteredSsblobs[bi] );
//        }
//
//        std::set< surf::ScaleSpaceBlob *> auxBot = findBlob(ssblobs, filteredSsblobs[i]->subject, filteredSsblobs[i]->index)->bottomBlobs;
//        for ( it = auxBot.begin() ; it != auxBot.end() ; it ++ ) {
//            int bi = findBlobIndex ( filteredSsblobs, (*it)->subject, (*it)->index ); 
//            if ( bi != -1 )            
//                filteredSsblobs[i]->bottomBlobs.insert(filteredSsblobs[bi]);
//        }
//    }
//}
