#include <cortical_surface/structuralanalysis/subjectdata.h>
#include <aims/data/data_g.h>
#include <aims/io/io_g.h>




void GroupData::readData( bool verbose ) {

    meshes = std::vector < AimsSurfaceTriangle > ( graphs.size() );
    textures = std::vector < TimeTexture<float> > ( graphs.size() );
    latitudes = std::vector < TimeTexture<float> > ( graphs.size() );
    longitudes = std::vector < TimeTexture<float> > ( graphs.size() );
    subjects_id = std::map< std::string, uint > ();

    std::string meshPath, subject_id, texPath, latPath, lonPath;
    std::map<std::string, SubjectData*>::iterator it;

    // Processing every subject...
    for (uint i = 0 ; i < graphs.size() ; i++) {

        Graph *graph = graphs[i];

        graph->getProperty("mesh", meshPath);
        graph->getProperty("subject", subject_id);
        graph->getProperty("texture", texPath);
        graph->getProperty("latitude", latPath);
        graph->getProperty("longitude", lonPath);

        aims::Reader<AimsSurfaceTriangle> meshRdr ( meshPath ) ;
        meshRdr.read( meshes[i] );
        aims::Reader<TimeTexture<float> > texRdr ( texPath ) ;
        texRdr.read( textures[i] );

        SubjectData *subject = new SubjectData ( subject_id, meshPath, texPath, latPath, lonPath );
        std::pair< std::string, SubjectData *> p;
        p.first = subject_id ;
        if ( verbose )
            std::cout << subject_id << std::endl;
        p.second = subject;
        subjects_id[subject_id] = i;
        this->insert ( p );

        (*this)[subject_id]->storeData ( & (meshes[i][0]), & (textures[i][0]), false, true, verbose );
        (*this)[subject_id]->graph = graph;

        if ( latPath != "" ) {
           aims::Reader<TimeTexture<float> > latRdr ( latPath );
           latRdr.read( latitudes[i] );

           if ( lonPath != "" ) {
               aims::Reader<TimeTexture<float> > lonRdr ( lonPath );
               lonRdr.read( longitudes[i] );
               (*this)[subject_id]->storeCoordinates ( & (latitudes[i][0]), & (longitudes[i][0]), verbose );
           }
           else
               (*this)[subject_id]->storeCoordinates ( & (latitudes[i][0]), NULL, verbose );
        }
    }
}

inline float calcule_distance(const Point3df &p, const Point3df &t){
    Point3df aux(p-t);
    return aux.dnorm();
}

inline float calcule_distance(const Point3df &p, const Point3d &t){
    Point3df aux(p);
    aux[0] -= t[0]; aux[1] -= t[1];  aux[2] -= t[2];
    return aux.dnorm();
}


