#include <cortical_surface/structuralanalysis/subjectdata.h>
#include <aims/data/data_g.h>
#include <aims/io/io_g.h>




void GroupData::readData( ) {

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
        std::cout << subject_id << std::endl;
        p.second = subject;
        subjects_id[subject_id] = i;
        this->insert ( p );

        (*this)[subject_id]->storeData ( & (meshes[i][0]), & (textures[i][0]), false, true );
        (*this)[subject_id]->graph = graph;

        if ( latPath != "" ) {
           aims::Reader<TimeTexture<float> > latRdr ( latPath );
           latRdr.read( latitudes[i] );

           if ( lonPath != "" ) {
               aims::Reader<TimeTexture<float> > lonRdr ( lonPath );
               lonRdr.read( longitudes[i] );
               (*this)[subject_id]->storeCoordinates ( & (latitudes[i][0]), & (longitudes[i][0]) );
           }
           else
               (*this)[subject_id]->storeCoordinates ( & (latitudes[i][0]) );
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

std::vector<std::map<uint,float> > SubjectData::GetSecondOrderNeighbours ( ) {
    assert( this->neighbours.size() > 0 );
    std::vector<std::set<uint> > voisins ( this->neighbours );
    std::set<uint>::iterator it, it2, it3, it4, it_min;
    std::set<uint> vois2;
    std::vector<std::map<uint,float> > voisins2( voisins.size() );

    std::cerr << "Computing the 2nd-order neighbors :" << std::endl;
    for ( uint i = 0 ; i < voisins.size() ; i++ ) {
        if ( i%1000 == 0 ) std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b" << i << "/" << voisins.size() << std::flush;
        for ( it = voisins[i].begin() ; it != voisins[i].end() ; it++ ) { // on parcourt les voisins de i
            voisins2[i][*it] = calcule_distance( mesh[0].vertex()[i], mesh[0].vertex()[*it] );
            vois2.clear();

            for ( it2 = voisins[*it].begin() ; it2 != voisins[*it].end() ; it2++ ) { // on parcourt les voisins de *it
                if ( voisins[i].find(*it2) == voisins[i].end() && *it2 != i ) { // il ne faut pas que *it2 soit un voisin de premier ordre de i ni i
                    vois2.insert(*it2);
                }
            }
            for ( it2 = vois2.begin() ; it2 != vois2.end() ; it2++ ) {
                float distance = 1000.0;
                for ( it3 = voisins[*it2].begin() ; it3 != voisins[*it2].end() ; it3++ ) {
                    it4 = voisins[i].find(*it3);
                    if ( it4 != voisins[i].end() ) { // on considère les voisins de *it2 qui sont aussi voisins de i
                        float aux = calcule_distance( mesh[0].vertex()[*it2],
                                        mesh[0].vertex()[*it4]) + calcule_distance ( mesh[0].vertex()[*it4], mesh[0].vertex()[i] );
                        if ( aux < distance )
                            distance = aux;
                    }
                }
                voisins2[i][*it2] = distance;
//                std::cout << distance << " " << std::flush;
            }
        }
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b" << voisins.size() << "/" << voisins.size() << std::endl;
    return voisins2;

}

