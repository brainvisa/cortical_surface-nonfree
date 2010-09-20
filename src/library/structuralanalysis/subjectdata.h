#ifndef SUBJECTDATA_H_
#define SUBJECTDATA_H_

#include <aims/getopt/getopt2.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/data/data_g.h>
#include <aims/mesh/surfaceOperation.h>





enum coordinatesMode{
    RAW_3D, LATLON_2D, LAT_1D
};

class SubjectPaths{
    public:

        std::string meshPath;
        std::string texPath;
        std::string latPath;
        std::string lonPath;
        std::string graphPath;

        SubjectPaths ( ) { meshPath = ""; texPath = ""; latPath = ""; lonPath = ""; }
        SubjectPaths ( SubjectPaths &s ) {
            meshPath = s.meshPath;
            texPath = s.texPath;
            latPath = s.latPath;
            lonPath = s.lonPath;
        }
        SubjectPaths ( std::string _meshPath, std::string _texPath, std::string _latPath = "", std::string _lonPath = "" )
            { meshPath = _meshPath; texPath = _texPath; latPath = _latPath; lonPath = _lonPath; }

};


class SubjectData{
    public :

        int coordinates;
        std::string subject_id;
        SubjectPaths paths;

        AimsSurface<3,Void> *mesh;

        Texture<float> *tex;
        Texture<float> *lat;
        Texture<float> *lon;
        std::map<unsigned, std::set< std::pair<unsigned,float> > > weightLapl;
        std::vector<std::set<uint> > neighbours;

        Graph *graph;

        SubjectData ( std::string id = "UNKNOWN_SUBJECT_ID" ) { subject_id = id; }
        SubjectData ( SubjectData &subject ) {
            subject_id = subject.subject_id;
            coordinates = subject.coordinates;
            paths = subject.paths;
            mesh = subject.mesh;
            tex = subject.tex;
            lat = subject.lat;
            lon = subject.lon;
            weightLapl = subject.weightLapl;
        }



        SubjectData ( std::string id, std::string meshPath, std::string texPath, std::string latPath = "", std::string lonPath = "" ) {
            subject_id = id;
            paths.meshPath = meshPath;
            paths.texPath = texPath;
            paths.latPath = latPath;
            paths.lonPath = lonPath;
            mesh = NULL; tex = NULL; lat = NULL; lon = NULL;
        }

//        void readData ( std::string _meshPath, std::string _texPath, AimsSurfaceTriangle) ;


        void storeData ( AimsSurface<3,Void> *mesh, Texture<float> *tex, bool computeWeights = true, bool computeNeighbours = true ) {
            this->mesh = mesh;
            this->tex = tex;

            std::cout << "  mesh : " << this->mesh->vertex().size() << " vertices" << std::endl;
            std::cout << "  tex : " << this->tex->nItem() << " values" << std::endl;

            this->coordinates = RAW_3D;
            if ( computeWeights ) {
                this->weightLapl = AimsMeshWeightFiniteElementLaplacian ( *(this->mesh), 0.98 );
                std::cout << "  weights : " << this->weightLapl.size() << " weights" << std::endl;
            }
            if ( computeNeighbours ) {
                AimsSurfaceTriangle mesh;
                mesh[0] = *(this->mesh);
                this->neighbours = aims::SurfaceManip::surfaceNeighbours( mesh );
                std::cout << "  neighbours : " << this->neighbours.size() << " sets" << std::endl;

            }

        }

        void storeData ( AimsSurface<3,Void> *mesh, Texture<float> *tex, std::map<unsigned, std::set< std::pair<unsigned,float> > > & weights ) {
            storeData ( mesh, tex, false );
            this->weightLapl = weights;
        }

        void storeCoordinates ( Texture<float> *lat, Texture<float> *lon = NULL ) {
            this->lat = lat;
            this->lon = lon;


            std::cout << "  lat : " << this->lat->nItem() << " values" << std::endl;


            if ( this->lon != NULL ) {
                std::cout << "  lon : " << this->lon->nItem() << " values" << std::endl;
                this->coordinates = LATLON_2D;
            }
            else {
                //bool isCoordinatesOrFilteringTexture = false;
                //std::set<float> possibleValues;
                //possibleValues.insert(-1.0);
                //possibleValues.insert(-2.0);
                //possibleValues.insert(0.0);
                //possibleValues.insert(1.0);
                //for ( uint i = 0 ; i < this->lat->nItem() && !isCoordinatesOrFiltering; i ++ ) {
                //    if (this->lat->item(i) != -1 && this->
                //
                //}
                this->coordinates = LAT_1D;
            }
        }
        std::vector<std::map<uint, float> > GetSecondOrderNeighbours();

};

class GroupData:
    public std::map< std::string, SubjectData *> {
    public:
      void readData( void );
      std::map< std::string, uint > subjects_id;
      std::vector < AimsSurfaceTriangle > meshes;
      std::vector < TimeTexture<float> > textures;
      std::vector < TimeTexture<float> > latitudes;
      std::vector < TimeTexture<float> > longitudes;
      std::vector < Graph *> graphs;

};





#endif /*SUBJECTDATA_H_*/
