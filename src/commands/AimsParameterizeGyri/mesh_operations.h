
#ifndef AIMS_PARAMETERIZEGYRI_MESH_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_MESH_OPERATIONS_H

using namespace aims;
using namespace std;


AimsSurfaceTriangle getGyrusMesh(AimsSurface<3, Void> &inMesh, const vector<uint> &gyrusVertices, vector<uint> &corres);

map<unsigned, set<pair<unsigned,float> > > getGyrusWeight (map<unsigned, set<pair<unsigned,float> > > &poids, vector<uint> &gyrusVertices, vector<uint> &corres);

AimsSurface<3,Void> getFlatMesh(const AimsSurface<3,Void> &mesh, const vector<uint> &vertices, const vector<uint> &corres,
      TimeTexture<float> &paramTex);

Texture<double> AimsMeshLaplacian( const Texture<double> &inittex,
                      const map<unsigned, set< pair<unsigned,float> > > &lapl);

Texture<double> diffusion(map<unsigned, set<pair<unsigned,float> > > &poidsGyrus,
      AimsSurface<3,Void> &mesh_base, vector<uint> &haut, vector<uint> &bas, const vector<pair<vector<uint>,short> > &constraints
      , short init, vector<uint> &gyrusvertices, vector<uint> &corr, const float criter, const float dt);



#endif

