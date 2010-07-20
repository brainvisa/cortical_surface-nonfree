
#ifndef AIMS_PARAMETERIZEGYRI_MESH_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_MESH_OPERATIONS_H

using namespace std;


AimsSurfaceTriangle getGyrusMesh(AimsSurface<3, Void> &inMesh, const std::vector<uint> &gyrusVertices, std::vector<uint> &corres);

std::map<unsigned, std::set<pair<unsigned,float> > > getGyrusWeight (std::map<unsigned, std::set<std::pair<unsigned,float> > > &poids, std::vector<uint> &gyrusVertices, std::vector<uint> &corres);

AimsSurface<3,Void> getFlatMesh(const AimsSurface<3,Void> &mesh, const std::vector<uint> &vertices, const std::vector<uint> &corres,
      TimeTexture<float> &paramTex);

Texture<double> AimsMeshLaplacian( const Texture<double> &inittex,
                      const std::map<unsigned, std::set< std::pair<unsigned,float> > > &lapl);

Texture<double> diffusion(std::map<unsigned, std::set<std::pair<unsigned,float> > > &poidsGyrus,
      AimsSurface<3,Void> &mesh_base, std::vector<uint> &haut, std::vector<uint> &bas, const std::vector<std::pair<std::vector<uint>,short> > &constraints
      , short init, std::vector<uint> &gyrusvertices, std::vector<uint> &corr, const float criter, const float dt);



#endif

