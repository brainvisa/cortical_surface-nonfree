
#ifndef AIMS_PARAMETERIZEGYRI_VERIF_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_VERIF_OPERATIONS_H


bool verifMeshTexMatch(const AimsSurface<3,Void> &inMesh, const Texture<short> &inTex);

bool verifTableTexMatch(const Texture<short> &inTex, const std::vector<uint> &corres);

bool verifVoisins(const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex, const std::vector<uint> &corres);

bool verifIntersectionsModel(const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex, const std::vector<uint> &corres);

bool verifGyriExistence(short g, const Point3d &gHaut, const Point3d &gBas, const std::vector<uint> &corres);



#endif

