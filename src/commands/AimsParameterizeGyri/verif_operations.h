
#ifndef AIMS_PARAMETERIZEGYRI_VERIF_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_VERIF_OPERATIONS_H

using namespace aims;
using namespace std;

bool verifMeshTexMatch(const AimsSurface<3,Void> &inMesh, const Texture<short> &inTex);

bool verifTableTexMatch(const Texture<short> &inTex, const vector<uint> &corres);

bool verifVoisins(const vector<set<uint> > &voisins, const Texture<short> &inTex, const vector<uint> &corres);

bool verifIntersectionsModel(const vector<set<uint> > &voisins, const Texture<short> &inTex, const vector<uint> &corres);

bool verifGyriExistence(short g, const Point3d &gHaut, const Point3d &gBas, const vector<uint> &corres);



#endif

