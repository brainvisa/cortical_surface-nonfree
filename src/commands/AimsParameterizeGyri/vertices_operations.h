
#ifndef AIMS_PARAMETERIZEGYRI_VERTICES_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_VERTICES_OPERATIONS_H

using namespace aims;
using namespace std;

vector<vector<uint> > sortVertices(uint gyruslabel, const vector<set<uint> > &voisins, const Texture<short> &inTex) ;

vector<uint> getIntersection(short gyruslabel, short voisin1, short voisin2, const vector<set<uint> > &voisins,
      const Texture<short> &inTex);

vector<uint> getRealIntersection(short gyruslabel, short voisin1, short voisin2, const vector<set<uint> > &voisins,
      const Texture<short> &inTex);

vector<uint> getIntersection(short gyruslabel, short voisin1, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<uint> getThirdPoints(uint a, uint b, const vector<set<uint> > &voisins);

vector<uint> getBorderNeighbours(uint v, const vector<uint> &borderVertices, const vector<set<uint> > &voisins);

double getdistance(Point3df &a, Point3df &b);

bool distanceCompare(Point3df &a, Point3df &b, Point3df &c);

vector<uint> isolineExtraction(double value, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &inTex);

vector<uint> lineExtraction2(uint start, uint end, const pair<vector<uint>, vector<uint> > &hautBas, AimsSurface<3,Void> &gyrusSurf);

Point3df cross(Point3df p1, Point3df p2);

vector<uint> lineExtraction(uint start, uint end, const pair<vector<uint>, vector<uint> > &hautBas, AimsSurface<3,Void> &gyrusSurf);

vector<uint> lineExtraction(uint start, uint end, AimsSurface<3,Void> &gyrusSurf);

uint getNearestPoint(const vector<uint> &vec, const Texture<double> &inTex, double value);



#endif

