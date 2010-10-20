
#ifndef AIMS_PARAMETERIZEGYRI_VERTICES_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_VERTICES_OPERATIONS_H

#include <aims/mesh/surface.h>


std::vector<std::vector<uint> > sortVertices(uint gyruslabel, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex) ;

std::vector<uint> getIntersection(short gyruslabel, short voisin1, short voisin2, const std::vector<std::set<uint> > &voisins,
      const Texture<short> &inTex);

std::vector<uint> getRealIntersection(short gyruslabel, short voisin1, short voisin2, const std::vector<std::set<uint> > &voisins,
      const Texture<short> &inTex);

std::vector<uint> getIntersection(short gyruslabel, short voisin1, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<uint> getThirdPoints(uint a, uint b, const std::vector<std::set<uint> > &voisins);

std::vector<uint> getBorderNeighbours(uint v, const std::vector<uint> &borderVertices, const std::vector<std::set<uint> > &voisins);

double getdistance(Point3df &a, Point3df &b);

bool distanceCompare(Point3df &a, Point3df &b, Point3df &c);

std::vector<uint> isolineExtraction(double value, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &inTex);

std::vector<uint> lineExtraction2(uint start, uint end, const std::pair<std::vector<uint>, std::vector<uint> > &hautBas, AimsSurface<3,Void> &gyrusSurf);

Point3df cross(Point3df p1, Point3df p2);

std::vector<uint> lineExtraction(uint start, uint end, const std::pair<std::vector<uint>, std::vector<uint> > &hautBas, AimsSurface<3,Void> &gyrusSurf);

std::vector<uint> lineExtraction(uint start, uint end, AimsSurface<3,Void> &gyrusSurf);

uint getNearestPoint(const std::vector<uint> &vec, const Texture<double> &inTex, double value);



#endif

