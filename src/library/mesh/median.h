 

#ifndef AIMS_MESH_MEDIAN_H
#define AIMS_MESH_MEDIAN_H

using namespace aims;
using namespace std;

set<uint> nearest_vertices(Point3df pf, AimsSurfaceTriangle &mesh, float rayon);

pair<int,float> plus_proche_point_normal(Point3df p, Point3df n, AimsSurfaceTriangle &mesh, set<uint> &vertices);

pair<int,float> plus_proche_point_normal(Point3df p, Point3df n, AimsSurfaceTriangle &mesh);

pair<Point3df,bool> isInsideTriangle(const Point3df &pt, const Point3df &vertex0, const Point3df &vertex1, const Point3df &vertex2);

pair<Point3df, int> plus_proche_point_sur_triangle(Point3df p, Point3df n, AimsSurfaceTriangle &mesh, set<uint> &vertices, vector<set<uint> > &voisins);

pair<Point3df, int> plus_proche_point_sur_triangle(Point3df p, Point3df n, AimsSurfaceTriangle &mesh, vector<set<uint> > &voisins);

vector<set<uint> > readVoisinsFromDisk(string path);

void writeVoisinsToDisk(const string &path, vector<set<uint> > &extvoisins4);

vector<set<uint> > compute_neighbours_order(AimsSurfaceTriangle extmesh, uint order);

pair<AimsSurfaceTriangle, TimeTexture<float> > build_median_surface(AimsSurfaceTriangle &intmesh, AimsSurfaceTriangle &extmesh, vector<set<uint> > &extvoisins4, int op);

#endif

