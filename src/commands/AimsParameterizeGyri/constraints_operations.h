
#ifndef AIMS_PARAMETERIZEGYRI_CONSTRAINTS_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_CONSTRAINTS_OPERATIONS_H

using namespace aims;
using namespace std;


pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr,
      const vector<set<uint> > &voisins, const Texture<short> &inTex);

pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr, const vector<set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &vertTex, const Texture<double> &horizTex);

void init_constvector(vector<uint> &haut, const vector<uint> &prohaut, const Texture<double> &horizTex);

double getRatioValue(uint valX, uint val1, uint val2, double valA, double valB);

void update_result(vector<pair<vector<uint>,short> > &result, const vector<uint> &haut,
      const vector<uint> &bas, const vector<uint> &prohaut, const vector<uint> &probas, const vector<uint> &corr,
      AimsSurface<3,Void> &flatMesh, const Texture<double> &horizTex);

pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr, const vector<set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, const Texture<double> &vertTex, const Texture<double> &horizTex);

pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr, const vector<set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, const Texture<double> &vertTex, const Texture<double> &horizTex,
      const Texture<float> &spmTex);

pair<vector<pair<vector<uint>,short> >,vector<pair<vector<uint>,short> > > getConstraints(short gyrus,
      const vector<vector<uint> > &constrMod, const pair<vector<uint>,vector<uint> > &hautBas,
      const pair<vector<uint>,vector<uint> > &gaucheDroite, const vector<uint> &corr, const vector<set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &vertTex,
      const Texture<double> &horizTex, const Texture<float> &spmTex);



#endif


