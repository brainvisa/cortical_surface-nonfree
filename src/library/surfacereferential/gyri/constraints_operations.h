
#ifndef AIMS_PARAMETERIZEGYRI_CONSTRAINTS_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_CONSTRAINTS_OPERATIONS_H



std::pair<std::vector<std::pair<std::vector<uint>,short> >,std::vector<std::pair<std::vector<uint>,short> > > getConstraints(short gyrus,
      const std::vector<std::vector<uint> > &constrMod, const std::pair<std::vector<uint>,std::vector<uint> > &hautBas,
      const std::pair<std::vector<uint>,std::vector<uint> > &gaucheDroite, const std::vector<uint> &corr,
      const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::pair<std::vector<std::pair<std::vector<uint>,short> >,std::vector<std::pair<std::vector<uint>,short> > > getConstraints(short gyrus,
      const std::vector<std::vector<uint> > &constrMod, const std::pair<std::vector<uint>,std::vector<uint> > &hautBas,
      const std::pair<std::vector<uint>,std::vector<uint> > &gaucheDroite, const std::vector<uint> &corr, const std::vector<std::set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &vertTex, const Texture<double> &horizTex);

void init_constvector(std::vector<uint> &haut, const std::vector<uint> &prohaut, const Texture<double> &horizTex);

double getRatioValue(uint valX, uint val1, uint val2, double valA, double valB);

void update_result(std::vector<std::pair<std::vector<uint>,short> > &result, const std::vector<uint> &haut,
      const std::vector<uint> &bas, const std::vector<uint> &prohaut, const std::vector<uint> &probas, const std::vector<uint> &corr,
      AimsSurface<3,Void> &flatMesh, const Texture<double> &horizTex);

std::pair<std::vector<std::pair<std::vector<uint>,short> >,std::vector<std::pair<std::vector<uint>,short> > > getConstraints(short gyrus,
      const std::vector<std::vector<uint> > &constrMod, const std::pair<std::vector<uint>,std::vector<uint> > &hautBas,
      const std::pair<std::vector<uint>,std::vector<uint> > &gaucheDroite, const std::vector<uint> &corr, const std::vector<std::set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, const Texture<double> &vertTex, const Texture<double> &horizTex);

std::pair<std::vector<std::pair<std::vector<uint>,short> >,std::vector<std::pair<std::vector<uint>,short> > > getConstraints(short gyrus,
      const std::vector<std::vector<uint> > &constrMod, const std::pair<std::vector<uint>,std::vector<uint> > &hautBas,
      const std::pair<std::vector<uint>,std::vector<uint> > &gaucheDroite, const std::vector<uint> &corr, const std::vector<std::set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, const Texture<double> &vertTex, const Texture<double> &horizTex,
      const Texture<float> &spmTex);

std::pair<std::vector<std::pair<std::vector<uint>,short> >,std::vector<std::pair<std::vector<uint>,short> > > getConstraints(short gyrus,
      const std::vector<std::vector<uint> > &constrMod, const std::pair<std::vector<uint>,std::vector<uint> > &hautBas,
      const std::pair<std::vector<uint>,std::vector<uint> > &gaucheDroite, const std::vector<uint> &corr, const std::vector<std::set<uint> > &voisins,
      const Texture<short> &inTex, AimsSurface<3,Void> &flatMesh, AimsSurfaceTriangle &gyrusMesh, const Texture<double> &vertTex,
      const Texture<double> &horizTex, const Texture<float> &spmTex);



#endif


