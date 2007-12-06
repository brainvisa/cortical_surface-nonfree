
#ifndef AIMS_PARAMETERIZEGYRI_GYRI_PARAMETERIZATION_H
#define AIMS_PARAMETERIZEGYRI_GYRI_PARAMETERIZATION_H

using namespace aims;
using namespace std;

TimeTexture<float> GyriParamTexture(AimsSurface<3,Void> &inMesh, Texture<short> &inTex, const Texture<float> &spmTex,
      const vector<uint> &corres, const pair<vector<vector<uint> >,vector<vector<uint> > > &diffMod, short option);

TimeTexture<float> GyriParamTexture(string meshfile, string intexfile, string gyrifile, string diffmodfile, string spmtexfile, uint constraint_method, float criter, float dt);

#endif

