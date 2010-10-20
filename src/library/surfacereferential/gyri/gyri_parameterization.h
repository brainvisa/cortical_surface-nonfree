
#ifndef AIMS_PARAMETERIZEGYRI_GYRI_PARAMETERIZATION_H
#define AIMS_PARAMETERIZEGYRI_GYRI_PARAMETERIZATION_H



TimeTexture<float> GyriParamTexture(AimsSurface<3,Void> &inMesh, Texture<short> &inTex, const Texture<float> &spmTex,
      const std::vector<uint> &corres, const std::pair<std::vector<std::vector<uint> >, std::vector< std::vector<uint> > > &diffMod, short option);

TimeTexture<float> GyriParamTexture(std::string meshfile, std::string intexfile, std::string gyrifile, std::string diffmodfile, std::string spmtexfile, uint constraint_method, float criter, float dt);

#endif

