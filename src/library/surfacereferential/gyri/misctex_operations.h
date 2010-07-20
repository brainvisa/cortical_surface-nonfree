
#ifndef AIMS_PARAMETERIZEGYRI_MISCTEX_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_MISCTEX_OPERATIONS_H



void makeGenericTexture(Texture<short> &inTex, const std::vector<uint> &corres);

void createAuxTexture(uint size,const std::vector<std::pair<std::pair<std::vector<uint>, std::vector<uint> >, std::pair<std::vector<uint>, std::vector<uint> > > > &points,
    char* path);




#endif

