
#ifndef AIMS_PARAMETERIZEGYRI_MISCTEX_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_MISCTEX_OPERATIONS_H

using namespace aims;
using namespace std;


void makeGenericTexture(Texture<short> &inTex, const vector<uint> &corres);

void createAuxTexture(uint size,const vector<pair<pair<vector<uint>, vector<uint> >, pair<vector<uint>, vector<uint> > > > &points,
    char* path);




#endif

