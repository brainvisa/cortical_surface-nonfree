#ifndef AIMS_FUNCTION_PROJECTION_H
#define AIMS_FUNCTION_PROJECTION_H

using namespace aims;
using namespace carto;
using namespace std;
using namespace aims::meshdistance;

vector<AimsData<float> > load_kernel(string path);

Texture<float> deconvolve(AimsData<float> inFuncData, vector<AimsData<float> > &kernel, AimsSurfaceTriangle mesh);

#endif

