
#ifndef AIMS_PARAMETERIZEGYRI_MODEL_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_MODEL_OPERATIONS_H

using namespace aims;
using namespace std;

vector<string> readGyriToTextureFile(const char *path);

pair<vector<vector<uint> >, vector<vector<uint> > > getDiffusionModel(const char *path);

vector<uint> gyrToTexCorres(const vector<string> &tab);

vector<uint> getGyriToTextureCorres(const char *path);

pair<Point3d, Point3d> getGyrusModel(uint gyrus, const vector<vector<uint> > &diffMod);

vector<vector<short> > getGyrusConstraints(uint gyrus, const vector<vector<uint> > &diffMod);


#endif

