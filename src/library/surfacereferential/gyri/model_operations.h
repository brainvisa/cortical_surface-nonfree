
#ifndef AIMS_PARAMETERIZEGYRI_MODEL_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_MODEL_OPERATIONS_H



std::vector<std::string> readGyriToTextureFile(const char *path);

std::pair<std::vector<std::vector<uint> >, std::vector<std::vector<uint> > > getDiffusionModel(const char *path);

std::vector<uint> gyrToTexCorres(const std::vector<std::string> &tab);

std::vector<uint> getGyriToTextureCorres(const char *path);

std::pair<Point3d, Point3d> getGyrusModel(uint gyrus, const std::vector<std::vector<uint> > &diffMod);

std::vector<std::vector<short> > getGyrusConstraints(uint gyrus, const std::vector<std::vector<uint> > &diffMod);


#endif

