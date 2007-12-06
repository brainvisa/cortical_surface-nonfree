
#ifndef AIMS_PARAMETERIZEGYRI_BASICTEX_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_BASICTEX_OPERATIONS_H

using namespace aims;
using namespace std;

TimeTexture<short> getTimeTex(const Texture<short> &inTex);
TimeTexture<float> getTimeTex(const Texture<double> &inTex);
TimeTexture<float> getTimeTex(const Texture<float> &inTex);
void writeTexture(const Texture<short> &inTex, char *path);
void writeTexture(const TimeTexture<short> &inTex, char *path);
void writeTexture(const Texture<double> &inTex, char *path);
void writeTexture(const Texture<float> &inTex, char *path);
void addToTexture(const Texture<double> &inTex, char *path);


#endif

