
#ifndef AIMS_PARAMETERIZEGYRI_GYRI_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_GYRI_OPERATIONS_H

using namespace std;


std::vector<short> getNeighbours(short gyrus1, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<short> getNeighbouringLabels(short gyrus1, short gyrus2, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<short> getNeighbouringLabels(uint n, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<short> getNeighbouringLabels(const std::vector<uint> &v, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<short> getCommonNeighbours(short gyrus1, short gyrus2, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<std::pair<short,short> > getInBetweenLabels(short gyrus1, short gyrus2, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);



#endif

