
#ifndef AIMS_PARAMETERIZEGYRI_GYRI_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_GYRI_OPERATIONS_H

using namespace aims;
using namespace std;


vector<short> getNeighbours(short gyrus1, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<short> getNeighbouringLabels(short gyrus1, short gyrus2, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<short> getNeighbouringLabels(uint n, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<short> getNeighbouringLabels(const vector<uint> &v, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<short> getCommonNeighbours(short gyrus1, short gyrus2, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<pair<short,short> > getInBetweenLabels(short gyrus1, short gyrus2, const vector<set<uint> > &voisins, const Texture<short> &inTex);



#endif

