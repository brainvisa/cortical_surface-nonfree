
#ifndef AIMS_PARAMETERIZEGYRI_INTERSEC_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_INTERSEC_OPERATIONS_H

using namespace aims;
using namespace std;



vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, short stop, const vector<set<uint> > &voisins,
   const Texture<short> &inTex);

vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, const pair<short,short> &stop, const vector<set<uint> > &voisins,
   const Texture<short> &inTex);

vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, const vector<pair<short,short> > &stop, const vector<set<uint> > &voisins,
   const Texture<short> &inTex);


uint lookUpIntersectionCase(short gyruslabel, short left, short central, short right, const vector<set<uint> > &voisins,
   const Texture<short> &inTex);

vector<uint> substractIntersections(const vector<uint> &intersection, short gyruslabel, const vector<short> &parcours, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<uint> getIntersection(uint code, short gyruslabel, short left, short central, short right, const vector<set<uint> > &voisins,
                                    const Texture<short> &inTex);



#endif

