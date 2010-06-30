
#ifndef AIMS_PARAMETERIZEGYRI_COMPCONN_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_COMPCONN_OPERATIONS_H

using namespace aims;
using namespace std;


vector<vector<uint> > getComposantesConnexes(const set<uint> &v, const vector<set<uint> > &voisins);

vector<vector<uint> > getComposantesConnexes2(const vector<uint> &v, const vector<set<uint> > &voisins);

vector<vector<uint> > getComposantesConnexes(short gyruslabel, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<uint> fusionComposantesConnexes(uint newpoint, const vector<vector<uint> > &compConn, const vector<set<uint> > &voisins);

vector<vector<uint> > fusionComposantesConnexes(uint comp1, uint comp2, const vector<vector<uint> > &compConn);

void raccomodage(vector<vector<uint> > &compConn, const vector<uint> &candidates, const vector<set<uint> > &voisins);

pair<vector<uint>, vector<uint> > sortRightLeft(AimsSurface<3,Void> &inMesh, const pair<vector<uint>, vector<uint> > &hautBas,
      const pair<vector<uint>, vector<uint> > &gaucheDroite, const vector<set<uint> > &voisins);

pair<vector<uint>, vector<uint> > getOppositeSides(pair<vector<uint>, vector<uint> > &hautBas,
               const vector<vector<uint> > &vertices, const vector<set<uint> > &voisins, const Texture<short> &inTex);

vector<vector<uint> > nettoyerTaches(Texture<short> &inTex, const vector<set<uint> > &voisins);



#endif

