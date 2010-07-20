
#ifndef AIMS_PARAMETERIZEGYRI_COMPCONN_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_COMPCONN_OPERATIONS_H




std::vector<std::vector<uint> > getComposantesConnexes(const std::set<uint> &v, const std::vector<std::set<uint> > &voisins);

std::vector<std::vector<uint> > getComposantesConnexes2(const std::vector<uint> &v, const std::vector<std::set<uint> > &voisins);

std::vector<std::vector<uint> > getComposantesConnexes(short gyruslabel, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<uint> fusionComposantesConnexes(uint newpoint, const std::vector<std::vector<uint> > &compConn, const std::vector<std::set<uint> > &voisins);

std::vector<std::vector<uint> > fusionComposantesConnexes(uint comp1, uint comp2, const std::vector<std::vector<uint> > &compConn);

void raccomodage(std::vector<std::vector<uint> > &compConn, const std::vector<uint> &candidates, const std::vector<std::set<uint> > &voisins);

std::pair<std::vector<uint>, std::vector<uint> > sortRightLeft(AimsSurface<3,Void> &inMesh, const std::pair<std::vector<uint>, std::vector<uint> > &hautBas,
      const std::pair<std::vector<uint>, std::vector<uint> > &gaucheDroite, const std::vector<std::set<uint> > &voisins);

std::pair<std::vector<uint>, std::vector<uint> > getOppositeSides(std::pair< std::vector<uint>, std::vector<uint> > &hautBas,
               const std::vector<std::vector<uint> > &vertices, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<std::vector<uint> > nettoyerTaches(Texture<short> &inTex, const std::vector<std::set<uint> > &voisins);



#endif

