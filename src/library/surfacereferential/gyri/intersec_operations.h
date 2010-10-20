
#ifndef AIMS_PARAMETERIZEGYRI_INTERSEC_OPERATIONS_H
#define AIMS_PARAMETERIZEGYRI_INTERSEC_OPERATIONS_H




std::vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, short stop, const std::vector<std::set<uint> > &voisins,
   const Texture<short> &inTex);

std::vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, const std::pair<short,short> &stop, const std::vector<std::set<uint> > &voisins,
   const Texture<short> &inTex);

std::vector<short> parcoursPerim(short gyruslabel, short start, short forbidden, const std::vector<std::pair<short,short> > &stop, const std::vector<std::set<uint> > &voisins,
   const Texture<short> &inTex);


uint lookUpIntersectionCase(short gyruslabel, short left, short central, short right, const std::vector<std::set<uint> > &voisins,
   const Texture<short> &inTex);

std::vector<uint> substractIntersections(const std::vector<uint> &intersection, short gyruslabel, const std::vector<short> &parcours, const std::vector<std::set<uint> > &voisins, const Texture<short> &inTex);

std::vector<uint> getIntersection(uint code, short gyruslabel, short left, short central, short right, const std::vector<std::set<uint> > &voisins,
                                    const Texture<short> &inTex);



#endif

