#ifndef AIMS_IOGRAPH_H
#define AIMS_IOGRAPH_H

using namespace aims;
using namespace carto;
using namespace std;

void LireGraphes(string graphFile, Graph &primal);
void SauvegarderGraphes(Graph &primal, string graphFile, string output );
set<string> RecupererSujets(Graph &primal);

#endif
