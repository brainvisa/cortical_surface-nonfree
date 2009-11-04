#ifndef AIMS_IOGRAPH_H
#define AIMS_IOGRAPH_H

using namespace aims;
using namespace carto;
using namespace std;

void ConcatenerGraphes( const vector<Graph*> in, Graph & out, const string & subjectatt );
vector<string> splitGraphFile(string graphFile);
void LireGraphes(string graphFile, Graph &primal);
void SauvegarderGraphes(Graph &primal, string graphFile, string output );
set<string> RecupererSujets(Graph &primal);
void RecupererGraphesIndividuels( const std::vector<Graph*> subjects, Graph & multi, const std::string & vertexatt );


#endif
