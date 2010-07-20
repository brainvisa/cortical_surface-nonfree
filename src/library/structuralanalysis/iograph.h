#ifndef AIMS_IOGRAPH_H
#define AIMS_IOGRAPH_H



void ConcatenerGraphes( const std::vector<Graph*> in, Graph & out, const std::string & subjectatt );
std::vector<std::string> splitGraphFile(std::string graphFile);
void LireGraphes(std::string graphFile, Graph &primal);
void SauvegarderGraphes(Graph &primal, std::string graphFile, std::string output );
std::set<std::string> RecupererSujets(Graph &primal);
void RecupererGraphesIndividuels( const std::vector<Graph*> subjects, Graph & multi, const std::string & vertexatt );


#endif
