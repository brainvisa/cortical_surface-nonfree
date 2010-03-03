#ifndef AIMS_GYRI_MODEL_H
#define AIMS_GYRI_MODEL_H

#include <aims/mesh/surfaceOperation.h>
#include <aims/io/io_g.h>

using namespace aims;
using namespace std;

class GyriRegularization
{
public:

     GyriRegularization(AimsSurfaceTriangle mesh, std::map<uint, std::set<int> > labelMap, TimeTexture<int> gyriTexture, TimeTexture<float> curvMap) :
     _mesh(mesh), _labelMap(labelMap), _gyriTexture(gyriTexture), _curvMap(curvMap) {_size=mesh.vertex().size(); computeNeighbours(); computeCliquesAndNodes(); computeGraphEnergy();}

     void computeGraphEnergy(); // global energy
     double computeLocalEnergyChange(uint node, int label);  // energy change when changing node label
     void updateLabelAndEnergy(uint node, int label); // actually changing node label and updating global energy

     double computeCliquePotential(uint cl, int l1, int l2);

     void runICM();

private:
     uint _size;
     AimsSurfaceTriangle _mesh;
     std::map<uint, std::set<int> > _labelMap; // the list of possible labels for each node;
     TimeTexture<int> _gyriTexture; // the current labels;
     TimeTexture<float> _curvMap;
     std::vector<std::set<uint> >  _neigh;
     std::vector<std::pair<uint, uint> > _cliques; // the list of cliques
     std::map<uint, std::set<uint> > _nodes; //map node index to the set of cliques it belongs to.
											 // these are only the nodes that have more than one possible label;

     double _currentE;

     inline void computeNeighbours() {_neigh= SurfaceManip::surfaceNeighbours(_mesh);}
     void computeCliquesAndNodes();
};



#endif
