#ifndef AIMS_GYRI_MODEL_H
#define AIMS_GYRI_MODEL_H

#include <aims/mesh/surfaceOperation.h>
#include <aims/io/io_g.h>

using namespace aims;
using namespace std;

class RingClique // cliques that include a node and all neighbors around
{
public:
	RingClique(uint node, std::set<uint> neigh) :
		_center(node), _ring(neigh) {};
	uint _center;
	std::set<uint> _ring;
	uint _max;
	uint _min;
};

class CurvClique // curv that includes a node and the two neighbours in the max and min curvature direction
{
public:
	CurvClique(uint node, uint min, uint min2, uint max, uint max2, float kmin, float kmax, float kmax2) :
		_center(node), _min(min), _min2(min2), _max(max), _max2(max2), _kmin(kmin), _kmax(kmax), _kmax2(kmax2) {};
	uint _center;
	uint _min;
	uint _min2;
	uint _max;
	uint _max2; // the direction opposite the one of max curvature
	float _kmin;
	float _kmax;
	float _kmax2;
};

class GyriRegularization
{
public:

     GyriRegularization(AimsSurfaceTriangle mesh, TimeTexture<int> gyriTexture, TimeTexture<float> curvMap, double weightData, int smooth) :
     _mesh(mesh), _gyriTexture(gyriTexture), _curvMap(curvMap), _weightData(weightData), _smooth(smooth) {_offset=0; _size=mesh.vertex().size(); computeNeighbours(); computeGyriProba(); /*computeCurvCliquesAndNodes();*/  compute2ndOrderCliquesAndNodes(); value2ndOrderCliques(); computeGraphEnergy(); initializeGyriEvolution();}

     void computeGraphEnergy(); // global energy
     double computeLocalEnergyChange(uint node, int label);  // energy change when changing node label
    // void updateLabelAndEnergy(uint node, int label); // actually changing node label and updating global energy
// function above is not used at the moment. Commented in the .cc as well

     double compute2ndOrderCliquePotential(uint cl, int l1, int l2);
     double computeCurvCliquePotential(uint cl, int l, int lmin, int lmin2, int lmax, int lmax2);
     double computeDataDrivenPotential(uint i, int l);

     void runICM();
     void runAnnealing(float T, float kT);
     void runICMdebug(uint node);
     void writeGyri(string fileOut);

     void debugCliques();

private:
     uint _size;
     AimsSurfaceTriangle _mesh;
     std::map<uint, std::set<int> > _labelMap; // the list of possible labels for each node;
     TimeTexture<int> _gyriTexture; // the current labels;
     TimeTexture<int> _gyriEvolution;
     TimeTexture<float> _curvMap;
     std::vector<std::set<uint> >  _neigh;
     std::vector<std::pair<uint, uint> > _2ndOrderCliques; // the vector of 2nd order cliques
     std::vector<RingClique> _ringCliques; // vector of 'ring' cliques (one node + direct neighbours)
     std::vector<CurvClique> _curvCliques; // vector of 'ring' cliques (one node + direct neighbours)
     std::vector<float> _cliqueValues; // if the cliques have to be valued, that is where we put the values
     std::map<uint, std::set<uint> > _nodes; //map node index to the set of cliques it belongs to.
											 // these are only the nodes that have more than one possible label;

     double _weightData;
     double _currentE;
     double _dataDrivenE;
     double _curvE;
     double _2ndE;
     int _offset; //just to know how many annealing steps there was before starting ICM
     int _smooth; // nb of iterations for smoothing the proba maps

     std::set<int> _seeds; // one node per gyri that keeps its original label;

     std::vector<double> _evolutionE;

     TimeTexture<float> _gyriProba;

     inline void computeNeighbours() {_neigh= SurfaceManip::surfaceNeighbours(_mesh);}
     void compute2ndOrderCliquesAndNodes();
     void computeCurvCliquesAndNodes();
     void computeRingCliquesAndNodes();
     void value2ndOrderCliques();
     void computeGyriProba();
     void initializeGyriEvolution();

     double computeLocalEnergyChangeDebug(uint node, int label);
};



#endif
