#include <aims/mesh/texture.h>
#include <aims/mesh/surface.h>
#include "gyriModel.h"

using namespace std;
using namespace aims;

void GyriRegularization::computeCliquesAndNodes()
{
	std::cout << "Building cliques and node list" << std::endl;
	_cliques=std::vector<std::pair<uint, uint> >();
	uint cliqueNb=0;
	for (uint i=0; i<_size; i++)
	{
		if (_labelMap[i].size() >= 2) // only nodes that have to or can be changed label
		{
			if (_nodes.find(i) == _nodes.end())
				_nodes[i]=std::set<uint>(); // if node does not exist yet it is created
			std::set<uint> neigh=_neigh[i];
			std::set<uint>::iterator neighIt=neigh.begin();
			for (; neighIt!=neigh.end(); neighIt++)
			{
				if ((*neighIt) > i) // clique has not been sen previously
				{
					_cliques.push_back(std::pair<uint, uint>(i, *neighIt) ); //clique created
					_nodes[i].insert(cliqueNb); // clique added to node
					if (_labelMap[*neighIt].size() >=2) // clique added to other node if necessary
					{
						_nodes[*neighIt]=std::set<uint>();
						_nodes[*neighIt].insert(cliqueNb);
					}
					cliqueNb++;
				}
			} // next neighbour
		}
	} // next node
//	// test
//	for (uint i=0; i<_cliques.size(); i++)
//	{
//		std::cout << "Clique " << i << " : (" << _cliques[i].first << ", " << _cliques[i].second << ")" << std::endl;
//	}
}

double GyriRegularization::computeCliquePotential(uint cl, int l1, int l2)
{
	if (l1==l2)
		return(0.0);
	else
	{
		float k1=(double) _curvMap[0].item(_cliques[cl].first);
		float k2=(double) _curvMap[0].item(_cliques[cl].second);
		return((k1+k2)/2.0);
	}
}

void GyriRegularization::computeGraphEnergy()
{

	uint nbCliques=_cliques.size();
	_currentE=0.0;
	for (uint i=0; i<nbCliques; i++)
		_currentE+=computeCliquePotential(i,_gyriTexture[0].item(_cliques[i].first),_gyriTexture[0].item(_cliques[i].second) );
}

double GyriRegularization::computeLocalEnergyChange(uint node, int label)
{
	int oldLabel=_gyriTexture[0].item(node);
	double oldLocalEnergy=0.0, localEnergy=0.0;
	std::set<uint> cliqueSet=_nodes[node];
	std::set<uint>::iterator cliqueIt=cliqueSet.begin();
	for (; cliqueIt!=cliqueSet.end(); cliqueIt++)
	{
		uint other;
		if (_cliques[*cliqueIt].first == node) other=_cliques[*cliqueIt].second;
		else other=_cliques[*cliqueIt].second;
		int label2=_gyriTexture[0].item(other);
		oldLocalEnergy+=computeCliquePotential(*cliqueIt, oldLabel, label2);
		localEnergy+=computeCliquePotential(*cliqueIt, label, label2);
	}
	return(localEnergy - oldLocalEnergy);
}

void GyriRegularization::updateLabelAndEnergy(uint node, int label)
{
	int oldLabel=_gyriTexture[0].item(node);
	double oldLocalEnergy=0.0, localEnergy=0.0;
	std::set<uint> cliqueSet=_nodes[node];
	std::set<uint>::iterator cliqueIt=cliqueSet.begin();
	for (; cliqueIt!=cliqueSet.end(); cliqueIt++)
	{
		uint other;
		if (_cliques[*cliqueIt].first == node) other=_cliques[*cliqueIt].second;
		else other=_cliques[*cliqueIt].second;
		int label2=_gyriTexture[0].item(other);
		oldLocalEnergy += computeCliquePotential(*cliqueIt, oldLabel, label2);
		localEnergy += computeCliquePotential(*cliqueIt, label, label2);
	}

	_currentE += (localEnergy - oldLocalEnergy);
	_gyriTexture[0].item(node)=label;
}


void GyriRegularization::runICM()
{
	TimeTexture<int> gyriEvolution(20, _size);
	uint i;
	int iter=1;
	for (uint j=0; j<20; j++)
		for (i=0; i<_size; i++)
			gyriEvolution[j].item(i)=_gyriTexture[j].item(i);
	std::map<uint, std::set<uint> >::iterator nodeIt;

	for ( ; iter<10000; iter++)      // iterations
	{
		std::cout << "Iteration " << iter << ": " << std::flush;
		uint nbCh=0;
		for (nodeIt=_nodes.begin(); nodeIt!=_nodes.end(); nodeIt++) // parcours des noeuds
		{
			i=(*nodeIt).first;
			double deMin=0, dE;
			int bestL=_gyriTexture[0].item(i), currentL=bestL;
			std::set<int> labelSet=_labelMap[i];
			std::set<int>::iterator labelIt=labelSet.begin();
			for (; labelIt!=labelSet.end(); labelIt++)
			{
				if ((*labelIt)!=currentL) dE=computeLocalEnergyChange(i, (*labelIt));
				if (dE<deMin)
				{
					deMin=dE;
					bestL=(*labelIt);
				}
			}
			if (bestL != currentL) nbCh++;
			_gyriTexture[0].item(i)=bestL; // this two line do what updateLabelAndEnergy()
			_currentE += deMin;			   // does, but so far it's more efficient not to call it.
			gyriEvolution[iter].item(i)=_gyriTexture[0].item(i);
		}
		std::cout << "nb changes : " << nbCh << "; Energy =" << _currentE << std::endl;
		if (nbCh==0)
			iter=20000;
	}
	std::cout << "Finished" << std::endl;
	std::cout << "Writing texture gyriEvolution..." << std::endl;
    Writer<TimeTexture<int> > texResultW( "gyriEvolution" );
    texResultW << gyriEvolution ;
}























