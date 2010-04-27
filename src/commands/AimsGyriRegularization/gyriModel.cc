#include <aims/mesh/texture.h>
#include <aims/mesh/surface.h>
#include <aims/math/random.h>
#include <aims/primalsketch/finiteElementSmoother_d.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshmorphomat.h>
#include <algorithm>
#include "gyriModel.h"
#include <algorithm>

using namespace std;
using namespace aims;
using namespace aims::meshdistance;


void GyriRegularization::compute2ndOrderCliquesAndNodes()
{
	std::cout << "Building cliques and node list" << std::endl;
	_2ndOrderCliques=std::vector<std::pair<uint, uint> >();
	uint cliqueNb=0;
	for (uint i=0; i<_size; i++)
	{
		if (_labelMap[i].size() >= 1) // ATTENTION C'est 2 pour :// only nodes that have to or can be changed label
		{
			if (_nodes.find(i) == _nodes.end())
				_nodes[i]=std::set<uint>(); // if node does not exist yet it is created
			std::set<uint> neigh=_neigh[i];
			std::set<uint>::iterator neighIt=neigh.begin();
			for (; neighIt!=neigh.end(); neighIt++)
			{
				if ((*neighIt) > i) // clique has not been sen previously
				{
					_2ndOrderCliques.push_back(std::pair<uint, uint>(i, *neighIt) ); //clique created
					_nodes[i].insert(cliqueNb); // clique added to node
					if (_labelMap[*neighIt].size() >=1) // ATTENTION 2 aussi // clique added to other node if necessary
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
	std::cout << "OK" << std::endl;
}

void GyriRegularization::computeRingCliquesAndNodes()
{
	std::cout << "Building cliques and node list" << std::endl;
//	_ringCliques=std::vector<std::set<uint> >();
//
//	uint cliqueNb=0;
//	for (uint i=0; i<_size; i++)
//	{
//		if (_labelMap[i].size() >= 1) // ATTENTION C'est 2 pour :// only nodes that have to or can be changed label
//		{
//			if (_nodes.find(i) == _nodes.end())
//				_nodes[i]=std::set<uint>(); // if node does not exist yet it is created
//			std::set<uint> neigh=_neigh[i];
//			std::set<uint>::iterator neighIt=neigh.begin();
//			for (; neighIt!=neigh.end(); neighIt++)
//			{
//				if ((*neighIt) > i) // clique has not been sen previously
//				{
//					_cliques.push_back(std::pair<uint, uint>(i, *neighIt) ); //clique created
//					_nodes[i].insert(cliqueNb); // clique added to node
//					if (_labelMap[*neighIt].size() >=1) // ATTENTION 2 aussi // clique added to other node if necessary
//					{
//						_nodes[*neighIt]=std::set<uint>();
//						_nodes[*neighIt].insert(cliqueNb);
//					}
//					cliqueNb++;
//				}
//			} // next neighbour
//		}
//	} // next node
////	// test
////	for (uint i=0; i<_cliques.size(); i++)
////	{
////		std::cout << "Clique " << i << " : (" << _cliques[i].first << ", " << _cliques[i].second << ")" << std::endl;
////	}
	std::cout << "OK" << std::endl;
}

void GyriRegularization::computeCurvCliquesAndNodes()
{
	std::cout << "Building cliques, node list, and everything necessary" << std::endl;
	_curvCliques=std::vector<CurvClique>();

	// construire le labelMap pour connaitre les relations de voisinage entre gyri
	// et les labels possibles pour chaque noeud.

	std::map<int, std::set<int> > labelVois;

	float likeliT=0.00001; // threshold on likelihood (P=-log(likeliT))
							// below this, nodes are not included in the optimization loop

	for (uint i=0; i<_size; i++)
	{
		std::set<uint> vois=_neigh[i];
		std::set<uint>::iterator voisIt;
		uint l=_gyriTexture[0].item(i);
		if (labelVois.find(l)==labelVois.end())
			labelVois[l]=std::set<int>();
		for (voisIt=vois.begin(); voisIt!=vois.end(); voisIt++)
		{
			if (_gyriTexture[0].item(*voisIt) != l)
				labelVois[l].insert(_gyriTexture[0].item(*voisIt));
		}
	}
	std::map<int, std::set<int> >::iterator lvIt;
	//looking at neighboring relationships
//	for (lvIt=labelVois.begin(); lvIt!=labelVois.end(); lvIt++)
//	{
//		std::cout << (*lvIt).first << " : " << std::flush;
//		std::set<int> voisin=(*lvIt).second;
//		set<int>::iterator vIt=voisin.begin();
//		for (; vIt != voisin.end(); vIt++)
//			std:cout << *vIt << " / " << std::flush;
//		std::cout << std::endl;
//	}

	//std::cout << labelVois << endl;

	for (uint i=0; i<_size; i++)
	{
		_labelMap[i]=std::set<int>();
		_labelMap[i].insert(_gyriTexture[0].item(i));
		std::set<int> setVois=labelVois[_gyriTexture[0].item(i)];
		std::set<int>::iterator setVoisIt=setVois.begin();
		for (; setVoisIt != setVois.end(); setVoisIt++)
			_labelMap[i].insert(*setVoisIt);
	}

	uint cliqueNb=0;
	uint count=0;
	for (uint i=0; i<_size; i++)
	{
		if ((_gyriTexture[0].item(i) > 1) && (_gyriProba[_gyriTexture[0].item(i)].item(i)>likeliT))
			// cingular pole excluded as well as points with too high a likelyhood.
		{
			if (_nodes.find(i) == _nodes.end())
				_nodes[i]=std::set<uint>(); // if node does not exist yet it is created
			std::set<uint> neigh=_neigh[i];
			std::set<uint>::iterator neighIt=neigh.begin();
			float k, kmin=10000.0, kmax=-10000.0;
			uint imax, imin, imax2, imin2;
			for (; neighIt!=neigh.end(); neighIt++)
			{
				float nx, ny, nz, vx, vy, vz;
													// Computing curvatures
				nx=((_mesh.normal())[i])[0];
				ny=((_mesh.normal())[i])[1];
				nz=((_mesh.normal())[i])[2];
				vx=((_mesh.vertex())[*neighIt])[0] - ((_mesh.vertex())[i])[0];
				vy=((_mesh.vertex())[*neighIt])[1] - ((_mesh.vertex())[i])[1];
				vz=((_mesh.vertex())[*neighIt])[2] - ((_mesh.vertex())[i])[2];
				float k=(float) 2*((nx*vx)+(ny*vy)+(nz*vz))/float(vx*vx + vy*vy + vz*vz);
				if (k<kmin) {kmin=k; imin=(*neighIt);}
				if (k>kmax) {kmax=k; imax=(*neighIt);}
			}

			// now trying to look for the neighbor opposite the direction of max curvature.
			float pmax=0.0; float pmax2=0.0;
			for (neighIt=neigh.begin(); neighIt!=neigh.end(); neighIt++)
			{
				float mx, my, mz, vx, vy, vz, sum, mx2, my2, mz2;
				int p=*neighIt;
				if ((p!=imax) && (p!=imin))
				{
					vx=((_mesh.vertex())[p])[0] - ((_mesh.vertex())[i])[0];
					vy=((_mesh.vertex())[p])[1] - ((_mesh.vertex())[i])[1];
					vz=((_mesh.vertex())[p])[2] - ((_mesh.vertex())[i])[2];
					sum = sqrt(vx*vx + vy*vy + vz*vz);
					vx/=sum; vy/=sum; vz/=sum;

					mx=((_mesh.vertex())[imax])[0] - ((_mesh.vertex())[i])[0];
					my=((_mesh.vertex())[imax])[1] - ((_mesh.vertex())[i])[1];
					mz=((_mesh.vertex())[imax])[2] - ((_mesh.vertex())[i])[2];
					sum = sqrt(mx*mx + my*my + mz*mz);
					mx/=sum; my/=sum; mz/=sum;

					mx2=((_mesh.vertex())[imin])[0] - ((_mesh.vertex())[i])[0];
					my2=((_mesh.vertex())[imin])[1] - ((_mesh.vertex())[i])[1];
					mz2=((_mesh.vertex())[imin])[2] - ((_mesh.vertex())[i])[2];
					sum = sqrt(mx2*mx2 + my2*my2 + mz2*mz2);
					mx2/=sum; my2/=sum; mz2/=sum;

					float ps=vx*mx + vy*my + vz*mz;
					float ps2=vx*mx2 + vy*my2 + vz*mz2;
					if (ps<pmax)
					{
						pmax=ps;
						imax2=p;
					}
					if (ps2<pmax2)
					{
						pmax2=ps2;
						imin2=p;
					}
				}
			}

			// computing curvature in the imax2 direction
			float nx, ny, nz, vx, vy, vz;
			nx=((_mesh.normal())[i])[0];
			ny=((_mesh.normal())[i])[1];
			nz=((_mesh.normal())[i])[2];
			vx=((_mesh.vertex())[imax2])[0] - ((_mesh.vertex())[i])[0];
			vy=((_mesh.vertex())[imax2])[1] - ((_mesh.vertex())[i])[1];
			vz=((_mesh.vertex())[imax2])[2] - ((_mesh.vertex())[i])[2];
			float kmax2=(float) 2*((nx*vx)+(ny*vy)+(nz*vz))/float(vx*vx + vy*vy + vz*vz);

			_curvCliques.push_back(CurvClique(i, imin, imin2, imax, imax2, kmin, kmax, kmax2));
			_nodes[i].insert(cliqueNb);
			if ((_gyriTexture[0].item(imin) > 1) && (_gyriProba[_gyriTexture[0].item(imin)].item(imin)>likeliT))
			{
				if (_nodes.find(imin)==_nodes.end())
					_nodes[imin]=std::set<uint>();
				_nodes[imin].insert(cliqueNb);
			}
			if ((_gyriTexture[0].item(imax) > 1) && (_gyriProba[_gyriTexture[0].item(imax)].item(imax)>likeliT))
			{
				if (_nodes.find(imax)==_nodes.end())
					_nodes[imax]=std::set<uint>();
				_nodes[imax].insert(cliqueNb);
			}
			if ((_gyriTexture[0].item(imax2) > 1) && (_gyriProba[_gyriTexture[0].item(imax2)].item(imax2)>likeliT))
			{
				if (_nodes.find(imax2)==_nodes.end())
					_nodes[imax2]=std::set<uint>();
				_nodes[imax2].insert(cliqueNb);
			}
			cliqueNb++;

		}
		else
			count++;
	}

	TimeTexture<int> debugCurv(1,_size);
	for (uint i=0; i< _size; i++)
		debugCurv[0].item(i)=0;

//	for (uint i=0; i<_curvCliques.size(); i=i+1000)
//	{
//		debugCurv[0].item(_curvCliques[i]._center)=100;
//		debugCurv[0].item(_curvCliques[i]._min)=1;
//		debugCurv[0].item(_curvCliques[i]._max2)=200;
//		debugCurv[0].item(_curvCliques[i]._max)=300;
//	}
	std::map<uint, std::set<uint> >::iterator nodeIt;
	for (nodeIt=_nodes.begin(); nodeIt!=_nodes.end(); nodeIt++)
	{
		debugCurv[0].item((*nodeIt).first)=1;
	}
	Writer<TimeTexture<int> > texDebugW( "debugNodes" );
    texDebugW << debugCurv ;

	std::cout << "Built " << _curvCliques.size() << " cliques" << std::endl;
	std::cout << "and " << _nodes.size() << " nodes" << std::endl;
	std::cout << count << " nodes were discarded over a total of " << _size << std::endl;
}

void GyriRegularization::value2ndOrderCliques()
{
	int sC=_2ndOrderCliques.size();
	std::cout << "Computing clique values" << std::endl;
	_mesh.updateNormals();
	float sig=0.3;
	float min=100.0, max=-100.0;
	for (int i=0; i<sC; i++)
	{
		uint i1=_2ndOrderCliques[i].first;
		uint i2=_2ndOrderCliques[i].second;

		//here we compute the (normal) curvature at the node in the direction of the neighbour
		// and at the neighbour in the direction of the node.
		// The average is then the clique value.

		float val;
		float nk=0.0;
		float nx1, ny1, nz1, vx1, vy1, vz1;
		float nx2, ny2, nz2, vx2, vy2, vz2;

		nx1=((_mesh.normal())[i1])[0];nx2=((_mesh.normal())[i2])[0];
		ny1=((_mesh.normal())[i1])[1];ny2=((_mesh.normal())[i2])[1];
		nz1=((_mesh.normal())[i1])[2];nz2=((_mesh.normal())[i2])[2];
		vx1=((_mesh.vertex())[i2])[0] - ((_mesh.vertex())[i1])[0];
		vy1=((_mesh.vertex())[i2])[1] - ((_mesh.vertex())[i1])[1];
		vz1=((_mesh.vertex())[i2])[2] - ((_mesh.vertex())[i1])[2];
		vx2=((_mesh.vertex())[i1])[0] - ((_mesh.vertex())[i2])[0];
		vy2=((_mesh.vertex())[i1])[1] - ((_mesh.vertex())[i2])[1];
		vz2=((_mesh.vertex())[i1])[2] - ((_mesh.vertex())[i2])[2];
		float k1=(float) 2*((nx1*vx1)+(ny1*vy1)+(nz1*vz1))/float(vx1*vx1 + vy1*vy1 + vz1*vz1);
		float k2=(float) 2*((nx2*vx2)+(ny2*vy2)+(nz2*vz2))/float(vx2*vx2 + vy2*vy2 + vz2*vz2);
		float k=(k1+k2)/2.0;
		if (k<min) min=k;
		if (k>max) max=k;
		if (k>0)
			val=exp(-(k*k)/(2*sig*sig));
		else val=1.0;
		_cliqueValues.push_back(val);
	}
	std::cout << "OK. Min=" << min << " and Max=" << max << std::endl;
}

double GyriRegularization::compute2ndOrderCliquePotential(uint cl, int l1, int l2)
{
	if (l1==l2)
		return(0.0);
	else
	{
		return(_cliqueValues[cl]);
//		double k1=(double) _curvMap[0].item(_cliques[cl].first);
//		double k2=(double) _curvMap[0].item(_cliques[cl].second);
//		if (k1>=k2) return(k1);
//		else  return(k2);
			//return((k1+k2)/2.0);
	}
}

double GyriRegularization::computeCurvCliquePotential(uint cl, int l, int lmin, int lmax, int lmax2)
{

	double epsilon=0.1;
	double T0=0.3;
	double a, c;

	a=(1.0 - epsilon)/(T0*T0);
	c=epsilon;

	// max curvature
	float kmax=_curvCliques[cl]._kmax;
	float kmax2=_curvCliques[cl]._kmax2;
	if (lmin!=l)
		return(1.0);
	else
	{
		if ((lmax!=l) && (lmax2!=l))
			return(1.0);
		else if ((l==lmax2) && (l!=lmax))
		{
			if (kmax<=0)
				return(1.0);
			else if (kmax<=T0)
				return( (double) a*(kmax-T0)*(kmax-T0) + epsilon);
			else
				return(epsilon);
		}
		else if ((l==lmax) && (l!=lmax2))
		{
			if (kmax2<=0)
				return(1.0);
			else if (kmax2<=T0)
				return( (double) a*(kmax2-T0)*(kmax2-T0) + epsilon);
			else
				return(epsilon);
		}
		else return(0.0);
	}
//	if (lmin!=l)
//		return(1.0);
//	else
//	{
//		if (lmax==l)
//			return(0.0);
//		else
//		{
//			if (kmax<=0)
//				return(1.0);
//			else if (kmax<=0.4)
//				return( (double) 6.25*(kmax-0.4)*(kmax-0.4) );
//			else
//				return(0.0);
//		}
//	}
}

double GyriRegularization::computeDataDrivenPotential(uint i, int l)
{
//	if (_gyriProba[l].item(i) < 0.01)
//		return(double(_weightData*(-log(0.01))));
//	else
//		return(double(_weightData*(-log(_gyriProba[l].item(i)))));
	return(double(_weightData*(_gyriProba[l].item(i))));
}

void GyriRegularization::computeGraphEnergy()
{

	std::cout << "Computing global energy" << std::endl;
	uint nbCliques=_curvCliques.size();
	_currentE=0.0;
	_dataDrivenE=0.0;
	_curvE=0.0;
	for (uint i=0; i<_size; i++)
	{
		_dataDrivenE+=computeDataDrivenPotential(i, _gyriTexture[0].item(i));
	}
	for (uint i=0; i<nbCliques; i++)
		_curvE+=computeCurvCliquePotential(i,_gyriTexture[0].item(_curvCliques[i]._center),_gyriTexture[0].item(_curvCliques[i]._min) ,_gyriTexture[0].item(_curvCliques[i]._max), _gyriTexture[0].item(_curvCliques[i]._max2) );

	_currentE = _dataDrivenE + _curvE;

	std::cout << "Found : E=" << _currentE << ", with dataDrivenE=" << _dataDrivenE << ", and curvE=" << _curvE << std::endl;
}

double GyriRegularization::computeLocalEnergyChange(uint node, int label)
{
	int oldLabel=_gyriTexture[0].item(node);
	double oldLocalEnergy=0.0, localEnergy=0.0;
	std::set<uint> cliqueSet=_nodes[node];
	std::set<uint>::iterator cliqueIt=cliqueSet.begin();

	// first the data-driven term

	oldLocalEnergy+=computeDataDrivenPotential(node, oldLabel);
	localEnergy+=computeDataDrivenPotential(node, label);

	// then the list of cliques the node is involved in

	for (; cliqueIt!=cliqueSet.end(); cliqueIt++)
	{
		int l=_gyriTexture[0].item(_curvCliques[*cliqueIt]._center);
		int l1=_gyriTexture[0].item(_curvCliques[*cliqueIt]._min);
		int l2=_gyriTexture[0].item(_curvCliques[*cliqueIt]._max);
		int l3=_gyriTexture[0].item(_curvCliques[*cliqueIt]._max2);

		if (node==_curvCliques[*cliqueIt]._center)
		{
			oldLocalEnergy+=computeCurvCliquePotential(*cliqueIt, oldLabel, l1, l2, l3);
			localEnergy+=computeCurvCliquePotential(*cliqueIt, label, l1, l2, l3);
		}
		else if (node==_curvCliques[*cliqueIt]._min)
		{
			oldLocalEnergy+=computeCurvCliquePotential(*cliqueIt, l, oldLabel, l2, l3);
			localEnergy+=computeCurvCliquePotential(*cliqueIt, l, label, l2, l3);
		}
		else if (node==_curvCliques[*cliqueIt]._max)
		{
			oldLocalEnergy+=computeCurvCliquePotential(*cliqueIt, l, l1, oldLabel, l3);
			localEnergy+=computeCurvCliquePotential(*cliqueIt, l, l1, label, l3);
		}
		else if (node==_curvCliques[*cliqueIt]._max2)
		{
			oldLocalEnergy+=computeCurvCliquePotential(*cliqueIt, l, l1, l2, oldLabel);
			localEnergy+=computeCurvCliquePotential(*cliqueIt, l, l1, l2, label);
		}
		else
		{
			std::cerr << "ComputeLocalEnergyChange(node=" << node << ", label=" << label << ") : node not in clique" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	return(localEnergy - oldLocalEnergy);
}


void GyriRegularization::computeGyriProba()
{

	std::cout << "Computing probability map for each gyrus" << std::endl;
	uint i=0;
	int j;
	int l, lmax=0;

	// max label
	for (i=0; i<_size; i++)
		if (_gyriTexture[0].item(i) > lmax)
			lmax=_gyriTexture[0].item(i);


	_gyriProba=TimeTexture<float>(lmax+1, _size);

	// proba equal 0 to start with
	TimeTexture<float> tempProba(lmax+1, _size), tempProbaDil(lmax+1, _size);
	for (j=0; j<=lmax; j++)
		for (i=0; i<_size; i++)
			tempProba[j].item(i)=0.0;

	// then proba(label)=100.0 for the region that has this label
	for (i=0; i<_size; i++)
			tempProba[_gyriTexture[0].item(i)].item(i)=100.0;

	// dilation of the gyri to extend their possible localization a bit
	std::cout << "Dilating likelihood maps before smoothing" << std::endl;
	for (j=0; j<=lmax; j++)
	{
		tempProbaDil[j]=MeshDilation<float>( _mesh[0], tempProba[j], 0.0, -1, 3.0, true);
	}

	Writer<TimeTexture<float > > tprobW("tempProba");
	tprobW << tempProba;
	Writer<TimeTexture<float > > tprobdW("tempProbaDil");
	tprobdW << tempProbaDil;

	//Smoothing of the textures
	FiniteElementSmoother<3, float> *smooth;
	smooth=new FiniteElementSmoother<3, float>(float(0.1), &_mesh[0]);
	std::cout << "Smoothing probability maps" << std::endl;
	for (j=0; j<=lmax; j++)
	{
		_gyriProba[j]=smooth->doSmoothing(tempProbaDil[j], _smooth, false);
	}


	for (i=0; i<_size; i++)
	{
		float sum=0.0;
		for (j=0; j<=lmax; j++)
			sum += _gyriProba[j].item(i);

		for (j=0; j<=lmax; j++)
		{
			float p=_gyriProba[j].item(i)/float(sum);
			if (p<0.0000000001)
				_gyriProba[j].item(i)=-log(0.0000000001);
			else
				_gyriProba[j].item(i)=-log(p);
		}
	}

	Writer<TimeTexture<float > > probW("gyriProbaMaps");
	probW << _gyriProba;
}



//double GyriRegularization::computeLocalEnergyChange(uint node, int label)
//{
//	int oldLabel=_gyriTexture[0].item(node);
//	double oldLocalEnergy=0.0, localEnergy=0.0;
//	std::set<uint> cliqueSet=_nodes[node];
//	std::set<uint>::iterator cliqueIt=cliqueSet.begin();
//	for (; cliqueIt!=cliqueSet.end(); cliqueIt++)
//	{
//		uint other;
//		if (_2ndOrderCliques[*cliqueIt].first == node) other=_2ndOrderCliques[*cliqueIt].second;
//		else other=_2ndOrderCliques[*cliqueIt].first;
//		int label2=_gyriTexture[0].item(other);
//		oldLocalEnergy+=compute2ndOrderCliquePotential(*cliqueIt, oldLabel, label2);
//		localEnergy+=compute2ndOrderCliquePotential(*cliqueIt, label, label2);
//	}
//	return(localEnergy - oldLocalEnergy);
//}


//void GyriRegularization::updateLabelAndEnergy(uint node, int label)
//{
//	int oldLabel=_gyriTexture[0].item(node);
//	double oldLocalEnergy=0.0, localEnergy=0.0;
//	std::set<uint> cliqueSet=_nodes[node];
//	std::set<uint>::iterator cliqueIt=cliqueSet.begin();
//	for (; cliqueIt!=cliqueSet.end(); cliqueIt++)
//	{
//		uint other;
//		if (_2ndOrderCliques[*cliqueIt].first == node) other=_2ndOrderCliques[*cliqueIt].second;
//		else other=_2ndOrderCliques[*cliqueIt].second;
//		int label2=_gyriTexture[0].item(other);
//		oldLocalEnergy += compute2ndOrderCliquePotential(*cliqueIt, oldLabel, label2);
//		localEnergy += compute2ndOrderCliquePotential(*cliqueIt, label, label2);
//	}
//
//	_currentE += (localEnergy - oldLocalEnergy);
//	_gyriTexture[0].item(node)=label;
//}



void GyriRegularization::initializeGyriEvolution()
{
	int nb=2000;
	_gyriEvolution=TimeTexture<int>(nb, _size);
	for (int j=0; j<nb; j++)
		for (int i=0; i<_size; i++)
		{
			_gyriEvolution[j].item(i)=_gyriTexture[0].item(i);
		}
}

void GyriRegularization::runICM()
{
	int nb=100;
	TimeTexture<float> energyEvolution(nb, _size);
	uint i;
	int iter=1;
	for (int j=0; j<nb; j++)
		for (i=0; i<_size; i++)
		{
			energyEvolution[j].item(i)=0.0;
		}
	uint size_cliques=_curvCliques.size();

	std::cout << "Running ICM..." << std::endl;

	for (uint j=0; j<size_cliques; j++)
	{
		int l=_gyriTexture[0].item(_curvCliques[j]._center);
		int l1=_gyriTexture[0].item(_curvCliques[j]._min);
		int l2=_gyriTexture[0].item(_curvCliques[j]._max);
		int l3=_gyriTexture[0].item(_curvCliques[j]._max2);
		energyEvolution[0].item(_curvCliques[j]._center)=(float) computeCurvCliquePotential(j, l, l1, l2, l3);
	}

	std::map<uint, std::set<uint> >::iterator nodeIt;

	std::vector<uint> reorder;				// This is for a random order during optimisation
	for (nodeIt=_nodes.begin(); nodeIt!=_nodes.end(); nodeIt++)
	{
		reorder.push_back((*nodeIt).first);
	}

	std::random_shuffle(reorder.begin(), reorder.end());
	uint size_nodes=reorder.size();

	for ( ; iter<100; iter++)      // iterations
	{
		std::cout << "Iteration " << iter << ": " << std::flush;
		uint nbCh=0;
		for (uint j=0; j<size_nodes; j++)  // parcours des noeuds concernŽs
		{
			i=reorder[j];
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
			if (bestL != currentL)
			{
				nbCh++;
				_gyriTexture[0].item(i)=bestL; // this two line do what updateLabelAndEnergy()
				_currentE += deMin; // does, but so far it's more efficient not to call it.
			}
			_gyriEvolution[iter+_offset].item(i)=_gyriTexture[0].item(i);
		}
		for (uint j=0; j<size_cliques; j++)
		{
			int l=_gyriTexture[0].item(_curvCliques[j]._center);
			int l1=_gyriTexture[0].item(_curvCliques[j]._min);
			int l2=_gyriTexture[0].item(_curvCliques[j]._max);
			int l3=_gyriTexture[0].item(_curvCliques[j]._max);
			energyEvolution[iter].item(_curvCliques[j]._center)=(float) computeCurvCliquePotential(j, l, l1, l2, l3);
		}
		std::cout << "nb changes : " << nbCh << "; Energy =" << _currentE << std::endl;
		if (nbCh==0)
			iter=20000;
	}
	std::cout << "Finished" << std::endl;
	std::cout << "Writing texture gyriEvolution..." << std::endl;
    Writer<TimeTexture<int> > texResultW( "gyriEvolution" );
    texResultW << _gyriEvolution ;
    std::cout << "and texture curvEvolution." << std::endl;
	Writer<TimeTexture<float> > texEnergyW( "energyEvolution" );
	texEnergyW << energyEvolution ;
}


void GyriRegularization::runAnnealing(float T, float kT)
{
	int nb=200;
	//TimeTexture<float> energyEvolution(nb, _size);
	uint i;
	int iter=1;
//	for (int j=0; j<nb; j++)
//		for (i=0; i<_size; i++)
//		{
//			//gyriEvolution[j].item(i)=0;
//			energyEvolution[j].item(i)=0.0;
//		}
	uint size_cliques=_curvCliques.size();

//	for (uint j=0; j<size_cliques; j++)
//	{
//		int l=_gyriTexture[0].item(_curvCliques[j]._center);
//		int l1=_gyriTexture[0].item(_curvCliques[j]._min);
//		int l2=_gyriTexture[0].item(_curvCliques[j]._max);
//		//energyEvolution[0].item(_curvCliques[j]._center)=(float) computeCurvCliquePotential(j, l, l1, l2);
//	}

	std::map<uint, std::set<uint> >::iterator nodeIt;

	std::cout << "Starting annealing" << endl;

	std::vector<uint> reorder;				// This is for a random order during optimisation
	for (nodeIt=_nodes.begin(); nodeIt!=_nodes.end(); nodeIt++)
	{
		reorder.push_back((*nodeIt).first);
	}

	std::random_shuffle(reorder.begin(), reorder.end());
	uint size_nodes=reorder.size();

	uint nbZero=0;
	for ( ; iter<2000; iter++)      // iterations
	{
		std::cout << "Iteration " << iter << "-> T=" << T << " : " << std::flush;
		uint nbCh=0;
		for (uint j=0; j<size_nodes; j++)  // parcours des noeuds concernŽs
		{
			i=reorder[j];
			double deMin=0, dE;
			std::map<int, double> mapLabels;
			int bestL=_gyriTexture[0].item(i);
			int currentL=bestL;
			std::set<int> labelSet=_labelMap[i];
			std::set<int>::iterator labelIt=labelSet.begin();
			for (; labelIt!=labelSet.end(); labelIt++)
			{
				dE=computeLocalEnergyChange(i, (*labelIt));
				mapLabels[*labelIt]=dE;
			}
			std::map<int, double>::iterator mlIt=mapLabels.begin();
			double sum=0.0;
			double P;
			// Non-normalised probability computation
			for (; mlIt!=mapLabels.end(); mlIt++)
			{
//				P=exp(-(_currentE+(*mlIt).second)/double(T));
				P=exp(-((*mlIt).second)/double(T));
				(*mlIt).second=P;
				sum+=P;
			}
			double tirage = UniformRandom()*sum;
			sum=0.0;
			int label;
			for (mlIt=mapLabels.begin(); mlIt!=mapLabels.end(); mlIt++) // random draw according to probabilities
			{
				label=(*mlIt).first;
				sum+=(*mlIt).second;
				if (sum>tirage)
				{
					bestL=label;
					mlIt=mapLabels.end();
					break;
				}
			}

			//bestL=label;

			if (bestL != currentL)
			{
				nbCh++;
				_currentE += computeLocalEnergyChange(i, bestL); // does, but so far it's more efficient not to call it.
				_gyriTexture[0].item(i)=bestL; // this two line do what updateLabelAndEnergy()
				_gyriEvolution[iter].item(i)=bestL;
			}
			else
				_gyriEvolution[iter].item(i)=currentL;
		}
//		for (uint j=0; j<size_cliques; j++)
//		{
//			int l=_gyriTexture[0].item(_curvCliques[j]._center);
//			int l1=_gyriTexture[0].item(_curvCliques[j]._min);
//			int l2=_gyriTexture[0].item(_curvCliques[j]._max);
//			energyEvolution[iter].item(_curvCliques[j]._center)=(float) computeCurvCliquePotential(j, l, l1, l2);
//		}
		std::cout << "nb changes : " << nbCh << "; Energy =" << _currentE << "\r" << std::endl;
		_evolutionE.push_back(_currentE);
		if (nbCh==0)
			nbZero++;
		if (nbZero == 10)
		{
			_offset=iter;
			iter=20000;
		}
		T=T*kT;
	}
	if (iter < 20000) _offset=iter;
	std::cout << "Offset : " << _offset << std::endl;
	runICM();
	std::cout << "Finished" << std::endl;


	std::cout << "Computing data-driven local energy" << std::endl;
	Writer<TimeTexture<int> > resultW( "resultGyri" );
	resultW << _gyriTexture ;
	TimeTexture<float> dataD(1, _size);
	for (i=0; i< _size; i++)
	{
		dataD[0].item(i)=(float) computeDataDrivenPotential(i, _gyriTexture[0].item(i));
	}
	Writer<TimeTexture<float> > dataDW( "dataDriven" );
	dataDW << dataD ;

    std::cout << "Writing energy.txt" << endl;
    FILE *eFile=fopen("energy.txt", "w");
    for (i=0; i<_evolutionE.size(); i++)
    {
    	fprintf(eFile, "%f\n", _evolutionE[i]);
    }
    fclose(eFile);
//	Writer<TimeTexture<float> > texEnergyW( "energyEvolution" );
//	texEnergyW << energyEvolution ;
}

void GyriRegularization::writeGyri(string fileOut)
{
	Writer<TimeTexture<int> > texResultW( fileOut );
	texResultW << _gyriTexture ;
}





















