
#include <cortical_surface/surfacereferential/corticalConstraints.h>
#include <cortical_surface/surfacereferential/corticalTools.h>
#include <aims/distancemap/meshvoronoi.h>
#include <cortical_surface/surfacereferential/sulcusCorticalSnake_energy.h>
#include <float.h>

using namespace aims;
using namespace aims::meshdistance;
using namespace std;

namespace aims{
//////////////////////////////////////////////////////////////
//Total Energy
//////////////////////////////////////////////////////////////

	float SulcusCorticalSnake_energy::total_energy()
	{
// 		std::cout<<"SulcusCorticalSnake_energy : total_energy"<<std::endl;
		float res=0;
		res=alpha1*tex_distance[0].item(list_points[index_courant]) + alpha2*curv[0].item(list_points[index_courant]) + alpha3*elastic_energy() ;
// 		std::cout<<"SulcusCorticalSnake_energy : total_energy : return "<<res<<")"<<std::endl;
		return(res);
	}



//////////////////////////////////////////////////////////////
//Distance Energy
//////////////////////////////////////////////////////////////

/*	float SulcusCorticalSnake_energy::distance_energy()
{
// 	std::cout<<"distance energy"<<endl;
	float D_total=0;
	if(alpha1==0)
		return 0;
	if(size_vector!=0)
	{
		D_total=distance_Pi_to_bucket(list_points[index_courant]);
	}
	return ( D_total );
}

float SulcusCorticalSnake_energy::distance_Pi_to_bucket(uint index_vector)
{
// 	std::cout<<"distance pi to bucket"<<endl;
	float d_min=10000;
	float d_max=0;
	float d_curr=0;

	//for now, weight = 1
	//weight will depend of sulci length
	float weight=1;

	//CETTE BOUCLE DOIT ETRE REMPLACEE PAR LA DISTANCE AUX BUCKETS!!!
	for(uint i=0; i<voxel_bucket.size(); i++)
	{
		d_curr=weighted_square_distance(index_vector, voxel_bucket[i], weight);
		
		if(d_curr>d_max)
		{
			d_max=d_curr;
		}
		if(d_curr<d_min)
		{
			d_min=d_curr;
		}
	}
	return ( d_min );
}

float SulcusCorticalSnake_energy::weighted_square_distance(uint i, Point3d j, float weight)
{
	return ( square_distance(i,j)*weight );
}
*/
//Distance euclidienne au carre entre le sommet i et le bucket j
// float SulcusCorticalSnake_energy::square_distance(uint p, Point3d j)
// {
// 	return pow( ( ( (float)vert[p][0] - (float)j[0] ) + ( (float)vert[p][1] - (float)j[1] ) + ( (float)vert[p][2] - (float)j[2] ) ), 2);
// }


//////////////////////////////////////////////////////////////
//Curvature Energy
//////////////////////////////////////////////////////////////

float SulcusCorticalSnake_energy::curvature_energy()
{
// 	std::cout<<"SulcusCorticalSnake_energy : curvature_energy"<<std::endl;
// 	std::cout<<"curvature energy"<<endl;
	float C_total=0;
	if(alpha2==0)
		return 0;
	if(size_vector!=0)
	{
		//Test avec le calcul juste a l'index changeant
		C_total=curv[0].item(list_points[index_courant]);
	}
// 	std::cout<<"SulcusCorticalSnake_energy : curvature_energy : return "<<C_total<<std::endl;
	return ( C_total );
}


//////////////////////////////////////////////////////////////
//Elastic Energy
//////////////////////////////////////////////////////////////
float SulcusCorticalSnake_energy::elastic_energy()
{
// 	std::cout<<"SulcusCorticalSnake_energy : elastic_energy "<<std::endl;
	if(alpha3==0)
		return 0;
	
	float G_total=0;
// 	float tmp=0;
	std::vector<float> tab;
	std::vector<float> vect;
	std::vector<uint> index;
	
	//Si on utilise que l'index changeant
	if(index_courant==0)
	{
		index.push_back(n1);
	}
	

	if( index_courant>0 )
	{
		index.push_back(list_points[index_courant-1]);
	}
	if( index_courant<(size_vector-1) )
	{
		index.push_back(list_points[index_courant+1]);
	}
	
	//de pn à n2
	if( index_courant==(size_vector-1) )
	{
		index.push_back(n2);
	}
	
	vect=geodesic_distance(list_points[index_courant],index[0],index[1]);
	G_total=pow( vect[0], 2) + pow(vect[1], 2);
	G_total=G_total/pow( ( max/size_vector ), 2);
	
	tab.clear();
	vect.clear();
	index.clear();
	
	return ( G_total );
}

vector<float> SulcusCorticalSnake_energy::geodesic_distance(uint origin, uint i, uint j)
{
// 	std::cout<<"SulcusCorticalSnake_energy : geodesic_distance "<<std::endl;
	
	vector<float> v;
	
	TimeTexture<float> origin_distance(1,size);
	init_texture_single(origin_distance);
	
	origin_distance[0].item( origin )=10;
	
	v=MeshDistance_adapt( origin_distance[0],i,j);
	
	origin_distance.erase();
	
	return ( v );
}


////////////////////////////////////////////////////////
//Calcul distance geodesique (distance_map modif)
////////////////////////////////////////////////////////

vector<float> SulcusCorticalSnake_energy::MeshDistance_adapt( const Texture<float> & inittex, uint ind1, uint ind2 )
{
// 	std::cout<<"SulcusCorticalSnake_energy : MeshDistance_adapt "<<std::endl;
	
	Texture<float> tex;
	unsigned i, n = vert.size();

	ASSERT( inittex.nItem() == n );
	tex.reserve( n );

  // neighbours map

	std::map<unsigned, std::set<unsigned> >	neighbours;
	unsigned v1, v2, v3;

	for( i=0; i<poly.size(); ++i )
	{
		v1 = poly[i][0];
		v2 = poly[i][1];
		v3 = poly[i][2];
		if(inittex.item(v1)!=MESHDISTANCE_FORBIDDEN 
				 && inittex.item(v2)!=MESHDISTANCE_FORBIDDEN)
		{
			neighbours[v1].insert( v2 );
			neighbours[v2].insert( v1 );
		}
		if(inittex.item(v1)!=MESHDISTANCE_FORBIDDEN 
				 && inittex.item(v3)!=MESHDISTANCE_FORBIDDEN) 
		{
			neighbours[v1].insert( v3 );
			neighbours[v3].insert( v1 );
		}
		if(inittex.item(v2)!=MESHDISTANCE_FORBIDDEN 
				 && inittex.item(v3)!=MESHDISTANCE_FORBIDDEN)  
		{
			neighbours[v2].insert( v3 );
			neighbours[v3].insert( v2 );
		}
	}

  // init texture

	for( i=0; i<n; ++i )
	{
		if( inittex.item(i) == 0 )
			tex.push_back( FLT_MAX );
		else if( inittex.item(i) == MESHDISTANCE_FORBIDDEN )
			tex.push_back( MESHDISTANCE_FORBIDDEN );
		else
			tex.push_back( 0 );
	}

	std::multimap<float,unsigned>	front1, front2;
	std::multimap<float,unsigned>	*cfront = &front1, *nfront = &front2, *tmpf;
	std::multimap<float,unsigned>::iterator	iv, fv;
	std::set<unsigned>				neigh;
	std::set<unsigned>::iterator		in, fn;
	float					d, d2, l;
	Point3df				pos;
	vector<float> vect;
//	init first front

	for( i=0; i<n; ++i )
		if( tex.item(i) == 0 )
			front1.insert( std::pair<float,unsigned>( 0, i ) );

  //	loop until current front is empty
// 	int test=0;
	
	vect.push_back( 0 );
	vect.push_back( 0 );
	
// 	cout<<"ind1="<<ind1<<" - ind2="<<ind2<<endl;
	while( ( vect[0] <1 ) || ( vect[1] <1 ) )
	{
		nfront->clear();
		neigh.clear();

		for( iv=cfront->begin(), fv=cfront->end(); iv!=fv; ++iv )
		{
			i = (*iv).second;
			d = (*iv).first;
			for( in=neighbours[i].begin(), fn=neighbours[i].end(); in!=fn; ++in )
			{
				d2 = tex.item( *in );
				pos = vert[i] - vert[*in];
				l = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );
				if( d2 > d + l )
				{
					tex.item( *in ) = d + l;
					neigh.insert( *in );
// 					result_lim[0].item( *in )=tex.item( *in );
					if( (*in)==ind1 )
						vect[0]++;
					if( (*in)==ind2 )
						vect[1]++;
				}
			}
		}

		for( in=neigh.begin(), fn=neigh.end(); in!=fn; ++in )
			nfront->insert( std::pair<float,unsigned>( tex.item( *in ), *in ) );

		tmpf = cfront;
		cfront = nfront;
		nfront = tmpf;
	}
	
	vect[0]=tex.item(ind1);
	vect[1]=tex.item(ind2);
	
	neighbours.clear();
	front1.clear();
	front2.clear();
	(*cfront).clear();
	(*nfront).clear();
	(*tmpf).clear();
	neigh.clear();
	
// 	std::cout<<"SulcusCorticalSnake_energy : MeshDistance_adapt : return "<<vect[0]<<" et "<<vect[1]<<std::endl;
	return( vect );
}

}
















