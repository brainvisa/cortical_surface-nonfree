
#include <cstdlib>
#include <aims/mesh/geometric.h>
#include <aims/mesh/texture.h>
#include <aims/bucket/bucket.h>
#include <cortical_surface/surfacereferential/corticalTools.h>
#include <aims/distancemap/meshdistance.h>
#include <graph/graph/graph.h>
#include <float.h>

#include <cortical_surface/surfacereferential/sulcusCorticalSnake.h>
#include <cortical_surface/surfacereferential/sulcusCorticalSnake_energy.h>

using namespace aims;
using namespace aims::meshdistance;
using namespace std;
using namespace carto;



void SulcusCorticalSnake::createDistanceTexture()
{
	
	std::cout<<"SulcusCorticalSnake : createDistanceTexture"<<std::endl;
	for(uint j=0; j<size; j++)
		tex_distance[0].push_back(0);

	TimeTexture<float> init(1,size);	
	
	for(uint j=0; j<size;j++)
	{
		init[0].item(j)=0;
// 		std::cout<<" "<<constraint[0].item(j);
		if(constraint[0].item(j)==value)
			init[0].item(j)=10;
	}

	tex_distance[0]= MeshDistance( mesh[0], init[0],true);

	for(uint j=0; j<size;j++)
	{
		tex_distance[0].item(j)=tex_distance[0].item(j)*tex_distance[0].item(j);
	}
}

//Compute the total snake
TimeTexture<float> SulcusCorticalSnake::compute_snake()
{
	std::cout<<"SulcusCorticalSnake : ComputeSnake"<<std::endl;
	TimeTexture<float> cartea(1,size);
	for(uint i=0; i<size; i++)
		if(i==n1)
			cartea[0].item(i)=10;
	else
		cartea[0].item(i)=0;
	
	compute_curv();
	createDistanceTexture();
	
	for(uint j=0; j<size; j++)
		result_total[0].item(j)=0;
	
	int test=1;
	test=define_extremities();
	if(test==0)
		return( result_total );
	
	uint n=0;
	
	if( are_they_neighbours(n1,n2) == 0 )
	{
		n=define_new_middle_point(n1,n2);
		list_points.push_back(n);
		cpt_points.push_back(0);
	}
	
	//defini la taille entre les 2 estremites (le max en gros)
	
	
	TimeTexture<float> carte(1,size);
	for(uint i=0; i<size; i++)
		if(i==n1)
			carte[0].item(i)=10;
	else
		carte[0].item(i)=0;
	
	
	max=MeshDistance_adapt( carte[0],true,n2);
	
	TimeTexture<float> result(100,size);
	for(int i=0; i<100; i++)
		for(uint j=0; j<size; j++)
			result[i].item(j)=0;
	
	while(stop_total()==0)
	{
// 		std::cout<<"SulcusCorticalSnake : ComputeSnake : while stoptotal==0, resolution="<<cpt_resolution<<std::endl;
// 		for(uint j=0; j<size; j++)
// 		{
// 			result[cpt_multi].item(j)=0;
// 			for(uint i=0; i<list_points.size(); i++)
// 			{
// 				if(list_points[i]==j)
// 				{
// 					result[cpt_multi].item(j)=i+2;
// 				}
// 				else	
// 				{
// 				}
// 			}
// 		}
// 		result[cpt_multi].item(n1)=1;
// 		result[cpt_multi].item(n2)=10;
// 		
// 		Writer<Texture1d> wr("test111T.tex");
// 		wr.write(result);
// 		cpt_multi++;

		
		//On relance
		
		avant_list = list_points;
		
		//Listing des points
/*		cout<<"LISTE AVANT:"<<endl;
		cout<<n1<<endl;
		for(uint i=0; i<list_points.size(); i++)
		{
			cout<<list_points[i]<<endl;
		}
		cout<<n2<<endl;*/
		
// 		compute_snake_at_1_resolution(result);
		compute_snake_at_1_resolution();
		
		for(uint i=0; i<list_points.size(); i++)
		{
			if( avant_list[i] == list_points[i] )
				cpt_points[i]++;
			if( avant_list[i] != list_points[i] )
				cpt_points[i] = 0;
		}

		
		//On calcule le snake (rajouter boucle!)
		refine_vector();
		//Listing des points
/*		cout<<"LISTE APRES:"<<endl;
		cout<<n1<<endl;
		for(uint i=0; i<list_points.size(); i++)
		{
			cout<<list_points[i]<<endl;
		}
		cout<<n2<<endl;*/
		//update old vector (list_points) with new_vector (CHECK THE '=' OPERATION!!)
		cpt_resolution++;
	}
	cout<<"STOP TOTAL!!"<<endl;
	
	
	//ECRITURE TEXTURE TEMPORAIRE!!
	for(uint j=0; j<size; j++)
	{
		result_total[0].item(j)=0;
		for(uint i=0; i<list_points.size(); i++)
		{
			if(list_points[i]==j)
			{
				result_total[0].item(j)=value;
			}
			else	
			{
			}
		}
	}
	result_total[0].item(n1)=value;
	result_total[0].item(n2)=value;
	
//	Writer<Texture1d> wer1("test111T.tex");
//	wer1.write(result);
	
	return result_total;
}


//Compute the snake for a given level
// void SulcusCorticalSnake::compute_snake_at_1_resolution(TimeTexture<float> & result_multi)
void SulcusCorticalSnake::compute_snake_at_1_resolution()
{
	std::map< int, std::map< int, int > > count;
// 	for(uint i=0; i<list_points.size(); i++)
// 	{
// 			std::map< int, int > loc;
// 			loc[ list_points[i] ]=0;
// 		count[i][ list_points[i] ]=0;
// 	}
	std::cout<<"compute_snake_at_1_resolution"<<endl;
	while( stop_condition_1_resolution()==0 )
	{
		
// 		std::cout<<"SulcusCorticalSnake : compute_snake_at_1_resolution : while stop_condition_1==0"<<std::endl;
		//updates the vertices vector
// 		cout<<"update vector!"<<endl;
		new_vector_res=list_points;
		
		
		//Process new vector
		process_list( count );
		
		//ecriture dans result_multi
// 		for(uint j=0; j<size; j++)
// 		{
// 			result_multi[cpt_multi].item(j)=0;
// 			for(uint i=0; i<list_points.size(); i++)
// 			{
// 				if(list_points[i]==j)
// 				{
// 					result_multi[cpt_multi].item(j)=i+2;
// 				}
// 				else	
// 				{
// 				}
// 			}
// 		}
// 		result_multi[cpt_multi].item(n1)=1;
// 		result_multi[cpt_multi].item(n2)=10;
// 		
// 		Writer<Texture1d> wr("test111T.tex");
// 		wr.write(result_multi);
// 		cpt_multi++;

		
	}
}



//Total stop condition computing (all points in the snake have 2 marked neighbours - comparer new_vector et list_points, apres l'insertion des nouveaux points)
int SulcusCorticalSnake::stop_total()
{
// 	std::cout<<"SulcusCorticalSnake : test_stop_total"<<std::endl;
	int stop=1;
	std::set<uint>::const_iterator itneigh;
	
	//si le truc n'a pas de point ou n'en a qu'un
	if( list_points.size()==0 )
	{
		if( are_they_neighbours(n1,n2)==1 )
			return(1);
		else{}
	}
	else
	{
		if( ( are_they_neighbours(n1, list_points[0])==1 ) && ( are_they_neighbours(n2, list_points[0]) == 1 ) )
			return(1);
	}
	
	for(uint i=0; i<list_points.size(); i++)
	{
		int cpt=0;
		
		//si on est sur l'une des 2 extremites, on fait comme i on a deja rencontre un voisin sur le snake
		if( (i==0) || (i==list_points.size()-1) ) 
		    cpt++;
		
		itneigh=neigh[ list_points[i] ].begin();
		for ( ; itneigh != neigh[ list_points[i] ].end(); ++itneigh )
		{
			if( (*itneigh)==list_points[i-1] || (*itneigh)==list_points[i+1] )
				cpt++;
		}
		if(cpt<2)
		{
			stop = 0;
		}
	}
	
	if( list_points.size() < 2 )
	{
		stop=0;
	}
	return stop;
}



//Stop condition computing (stops when points don't move anymore)
int SulcusCorticalSnake::stop_condition_1_resolution()
{
// 	std::cout<<"SulcusCorticalSnake : test_stop_condition_1_resolution"<<std::endl;
	int stop=1;
	if( new_vector_res.size() != list_points.size() )
	{
		return 0;
	}
	
	for( uint i=0; i<new_vector_res.size(); i++)
	{
		if( list_points[i]!=new_vector_res[i] )
		{
			return 0;
		}
	}
	
	return stop;
}



//Defines the 2 extremities n1 & n2 of the snake
//POUR LE MOMENT? FAIT AVEC DES TEXTURES. A METTRE A JOUR AVEC LES BUCKETS
int SulcusCorticalSnake::define_extremities()
{
	std::cout<<"SulcusCorticalSnake : define_extremities"<<std::endl;
	TimeTexture<float> tex_ext(1,size);
	init_texture_single(tex_ext);
	TimeTexture<float> tex_temp(1,size);
	init_texture_single(tex_temp);
	
	std::set<uint>::const_iterator itvois;
	
	std::vector<uint> pts;
	
		//On enleve les points non susceptibles d'etre des extremites
	
	for(uint j=0; j<size;j++)
	{
		int cpt=0;
		if(constraint[0].item(j)==value)
		{
			int mem_vois=0;
			std::set<uint>::const_iterator itvoisin;
			itvoisin=neigh[j].begin();
			for(; itvoisin!=neigh[j].end(); itvoisin++)
			{
				if(constraint[0].item(*itvoisin)==value)
				{
					cpt++;
					mem_vois=(*itvois);
				}
			}
			if( cpt==1 || cpt==3 || cpt==4 )
			{
				tex_ext[0].item(j)=value;
			}
			else
				tex_ext[0].item(j)=0;
		}
	}
	
	int cpt=0;
	for(uint j=0; j<size;j++)
	{
		if(tex_ext[0].item(j)!=0)
		{
			cpt++;
			pts.push_back(j);
		}
	}
	
	if(cpt==0)
		return 0;
	
	float path_max=0;
	float path_temp=0;
	TimeTexture<float> tex_path(1,size);
	TimeTexture<float> tex_path_temp(1,size);
	init_texture_single(tex_path);
	init_texture_single(tex_path_temp);
	
	// on prend les 2 plus loin
	
	for(uint j=0; j<pts.size() - 1; j++)
	{
			
		TimeTexture<float>tmp(1,size);
		init_texture_single(tmp);
		TimeTexture<float>rez_tmp(1,size);
		init_texture_single(rez_tmp);
		
		//En fait on prend les 2 points les plus eloignes en terme de distance geodesique
		
		tmp[0].item(pts[j])=10;
		rez_tmp[0]=MeshDistance( mesh[0], tmp[0],true);
		
		for(uint k=j+1; k<pts.size(); k++)
		{
			path_temp=rez_tmp[0].item(pts[k]);
			
			if(path_temp>path_max)
			{
				path_max=path_temp;
				n1=pts[j];
				n2=pts[k];
			}
		}
	}
	return 1;
}


//test if a vertex is already in a vector
int SulcusCorticalSnake::is_it_in_the_vector(uint vertex)
{
// 	std::cout<<"SulcusCorticalSnake : is_it_in_the_vector"<<std::endl;
// 	cout<<"in is it in the vector"<<endl;
	for(uint i=0; i<list_points.size(); i++)
	{
		if( list_points[i]==vertex)
		{
			return 1;
		}
	}
	
	if( n1==vertex || n2==vertex )
		return 1;
	
	return 0;
}



//defines new point in the snake for snake refinement
uint SulcusCorticalSnake::define_new_middle_point(uint origine, uint destination)
{
// 	std::cout<<"SulcusCorticalSnake : define_new_middle_point"<<std::endl;
	std::map<uint, float> ord_list;

	uint size_vect=0;
	
	TimeTexture<float> dist(1,size);
	TimeTexture<float> dist_retour(1,size);
	TimeTexture<float> dist_changed(1,size);
	TimeTexture<float> dist_contraint_limit(1,size);
	TimeTexture<float> carte(1,size);
	TimeTexture<float> carte_retour(1,size);
	TimeTexture<float> carte_limit(1,size);
	
	for(uint i=0; i<size; i++)
	{
		if(i==destination)
			carte[0].item(i)=10;
		else
			carte[0].item(i)=0;
		
		if(i==origine)
			carte_retour[0].item(i)=10;
		else
			carte_retour[0].item(i)=0;
	}
	
	dist[0]=MeshDistance_adapt_tex( carte[0],true,origine);
	dist_retour[0]=MeshDistance_adapt_tex( carte_retour[0],true,destination);
	
	int indicateur=0;
	if( cpt_resolution<=2 )
	{
		for(uint i=0; i<size; i++)
		{
			if( (dist[0].item(i)==FLT_MAX) || (dist_retour[0].item(i)==FLT_MAX) )
				dist_changed[0].item(i)=0;
			else
			{
				dist_changed[0].item(i)=pow( dist[0].item(i), 2) + pow(dist_retour[0].item(i) ,2);
			}
			if(dist_changed[0].item(i)!=0 && constraint[0].item(i)!=0)
				indicateur=1;
		}
	}
		
	if( indicateur==1 )
	{
// 		std::cout<<"METHODE 1"<<endl;
		dist_contraint_limit[0]=MeshDistance( mesh[0], constraint[0],true);
	
		for(uint i=0; i<size; i++)
		{
			if( dist_changed[0].item(i)==0 )
				dist_contraint_limit[0].item(i)=0;
		}
		
		//choix du point resultant
// 		std::cout<<"choix du point resultant"<<endl;
		float dist_test=10000;
		uint t=0;
		int indic_move=0;
		for(uint i=0; i<size; i++)
		{
			if( dist_changed[0].item(i)!=0 && constraint[0].item(i)!=0 )
			{
				if(dist_changed[0].item(i)<dist_test && (is_it_in_the_vector( i ) == 0) )
				{
					t=i;
					dist_test=dist_changed[0].item(i);
					indic_move=1;
				}
			}
		}
// 		std::cout<<"choix du point resultant FINI : "<<t<<endl;
		if(indic_move==1)
			return(t);
		else
		{}
	}
	
	
// 	std::cout<<"METHODE 2"<<endl;
	std::set<uint>::const_iterator itneigh;
	uint result=origine;
	float distm, temp;

	temp=dist[0].item(origine);

	int indic=0;
	float longueur=0;
	
// 	std::cout<<"SulcusCorticalSnake : define_new_middle_point : avant while_result!=destination"<<std::endl;
	while(result!=destination)
	{
// 		std::cout<<"Noeud "<<result<<std::endl;
		int result_temp=result;
		itneigh=neigh[result].begin();
		indic ++;
		for ( ; itneigh != neigh[result].end(); ++itneigh )
		{
			//Teste si origine est voisin de destination
			distm=dist[0].item(*itneigh);
			if( ( indic==1) && ( (*itneigh)==destination ) )
			{
// 				cout<<"ORIGINE voisin de DESTINATION"<<endl;
				return(origine);
			}
			//et retourne "origine" dans ce cas
			
			if( ( distm < temp ) )//&& ( tex_origine[0].item(*itneigh)!=(float)flag ) )
			{
				result_temp=(*itneigh);
				//std::cout<<"\tVoisin "<<result <<" plus proche avec dist="<<dist<< " et dist prec.="<<temp<<std::endl;
				temp=distm;
			}
			else {}
		}

		longueur+= sqrt( (vert[result][0]-vert[result_temp][0])*(vert[result][0]-vert[result_temp][0])
				+ (vert[result][1]-vert[result_temp][1])*(vert[result][1]-vert[result_temp][1])
				+ (vert[result][2]-vert[result_temp][2])*(vert[result][2]-vert[result_temp][2]) );
		result=result_temp;
		if(result!=destination)
		{
			ord_list[result]=longueur;
			size_vect++;
		}
	}
// 	std::cout<<"SulcusCorticalSnake : define_new_middle_point : apres while_result!=destination"<<std::endl;
	
	
	std::map<uint, float>::iterator ord_pt=ord_list.begin();
	float prout=1000.0;
	uint resultat=0;
	for ( ; ord_pt != ord_list.end(); ++ord_pt)
	{
		if (fabs((*ord_pt).second - (longueur/2.0)) < prout)
		{
			prout=fabs((*ord_pt).second - (longueur/2.0));
			resultat=(*ord_pt).first;
		}
	}
	
	if( is_it_in_the_vector( resultat ) == 1 )
	{
// 		std::cout<<"SulcusCorticalSnake : define_new_middle_point : is_it_in_the_vector( resultat ) == 1"<<std::endl;
		uint tmp=resultat;
		std::map<uint, float>::iterator ord_pt=ord_list.begin();
		for ( ; ord_pt != ord_list.end(); ++ord_pt)
		{
			std::set<uint>::const_iterator itvoisin;
			itvoisin = neigh[resultat].find( (*ord_pt).first );
			if( itvoisin!=neigh[resultat].end() )
			tmp=(*itvoisin);
		}
		resultat=tmp;
	}
	
	return(resultat);
	
}


int SulcusCorticalSnake::are_they_neighbours(uint one, uint two)
{

// 	std::cout<<"SulcusCorticalSnake : are_they_neighbours"<<std::endl;
	std::set<uint>::const_iterator itvoisin;
	itvoisin=neigh[one].begin();
	for(; itvoisin!=neigh[one].end(); itvoisin++)
	{
		if( (*itvoisin)==two )
		{
			return(1);
		}
	}
	return(0);

}


//refines vector list
void SulcusCorticalSnake::refine_vector()
{
// 	std::cout<<"SulcusCorticalSnake : refine_vector"<<std::endl;
	vector<uint> new_vector_temp;
	
	add_pts=0;
	
	int pos=0;
	
	uint s=list_points.size();
	
	uint new_point;
	
	//insert new point between n1 and p0
	
	if( !( neigh[n1].find(list_points[0])!=neigh[n1].end() ) )
	{
		new_point = define_new_middle_point( n1, list_points[0] );
	
	//on teste si on va bien rajouter un nouveau point (s'il est deja marque, alors on ne fait rien)
		
		list_points.insert( list_points.begin( ) , new_point );
		cpt_points.insert( cpt_points.begin( ) , 0 );
		pos=1;
		add_pts++;
	}
	
// 	std::cout<<"SulcusCorticalSnake : refine_vector : loop to insert new point between pn and pn+1"<<std::endl;
	//loop that inserts new points between pn and pn+1
	for(uint i=0; i<s-1; i++)
	{
		//adds in new vector point already existing
		
		//creates new point
		if( !( neigh[ list_points[pos] ].find(list_points[pos+1])!=neigh[ list_points[pos] ].end() ) )
		{
			new_point = define_new_middle_point( list_points[pos], list_points[pos+1] ) ;
			pos++;
			list_points.insert( list_points.begin( ) + pos , new_point );
			cpt_points.insert( cpt_points.begin( ) + pos , 0 );
			add_pts++;
		}
		pos++;
	}
	
	//adds in new vector point already existing (last point of the old vector
	
	//insert new point between pn and n2 (carefull with conditions (i.e. if list has only 1 element)

	uint tr=list_points.size();
	if( !( neigh[ list_points[tr-1] ].find(n2)!=neigh[ list_points[tr-1] ].end() ) )
	{
		new_point = define_new_middle_point( list_points[tr-1], n2 ) ;
		list_points.push_back(new_point);
		cpt_points.push_back(0);
		add_pts++;
	}

}



//computes mesh mean curvature (finite element method)
void SulcusCorticalSnake::compute_curv()
{
	std::cout<<"SulcusCorticalSnake : compute_curv"<<std::endl;
	CurvatureFactory CF;
	Curvature * curvat = CF.createCurvature(mesh,"fem");
	curv[0] = curvat->doIt();
	curvat->regularize(curv[0],1);
	curvat->getTextureProperties(curv[0]);
	delete curvat;
	
	h_min=10000;
	h_max=0;
	
	for(uint j=0; j<size; j++)
	{
		if( curv[0].item(j) > h_max )
			h_max=curv[0].item(j);

		if( curv[0].item(j) < h_min )
			h_min=curv[0].item(j);
	}
		
	for(uint j=0; j<size; j++)
	{
		curv[0].item(j) = (curv[0].item(j) -h_min) / (h_max-h_min);
	}
}



//Compute energy with new parameters
float SulcusCorticalSnake::compute_energy(uint index_courant)
{
// 	std::cout<<"SulcusCorticalSnake : compute_energy"<<std::endl;
	float result;
	
	uint origin=list_points[index_courant];
	
	SulcusCorticalSnake_energy *NRJ = new SulcusCorticalSnake_energy( list_points, index_courant, n1, n2, alpha1, alpha2, alpha3, mesh, size, h_min, h_max, curv, max, tex_distance );
	
// 	std::cout<<"SulcusCorticalSnake : compute_energy : appel de total_energy"<<std::endl;
	result=(*NRJ).total_energy();
// 	std::cout<<"SulcusCorticalSnake : compute_energy : total_energy="<<result<<std::endl;
	
	delete NRJ;
	
	return( result );
}



//Process one point of the list (index=index of the point IN THE VECTOR): first order neighborhood
void SulcusCorticalSnake::treat_list_point(uint index, std::map< int, std::map< int, int > > & count )
{
// 	std::cout<<"SulcusCorticalSnake : treat_list_point"<<std::endl;
	uint vertex_origin = list_points[index];
	uint new_index = vertex_origin; //index of the new point ON THE MESH!!
	
	//Calcul de l'energie avec le point d'origine
	float energy_min = compute_energy(index);
	
	std::set<uint>::const_iterator itneigh;
	itneigh=neigh[vertex_origin].begin();
	
// 	cout<<"voisins: "<<endl;
	//for each neighbour, computes the total energys
// 	std::cout<<"SulcusCorticalSnake : treat_list_point : for each neighbour, computes the total energy"<<std::endl;
	for ( ; itneigh != neigh[vertex_origin].end(); ++itneigh )
	{
// 		cout<<"test IS IT IN THE VECTOR"<<endl;
		if( is_it_in_the_vector( (*itneigh) )== 0 )
		{
			//Replace "index" vertex by new neighbor
			list_points[index]=(*itneigh);
			
// 			cout<<list_points[index]<<endl;
			
			//compute new energy with new point in the list
			float tmp = compute_energy(index);
		
			if( tmp < energy_min )
			{
				new_index = (*itneigh);
				energy_min = tmp;
			}
		}
	}
	
	//updates the vector with the point that minimize the energy
	list_points[index]=new_index;
	count[index][new_index]++;
// 	std::cout<<"nouveau count["<<index<<"]["<<new_index<<"]]="<<count[index][new_index]<<std::endl;
}


//randomly process the points of the vector
void SulcusCorticalSnake::process_list( std::map< int, std::map< int, int > > & count )
{
// 	std::cout<<"SulcusCorticalSnake : process_list"<<std::endl;
// 	std::cout<<"process_list"<<endl;
	srand(time(NULL));
	
	uint cpt=0;
	uint s=list_points.size();
	int tab[s];
	
	//tableau pour voir si le point du vecteur a ete traite ou pas
	for(uint i=0; i<s; i++)
		tab[i]=0;
	
// 	std::cout<<"SulcusCorticalSnake : process_list : tirage au sort des points"<<std::endl;
	while( cpt < s )
	{
		//On tire au sort un index du vecteur qui n'est pas encore sorti
		int ind=0;
		do
		{
			float r=rand();
			ind=(int)( ( s * r )/( RAND_MAX ) );
		}
		while( tab[ind]!=0 );
		
		//updates indicators
		tab[ind]++;
		cpt++;
		
		//run point "ind" update
// 		cout<<"Point traite="<< list_points[ind] <<endl;
// 		if( cpt_points[ind] < 3 )
// 		std::cout<<"SulcusCorticalSnake : process_list : treat_list_point("<<ind<<")"<<std::endl;
// 		std::cout<<"Traitement de count["<<ind<<"]["<<list_points[ind]<<"]]="<<count[ind][ list_points[ind] ]<<std::endl;
		if( count[ind][ list_points[ind] ]<=3 )
			treat_list_point(ind, count);
	}
}


////////////////////////////////////////////////////////
//Calcul distance geodesique (distance_map modif)
////////////////////////////////////////////////////////

float SulcusCorticalSnake::MeshDistance_adapt( const Texture<float> & inittex, bool allowUnreached, uint ind )
{
	Texture<float> tex;

	TimeTexture<float> result_lim(1,size);
	init_texture_single(result_lim);

	unsigned				i, n = vert.size();

	ASSERT( inittex.nItem() == n );
	tex.reserve( n );

  // neighbours map

	allowUnreached=true;

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
	std::set<unsigned>				neigh_local;
	std::set<unsigned>::iterator		in, fn;
	float					d, d2, l;
	Point3df				pos;
	int vect=0;
//	init first front

	for( i=0; i<n; ++i )
		if( tex.item(i) == 0 )
			front1.insert( std::pair<float,unsigned>( 0, i ) );

  //	loop until current front is empty
// 	int test=0;
	
	while( vect <1 )
	{
		nfront->clear();
		neigh_local.clear();

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
// 					result_lim[0].item( *in )=tex.item( *in );
					neigh_local.insert( *in );
					if( (*in)==ind )
						vect++;
				}
			}
		}

		for( in=neigh_local.begin(), fn=neigh_local.end(); in!=fn; ++in )
			nfront->insert( std::pair<float,unsigned>( tex.item( *in ), *in ) );

		tmpf = cfront;
		cfront = nfront;
		nfront = tmpf;
	}
	
	neighbours.clear();
	front1.clear();
	front2.clear();
	(*cfront).clear();
	(*nfront).clear();
	(*tmpf).clear();
	
	return( tex.item(ind) );
}

Texture<float> SulcusCorticalSnake::MeshDistance_adapt_tex( const Texture<float> & inittex, bool allowUnreached, uint ind )
{
	Texture<float> tex;
	
	TimeTexture<float> result_lim(1,size);
	init_texture_single(result_lim);

	unsigned				i, n = vert.size();

	ASSERT( inittex.nItem() == n );
	tex.reserve( n );

  // neighbours map

	allowUnreached=true;

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
	std::set<unsigned>				neigh_local;
	std::set<unsigned>::iterator		in, fn;
	float					d, d2, l;
	Point3df				pos;
	int vect=0;
//	init first front

	for( i=0; i<n; ++i )
		if( tex.item(i) == 0 )
			front1.insert( std::pair<float,unsigned>( 0, i ) );

  //	loop until current front is empty
// 	int test=0;
	
// 	cout<<"ind1="<<ind1<<" - ind2="<<ind2<<endl;
	while( vect <1 )
	{
		nfront->clear();
		neigh_local.clear();

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
// 					result_lim[0].item( *in )=tex.item( *in );
					neigh_local.insert( *in );
					if( (*in)==ind )
						vect++;
				}
			}
		}

		for( in=neigh_local.begin(), fn=neigh_local.end(); in!=fn; ++in )
			nfront->insert( std::pair<float,unsigned>( tex.item( *in ), *in ) );

		tmpf = cfront;
		cfront = nfront;
		nfront = tmpf;
	}
	
	
	neighbours.clear();
	front1.clear();
	front2.clear();
	(*cfront).clear();
	(*nfront).clear();
	(*tmpf).clear();
	
	return( tex );
}

std::map< uint, float> SulcusCorticalSnake::MeshDistance_adapt_local( const Texture<float> & inittex, uint ind1, uint ind2 )
{
	Texture<float> tex;
/*	TimeTexture<float> result_lim(1,size);
	init_texture_single(result_lim);*/
	unsigned i, n = vert.size();

	ASSERT( inittex.nItem() == n );
	tex.reserve( n );

	std::map< uint, float> map_temp_l;
	
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
	std::set<unsigned>				neigh_local;
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
		neigh_local.clear();

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
					neigh_local.insert( *in );
// 					result_lim[0].item( *in )=tex.item( *in );
					map_temp_l[ *in ] = tex.item( *in );
					if( (*in)==ind1 )
						vect[0]++;
					if( (*in)==ind2 )
						vect[1]++;
				}
			}
		}

		for( in=neigh_local.begin(), fn=neigh_local.end(); in!=fn; ++in )
			nfront->insert( std::pair<float,unsigned>( tex.item( *in ), *in ) );

		tmpf = cfront;
		cfront = nfront;
		nfront = tmpf;
	}
	
// 	cout<<"tex.item(ind1)="<<tex.item(ind1)<<" - tex.item(ind2)="<<tex.item(ind2)<<endl;
	
/*	Writer<Texture1d> wer1("distance_limitee.tex");
	wer1.write(result_lim);*/
	
	neighbours.clear();
	front1.clear();
	front2.clear();
	(*cfront).clear();
	(*nfront).clear();
	(*tmpf).clear();
	
	return( map_temp_l );
}

