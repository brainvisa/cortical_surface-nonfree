
#include <cstdlib>
#include <cortical_surface/surfacereferential/corticalConstraints.h>
#include <cortical_surface/surfacereferential/shortestPath.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshmorphomat.h>
#include <cortical_surface/surfacereferential/sulcusCorticalSnake.h>

/****************CONTRAINTS METHODS****************/

using namespace aims;
using namespace aims::meshdistance;
using namespace std;

TimeTexture<float> aims::rescaleConstraints(TimeTexture<float> tex, std::map< int, std::map<int, std::string> > map_global)
{
	int size=tex[0].nItem();
	TimeTexture<float> result(1,size);
	init_texture_single(result);
	std::map< int, std::map<int, std::string> >::iterator itMap;
	for(int i=0;i<size; i++)
	{
		itMap = map_global.begin();
		for ( ; itMap != map_global.end(); ++itMap)
		{
			std::map<int, std::string>::iterator itMap2;
			itMap2 = (*itMap).second.begin();
			for ( ; itMap2 != (*itMap).second.end(); ++itMap2)
			{
				if( (int(rint(tex[0].item(i))))==(*itMap2).first )
					result[0].item(i) = (float)((*itMap).first);
			}
		}
	}
	
	return result;
}


TimeTexture<float> aims::rescaleLongitude(TimeTexture<float> tex)
{

	int size=tex[0].nItem();
	//std::cout<<"Taille="<<size<<std::endl;
	TimeTexture<float> result(1,size);
	init_texture_single(result);
	for(int i=0;i<size; i++)
	{
		switch(int(rint(tex[0].item(i))))
		{
		//Pour l'hemisphere gauche
			case(25):  //S.C.
			case(27):  //S.C.inf
			case(29):  //S.C.sup
				result[0].item(i) = 360;
				break;
			case(61):  //S.Pe.C.sup
			case(57):  //S.Pe.C.inf
			//case(59):  //S.Pe.C.median
				result[0].item(i) = 16;
				break;
			case(1):  //F.C.L.r.asc.
				result[0].item(i) = 43;
				break;
			case(41):  //S.F.Marginal
			case(43):  //S.F.Orbitaire
				result[0].item(i) = 58;
				break;
			case(9):  //F.Cal.ant.-Sc.Cal
				result[0].item(i) = 287;
				break;
			case(19):  //F.P.O
// 			case(21):  //F.P.O.inf
				result[0].item(i) = 303;
				break;
			case(63):  //S.Po.C.Sup
			case(15):  //F.I.P.Po.C.inf
			case(17):  //F.I.P.Po.C.sup
				result[0].item(i) = 340;
				break;
		
		//Pour l'hemisphere droite
			
			case(26):  //S.C.
			case(28):  //S.C.inf
			case(30):  //S.C.sup
				result[0].item(i) = 360;
				break;
			case(62):  //S.Pe.C.sup
			case(58):  //S.Pe.C.inf
			case(60):  //S.Pe.C.median
				result[0].item(i) = 16;
				break;
			case(2):  //F.C.L.r.asc.
				result[0].item(i) = 43;
				break;
			case(42):  //S.F.Marginal
			case(44):  //S.F.Orbitaire
				result[0].item(i) = 58;
				break;
			case(10):  //F.Cal.ant.-Sc.Cal
				result[0].item(i) = 287;
				break;
			case(20):  //F.P.O
			case(22):  //F.P.O.inf
				result[0].item(i) = 303;
				break;
			case(64):  //S.Po.C.Sup
			case(16):  //F.I.P.Po.C.inf
			case(18):  //F.I.P.Po.C.sup
				result[0].item(i) = 340;
				break;

		
		}
	}
	return result;

}


TimeTexture<float> aims::rescaleLatitude(TimeTexture<float> tex)

{	int size=tex[0].nItem();
	TimeTexture<float> constraint(1,size);
	init_texture_single(constraint);
	//Contraintes
	for(int i=0; i<size; i++)
	{
		switch(int(rint(tex[0].item(i))))
		{
			//case(14):  //Pole calleux
			//	constraint[0].item(i)=1;
			//	break;
			
		//Pour l'hemisphere gauche (mise a l'echelle pour l'histoire des poles,
		//et +1, car diffusion entre 1 et 181
			
			case(5):  //F.C.M.ant
// 			case(7):  //F.C.M.asc
// 			case(3):  //F.C.M.AMS
			case(77):  //S.intraCing
			case(79):  //S.s.P
			case(11):  //F.Coll
				constraint[0].item(i)=51;
				break;
			case(45):  //S.F.Sup
// 			case(47):  //S.F.Sup.ant
// 			case(49):  //S.F.Sup.moy
// 			case(51):  //S.F.Sup.post
			case(53):  //S.O.T.lat.post
			case(55):  //S.Olf
				constraint[0].item(i)=84;
				break;
			case(39):  //S.F.inter
			case(65):  //S.T.i.ant
			case(75):  //S.T.s.ter.asc.post
// 			case(69):  //S.T.post
			case(67):  //S.T.i.post
				constraint[0].item(i)=96;
				break;
			case(33):  //S.F.inf
// 			case(35):  //S.F.inf.moy
// 			case(37):  //S.F.inf.post
//			case(13):  //F.I.P.Horiz
			case(73):  //S.T.s.ter.asc.ant
			case(71):  //S.T.s
				constraint[0].item(i)=114;
				break;
		

		//Pour l'hemisphere droite
				
		
			case(6):  //F.C.M.ant
			case(8):  //F.C.M.asc
			case(4):  //F.C.M.AMS
			case(78):  //S.intraCing
			case(80):  //S.s.P
			case(12):  //F.Coll
				constraint[0].item(i)=51;
				break;
			case(46):  //S.F.Sup
			case(48):  //S.F.Sup.ant
			case(50):  //S.F.Sup.moy
			case(52):  //S.F.Sup.post
			case(54):  //S.O.T.lat.post
			case(56):  //S.Olf
				constraint[0].item(i)=84;
				break;
			case(40):  //S.F.inter
			case(66):  //S.T.i.ant
			case(76):  //S.T.s.ter.asc.post
			case(70):  //S.T.post
			case(68):  //S.T.i.post
				constraint[0].item(i)=96;
				break;
			case(34):  //S.F.inf
			case(36):  //S.F.inf.moy
			case(38):  //S.F.inf.post
			case(14):  //F.I.P.Horiz
			case(74):  //S.T.s.ter.asc.ant
			case(72):  //S.T.s
				constraint[0].item(i)=114;
				break;

			default:
				constraint[0].item(i)=0;
				break;

		}
	}
	return constraint;
}


std::map< int, std::map<int, std::string> > aims::createCorrespMap( std::string & _adr_cor, std::string & _adr_file, std::string & side )
{
	cout << "Parsing constraint values"<< endl;
			
	std::map< int, std::map<int, std::string> > map_global;
	
	const char * adr_cor= _adr_cor.c_str();
	const char * adr_file= _adr_file.c_str();
	
	FILE *corres; //correspondance between contraints names and real values
	FILE *file; //correspondance between contraints names and projected values
	cout << "Opening files" << endl;
	cout << "Parsing constraint values from " << adr_cor <<" and "<<adr_file<< endl;
	if ((corres=fopen(adr_cor, "r")) == NULL)
	{
		cerr << "Cannot open file " << adr_cor << endl;
		exit(EXIT_FAILURE);
	}
	if ((file=fopen(adr_file, "r")) == NULL)
	{
		cerr << "Cannot open file " << adr_file << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Opening files OK" << endl;
	
	std::string arg1;
	std::string arg2;
/*	char *  arg1;
	char *  arg2;*/
	int val_contraint;
	int val_projection;
	
// 	cout << "Reading files" << endl;
	
	//LECTURE
	map< int, std::string> map_local2;
	
// 	std::vector< std::map<int, std::string > > vect;
/*	while (!feof(file))
	{
		char c;
		fscanf(file, "%c", &c);
		cout<<c;
	}*/
		
		
	while (!feof(file))
	{
		char c[40];
		fscanf(file, "%s %i\n", c, &val_projection);
		arg2=c;
		map_local2[val_projection]=arg2;
	}
	
	map< int, std::vector< std::string > > map_local1;
	while (!feof(corres))
	{
		char c[40];
		char a[5];
		fscanf(corres, "%s %s %i\n", a, c, &val_contraint);
		arg1=c + side;
		map_local1[val_contraint].push_back(arg1);
	}
	
// 	std::string, std::pair<float, float> p;
	
/*	cout<<"LISTE!!"<<endl;
	std::map<int, std::vector< string > >::iterator ord_ptz=map_local1.begin();
	for ( ; ord_ptz != map_local1.end(); ++ord_ptz)
	{
		for ( uint i=0; i<(*ord_ptz).second.size(); i++ )
			cout<<(*ord_ptz).first<<" - "<<(*ord_ptz).second[i]<<endl;
	}
	
	cout<<"LISTE2!!"<<endl;
	std::map<int, std::string >::iterator ord_pt2z=map_local2.begin();
	for ( ; ord_pt2z != map_local2.end(); ++ord_pt2z)
	{
		cout<<(*ord_pt2z).first<<" - "<<(*ord_pt2z).second<<endl;
	}*/
	
	std::map<int, std::vector< std::string > >::iterator ord_pt;
	std::map<int, std::string >::iterator ord_pt2;
	
	ord_pt=map_local1.begin();
	for ( ; ord_pt != map_local1.end(); ++ord_pt)
	{
		for ( uint i=0; i<(*ord_pt).second.size(); i++ )
		{
			ord_pt2=map_local2.begin();
			for ( ; ord_pt2 != map_local2.end(); ++ord_pt2)
			{
				if(strcmp( ( (*ord_pt).second[i] ).c_str(), ( (*ord_pt2 ).second ).c_str() )==0 )
				{
	//				tmp=(*ord_pt2).first;
					map_global[(*ord_pt).first][(*ord_pt2 ).first]=(*ord_pt2 ).second;
				}
			}
		}
		/*if(tmp!=-1)
		{
			int fl=(*ord_pt).first;
			p.first=tmp;
			map_global[fl]=p;
		}*/
	}
	
	cout << "listing" << endl;
	
	std::map< int, std::map<int, std::string> >::iterator itMap = map_global.begin();
	for ( ; itMap != map_global.end(); ++itMap)
	{
		std::map<int, std::string >::iterator ord_pt2;
		ord_pt2=(*itMap).second.begin();
		for ( ; ord_pt2 != (*itMap).second.end(); ++ord_pt2)
		{
			cout<<(*itMap).first<<" - "<<(*ord_pt2).first<<" - "<<(*ord_pt2).second<<endl;
		}
	}
	
	cout << "done" << endl;
	
	
	fclose(corres);
	fclose(file);
	
	return map_global;
}

TimeTexture<float> aims::constraintCleaner(TimeTexture<float> & texture, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh, float contr, float curvature, float elasticity )
{
	int size=texture[0].nItem();
	Texture<float> tex;
	TimeTexture<float> tex_process(1,size);
	init_texture_single(tex_process);
			
	tex=texture[0];

	//label des poles �ne pas fermer topo.
// 	std::set<float> set_poles;
// 	//liste des valeurs des poles
// 	set_poles.insert(1);
// 	set_poles.insert(180);

	//On fait � pour chaque valeur
	
	
	std::vector<float> vals;
	std::vector<float>::iterator it;
	
	
	//Pour tester, a enlever apres///////////////////
/*	for(int j=0; j<size; j++)
	{
		if(texture[0].item(j)!=33)
			texture[0].item(j)=0;
	}*/
	/////////////////////////////////////////////////
	
	
	for(int i=0; i<size; i++)
	{
		if(texture[0].item(i)!=0)
		{
			bool present=false;
			it=vals.begin();
			for( ; it!=vals.end(); ++it) 
				if( (*it)==texture[0].item(i) )
					present=true;
			if(present==false)
				vals.push_back( texture[0].item(i) );
		}
	}
	
	
	//pour chaque valeur
	
	it=vals.begin();
	
	for(int j=0; j<size; j++)
	{
		tex_process[0].item(j)=0;
	}
	
	
// 	float alpha1=20, alpha2=500, alpha3=200;
// 	float alpha1=20, alpha2=20, alpha3=50;
	
	for(; it!=vals.end(); ++it)
	{
		TimeTexture<float> thin(1,size);
		std::cout<<"Cleaner pour Valeur="<<(*it)<<std::endl;
		//if( (*it)==71 )
// 		thin=EconstraintCleaner(texture,neigh,mesh, (*it) );
		SulcusCorticalSnake *k = new SulcusCorticalSnake( texture, (*it), contr, curvature, elasticity, mesh );
		thin = (*k).compute_snake();
		
		for(int j=0; j<size; j++)
		{
			if(thin[0].item(j)!=0)
				tex_process[0].item(j)=(*it);
		}
		delete k;
	}
// 	std::cout<<"Done Cleaner"<<std::endl;
	
	return tex_process;
}

TimeTexture<float> aims::EconstraintCleaner(TimeTexture<float> & texture, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh, float valeur)
{
	//Separation des contraintes non connexes
	
	int size=texture[0].nItem();
	TimeTexture<float> result(1,size);
	std::set<uint>::const_iterator itvois;
	
	//Texture qui sera le resultat de la procedure
	init_texture_single(result);
	
	//Textures temporaires
	TimeTexture<float> tex_result(1,size);
	init_texture_single(tex_result);
	TimeTexture<float> tex_rem(1,size);
	init_texture_single(tex_rem);
	TimeTexture<float> tex_single_cc(1,size);
	init_texture_single(tex_single_cc);
	float val=0;
	
	//identification des composantes connexes
//          cout<<"EConstraints: Connected parts Id"<<endl;
	for(int i=0; i<size; i++)
	{
		tex_single_cc[0].item(i) = texture[0].item(i);
		if( (texture[0].item(i) == valeur)  && (tex_result[0].item(i)==0) )
		{
			val++;
			cout<<"Valeur="<<valeur<<endl;
			recurs_proc(i, texture, tex_result, neigh, val);
		}
		
	}
// 		Writer<Texture1d> ws2res("resultConstraints.tex");
// 		ws2res.write(tex_result);
	
	TimeTexture<float> res(1,size);
	init_texture_single(res);
	//On traite chaque morceau de la contrainte pour ebarbuler
	float value;
	for( value=1; value<=val; value++)
	{
		TimeTexture<float> tex_ext(1,size);
		init_texture_single(tex_ext);
		TimeTexture<float> tex_temp(1,size);
		init_texture_single(tex_temp);
	
		int compteur=0;
		//On enl�e les points non susceptibles d'etre des extremites
		for(int j=0; j<size;j++)
		{
			int cpt=0;
			if(tex_result[0].item(j)==value)
			{
				int mem_vois=0;
				std::set<uint>::const_iterator itvoisin;
				itvoisin=neigh[j].begin();
				for(; itvoisin!=neigh[j].end(); itvoisin++)
				{
					if(tex_result[0].item(*itvoisin)==value)
					{
						cpt++;
						mem_vois=(*itvoisin);
					}
				}
				if( cpt==1 || cpt==3 || cpt==4 )
				{
					compteur++;
					tex_ext[0].item(j)=value;
				}
				else
					tex_ext[0].item(j)=0;
			}
		}
		int tab[compteur];
		int ind=0;
		for(int j=0; j<size;j++)
		{
			if(tex_ext[0].item(j)==value)
			{
				tab[ind]=j;
				ind++;
			}
		}
		
		//On passe un coup de PCC pour reperer les extremites les plus eloignees
		int path_max=0;
		int path_temp=0;
		int ext1=0, ext2=0;
		TimeTexture<float> tex_path(1,size);
		TimeTexture<float> tex_path_temp(1,size);
		init_texture_single(tex_path);
		init_texture_single(tex_path_temp);
		for(int j=0; j<compteur-1; j++)
		{
			for(int k=j+1; k<compteur; k++)
			{
				path_temp=0;
				
				GraphPath<float> gr;
				tex_path=gr.process(tex_result,mesh,value,tab[j],tab[k]);
				//delete &gr;
				
				for(int x=0; x<size;x++)
				{
					if(tex_path[0].item(x)==value)
						path_temp++;
				}
				if(path_temp>path_max)
				{
					path_max=path_temp;
					ext1=tab[j];
					ext2=tab[k];
				}
			}
		}
		
		//on recupere le plus court chemin
		if(path_max!=0)
		{
			GraphPath<float> gr1;
			tex_path_temp=gr1.process(tex_result,mesh,value,ext1,ext2);
			//delete &gr1;
		}
		
		for(int j=0; j<size; j++)
		{
			if(tex_path_temp[0].item(j)!=0)
			{
				res[0].item(j)=value;
				tex_rem[0].item(j)=tex_path_temp[0].item(j);
			}
		}
	}
	
	for(int j=0; j<size; j++)
	{
		tex_result[0].item(j)=res[0].item(j);
		tex_single_cc[0].item(j)=0;
		if(tex_result[0].item(j)!=0)
			tex_single_cc[0].item(j)=valeur;
	}
	
	Writer<Texture1d> wresultrr("tex_result_PRE.tex");
	wresultrr.write(tex_result);
	//On traite chaque morceau de la contrainte
// 	cout<<"chaque morceau de "<<valeur<<endl;
	
	std::vector< std::pair < int, float > > vect_extr;
	
	TimeTexture<float> tex_t(1,size);
	init_texture_single(tex_t);
	
	for( value=1; value<=val; value++)
	{
//                 cout<<"Value="<<value<<endl;
		TimeTexture<float> tex_ext(1,size);
		init_texture_single(tex_ext);
		TimeTexture<float> tex_temp(1,size);
		init_texture_single(tex_temp);
	
		int compteur=0;
		//On enleve les points non susceptibles d'etre des extremites
   //              cout<<"Removing non Xt points"<<endl;
		for(int j=0; j<size;j++)
		{
			int cpt=0;
			if(tex_result[0].item(j)==value)
			{
				int mem_vois=0;
				std::set<uint>::const_iterator itvoisin;
				itvoisin=neigh[j].begin();
				for(; itvoisin!=neigh[j].end(); itvoisin++)
				{
					if(tex_result[0].item(*itvoisin)==value)
					{
						cpt++;
						mem_vois=(*itvoisin);
					}
				}
				if( cpt==1 )
				{
					compteur++;
					tex_ext[0].item(j)=value;
				}
				else
					tex_ext[0].item(j)=0;
			}
		}
		int tab[compteur];
		int ind=0;
		for(int j=0; j<size;j++)
		{
			if(tex_ext[0].item(j)==value)
			{
				tab[ind]=j;
				ind++;
				
				//Stockage des extremites dans le vector de pair
				std::pair< int, float > p(j, value);
				vect_extr.push_back(p);
				tex_t[0].item(j)=value;
			}
		}
		
	}
		
// 	cout<<"voila"<<endl;
		
// 	cout<<"relie les CC"<<endl;
	///////////////////////////////////////////////////////////////
	// IL FAUT RAJOUTER: ON RELIE TOUTES LES COMPOSANTES CONNEXES
	// POUR NE FAIRE PLUS QU'UNE SEULE COMPOSANTE CONNEXE
	///////////////////////////////////////////////////////////////
	
		//On passe un coup de PCC pour reperer les extremites les plus eloignees
	
	float path_max=0;
	float path_temp=0;
	int ext1=0, ext2=0;
	TimeTexture<float> tex_path(1,size);
	TimeTexture<float> tex_path_temp(1,size);
	init_texture_single(tex_path);
	init_texture_single(tex_path_temp);
	TimeTexture<float> tex_path_extr(1,size);
	init_texture_single(tex_path_extr);
// 	cout<<"taille du vect_extr="<<vect_extr.size()<<endl;
	for(uint j=0; j<vect_extr.size() - 1; j++)
	{
// 		cout<<"j="<<j<<endl;
			
		TimeTexture<float>tmp(1,size);
		init_texture_single(tmp);
		TimeTexture<float>rez_tmp(1,size);
		init_texture_single(rez_tmp);
		
		//En fait on prend les 2 points les plus eloignes en terme de distance geodesique
		
		tmp[0].item(vect_extr[j].first)=10;
		rez_tmp[0]=MeshDistance( mesh[0], tmp[0],true);
		
		for(uint k=j+1; k<vect_extr.size(); k++)
		{
// 			cout<<"indice j="<<vect_extr[j].first<<"indice k="<<vect_extr[k].first<<endl;
// 			path_temp=0;
			
//				Writer< Texture1d > wrtr("test_pcc.tex");
//				wrtr.write(tex_result);
//				cout<<"PCC..."<<endl;
// 			GraphPath<float> gr;
// 			tex_path=gr.process(tex_single_cc,mesh,valeur,vect_extr[j].first,vect_extr[k].first);
/*			for(int x=0; x<size;x++)
			{
			if(tex_path[0].item(x)==valeur)
			path_temp++;
		}*/
			
			path_temp=rez_tmp[0].item(vect_extr[k].first);
			
			if(path_temp>path_max)
			{
// 				cout<<"path_max="<<path_max<<" et path_temp="<<path_temp<<endl;
				path_max=path_temp;
				ext1=vect_extr[j].first;
				ext2=vect_extr[k].first;
			}
		}
		
		
	}
	
	//On relie tout ca
	
	TimeTexture<float> link_extr(1,size);
	init_texture_single(link_extr);
	TimeTexture<float> mark_extr(1,size);
	init_texture_single(mark_extr);
	TimeTexture<float> pass(1,size);
	init_texture_single(pass);
	
	pass[0].item(ext1)=10;
	pass[0].item(ext2)=10;
	
// 	cout<<"extrs="<<ext1<<" et "<<ext2<<endl;
	//pour chaque extremite
	for(uint i=0; i<vect_extr.size(); i++)
	{
		//On récupère l'extremite la plus proche ne faisant pas partie de la meme CC
		int index_proche=0;
		TimeTexture<float>tmp(1,size);
		init_texture_single(tmp);
		TimeTexture<float>rez_tmp(1,size);
		init_texture_single(rez_tmp);
		
		tmp[0].item(vect_extr[i].first)=10;
		rez_tmp[0]=MeshDistance( mesh[0], tmp[0],true);
		
		float dist_min=1111111;
		int trigger=0;
		//for(uint j=0; j<vect_extr.size(); j++)
		if( vect_extr[i].first!=ext1 && vect_extr[i].first!=ext2 )//&& pass[0].item(vect_extr[i].first)==0) //Si on n'est pas sur un point d'extremite globale
		{
// 			cout<<"presque presque Trigger"<<endl;
// 			cout<<"vect(i)="<<vect_extr[i].first<<" - val="<<vect_extr[i].second<<endl;
			//for(uint j=0; j<size; j++)
			for(uint j=0; j<vect_extr.size(); j++)
			{
// 				cout<<"vect(j)="<<vect_extr[j].first<<" - val="<<vect_extr[j].second<<" - pass="<<pass[0].item(vect_extr[j].first)<<endl;
				//if( tex_single_cc[0].item( j ) == valeur  && j!=ext1 && j!=ext2 )
				if( vect_extr[j].first!=ext1 && vect_extr[j].first!=ext2 && pass[0].item(vect_extr[j].first)==0 )
				{
// 					cout<<"presque Trigger"<<endl;
					if( ( rez_tmp[0].item( vect_extr[j].first ) < dist_min) && ( vect_extr[j].second != vect_extr[i].second ) )
					{
						dist_min=rez_tmp[0].item( vect_extr[j].first );
						index_proche=vect_extr[j].first;
						pass[0].item(index_proche)=10;
						pass[0].item(vect_extr[i].first)=10;
			//			dist_min=rez_tmp[0].item( j );
			//			index_proche=j;
						trigger=1;
// 						cout<<"Trigger"<<endl;
					}
				}
			}
		}
		//ok
		
		//Si la distance du sommet en question n'est pas trop grande
		//i.e. superieure à sa plus grande distance intra elle-même
		
		float max_intra=0;
		for(uint j=0; j<vect_extr.size(); j++)
		{
			if( vect_extr[j].second == vect_extr[i].second )
				if( rez_tmp[0].item( vect_extr[j].first ) > max_intra )
					max_intra=rez_tmp[0].item( vect_extr[j].first );
		}
		
		
		//Now, on relie les 2 extremites en question SI la condition d'avant
		//if( dist_min < max_intra && trigger==1 && mark_extr[0].item(index_proche)==0 )
		//if( dist_min < max_intra && trigger==1 )
		if( trigger==1 )
		{
//			cout<<"ok boy"<<endl;
			findNearNeigh(vect_extr[i].first, index_proche, tex_single_cc, (int)valeur, mesh, neigh );
			mark_extr[0].item(index_proche)=10;
		}
	}
	
//	cout<<"voila"<<endl;
		
	cout<<"PCC"<<endl;
	
	////////////////////////////////////////////////////////////////
	// VOILA
	////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////
	// Maintenant on fait un pcc pour garder un seul chemin
	////////////////////////////////////////////////////////////////
		
	
		
	//On garde la plus grande des CC qui reste
	
	TimeTexture<float> t(1,size);
	init_texture_single(t);
	float vale=0;
	//on compte le nb de CC
	for(int i=0; i<size; i++)
	{
		if( (tex_single_cc[0].item(i) == valeur)  && (t[0].item(i)==0) )
		{
			vale++;
			recurs_proc(i, tex_single_cc, t, neigh, vale);
		}
	}
	
	//on compte le nb de sommets par CC
	
	int ta[(int)vale];
	for(int i=0; i<(int)vale; i++)
		ta[i] =0;
	
	for(int i=0; i<size; i++)
	{
		if( t[0].item(i)!=0 )
		{
			float ij=t[0].item(i)-1;
			ta[ (int)ij ] ++;
			//cout<<"+-ij="<<ij<<"-ta[ij-1]="<<ta[ (int)ij ];
		}
	}
	
	int ind_max=0;
	int nb_max=0;
	for(int i=0; i<(int)vale; i++)
	{
		//cout<<"taille de CC "<<i+1<<"="<<ta[i]<<endl;
		if( ta[i] > nb_max )
		{
			ind_max=i+1;
			nb_max=ta[i];
		}
	}
	
	//cout<<"nb de cc="<<vale<<endl;
	//cout<<"CC max="<<ind_max<<endl;
	//on vire de l'image les sommets appartenant aux CC les plus petites
	for(int i=0; i<size; i++)
	{
		if(t[0].item(i) != ind_max)
			tex_single_cc[0].item(i) = 0;
	}
	
	
	//if(valeur==71 || valeur==67 || valeur==43)
	if(valeur==71)
	{
		Writer<Texture1d> wresult("tex_result_lat_STs.tex");
		wresult.write(tex_single_cc);
		Writer<Texture1d> wresultt("tex_result_lat_STs_extr.tex");
		wresultt.write(tex_t);
	}
	
	//on recupere les pts les plus eloignes s'ils ont ete effaces (on connait ind1)
	int ind1=0, ind2=0, fl=0;
	if( tex_single_cc[0].item(ext1) == 0 )
	{
		ind1=ext2;
		fl=1;
		//cout<<"ext1=0"<<endl;
	}
	if( tex_single_cc[0].item(ext2) == 0 )
	{
		ind1=ext1;
		fl=2;
		//cout<<"ext2=0"<<endl;
	}
	
	if(fl!=0)
	{
		cout<<"Remplacement xtremite"<<endl;
		TimeTexture<float> tmpt(1,size);
		TimeTexture<float> rez_tmpt(1,size);
		float val_maxy=0;
		tmpt[0].item(ind1)=10;
		rez_tmpt[0]=MeshDistance( mesh[0], tmpt[0],true);
		for(int i=0; i<size; i++)
		{
			if( tex_single_cc[0].item(i)!=0 && rez_tmpt[0].item(i)>val_maxy)
			{
				ind2=i;
				val_maxy=rez_tmpt[0].item(i);
			}
		}
		if(fl==1)
		{
			ext1=ind2;
		//	cout<<"ext1="<<ext1<<endl;
		}
		else
		{
			ext2=ind2;
		//	cout<<"ext2="<<ext2<<endl;
		}
	}
	
	
	//cout<<"longest pcc"<<endl;
	
	
	
	
	//PCC entre les 2 points les plus eloignes
	if(path_max!=0)
	{
//		GraphPath<float> gr1;
//		cout<<"longest pcc"<<endl;
		GraphPath<float> gr;
		tex_path_temp=gr.process(tex_single_cc,mesh,valeur,ext1,ext2);
		//delete &gr;
	}
	tex_path_extr[0].item(ext1)=valeur;
	tex_path_extr[0].item(ext2)=valeur;
//                  cout<<"done longest pcc"<<endl;
	for(int j=0; j<size; j++)
	{
		if(tex_path_temp[0].item(j)!=0)
		{
			result[0].item(j)=valeur;
			tex_rem[0].item(j)=tex_path_temp[0].item(j);
		}
		tex_path_temp[0].item(j)=0;
	}
	
	
	//2° passe pour virer les boucles à la con
	GraphPath<float> grt;
	tex_path_temp=grt.process(result,mesh,valeur,ext1,ext2);
	
	//if(valeur==71 || valeur==67 || valeur==43 )
	if(valeur==71)
	{
		Writer<Texture1d> wresult("result_lat.tex");
		wresult.write(tex_path_temp);
		Writer<Texture1d> wresulte("tex_rem.tex");
		wresulte.write(tex_rem);
		Writer<Texture1d> wresultex("tex_path_extr.tex");
		wresultex.write(tex_path_extr);
	}
	
	return(tex_path_temp);
	//return(result);
	
	
// 		GraphPath<float> gr;
// 		
// 		//On passe un coup de PCC pour reperer les extremites les plus eloignees
//              cout<<"PCC constraints"<<endl;
// 		int path_max=0;
// 		int path_temp=0;
// 		int ext1=0, ext2=0;
// 		TimeTexture<float> tex_path(1,size);
// 		TimeTexture<float> tex_path_temp(1,size);
// 		init_texture_single(tex_path);
// 		init_texture_single(tex_path_temp);
// 		for(int j=0; j<compteur-1; j++)
// 		{
// 			for(int k=j+1; k<compteur; k++)
// 			{
// 				path_temp=0;
// //				Writer< Texture1d > wrtr("test_pcc.tex");
// //				wrtr.write(tex_result);
// //				cout<<"PCC..."<<endl;
// 				tex_path=gr.process(tex_result,mesh,value,tab[j],tab[k]);
// //				cout<<"OK!"<<endl;
// 				for(int x=0; x<size;x++)
// 				{
// 					if(tex_path[0].item(x)==value)
// 						path_temp++;
// 				}
// 				if(path_temp>path_max)
// 				{
// 					path_max=path_temp;
// 					ext1=tab[j];
// 					ext2=tab[k];
// 				}
// 			}
// 		}
// //                 cout<<"longest pcc"<<endl;
// 		
// 		//on recupere le plus long des pcc
// 		if(path_max!=0)
// 		{
// // 			GraphPath<float> gr1;
// 			tex_path_temp=gr.process(tex_result,mesh,value,ext1,ext2);
// 		}
// 		
// //                 cout<<"done longest pcc"<<endl;
// 		for(int j=0; j<size; j++)
// 		{
// 			if(tex_path_temp[0].item(j)!=0)
// 			{
// 				result[0].item(j)=valeur;
// 				tex_rem[0].item(j)=tex_path_temp[0].item(j);
// 			}
// 		}
// 		
// 	}
	

//         cout<<"Removing small parts"<<endl;
	//Removing small pieces
// 	int tailles[(int)val];
// 	int moy=0;
// 	std::cout<<"Contrainte="<<valeur<<std::endl;
// 	for(int i=0; i<(int)val; i++)
// 	{
// 		tailles[i]=nb_vertex(tex_rem, i+1, size);
// 		moy+=tailles[i];
// 		//Affichage de la taille du morceau
// 		std::cout<<tailles[i]<<std::endl;
// 	}
// 	moy/=(int)val;
// 	
// 	std::cout<<"moyenne="<<moy<<std::endl;
// 	
// 	int seuil=moy/2;
// 	std::cout<<"seuil="<<seuil<<std::endl;
// 	for(int i=0; i<(int)val; i++)
// 	{
// 		if(tailles[i]<seuil)
// 		{
// 			for(int j=0; j<size; j++)
// 				if(tex_rem[0].item(j)==(i+1))
// 					result[0].item(j)=0;
// 		}
// 	}
	
// 	return(result);
}


void aims::recurs_proc(int index, TimeTexture<float> & tex, TimeTexture<float> & result, std::vector<std::set<uint> > & neigh, float & value)
{
	std::set<uint>::const_iterator itvois;
	int cpt=0;

	//Test d'arret
	itvois=neigh[index].begin();

	for(; itvois!=neigh[index].end(); itvois++)
	{
		if( (tex[0].item(*itvois)!=0) && (result[0].item(*itvois)==0) )
		{
			cpt++;
		}
	}
	if(cpt==0)  //on n'a plus de voisin a traiter
	{
		return ;
	
	}
	//sinon on traite les voisins non marques (enfin, on essaie...)
	else
	{
		itvois=neigh[index].begin();
		for(; itvois!=neigh[index].end(); itvois++)
		{
			if( (tex[0].item(*itvois)!=0) && (result[0].item(*itvois)==0) )
			{
				cpt--;
				result[0].item(*itvois)=value;
				recurs_proc( (int)(*itvois), tex, result, neigh, value);
			}
		}
		return ;
	}
}


int aims::nb_vertex(TimeTexture<float>tex, float value, int size)
{
	int cpt=0;
	for(int i=0; i<size; i++)
	{
		if(tex[0].item(i)==value)
			cpt++;
	}
	return cpt;
}



TimeTexture<float> aims::origin_meridian(TimeTexture<float> & tex, int nord, int sud, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh, TimeTexture<float> &  poles)
{
	int size=tex[0].nItem();
	
	std::vector<unsigned> vect_side_meridian;

	TimeTexture<float> pieces_meridian(1,size);
	TimeTexture<float> constraint_meridian(1,size);
	TimeTexture<float> side_meridian(1,size);
	TimeTexture<float> thin(1,size);
	TimeTexture<float> testmor(1,size);

	init_texture_single(pieces_meridian);
	init_texture_single(constraint_meridian);
	init_texture_single(side_meridian);
	init_texture_single(thin);
	init_texture_single(testmor);

	for(int i=0;i<size; i++) 
	{

		if(tex[0].item(i) == 360)
		{
			pieces_meridian[0].item(i) = 360;
// 			cout<<"PIECEMERIDIAN"<<endl;
		}
		else
			pieces_meridian[0].item(i) = 0;

	}
	for(int i=0;i<size; i++) 
	{

		if(tex[0].item(i) == 360)
			testmor[0].item(i) = 360;
		else
			testmor[0].item(i) = 0;

	}
	testmor[0].item(nord) = 360;
	testmor[0].item(sud) = 360;
/* 	Writer<Texture1d> wTp1("test_morceaux_meridiens.tex");
 	wTp1.write(pieces_meridian);
 	Writer<Texture1d> wTp2("test_morceaux_avec_poles.tex");
 	wTp2.write(testmor);*/
	//Relie les 2 poles en passant par le meridien d'origine
        cout<<"MeridianLink"<<endl;
	meridianLink( pieces_meridian, constraint_meridian, 360, sud, nord, neigh, mesh, poles);//, constraint_meridian );

	
	cout<<"Fin MeridianLink"<<endl;
	
// 	Writer<Texture1d> wT3ff("pieces_meridian.tex");
// 	wT3ff.write(pieces_meridian);
	
// 	Writer<Texture1d> wT3ffd("constraint_meridian.tex");
// 	wT3ffd.write(constraint_meridian);
	
	for(int i=0; i<size; i++)
	{
		if( constraint_meridian[0].item(i)!=0 || pieces_meridian[0].item(i)!=0 )
			constraint_meridian[0].item(i)=360;
	}
	
	std::set<uint>::const_iterator itneigh;
	
	int cptVois=0;
	
	for(int i=0; i<size; i++)
	{
		if( constraint_meridian[0].item(i)!=0 && i!=nord && i!=sud )
		{
			cptVois=0;
			itneigh=neigh[i].begin();
			for ( ; itneigh != neigh[i].end(); ++itneigh )
			{
				if( constraint_meridian[0].item(*itneigh)!=0 )
					cptVois++;
			}
			if( cptVois<=1 )
			{
				itneigh=neigh[i].begin();
				for ( ; itneigh != neigh[i].end(); ++itneigh )
				{
					constraint_meridian[0].item(*itneigh) = constraint_meridian[0].item(i);
				}
			}
		}
	}
//  	Writer<Texture1d> wT3f("mer_fermeture.tex");
//   	wT3f.write(constraint_meridian);
	
	std::cout<<"Debut graph"<<std::endl;
	GraphPath<float> gr;
	thin=gr.process(constraint_meridian, mesh, 360, nord, sud);

	std::cout<<"Fin graph"<<std::endl;
	
// 	Writer<Texture1d> wT3h("test_long_thin.tex");
// 	wT3h.write(thin);
	
	//Fabriquer le vecteur du meridien d'origine
	//vect_side_meridian=buildOriginVector( thin, nord, sud, neigh );

	//Marque les 2 cotes du meridien d'origine
	side_meridian=originNeighbourgs(thin, nord, sud, mesh, neigh, poles);

	//On remet la texture de contrainte a jour avec le nouveau meridien d'origine
	for(int i=0;i<size; i++) 
	{
		if(tex[0].item(i) == 360)
			tex[0].item(i) = 0;
			
		if(thin[0].item(i) == 360)
			tex[0].item(i) = 360;
		else {}
	}

	return side_meridian;
}


void aims::meridianLink(TimeTexture<float> & origine, TimeTexture<float> & finish, int flag, int nord, int sud, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh, TimeTexture<float> & poles)
{
	int size=origine[0].nItem();

	TimeTexture<float> result_link(1,size);
	TimeTexture<float> result_pairs(1,size);
	TimeTexture<float> temp(1,size);
	TimeTexture<float> temp1(1,size);
	std::vector< std::vector< int > > meridian_parts;
	std::vector< std::vector< int > > meridian_pairs;


	//iterators
	std::set<uint>::const_iterator itvois;
	int index_mem1=0,cpt=0;
	bool ind_final=false;

	for(int i=0;i<size;i++) 
	{
		result_pairs[0].item(i)=0;
		result_link[0].item(i)=0;
		temp[0].item(i)=0;
		temp1[0].item(i)=0;
	}
	
	int avance=0;
        int c=1;
	
	//**************************************************************//
	//**********Identification des composantes connexes*************//
	//**************************************************************//
	
	//Textures temporaires
	TimeTexture<float> tex_result(1,size);
	init_texture_single(tex_result);
	TimeTexture<float> texture(1,size);
	init_texture_single(texture);
	float val=0;
	
	//identification des composantes connexes
//         cout<<"EConstraints: Connected parts Id"<<endl;
	int ct=0;
	for(int i=0; i<size; i++)
	{
		if( (origine[0].item(i) != 0)  && (texture[0].item(i)==0) )
		{
			val++;
// 			cout<<"Valeur="<<val<<" avec texture="<<origine[0].item(i)<<endl;
			recurs_proc(i, origine, tex_result, neigh, val);
			for(int i=0; i<size; i++)
			{
				if(tex_result[0].item(i)!=0)
					texture[0].item(i)=tex_result[0].item(i);
			}
		}
	}
// 	Writer<Texture1d> ws2res("result.tex");
// 	ws2res.write(tex_result);
// 	Writer<Texture1d> ws2resor("texture.tex");
// 	ws2resor.write(texture);

	//On traite chaque morceau de la contrainte
	float value;
	for( value=1; value<=val; value++)
	{
		std::vector< int > vect_temp;
		std::vector< int > pairs_temp;
                cout<<"Value="<<value<<endl;
		TimeTexture<float> tex_ext(1,size);
		init_texture_single(tex_ext);
		TimeTexture<float> tex_temp(1,size);
		init_texture_single(tex_temp);
	
		int compteur=0;
		//On enl�e les points non susceptibles d'etre des extremites
                cout<<"Removing non Xt points"<<endl;
		for(int j=0; j<size;j++)
		{
			int cpt=0;
			if(tex_result[0].item(j)==value)
			{
				vect_temp.push_back(j);
				int mem_vois=0;
				std::set<uint>::const_iterator itvoisin;
				itvoisin=neigh[j].begin();
				for(; itvoisin!=neigh[j].end(); itvoisin++)
				{
					if(tex_result[0].item(*itvoisin)==value)
					{
						cpt++;
						mem_vois=(*itvoisin);
					}
				}
				if( cpt==1 || cpt==3 || cpt==4 )
				{
					compteur++;
					tex_ext[0].item(j)=value;
				}
				else
					tex_ext[0].item(j)=0;
			}
		}
		int tab[compteur];
		int ind=0;
		for(int j=0; j<size;j++)
		{
			if(tex_ext[0].item(j)==value)
			{
				tab[ind]=j;
				ind++;
			}
		}
		
		
		GraphPath<float> gr;
		
		//On passe un coup de PCC pour reperer les extremites les plus eloignees
		int path_max=0;
		int path_temp=0;
		int ext1=0, ext2=0;
		TimeTexture<float> tex_path(1,size);
		TimeTexture<float> tex_path_temp(1,size);
		init_texture_single(tex_path);
		init_texture_single(tex_path_temp);
		for(int j=0; j<compteur-1; j++)
		{
			for(int k=j+1; k<compteur; k++)
			{
				path_temp=0;
				tex_path=gr.process(tex_result,mesh,value,tab[j],tab[k]);
				for(int x=0; x<size;x++)
				{
					if(tex_path[0].item(x)==value)
						path_temp++;
				}
				if(path_temp>path_max)
				{
					path_max=path_temp;
					ext1=tab[j];
					ext2=tab[k];
				}
			}
		}
                cout<<"longest pcc"<<endl;
		
		//on recupere le plus long des pcc
		if(path_max!=0)
		{
// 			GraphPath<float> gr1;
			tex_path_temp=gr.process(tex_result,mesh,value,ext1,ext2);
		}
		
                cout<<"done longest pcc"<<endl;
		
		if(ext1==ext2)
		{
// 			pairs_temp.push_back(ext1);
// 			pairs_temp.push_back(ext1);
			//temp[0].item(*itvois)=flag;
		}
		else
		{
			
			pairs_temp.push_back(ext1);
			pairs_temp.push_back(ext2);
			std::cout<<"Push "<<ext1<<" et "<<ext2<<std::endl;
			meridian_pairs.push_back(pairs_temp);
		}
		
// 		Writer<Texture1d> ws2y("or.tex");
// 		ws2y.write(tex_path_temp);
		
	}

	//****Rep?rage des extr?mit?s des morceaux****
//          cout<<"reperage des extremites et taille de meridian_pairs:"<<meridian_pairs.size()<<endl;

	for(unsigned i=0; i<meridian_pairs.size(); i++)
	{
		for(unsigned j=0; j<meridian_pairs[i].size(); j++)
		{
			//result_pairs[0].item(meridian_pairs[i][j]) = flag*(i+1)*2;
// 			cout<<"taille de meridian_pairs["<<i<<"]:"<<meridian_pairs[i].size()<<endl;
			result_pairs[0].item(meridian_pairs[i][j]) = 1;
// 			std::cout<<"Paint "<<meridian_pairs[i][j]<<std::endl;
		}
	}
//	result_pairs[0].item(0)=50;
// 	Writer<Texture1d> ws2e("ext.tex");
// 	ws2e.write(result_pairs);
	// 	
// 	cout<<"Top 1"<<endl;
	for(unsigned i=0; i<meridian_parts.size(); i++)
	{
		for(unsigned j=0; j<meridian_parts[i].size(); j++)
		{
			result_link[0].item(meridian_parts[i][j]) = flag;
		}
	}
// 	cout<<"Top 2"<<endl;
/*	Writer<Texture1d> wTp2parts("meridian_parts.tex");
	wTp2parts.write(result_link);*/
	TimeTexture<float> dist_remp(1,size);
	TimeTexture<float> distance_poles(1,size);
	TimeTexture<float> distance_pole_calleux_entier(1,size);
	TimeTexture<float> distance_pole_insula_entier(1,size);
	init_texture_single(distance_pole_calleux_entier);
	init_texture_single(distance_pole_insula_entier);
	
	for(int i=0; i<size; i++)
	{
		if(i==nord)
			dist_remp[0].item(i)=10;
		else
			dist_remp[0].item(i)=0;
	}

// 	cout<<"Top 3"<<endl;

	distance_poles[0]=MeshDistance( mesh[0], dist_remp[0],true);

// 	cout<<"Top 4"<<endl;
	TimeTexture<float> yeah(1,size);
	init_texture_single(yeah);
	int cpty=1;
// 	cout<<"meridian_pairs.size()="<<meridian_pairs.size()<<endl;
/*	for(unsigned i=0; i<meridian_pairs.size(); i++)
	{
		cout<<"meridian_pairs["<<i<<"][0]"<<meridian_pairs[i][0]<<endl;
		cout<<"meridian_pairs["<<i<<"][1]"<<meridian_pairs[i][1]<<endl;
	}*/
	for(unsigned i=0; i<meridian_pairs.size(); i++)
	{
// 		cout<<"iter="<<i<<endl;
		
// 		cout<<"avant Test: [0]="<<distance_poles[0].item( meridian_pairs[i][0])<<" et [1]="<<distance_poles[0].item( meridian_pairs[i][1] )<<endl;
		//L'element [0] doit �re le plus proche du pole sur la carte de distance
		if( distance_poles[0].item( meridian_pairs[i][0])  > distance_poles[0].item( meridian_pairs[i][1] ) )
		{
			int pt= meridian_pairs[i][1];
			cout<<"pt="<<pt<<endl;
			meridian_pairs[i][1]=meridian_pairs[i][0];
			meridian_pairs[i][0]=pt;
		}
		else {}
		
// 		cout<<"meridian_pairs["<<i<<"][0]"<<meridian_pairs[i][0]<<endl;
// 		cout<<"meridian_pairs["<<i<<"][1]"<<meridian_pairs[i][1]<<endl;
		yeah[0].item(meridian_pairs[i][0])=cpty++;
		yeah[0].item(meridian_pairs[i][1])=cpty++; 

	}
// 	Writer<Texture1d> wTp2yeah("yeah.tex");
// 	wTp2yeah.write(yeah);
	
// 
// 	cout<<"Top 5"<<endl;
	std::vector< std::vector< int > > mark = meridian_pairs;
	int tab[mark.size()];
	
	// Classement
	for(unsigned int i=0; i<mark.size(); i++)
	{
		tab[i]=0;
	}


	TimeTexture<float> dist(1,size);
	TimeTexture<float> carte(1,size);
// 	int cp=0;

	//Mise �jour de la carte de distance reference
	
	std::cout<<"Mise a jour carte"<<std::endl;
	for(int k=0; k<size; k++)
		if(k==nord)
		{
			carte[0].item(k)=10;
		}
	else
		carte[0].item(k)=0;
	

// 	cout<<"Top 6"<<endl;
	for(unsigned int i=0; i<mark.size(); i++)
	{
		float top;
		int ind=i;
		dist[0]=MeshDistance( mesh[0], carte[0],true );
		top=dist[0].item(meridian_pairs[0][0]);
		for(unsigned int j=0; j<mark.size(); j++)
		{
			if(tab[j]==0)
			{
				if(meridian_pairs[j][0]>=size) std::cout<<"OVERFLOW!!!!"<<std::endl;
				if( dist[0].item(meridian_pairs[j][0]) <= top)
				{
					for(unsigned h=0; h<meridian_pairs[j].size(); h++)
					{
						mark[i][h] = meridian_pairs[j][h];
						mark[i][h] = meridian_pairs[j][h];
					}
					top=dist[0].item(mark[i][0]);
					ind=j;
				}
				else{}
			}
		}
		tab[ind]=1;
		//Remise �jour de la carte de distance avec une nouvelle branche
		if(i<mark.size()-1)
		{
			std::cout<<"Remise a jour carte"<<std::endl;
			for(int k=0; k<size; k++)
				if(k==meridian_pairs[i+1][0])
				{
					carte[0].item(k)=10;
				}
				else
					carte[0].item(k)=0;
		}
	}
	

// 	cout<<"Top 7"<<endl;
	for(int i=0; i<size; i++)
	{
		finish[0].item(i) = result_link[0].item(i);
	}

	TimeTexture<float> tex_mark(1, size);
	TimeTexture<float> tex_meridian_pairs(1, size);
	init_texture_single(tex_mark);
	init_texture_single(tex_meridian_pairs);
	
	int cm=1;
	int cmp=1;
	for(unsigned int i=0; i<mark.size(); i++)
	{
		for(unsigned int j=0; j<mark[i].size(); j++)
		{
			tex_mark[0].item(mark[i][j])=cm++;
		}
	}

	for(unsigned int i=0; i<meridian_pairs.size(); i++)
	{
		for(unsigned int j=0; j<meridian_pairs[i].size(); j++)
			tex_meridian_pairs[0].item(meridian_pairs[i][j])=cmp++;
	}
	

// 	cout<<"Top 8"<<endl;
// 	Writer<Texture1d> wTptyeahzf("tex_mark.tex");
// 	wTptyeahzf.write(tex_mark);
// 	
// 	Writer<Texture1d> wTptyeahzzt("tex_meridian_pairs.tex");
// 	wTptyeahzzt.write(tex_meridian_pairs);
	
	
	TimeTexture<float> carte_p(1, size);
	TimeTexture<float> dist_p(1, size);
	for(int i=0; i<size; i++)
	{
		if(i==meridian_pairs[mark.size()-1][1])
			carte[0].item(i)=10;
		else
			carte[0].item(i)=0;
		
// 		if(poles[0].item(i)==180)
// 			carte_p[0].item(i)=10;
// 		else
// 			carte_p[0].item(i)=0;
	}
	dist[0]=MeshDistance( mesh[0], carte[0],true);
// 	dist_p[0]=MeshDistance( mesh[0], carte_p[0],true);
	
	TimeTexture<float> dist_temp(mark.size()+10,size);
	//std::cout<<"TAILLE="<<size<<std::endl;
//	TimeTexture<float> dist_temp(1,size);
	int cp=0;

	findNearNeigh(nord, mark[0][0], finish, flag, mesh, neigh);

	for(int i=0; i<size; i++)
	{
		dist_temp[cp].item(i)=finish[0].item(i);
	}
	cp++;
// 	Writer<Texture1d> wTptyeahzz("dist_liaison_nord.tex");
// 	wTptyeahzz.write(finish);
	

// 	cout<<"Top 9"<<endl;
	fflush(stdout);
	
 	for(unsigned int i=0; i<mark.size()-1; i++)
	{
// 		std::cout<<"Boucle lele"<<std::endl;
		findNearNeigh(mark[i][1], mark[i+1][0], finish, flag, mesh, neigh);
		for(int i=0; i<size; i++)
		{
			dist_temp[cp].item(i)=finish[0].item(i);
		}
		cp++;
	}

 /*	Writer<Texture1d> wTptyeahee("dist_liaison_suite.tex");
 	wTptyeahee.write(finish);*/
	

// 	cout<<"Top 10"<<endl;
	findNearNeigh(mark[mark.size()-1][1], sud, finish, flag, mesh, neigh);
	for(int i=0; i<size; i++)
	{
		dist_temp[cp].item(i)=finish[0].item(i);
	}
	

// 	cout<<"Top 11"<<endl;

//  	Writer<Texture1d> wTptyerahee("construction_mer_origin.tex");
//  	wTptyerahee.write(dist_temp);
}


std::vector<unsigned> aims::buildOriginVector(TimeTexture<float> & tex, int nord, int sud, std::vector<std::set<uint> > neigh)
{
	int size=tex[0].nItem();
	std::vector<unsigned> meridian_nodes;
	TimeTexture<float> marquage(1, size);
	
	std::set<uint> set_vois=neigh[nord];
	std::set<uint>::const_iterator itvois;
	
	init_texture_single(marquage);
	
	meridian_nodes.push_back(nord);
	int courant=nord;
	marquage[0].item(nord)=1;  //noeud du meridien marque
	/*******FABRICATION DU VECTEUR DE NOEUDS DU MERIDIEN D'ORIGINE***********/
	/****INUTILE POUR L'INSTANT*******/
// 	std::cout<<"Apres liaison7!"<<std::endl;

	while(courant!=sud)
	{
// 		std::cout<<"while avec courant="<<courant<<std::endl;
		set_vois=neigh[courant];
		itvois=set_vois.begin();
		int max;
		for ( ; itvois != set_vois.end(); itvois++ )
		{
			max=0;
			if( (marquage[0].item(*itvois)!=1) && (tex[0].item(*itvois)==360) )
			{
				int compteur=0;
				std::set<uint>::const_iterator ittest;
				ittest=neigh[*itvois].begin();
				for ( ; ittest != neigh[*itvois].end(); ittest++)
					if (tex[0].item(*ittest)==360)
						compteur++;
				if (compteur > max)
					courant=*itvois;
			}
			marquage[0].item(courant)=1;
			meridian_nodes.push_back(courant);
		}
	}
	return meridian_nodes;
}


void aims::findNearNeigh(int origine, int destination, TimeTexture<float> & tex_origine, int flag, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh )
{
	//TimeTexture<float> link(1, size);

	int size=tex_origine[0].nItem();

// 	std::cout<<"origine = "<<origine<<" et destination = "<<destination<<std::endl;
	TimeTexture<float> dist(1,size);
	TimeTexture<float> carte(1,size);
	for(int i=0; i<size; i++)
		if(i==destination)
			carte[0].item(i)=10;
		else
			carte[0].item(i)=0;

	dist[0]=MeshDistance( mesh[0], carte[0],true);

	//for(int i=0; i<size; i++)
		//link[0].item(i)=tex_origine[0].item(i);
		//link[0].item(i)=0;

	std::set<uint>::const_iterator itneigh;
	int result=origine;
	float distm, temp;

	temp=dist[0].item(origine);

	while(result!=destination)
	{
// 		std::cout<<"Noeud "<<result<<std::endl;
		int result_temp=result;
		itneigh=neigh[result].begin();
		for ( ; itneigh != neigh[result].end(); ++itneigh )
		{
			distm=dist[0].item(*itneigh);
			if( ( distm < temp ) )//&& ( tex_origine[0].item(*itneigh)!=(float)flag ) )
			{
				result_temp=(*itneigh);
				//std::cout<<"\tVoisin "<<result <<" plus proche avec dist="<<dist<< " et dist prec.="<<temp<<std::endl;
				temp=distm;
			}
			else {}
		}
		result=result_temp;
// 		std::cout<<"dehors du for"<<std::endl;
		tex_origine[0].item(result)=flag;
	}

}

void aims::findNearNeighPoles(int origine, int destination, TimeTexture<float> & tex_origine, int flag, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh, TimeTexture<float> & pol )
{
	int size=tex_origine[0].nItem();

	TimeTexture<float> dist(1,size);
	TimeTexture<float> carte(1,size);
	for(int i=0; i<size; i++)
		if(i==destination)
			carte[0].item(i)=10;
		else
			carte[0].item(i)=0;

	dist[0]=MeshDistance( mesh[0], carte[0],true);

	std::set<uint>::const_iterator itneigh;
	int result=origine;
	float distm, distp, temp, temp_p;

	temp=dist[0].item(origine);
	temp_p=pol[0].item(origine);

/*		Writer<Texture1d> wTptyre("res_orig.tex");
		wTptyre.write(tex_origine);*/
	while(result!=destination)
	{
		int result_temp=result;
		itneigh=neigh[result].begin();
		for ( ; itneigh != neigh[result].end(); ++itneigh )
		{
			distm=dist[0].item(*itneigh);
			distp=pol[0].item(*itneigh);
			if( distm < temp )
			{
				std::cout<<"\tVoisin "<<result <<" plus proche avec dist="<<distm<< " et dist pol.="<<distp<<"et dist_prec.="<<temp<<std::endl;
				if( distp < temp_p )
				{
					temp_p=distp;
					result_temp=(*itneigh);
					temp=distm;
				}
			}
			else {}
		}
		std::cout<<"dehors du for"<<std::endl;
		
// 		Writer<Texture1d> wTptye("res.tex");
// 		wTptye.write(tex_origine);

		result=result_temp;
		tex_origine[0].item(result)=flag;
	}

}


TimeTexture<float> aims::originNeighbourgs(TimeTexture<float> originMeridian, int nord, int sud, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh , TimeTexture<float> & poles)

{
	int size=originMeridian.nItem();
	TimeTexture<float> bothSides(1,size);
	TimeTexture<float> bothSides_origin(1,size);
	TimeTexture<float> side1(1,size);
	TimeTexture<float> side2(1,size);
	TimeTexture<float> points(1,size);
	TimeTexture<float> init(1,size);
	
	init_texture_single(bothSides);
	init_texture_single(bothSides_origin);
	init_texture_single(side1);
	init_texture_single(side2);
	init_texture_single(points);
	init_texture_single(init);

	
	//redefinie des poles temporaires
	int nord_temp;
	int sud_temp;
	int cptZero;
	int cptPole1;
	int cptPole2;
	
	nord_temp=nord;
	sud_temp=sud;
	
	std::set<uint>::const_iterator itneigh;

     std::cerr << "entering origin neighbours" << std::endl;
	
	for(int i=0; i<size; i++)
	{
		
		if( (originMeridian[0].item(i)!=0) && (poles[0].item(i)!=0) )
		{
			cptZero=0;
			cptPole1=0;
			cptPole2=0;
			itneigh=neigh[i].begin();
			for ( ; itneigh != neigh[i].end(); ++itneigh )
			{
				if( poles[0].item(*itneigh)==0 )
					cptZero++;
				if( poles[0].item(*itneigh)==1 )
					cptPole1++;
				if( poles[0].item(*itneigh)==180 )
					cptPole2++;
			}
			//determination des poles
			//EN COMMENTAIRE POUR RE-TESTER L'ANCIENNE VERSION!!
/*			if( cptZero!=0 && cptPole1!=0)
				nord_temp=i;
			if( cptZero!=0 && cptPole2!=0)
				sud_temp=i;*/
		}
	}
	

	points[0].item(nord_temp)=10;
	points[0].item(sud_temp)=20;
	
	for(int i=0; i<size; i++)
	{
		if( (originMeridian[0].item(i)!=0) && (poles[0].item(i)==0) )
			init[0].item(i)=originMeridian[0].item(i);
	}
		
	itneigh=neigh[nord_temp].begin();
	for ( ; itneigh != neigh[nord_temp].end(); ++itneigh )
	{
		init[0].item(*itneigh)=originMeridian[0].item(*itneigh);
	}
	
	itneigh=neigh[sud_temp].begin();
	for ( ; itneigh != neigh[sud_temp].end(); ++itneigh )
	{
		init[0].item(*itneigh)=originMeridian[0].item(*itneigh);
	}

	init[0].item(nord_temp)=originMeridian[0].item(nord_temp);
	init[0].item(sud_temp)=originMeridian[0].item(sud_temp);
	
// 	Writer<Texture1d> wsp("points.tex");
// 	wsp.write(points);
	
// 	Writer<Texture1d> wsi("init.tex");
// 	wsi.write(init);
	//Cree une texture avec les 2 cotes du meridien d'origine
	
	//dans cette boucle, init remplace originMeridian dans la version modifiee.
	for(int i=0; i<size; i++)
	{
		itneigh=neigh[i].begin();
		if( (originMeridian[0].item(i)!=0) && (i!=nord_temp) && (i!=sud_temp) )
		{
			for ( ; itneigh != neigh[i].end(); ++itneigh )
			{
				if( originMeridian[0].item(*itneigh)==0 && (*itneigh)!=(uint)nord_temp && (*itneigh)!=(uint)sud_temp )
				{
					bothSides[0].item(*itneigh)=1;
					bothSides_origin[0].item(*itneigh)=1;
				}
			}
		}
	}
	
	bothSides[0].item(nord_temp)=1;
	bothSides[0].item(sud_temp)=1;
	
// 	Writer<Texture1d> wsb("bothSides.tex");
// 	wsb.write(bothSides);
	
	GraphPath<float> gr;

     // ATTENTION : CECI EST UNE CORRECTION DE BUG  (OLIVIER)
     // LE PLUS COURT CHEMIN VIRE DES TRIANGLES QUI DEVRAIENT RESTER
     // LES CONSEQUENCES SUR LES COORDONNEES SONT ENORMES
     // SPARADRAP DE REPARATION : debug Sides et boucle en fin de fonction
     TimeTexture<float> debugSides(bothSides);
	
	//repere le premier cote par un pcc
     std::cerr << "\t origin neighbours : 1er graphe" << std::endl;

	side1=gr.process(bothSides, mesh, 1, nord_temp, sud_temp);
	/*
	for(int i=0; i<size; i++)
		if( bothSides[0].item(i)==1 && side1[0].item(i)==0 )
			side2[0].item(i)=1;
	
	for(int i=0; i<size; i++)
		if( bothSides[0].item(i)==1 && side2[0].item(i)==0 )
			side1[0].item(i)=1;
		
*/	

// 	Writer<Texture1d> ws1("side1.tex");
// 	ws1.write(side1);
	
	for(int i=0; i<size; i++)
		if( bothSides[0].item(i)==1 && side1[0].item(i)!=0 )
			bothSides[0].item(i)=0;

	bothSides[0].item(nord_temp)=1;
	bothSides[0].item(sud_temp)=1;
	
// 	Writer<Texture1d> wsI("sideInter.tex");
// 	wsI.write(bothSides);
     std::cerr << "\t origin neighbours : 2eme graphe" << std::endl;

	side2=gr.process(bothSides, mesh, 1, nord_temp, sud_temp);
	
// 	Writer<Texture1d> ws2("side2.tex");
// 	ws2.write(side2);
	
	for(int i=0; i<size; i++)
		if( side1[0].item(i)!=0 )
			bothSides[0].item(i)=4;
		else
			if(side2[0].item(i)!=0)
				bothSides[0].item(i)=2;
			else
				bothSides[0].item(i)=0;
	
	for(int i=0; i<size; i++)
	{
		if( bothSides_origin[0].item(i)==1 && bothSides[0].item(i)!=4 && bothSides[0].item(i)!=2 )
		{
			itneigh=neigh[i].begin();
			float ind=0;
			for ( ; itneigh != neigh[i].end(); ++itneigh )
			{
				if( bothSides[0].item(*itneigh)==2)
					ind=2;
				if( bothSides[0].item(*itneigh)==4)
					ind=4;
			}
			if( ind!=0 )
				bothSides[0].item(i)=ind;
		}
	}

     // LA BOUCLE DE CORRECTION DE BUGS

     for (int i=0; i<size; i++)
     {
          if ((debugSides[0].item(i) != 0) && (bothSides[0].item(i) == 0))
          {
               for (itneigh=neigh[i].begin(); itneigh!=neigh[i].end(); ++itneigh)
               {
                    if (bothSides[0].item(*itneigh) != 0)
                    {
                         bothSides[0].item(i) = bothSides[0].item(*itneigh);
                    }
               }
          }
     }
	
     std::cerr << "\t origin neighbours : sortie" << std::endl;

	return bothSides;
	
}

TimeTexture<float> aims::defineSides( TimeTexture<float> & sides, TimeTexture<float> & constraints, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh )
{
	//bool stop=false;
	int depart=0;
	int arrivee=0;
	int size=sides.nItem();
	TimeTexture<float> p( 1,sides.nItem() );
	init_texture_single(p);
	
	int changeTab[11];
	int departInd=0;
	int arriveeInd=0;
	
	for(int x=0;x<11;x++)
	{
		init_texture_single(p);
		for(int i=departInd+1; i<size; i++)
			if( sides[0].item(i)==2 )
			{
				depart=i;
				departInd=i;
				p[0].item(i)=10;
// 				cout<<"DepartInd="<<departInd;
				break;
			}
		
		for(int i=arriveeInd; i<size; i++)
			if( constraints[0].item(i)>10 && constraints[0].item(i)<80 )  //constraints[0].item(i)==16 || constraints[0].item(i)==43 || constraints[0].item(i)==58 )
			{
				arrivee=i;
				p[0].item(i)=10;
				break;
			}
			
		findNearNeigh(depart, arrivee, p, 10, mesh, neigh);
	
		changeTab[x]=0;
		for(int i=0; i<size; i++)
			if( sides[0].item(i)==4 && p[0].item(i)==10 )
				changeTab[x]=1;
			
		
	}
	
	
	bool change=false;
	int sum=0;
	for(int i=0;i<11;i++)
	{
		sum+=changeTab[i];
// 		cout<<"tab["<<i<<"] = "<<changeTab[i]<<"\t";
	}
	if(sum<6)
		change=true;
// 	cout<<"sum="<<sum<<endl;
// 	cout<<"change="<<change<<endl;
		
/*	for(int i=0; i<size; i++)
		if( sides[0].item(i)==4 && p[0].item(i)==10 )
		{
			change=true;
			break;
		}*/
	
	//std::cout<<"CHANGE="<<change<<std::endl;
	if(change==true)
	{
// 		cout<<"change=true alors on change true haha"<<endl;
		for(int i=0; i<size; i++)
			if(sides[0].item(i)==4)
				sides[0].item(i)=2;
			else
				if(sides[0].item(i)==2)
					sides[0].item(i)=4;
				else{}
	}
	
	//return sides;
// 	Writer<Texture1d> wT1c("sides.tex");
// 	wT1c.write(sides);
	return sides;
	
}



TimeTexture<float> aims::defineSidesPoles( TimeTexture<float> & sides, TimeTexture<float> & constraints, AimsSurfaceTriangle mesh, std::vector<std::set<uint> > neigh )
{
	//bool stop=false;
	int depart=0;
	int arrivee=0;
	int size=sides.nItem();
	TimeTexture<float> p( 1,sides.nItem() );
	init_texture_single(p);
	
	int changeTab[11];
	int departInd=0;
	int arriveeInd=0;
	
	for(int x=0;x<11;x++)
	{
		init_texture_single(p);
		for(int i=departInd+1; i<size; i++)
			if( sides[0].item(i)==2 )
		{
			depart=i;
			departInd=i;
			p[0].item(i)=10;
// 				cout<<"DepartInd="<<departInd;
			break;
		}
		
		for(int i=arriveeInd; i<size; i++)
			if( constraints[0].item(i)>35 && constraints[0].item(i)<65 )
		{
			arrivee=i;
			p[0].item(i)=10;
			break;
		}
			
		findNearNeigh(depart, arrivee, p, 10, mesh, neigh);
	
		changeTab[x]=0;
		for(int i=0; i<size; i++)
			if( sides[0].item(i)==4 && p[0].item(i)==10 )
				changeTab[x]=1;
			
		
	}
	
	
	bool change=false;
	int sum=0;
	for(int i=0;i<11;i++)
	{
		sum+=changeTab[i];
// 		cout<<"tab["<<i<<"] = "<<changeTab[i]<<"\t";
	}
	if(sum<6)
		change=true;
// 	cout<<"sum="<<sum<<endl;
// 	cout<<"change="<<change<<endl;
		
/*	for(int i=0; i<size; i++)
	if( sides[0].item(i)==4 && p[0].item(i)==10 )
	{
	change=true;
	break;
}*/
	
	//std::cout<<"CHANGE="<<change<<std::endl;
	if(change==true)
	{
// 		cout<<"change=true alors on change true haha"<<endl;
		for(int i=0; i<size; i++)
			if(sides[0].item(i)==4)
				sides[0].item(i)=2;
		else
			if(sides[0].item(i)==2)
				sides[0].item(i)=4;
		else{}
	}
	
	//return sides;
// 	Writer<Texture1d> wT1c("sides.tex");
// 	wT1c.write(sides);
	return sides;
	
}



