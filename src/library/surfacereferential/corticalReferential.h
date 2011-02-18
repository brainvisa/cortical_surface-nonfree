/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *  	CEA/DSV/SHFJ
 *  	4 place du General Leclerc
 *  	91401 Orsay cedex
 *  	France
 */


/**************** BRAND NEW DEV VERSION YEAH!!!!****************/


#ifndef AIMS_CORTICALREFERENTIAL_H
#define AIMS_CORTICALREFERENTIAL_H


#include <aims/mesh/texture.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <cortical_surface/surfacereferential/corticalTools.h>
#include <cortical_surface/surfacereferential/corticalConstraints.h>
#include <cortical_surface/surfacereferential/sulcusCleaner.h>

namespace aims
{

	class CorticalReferential
	{
	
		public:
	
			AimsSurfaceTriangle mesh;
			AimsSurface<3, Void> mesh_base;
			std::map<unsigned, std::set< std::pair<unsigned,float> > >  weightLapl;
		//Liste de voisinage
			std::vector<std::set<uint> > neigh;
		//Surface vector
			std::vector<Point3df> vert;
	
		//Poles Index and size (=vertex numbers)
			int size, nord, sud, ind_min;
			std::pair<int,int> poles_points;
	
		//Parametres pour l'execution
			std::string adr;
			std::string adr_par;
			std::string adr_mer;
			std::string adr_call;
			std::string adr_poles;
		
		//Condition d'arret des iterations (difference avec la moy du laplacien 1000 iter. avant
			float criterium;
	
		//Reglage du pas pour le calcul iteratif
			float _dt;
		
		//Si on doit faire les cotes du sillon central inverses
			int context;
		
		//Ce que l'on effectue
			int choice_process;			
	
		//Reglage du terme d'attache aux donn�s
			float _Beta;
                
		//Type d'attache aux données
			int typeBeta;
		
        //Adresse des outputs
			std::string adr_lat;
			std::string adr_lon;
	
		//Textures de contraintes
		
			TimeTexture<float> constraint_lat;
			TimeTexture<float> constraint_lat_cleaned;
			TimeTexture<float> constraint_lat_read;
			TimeTexture<float> constraint_long;
			TimeTexture<float> constraint_long_cleaned;
			TimeTexture<float> constraint_long_read;
			TimeTexture<short> cercle_polaire;
	
			TimeTexture<float> pole_call;
			TimeTexture<float> pole_insula;
			TimeTexture<float> poles;
		
			std::vector<unsigned> forb_list;
	
		// indique si on a deja trouve les poles
			bool poles_found;
			
		// indique si on doit parametriser l'insula ou non
			bool doInsulaParameterization;
	
		//A voir pour le reste
			TimeTexture<float> diff_meridian_origine;
			TimeTexture<float> distance_poles;
			float diametre;
	
		//Constructor
			CorticalReferential ( ) {}
			CorticalReferential ( std::string & adr_param, 
			                      std::string & adr_parallele, 
			                      std::string & adr_meridien, 
			                      std::string & adr_calleux, 
			                      std::string & _adr_poles, 
			                      float criter, 
			                      float dt, 
			                      int c, 
			                      int choice, 
			                      float Beta, 
			                      int tBeta, 
			                      std::string & adrlat, 
			                      std::string & adrlon,
			                      bool _doInsulaParam = true ) : 
			                          adr(adr_param), 
			                          adr_par(adr_parallele), 
			                          adr_mer(adr_meridien), 
			                          adr_call(adr_calleux), 
			                          adr_poles(_adr_poles), 
			                          criterium(criter), 
			                          _dt(dt), 
			                          context(c), 
			                          choice_process(choice), 
			                          _Beta(Beta), 
			                          typeBeta(tBeta), 
			                          adr_lat(adrlat), 
			                          adr_lon(adrlon),
			                          doInsulaParameterization(_doInsulaParam)
			{
			    
				std::cout << "constructeur" << std::endl;
				std::cout << "mesh adress : " << adr << std::endl;
				std::cout << "adr_par : " << adr_par << std::endl;
				std::cout << "adr_mer : " << adr_mer << std::endl;
				std::cout << "adr_call : " << adr_call << std::endl;
				std::cout << "adr_poles : " << adr_poles << std::endl;
				std::cout << "_dt : "<< _dt  << std::endl;
				std::cout << "criterium : " << criterium  << std::endl;
				std::cout << "context : " << context  << std::endl;
				std::cout << "choice : " << choice_process  << std::endl;
				std::cout << "Beta : " << _Beta  << std::endl;
	
			//Opening brain mesh
				Reader < AimsSurfaceTriangle > r(adr);
				r.read( mesh );
				std::cout << "reader mesh ok" << std::endl;
			
			//Mesh version AimsSurface
				mesh_base = mesh[0];
				std::cout << "mesh_base ok" << std::endl;
			//Surface vector
				vert = mesh_base.vertex();
				std::cout << "vert ok!" << std::endl;
	
			//Computing neighbourhood
                neigh = SurfaceManip::surfaceNeighbours ( mesh );
	
				std::cout << "neigh ok" << std::endl;
				weightLapl = AimsMeshWeightFiniteElementLaplacian ( mesh_base, 0.98 );
	
				std::cout << "weightlaplacian ok" << std::endl;
				size = vert.size();
				std::cout << "size :" << size << std::endl;
	
			//LECTURE DES POLES (A INTEGRER PAR LA SUITE AUX CONTRAINTES)
	
				Reader < TimeTexture<float> > lp1 ( adr_call );
				lp1.read( pole_call );
	
				Reader < TimeTexture<float> > lp2 ( adr_poles );
				lp2.read( poles );
	
				poles_found = false;
				std::cout << "Constructeur fini!!" << std::endl;
			}
			
	
		//Main methods
			void process();
			void constraintPreprocess();
			void latitudePropagation();
			void longitudePropagation();
		
		
		//Diffusion methods
			void cingularPoint();
			TimeTexture<float> diffusionLatitudeRelax ( TimeTexture<float> & tex );
			TimeTexture<float> diffusionLatitude ( TimeTexture<float> &);
			TimeTexture<float> diffusionLongitudeRelax ( TimeTexture<float> & tex, TimeTexture<float> & side);
			TimeTexture<float> diffusionLongitude ( TimeTexture<float> & tex, TimeTexture<float> & side, TimeTexture<float> & poleSave, int ind);
			Texture<float> AimsMeshLaplacian_meridian ( const Texture<float> &smooth, TimeTexture<float> &source, const std::map<unsigned, std::set< std::pair<unsigned,float> > > &lapl, TimeTexture<float> & side);
		
	};



}

#endif


