/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *  	CEA/DSV/SHFJ
 *  	4 place du General Leclerc
 *  	91401 Orsay cedex
 *  	France
 */



#ifndef AIMS_SULCUS_CLEANER_H
#define AIMS_SULCUS_CLEANER_H


#include <aims/mesh/texture.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <cortical_surface/surfacereferential/corticalTools.h>
#include <cortical_surface/surfacereferential/corticalConstraints.h>

namespace aims
{

	class SulcusCleaner
	{
	
		public:
	
		//Parametres pour l'execution
			std::string adr_mesh;
			std::string adr_calleux;
			std::string adr_poles;
			std::string adr_mer;
			std::string adr_par;
			std::string adr_long_cleaned;
			std::string adr_lat_cleaned;
			std::string & adr_cor;
			std::string & adr_file;
					
			int context;
			std::string & side;
			
			float contr;
			float curvature;
			float elasticity;
			
			AimsSurfaceTriangle mesh;
		//Liste de voisinage
			std::vector<std::set<uint> > neigh;
	
		//Textures de contraintes
		
			TimeTexture<short> constraint_lat_read_short;
			TimeTexture<short> constraint_long_read_short;
	
			TimeTexture<float> constraint_long_read;
			TimeTexture<float> constraint_lat_read;
	
			TimeTexture<float> pole_call;
			
		//structure avec les textures de distance depuis tous les commets
			
			
			//Map with values for each constraint's name
			std::map< int, std::map<int, std::string> > corres_map;
			
			TimeTexture<float> poles;
			TimeTexture<float> constraint_long_scaled;
			TimeTexture<float> constraint_lat_scaled;
			uint size;
			
			float value_pole1, value_pole2;
	
		//Constructor
			SulcusCleaner( std::string & _adr_mesh, std::string & _adr_calleux, std::string & _adr_poles, std::string & _adr_mer, std::string & _adr_par, std::string & _adr_long_cleaned, std::string & _adr_lat_cleaned, std::string & _adr_cor, std::string & _adr_file, int _context, std::string & _side, float _contr, float _curvature, float _elasticity ) : adr_mesh(_adr_mesh), adr_calleux(_adr_calleux), adr_poles(_adr_poles), adr_mer(_adr_mer), adr_par(_adr_par), adr_long_cleaned(_adr_long_cleaned), adr_lat_cleaned(_adr_lat_cleaned), adr_cor(_adr_cor), adr_file(_adr_file), context(_context), side(_side), contr(_contr), curvature(_curvature), elasticity(_elasticity)
			{
	
			//std::cout<<"SIZE DEBUT="<<size<<std::endl;
	
			//LECTURE DES POLES (A INTEGRER PAR LA SUITE AUX CONTRAINTES)
	
			//Opening brain mesh
				Reader < AimsSurfaceTriangle > r(adr_mesh);
				r.read( mesh );
	
				Reader < TimeTexture<float> > lp1(adr_calleux);
				lp1.read( pole_call );

				size=pole_call[0].nItem();
				neigh = SurfaceManip::surfaceNeighbours( mesh );
			}
	
		//Main methods
			
			void processConstraints();
			void preprocess();
			void constraintCleaning( int ct );
			


	};



}

#endif


