/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *  	CEA/DSV/SHFJ
 *  	4 place du General Leclerc
 *  	91401 Orsay cedex
 *  	France
 */


#include <aims/mesh/texture.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>

#include <aims/mesh/surfaceOperation.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshvoronoi.h>
#include <aims/scalespace/meshDiffuse.h>
#include <aims/distancemap/meshdistance_d.h>
#include <aims/distancemap/meshmorphomat_d.h>

#include <aims/surfacereferential/corticalTools.h>
#include <aims/surfacereferential/corticalConstraints.h>
#include <aims/surfacereferential/sulcusCleaner.h>

namespace aims
{

	void SulcusCleaner::processConstraints()
	{
		preprocess();
		constraintCleaning(context);  //0= longitude, 1=latitude, 2=both
	}
	
	
	void SulcusCleaner::preprocess()
	{
	//Preprocessing 
	
		std::cout<<"lecture"<<std::endl;
		Reader < TimeTexture<short> > rtp(adr_par);
		rtp.read( constraint_lat_read_short );
	
		Reader < TimeTexture<short> > rtm(adr_mer);
		rtm.read( constraint_long_read_short );
	
		Reader < TimeTexture<float> > rpc(adr_calleux);
		rpc.read( pole_call );
	
		std::cout<<"transfert"<<std::endl;
		std::cout<<"size="<<size<<std::endl;
		for(uint i=0;i<size;i++)
		{
			constraint_lat_read[0].push_back( (float)(constraint_lat_read_short[0].item(i)) );
			constraint_long_read[0].push_back( (float)(constraint_long_read_short[0].item(i)) );
		}
	
		corres_map=createCorrespMap( adr_cor, adr_file, side);
		
		std::map< int, std::map<int, std::string> >::iterator itMap;
		for(int i=0;i<size; i++)
		{
			itMap = corres_map.begin();
			for ( ; itMap != corres_map.end(); ++itMap)
			{
				std::map<int, std::string>::iterator itMap2;
				itMap2 = (*itMap).second.begin();
				for ( ; itMap2 != (*itMap).second.end(); ++itMap2)
				{
					if(strcmp( ( (*itMap2).second ).c_str(), ( "INSULA" +side).c_str() )==0 )
					{
						value_pole1=(*itMap2).first;
// 						std::cout<<"Pole Insula = "<<value_pole1<<std::endl;
					}
					if(strcmp( ( (*itMap2).second ).c_str(), ( "S.Call." +side).c_str() )==0 )
					{
						value_pole2=(*itMap2).first;
// 						std::cout<<"Corps Calleux = "<<value_pole2<<std::endl;
					}
				}
			}
		}
		
	//Rescale des textures : assigne des coordonn�s aux sillons projet�
		std::cout<<"preprocessing datas..........";
		fflush(stdout);
	
	//Fermeture
		TimeTexture<short> pole1( 1, size );
		TimeTexture<short> pole1_dilate( 1, size );
		TimeTexture<short> pole1_erod( 1, size );
		TimeTexture<float> pole1_erod_float( 1, size );
		TimeTexture<float> pole1_float( 1, size );
	
		TimeTexture<short> pole2( 1, size );
		TimeTexture<short> pole2_dilate( 1, size );
		TimeTexture<short> pole2_erod( 1, size );
		TimeTexture<float> pole2_erod_float( 1, size );
		TimeTexture<float> pole2_float( 1, size );
	
		for(uint i=0; i<size; i++)
		{
// 			if ( (constraint_lat_read[0].item(i)==23) || (constraint_lat_read[0].item(i)==24))
			if (constraint_lat_read[0].item(i)==value_pole1)
			{
				pole1[0].item(i)=100;
				constraint_lat_read[0].item(i)=0;
			}
			else
				pole1[0].item(i)=0;
			
// 			if ( (pole_call[0].item(i)!=0) || (constraint_lat_read[0].item(i)==31) || (constraint_lat_read[0].item(i)==32) )
			if ( (pole_call[0].item(i)!=0) || (constraint_lat_read[0].item(i)==value_pole2) )
			{
				pole2[0].item(i)=100;
				constraint_lat_read[0].item(i)=0;
			}
			else
				pole2[0].item(i)=0;
		}
		std::cout<<"dilat/erod start"<<std::endl;

	//Fermeture des poles
		pole1_dilate[0]=MeshDilation<short>( mesh[0], pole1[0], short(0), -1, 15.0, true);
		pole1_erod[0]=MeshErosion<short> ( mesh[0], pole1_dilate[0], short(0), -1, 15.0, true);
	
		pole2_dilate[0]=MeshDilation<short>( mesh[0], pole2[0], short(0), -1, 15.0, true);
		pole2_erod[0]=MeshErosion<short> ( mesh[0], pole2_dilate[0], short(0), -1, 15.0, true);
	
	//puis ouverture des poles
	
		pole1_dilate[0]=MeshErosion<short>( mesh[0], pole1_erod[0], short(0), -1, 5.0, true);
		pole1_erod[0]=MeshDilation<short> ( mesh[0], pole1_dilate[0], short(0), -1, 5.0, true);
	
		pole2_dilate[0]=MeshErosion<short>( mesh[0], pole2_erod[0], short(0), -1, 5.0, true);
		pole2_erod[0]=MeshDilation<short> ( mesh[0], pole2_dilate[0], short(0), -1, 5.0, true);
		
		std::cout<<"dilat/erod OK"<<std::endl;
	
	//Enleve les parties ou les 2 poles se touchent
        
// 		std::cout<<"iter:"<<std::endl;
		for(uint i=0; i<size; i++)
		{
// 			std::cout<<" "<<i;
			if( pole1_erod[0].item(i)!=0 && pole2_erod[0].item(i)!=0 )
			{
// 				std::cout<<std::endl;
// 				std::cout<<"Recouvrement..."<<std::endl;
				pole1_erod[0].item(i)=0;
				pole2_erod[0].item(i)=0;
				std::set<uint>::const_iterator itvois;
			
// 				std::cout<<"i="<<i<<std::endl;
				itvois=neigh[i].begin();
// 				std::cout<<"premier voisin = "<<(*itvois)<<std::endl;
				for(; itvois!=neigh[i].end(); ++itvois)
				{
// 					std::cout<<"voisin="<<(*itvois)<<std::endl;
					pole1_erod[0].item(*itvois)=0;
					pole2_erod[0].item(*itvois)=0;
				}
// 				std::cout<<"Recouvrement ok"<<std::endl;
			}
		}
        
		std::cout<<"reconstruct OK"<<std::endl;
// 		std::cout<<"poles OK"<<std::endl;
	//convert to float textures
		for(uint i=0; i<size; i++)
		{
			if(pole1_erod[0].item(i)!=0)
				pole1_erod_float[0].item(i)=100;
			else
				pole1_erod_float[0].item(i)=0;
			
			if(pole2_erod[0].item(i)!=0)
				pole2_erod_float[0].item(i)=100;
			else
				pole2_erod_float[0].item(i)=0;
		}
	
		for(uint i=0; i<size; i++)
		{
			if( pole2_erod_float[0].item(i)!=0 ) 
				poles[0].push_back(1);
			else
				if( pole1_erod_float[0].item(i)!=0 )
					poles[0].push_back(180);
			else
				poles[0].push_back(0);
		}
		
		std::cout<<"writing texture..."<<std::endl;
		Writer<Texture1d> wT1cclfr(adr_poles);
		wT1cclfr.write(poles);
		
		std::cout<<"done"<<std::endl;
	}
	
	
	
	void SulcusCleaner::constraintCleaning( int ct )
	{
		
		TimeTexture<float> contraint_tmp(1,size);
		std::map<int, int> value_sc;
		std::map< int, std::map<int, std::string> >::iterator itMap;
		for(int i=0;i<size; i++)
		{
			//repere la valeur du SC et le degage du traitement plus tard
			itMap = corres_map.begin();
			for ( ; itMap != corres_map.end(); ++itMap)
			{
				
				std::map<int, std::string>::iterator itMap2;
				itMap2 = (*itMap).second.begin();
				for ( ; itMap2 != (*itMap).second.end(); ++itMap2)
				{
					
					if(strcmp( ( (*itMap2).second ).c_str(), ( "S.C."+side ).c_str() )==0 )
					{
						value_sc[ (*itMap2).first ] = 1;
	// 					std::cout<<"SC!!"<<std::endl;
					}
					if(strcmp( ( (*itMap2).second ).c_str(), ( "S.C.inf."+side ).c_str() )==0 )
					{
						value_sc[ (*itMap2).first ] = 1;
	// 					std::cout<<"SCinf!!"<<std::endl;
					}
					if(strcmp( ( (*itMap2).second ).c_str(), ( "S.C.sup."+side ).c_str() )==0 )
					{
						value_sc[ (*itMap2).first ] = 1;
	// 					std::cout<<"SCsup!!"<<std::endl;
					}
				}
			}
		}
		
		std::cout<<std::endl;
		
		if(ct==0)
		{
			std::cout<<"LONGITUDE"<<std::endl;
			
			TimeTexture<float> constraint_long_cleaned(1,size);
			
			for(int i=0;i<size; i++)
			{
				constraint_long_cleaned[0].item(i)=0;
			}
			
/*			for(int i=0;i<size; i++)
			{
				contraint_tmp[0].item(i)=constraint_long_read[0].item(i);
				std::map<int, int>::iterator it;
				it=value_sc.begin();
				for( ; it!=value_sc.end(); ++it)
				{
					if( constraint_lat_read[0].item(i)==(*it).first )
						contraint_tmp[0].item(i)=0;
				}
			}*/
			
			std::cout<<"Cleaning constraints.........";
			fflush(stdout);
			
			constraint_long_cleaned=constraintCleaner(constraint_long_read, neigh, mesh, contr, curvature, elasticity);
			std::cout<<"done Cleaned"<<std::endl;
			
			constraint_long_scaled=rescaleConstraints(constraint_long_cleaned, corres_map);
		
			
			Writer<Texture1d> wT1cclu(adr_long_cleaned);
			wT1cclu.write(constraint_long_scaled);
		}
		
		if(ct==1)
		{
			std::cout<<"LATITUDE"<<std::endl;
			
			TimeTexture<float> constraint_lat_cleaned(1,size);
			TimeTexture<float> contraint_tmp(1,size);
			
			for(int i=0;i<size; i++)
			{
				contraint_tmp[0].item(i)=constraint_lat_read[0].item(i);
				constraint_lat_cleaned[0].item(i)=0;
			
				if( (constraint_lat_read[0].item(i) == value_pole1) || (constraint_lat_read[0].item(i) == value_pole2) )
					contraint_tmp[0].item(i)=0;
			}
			
			std::cout<<"Cleaning constraints.........";
			fflush(stdout);
			
			constraint_lat_cleaned=constraintCleaner(contraint_tmp, neigh, mesh, contr, curvature, elasticity);
			std::cout<<"done Cleaned"<<std::endl;
			
			constraint_lat_scaled=rescaleConstraints(constraint_lat_cleaned, corres_map);
			
			for(int i=0;i<size; i++)
			{
				if( (constraint_lat_read[0].item(i) == value_pole1) || (constraint_lat_read[0].item(i) == value_pole2) )
					constraint_lat_cleaned[0].item(i)=poles[0].item(i);
			}
			
			Writer<Texture1d> wT1cclu(adr_lat_cleaned);
			wT1cclu.write(constraint_lat_scaled);
		}
		
		if(ct==2)
		{
			std::cout<<"LONGITUDE"<<std::endl;
			
			TimeTexture<float> constraint_long_cleaned(1,size);
			
			for(uint i=0; i<size; i++)
			{
				constraint_long_cleaned[0].item(i)=0;
			}
			
			std::cout<<"Cleaning constraints.........";
			fflush(stdout);
			
			constraint_long_cleaned=constraintCleaner(constraint_long_read, neigh, mesh, contr, curvature, elasticity);
			std::cout<<"done Cleaned"<<std::endl;
			
			constraint_long_scaled=rescaleConstraints(constraint_long_cleaned, corres_map);
			
			Writer<Texture1d> wlong(adr_long_cleaned);
			wlong.write(constraint_long_scaled);
			
			/////////////////////////////////////////////////////////////////////
			
			std::cout<<"LATITUDE"<<std::endl;
			
			TimeTexture<float> constraint_lat_cleaned(1,size);
			TimeTexture<float> contraint_tmp(1,size);
			
			for(int i=0;i<size; i++)
			{
				contraint_tmp[0].item(i)=constraint_lat_read[0].item(i);
				constraint_lat_cleaned[0].item(i)=0;
			
				if( (constraint_lat_read[0].item(i) == value_pole1) || (constraint_lat_read[0].item(i) == value_pole2) )
					contraint_tmp[0].item(i)=0;
			}
			
			std::cout<<"Cleaning constraints.........";
			fflush(stdout);
			
			constraint_lat_cleaned=constraintCleaner(contraint_tmp, neigh, mesh, contr, curvature, elasticity);
			std::cout<<"done Cleaned"<<std::endl;
			
			constraint_lat_scaled=rescaleConstraints(constraint_lat_cleaned, corres_map);
			
			for(int i=0;i<size; i++)
			{
				if( (constraint_lat_read[0].item(i) == value_pole1) || (constraint_lat_read[0].item(i) == value_pole2) )
					constraint_lat_cleaned[0].item(i)=poles[0].item(i);
			}
			
			Writer<Texture1d> wlat(adr_lat_cleaned);
			wlat.write(constraint_lat_scaled);
		}
	}


}

