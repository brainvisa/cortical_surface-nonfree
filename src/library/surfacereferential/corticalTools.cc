/*
 *  Copyright (C) 2000-2004 CEA
 *
 *  This software and supporting documentation were developed by
 *  	CEA/DSV/SHFJ
 *  	4 place du General Leclerc
 *  	91401 Orsay cedex
 *  	France
 */


#include <cstdlib>
#include <cortical_surface/surfacereferential/corticalTools.h>
#include <cortical_surface/surfacereferential/corticalConstraints.h>

using namespace std;
/****************TOOLS****************/

//initialization of single textures
void aims::init_texture_single(TimeTexture<float> & texture)
{
	int size=texture[0].nItem();
	for (int i=0;i<size;i++)
		texture[0].item(i)=0;
}

TimeTexture<float> aims::dilate_texture(TimeTexture<float> & texture, float val, std::vector<std::set<uint> > neigh, AimsSurfaceTriangle mesh)
{
	//Val = coef multiplicateur pour la valeur des voisins (1 pour C, 0.5 pour Beta)
	
	int size=texture[0].nItem();
	TimeTexture<float> result(1,size);
	init_texture_single(result);
	//std::cout<<"dilate1"<<std::endl;
	for(int i=0; i<size; i++)
	{
		if(texture[0].item(i)!=0)
		{
			//std::cout<<"test OK"<<std::endl;
			result[0].item(i)=texture[0].item(i);
			std::set<uint>::const_iterator itvoisin;
			itvoisin=neigh[i].begin();
			for(; itvoisin!=neigh[i].end(); itvoisin++)
			{
				//std::cout<<"vois"<<std::endl;
				if(texture[0].item(*itvoisin)==0)
				      result[0].item(*itvoisin)=texture[0].item(i)*val;
			}
		}
	}
	return result;
}

float aims::compute2D_distance(float lon1, float lon2, float lat1, float lat2)
{
	float d, t1, t2, t;
	t1=pow(lon1-lon2,2);
	if(lon1>lon2)
		t2=pow(lon2-lon1-360,2);
	else
		t2=pow(lon1-lon2-360,2);
	if(t1<t2)	t=t1;
	else		t=t2;

	d=sqrt( t + pow( lat1-lat2,2) );

	return d;
}

//Poles points determination (to be improved)
std::pair<int,int> aims::find_poles ( TimeTexture<float> & tex, AimsSurfaceTriangle mesh )
{
	int size = tex[0].nItem();
	std::cout << "SIZE = " << size << std::endl;
	TimeTexture<float> distance_nord ( 1, size );
	TimeTexture<float> distance_sud ( 1, size );
	TimeTexture<float> pole1 ( 1, size );
	TimeTexture<float> pole2 ( 1, size );
	std::pair<int,int> poles_points;
	//Textures used in meshdistance creation
	for ( int i = 0 ; i < size ; i++ ) 
	{
		if ( tex[0].item(i) == 1 ) 
		{
			pole1[0].item(i) = 0;
			pole2[0].item(i) = 360;
		}
		else
			if ( tex[0].item(i) == 180 ) 
			{
				pole1[0].item(i) = 360;
				pole2[0].item(i) = 0;
			}
			else 
			{
				pole1[0].item(i) = 360;
				pole2[0].item(i) = 360;
			}
	}

	//std::cout<<"distance map for pole1 (corpus callosum)"<<std::endl;	
	distance_nord[0] = meshdistance::MeshDistance( mesh[0], pole1[0],true);
	//std::cout<<"distance map for pole1 (insula)"<<std::endl;	
	distance_sud[0] = meshdistance::MeshDistance( mesh[0], pole2[0],true);

// 	Writer<Texture1d> w("/home/go224932/test_pole_sud.tex");
// 	w.write(distance_sud);
// 	Writer<Texture1d> w1("/home/go224932/test_pole_nord.tex");
// 	w1.write(distance_nord);
	
	//CALCUL POLES NORD ET SUD DANS LE POLE CALLEUX ET LE POLE INSULAIRE

	int max_nord = 0, max_sud = 0;
	float max_value_nord = 0, max_value_sud = 0;

	for ( int a = 0; a < size ; a++)
	{
		if ( distance_nord[0].item(a) > max_value_nord )
		{
			max_nord = a;
			max_value_nord = distance_nord[0].item(a);
		}

		if ( distance_sud[0].item(a) > max_value_sud )
		{
			max_sud = a;
			max_value_sud = distance_sud[0].item(a);
		}
	}
	
	std::cout << "nord=" << max_nord << " et sud =" << max_sud << std::endl;
	
	//TEST autre methode de determination du pole nord

	float Xmin = 10000,
	      Zmin = 10000,
	      Xmax = 0,
	      Zmax = 0;
	TimeTexture<float> poleN ( 1, size );
	for ( int i = 0 ; i < size ; i++ )
	{
		float x = 0,
		      z = 0;
		if ( distance_nord[0].item(i) != 0 )
		{
			x = mesh[0].vertex()[i][1];
			z = mesh[0].vertex()[i][2];
			if ( z > Zmax )
				Zmax = z;
			if ( z < Zmin )
				Zmin = z;
			if ( x > Xmax )
				Xmax = x;
			if ( x < Xmin )
				Xmin = x;
		}
	}
	
	float xT, zT, Delta;
	xT = ( Xmin + Xmax ) / 2;
	zT = ( Zmin + Zmax ) / 2;
	Delta = 2;
	for ( int i = 0 ; i < size ; i++ )
	{
		if ( distance_nord[0].item(i) != 0 )
		{
			float xtemp, ztemp;
			xtemp = mesh[0].vertex()[i][1] - xT;
			ztemp = mesh[0].vertex()[i][2] - zT;
			if ( ( sqrt(xtemp*xtemp) <= Delta ) && ( sqrt(ztemp*ztemp) <= Delta ) )
				max_nord=i;
		}
	}
	
	//TEST autre methode de determination du pole sud
	Xmin = 10000;
	Zmin = 10000;
	Xmax = 0;
	Zmax = 0;
	TimeTexture<float> poleS (1, size );
	for ( int i = 0 ; i < size ; i++ )
	{
		float x = 0,
		      z = 0;
		if ( distance_sud[0].item(i) != 0 )
		{
			x = mesh[0].vertex()[i][1];
			z = mesh[0].vertex()[i][2];
			if ( z > Zmax )
				Zmax = z;
			if ( z < Zmin )
				Zmin = z;
			if ( x > Xmax )
				Xmax = x;
			if ( x < Xmin )
				Xmin = x;
		}
	}
	
	xT = ( Xmin + Xmax ) / 2;
	zT = ( Zmin + Zmax ) / 2;
	Delta = 2;
	for ( int i = 0 ; i < size ; i++ )
	{
		if ( distance_sud[0].item(i)!=0 )
		{
			float xtemp, ztemp;
			xtemp = mesh[0].vertex()[i][1] - xT;
			ztemp = mesh[0].vertex()[i][2] - zT;
			if( ( sqrt(xtemp*xtemp) <= Delta ) && ( sqrt(ztemp*ztemp) <= Delta ) )
				max_sud=i;
		}
	}
	
	//Afectation valeurs
	std::cout<<"nord="<<max_nord<<" et sud ="<<max_sud<<std::endl;
	poles_points.first = max_nord;
	poles_points.second = max_sud;

	return poles_points;
}
        
void aims::drawGyri(std::string & adLong, std::string & adLat, std::string & adOut, std::string & adr_corl, std::string & side)
{
	
	TimeTexture<float> longitude, latitude;
	int size;
	
	Reader< TimeTexture<float> > longR ( adLong );
	longR.read(longitude);
	Reader< TimeTexture<float> > latR ( adLat );
	latR.read(latitude);
	
	size = longitude.nItem();
	
	TimeTexture<float> gyri ( 1,size );
	
	std::map< std::string, std::set<float> > map_global;

	const char *adr_cor= adr_corl.c_str();

	FILE *corres; //correspondance between contraints names and real values
	
	if ((corres=fopen(adr_cor, "r")) == NULL)
	{
		cerr << "Cannot open file " << adr_cor << endl;
		exit(EXIT_FAILURE);
	}

	std::string arg1;
		
	int val_contraint;
	int val_projection;
	
	while ( ! feof(corres) )
	{
		char c[40];
		char a[5];
		if( fscanf( corres, "%s %s %i\n", a, c, &val_projection )
                    != 3 )
                  cout << "warning in drawGyri: val_projection not correctly "
                          "read\n";
		arg1 = a;
		map_global[a].insert ( val_projection );
	}
		
	float u = 0,
	      v = 0;
	
	std::set<float>::const_iterator it_lon, it_lat, it_lon2, it_lat2;
		
	for ( int i = 0 ; i < size ; i++ )
	{
		gyri[0].item(i) = 0;
		int cpt = 2;
		
		u = latitude[0].item(i);
		v = longitude[0].item(i);
		
		it_lat = map_global["lat"].begin();
		for ( ; it_lat != map_global["lat"].end() ; ++it_lat )
		{
			it_lat2 = it_lat;
			it_lat2++;
			if ( it_lat2 != map_global["lat"].end() )
			{
				it_lon = map_global["lon"].begin();
				for ( ; it_lon != map_global["lon"].end() ; ++it_lon )
				{
					it_lon2 = it_lon;
					it_lon2++;
					if ( it_lon2 != map_global["lon"].end() )
					{
						if ( u > (*it_lat) && u <= (*it_lat2) && v > (*it_lon) && v <= (*it_lon2) )
						{
							gyri[0].item(i)=cpt;
						}
						cpt++;
					}
				}
			}
		}
		
		//Pole cingulaire
		
		if ( u <= 30 )
			gyri[0].item(i) = 1;
		
		/*
		for(; it_lat!=map_global["lat"].end(); ++it_lat)
		{
			it_lat2=it_lat;
			it_lat2++;
			it_lon=map_global["lon"].begin();
// 				std::cout<<(*it_lat)<<" - "<<(*it_lat2)<<std::endl;
				if( it_lat2!=map_global["lat"].end() )
				{
					for(; it_lon!=map_global["lon"].end(); ++it_lon)
					{
						it_lon2=it_lon;
						it_lon2++;
						if( it_lon2!=map_global["lon"].end() )
						{
							if( u>(*it_lat) && u<=(*it_lat2) && v>(*it_lon) && v<=(*it_lon2) )
							{
								gyri[0].item(i)=cpt;
								cpt++;
							}
						}
					}
				}
			}*/
			
			/*//Pre-central
			if( u<=16 && v>2 && v<= 30)
				gyri[0].item(i)=1;
				
			if( u<=16 && v>30 && v<= 79)
				gyri[0].item(i)=2;
				
			if( u<=16 && v>79 && v<= 98)
				gyri[0].item(i)=3;
				
			if( u<=16 && v>98 && v<= 125)
				gyri[0].item(i)=4;
				
			if( u<=16 && v>125 && v<= 179)
				gyri[0].item(i)=5;
				
				
			//Frontal inferior
			if( u>16 && u<=43 && v>2 && v<= 30)
				gyri[0].item(i)=6;
				
			if( u>16 && u<=43 && v>30 && v<= 79)
				gyri[0].item(i)=7;
				
			if( u>16 && u<=43 && v>79 && v<= 98)
				gyri[0].item(i)=8;
				
			if( u>16 && u<=43 && v>98 && v<= 125)
				gyri[0].item(i)=9;
				
			if( u>16 && u<=43 && v>125 && v<= 179)
				gyri[0].item(i)=10;
			
			
			//bande 3
			if( u>43 && u<=58 && v>2 && v<= 30)
				gyri[0].item(i)=11;
				
			if( u>43 && u<=58 && v>30 && v<= 79)
				gyri[0].item(i)=12;
				
			if( u>43 && u<=58 && v>79 && v<= 98)
				gyri[0].item(i)=13;
				
			if( u>43 && u<=58 && v>98 && v<= 125)
				gyri[0].item(i)=14;
				
			if( u>43 && u<=58 && v>125 && v<= 179)
				gyri[0].item(i)=15;

				
			//bande 4
			if( u>58 && u<=100 && v>2 && v<= 30)
				gyri[0].item(i)=16;
				
			if( u>58 && u<=100 && v>30 && v<= 79)
				gyri[0].item(i)=17;
				
			if( u>58 && u<=100 && v>79 && v<= 98)
				gyri[0].item(i)=18;
				
			if( u>58 && u<=100 && v>98 && v<= 125)
				gyri[0].item(i)=19;
				
			if( u>58 && u<=100 && v>125 && v<= 179)
				gyri[0].item(i)=20;
				
			//bande 5
			if( u>100 && u<=287 && v>2 && v<= 30)
				gyri[0].item(i)=21;
				
			if( u>100 && u<=287 && v>30 && v<= 79)
				gyri[0].item(i)=22;
				
			if( u>100 && u<=287 && v>79 && v<= 98)
				gyri[0].item(i)=23;
				
			if( u>100 && u<=287 && v>98 && v<= 125)
				gyri[0].item(i)=24;
				
			if( u>100 && u<=287 && v>125 && v<= 179)
				gyri[0].item(i)=25;
				
			//bande 6
			if( u>287 && u<=303 && v>2 && v<= 30)
				gyri[0].item(i)=26;
				
			if( u>287 && u<=303 && v>30 && v<= 79)
				gyri[0].item(i)=27;
				
			if( u>287 && u<=303 && v>79 && v<= 98)
				gyri[0].item(i)=28;
				
			if( u>287 && u<=303 && v>98 && v<= 125)
				gyri[0].item(i)=29;
				
			if( u>287 && u<=303 && v>125 && v<= 179)
				gyri[0].item(i)=30;
				
			//bande 7
			if( u>303 && u<=340 && v>2 && v<= 30)
				gyri[0].item(i)=31;
				
			if( u>303 && u<=340 && v>30 && v<= 79)
				gyri[0].item(i)=32;
				
			if( u>303 && u<=340 && v>79 && v<= 98)
				gyri[0].item(i)=33;
				
			if( u>303 && u<=340 && v>98 && v<= 125)
				gyri[0].item(i)=34;
				
			if( u>303 && u<=340 && v>125 && v<= 179)
				gyri[0].item(i)=35;
				
			//bande 8
			if( u>340 && u<=360 && v>2 && v<= 30)
				gyri[0].item(i)=36;
				
			if( u>340 && u<=360 && v>30 && v<= 79)
				gyri[0].item(i)=37;
				
			if( u>340 && u<=360 && v>79 && v<= 98)
				gyri[0].item(i)=38;
				
			if( u>340 && u<=360 && v>98 && v<= 125)
				gyri[0].item(i)=39;
				
			if( u>340 && u<=360 && v>125 && v<= 179)
				gyri[0].item(i)=40;*/
	}
		
	Writer<Texture1d> gW(adOut);
	gW.write(gyri);
}


