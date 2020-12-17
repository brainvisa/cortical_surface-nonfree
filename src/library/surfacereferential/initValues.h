
#include <iostream>
#include <aims/mesh/texture.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <cortical_surface/mesh/isoLine.h>
#include <aims/utility/utility_g.h>
#include <aims/mesh/mesh_g.h>
#include <aims/io/io_g.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>


namespace aims
{
	
	class InitValues
	{
	
	public:
		
		std::map<float, std::vector<float> > values_meridian;
		std::map<float, std::vector<float> > values_parallel;
		
		std::map<float, int> cpt_meridian;
		std::map<float, int> cpt_parallel;
		
		std::vector<float> vals_mer;
		std::vector<float> vals_par;
		
		AimsSurfaceTriangle atlas;
		AimsSurfaceTriangle mesh;
		TimeTexture<float> longi, lati, texLongi, texLati;
		int nv;
		
		TimeTexture<float> V1p;
		TimeTexture<float> V3p;
		TimeTexture<float> V4p;
		TimeTexture<float> V6p;
		TimeTexture<float> V7p;
		TimeTexture<float> V8p;
		TimeTexture<float> V10p;
		TimeTexture<float> V11p;
		
		TimeTexture<float> V1m;
		TimeTexture<float> V3m;
		TimeTexture<float> V4m;
		TimeTexture<float> V6m;
		TimeTexture<float> V7m;
		TimeTexture<float> V8m;
		TimeTexture<float> V10m;
		TimeTexture<float> V11m;
	
		TimeTexture<float> P1p;
		TimeTexture<float> P3p;
		TimeTexture<float> P4p;
		TimeTexture<float> P6p;
		TimeTexture<float> P7p;
		TimeTexture<float> P8p;
		TimeTexture<float> P10p;
		TimeTexture<float> P11p;
		
		
		TimeTexture<float> P1m;
		TimeTexture<float> P3m;
		TimeTexture<float> P4m;
		TimeTexture<float> P6m;
		TimeTexture<float> P7m;
		TimeTexture<float> P8m;
		TimeTexture<float> P10m;
		TimeTexture<float> P11m;
		
			
		int size1;
		int size3;
		int size4;
		int size6;
		int size7;
		int size8;
		int size10;
		int size11;

		float result1;
		float result5;
		float result7;
		float result21;
		float result25;
		float result27;
		
		float result3;
		float result17;
		float result19;
		float result23;

		int cpt1;
		int cpt5;
		int cpt7;
		int cpt21;
		int cpt25;
		int cpt27;
		
		int cpt3;
		int cpt17;
		int cpt19;
		int cpt23;
		
		std::string adress;
		std::string adressMesh;
		
	InitValues(std::string & adressIn, std::string & adressMeshIn) : adress(adressIn),adressMesh(adressMeshIn)
	{
		Reader < AimsSurfaceTriangle > rm(adressMesh);
		rm.read( mesh );
	}
	
	void init();
	void CountMeridian();
	void CountParallels();
	void plotPar(float u, float v, float value, int n);
	void findLimitsGyrus();
	void createGyrus();
	void drawGyri();
	void drawGyriNew();
	void frontal();
	float calcAreas(TimeTexture<float> gyrus, float value);
	float areaTotal(TimeTexture<float> gyrus);
	float triangleAreas(Point3df & A, Point3df & B, Point3df & C );
	
	};
	
inline void InitValues::init(  )
	{
			//std::cout<<"CONSTRUCTOR DE LA MORT!!"<<std::endl;  
		std::cout << "reading atlas mesh   : " << std::flush;
		
// 		Reader<AimsSurfaceTriangle> triR( "/home/ced/data/atlas/atlas.mesh" );
		Reader<AimsSurfaceTriangle> triR( "V3/tri/V3_Lwhite_inflated.mesh" );
		triR >> atlas;
		std::cout << "done" << std::endl;
		
		std::cout << "reading atlas coordinates   : " << std::flush;
		
// 		Reader<TimeTexture<float> > lonR( "/home/ced/data/atlas/longitude.tex");
		Reader<TimeTexture<float> > lonR( "V3/cortical_referential_Left/longitude.tex");
		lonR >> longi;
// 		Reader<TimeTexture<float> > latR( "/home/ced/data/atlas/latitude.tex");
		Reader<TimeTexture<float> > latR( "V3/cortical_referential_Left/latitude.tex");
		latR >> lati;
		std::cout << "done" << std::endl;
		
		nv=longi[0].nItem();
		
		Reader<TimeTexture<float> > lonRt( "/home/ced/data/protocols/atlas/longitude.tex");
		lonRt >> texLongi;
		Reader<TimeTexture<float> > latRt( "/home/ced/data/protocols/atlas/latitude.tex");
		latRt >> texLati;
		
		for(int i=0; i<nv; i++)
		{
			texLongi[0].item(i)=0;
			texLati[0].item(i)=0;
		}
		
				/****Lecture des textures******/
		
		Reader<Texture1d> rtm1("V1/cortical_referential_Left/longitude.tex");
		rtm1.read(V1m);
		Reader<Texture1d> rtm3("V3/cortical_referential_Left/longitude.tex");
		rtm3.read(V3m);
		Reader<Texture1d> rtm4("V4/cortical_referential_Left/longitude.tex");
		rtm4.read(V4m);
		Reader<Texture1d> rtm6("V6/cortical_referential_Left/longitude.tex");
		rtm6.read(V6m);
		Reader<Texture1d> rtm7("V7/cortical_referential_Left/longitude.tex");
		rtm7.read(V7m);
		Reader<Texture1d> rtm8("V8/cortical_referential_Left/longitude.tex");
		rtm8.read(V8m);
		Reader<Texture1d> rtm10("V10/cortical_referential_Left/longitude.tex");
		rtm10.read(V10m);
		Reader<Texture1d> rtm11("V11/cortical_referential_Left/longitude.tex");
		rtm11.read(V11m);
		
		Reader<Texture1d> rt1("V1/cortical_referential_Left/latitude.tex");
		rt1.read(V1p);
		Reader<Texture1d> rt3("V3/cortical_referential_Left/latitude.tex");
		rt3.read(V3p);
		Reader<Texture1d> rt4("V4/cortical_referential_Left/latitude.tex");
		rt4.read(V4p);
		Reader<Texture1d> rt6("V6/cortical_referential_Left/latitude.tex");
		rt6.read(V6p);
		Reader<Texture1d> rt7("V7/cortical_referential_Left/latitude.tex");
		rt7.read(V7p);
		Reader<Texture1d> rt8("V8/cortical_referential_Left/latitude.tex");
		rt8.read(V8p);
		Reader<Texture1d> rt10("V10/cortical_referential_Left/latitude.tex");
		rt10.read(V10p);
		Reader<Texture1d> rt11("V11/cortical_referential_Left/latitude.tex");
		rt11.read(V11p);
		
		
		Reader<Texture1d> rpm1("V1/segment/V1_Lwhite_sulci_Mer_float.tex");
		rpm1.read(P1m);
		Reader<Texture1d> rpm3("V3/segment/V3_Lwhite_sulci_Mer_float.tex");
		rpm3.read(P3m);
		Reader<Texture1d> rpm4("V4/segment/V4_Lwhite_sulci_Mer_float.tex");
		rpm4.read(P4m);
		Reader<Texture1d> rpm6("V6/segment/V6_Lwhite_sulci_Mer_float.tex");
		rpm6.read(P6m);
		Reader<Texture1d> rpm7("V7/segment/V7_Lwhite_sulci_Mer_float.tex");
		rpm7.read(P7m);
		Reader<Texture1d> rpm8("V8/segment/V8_Lwhite_sulci_Mer_float.tex");
		rpm8.read(P8m);
		Reader<Texture1d> rpm10("V10/segment/V10_Lwhite_sulci_Mer_float.tex");
		rpm10.read(P10m);
		Reader<Texture1d> rpm11("V11/segment/V11_Lwhite_sulci_Mer_float.tex");
		rpm11.read(P11m);
		
		
		Reader<Texture1d> rp1("V1/segment/V1_Lwhite_sulci_Par_float.tex");
		rp1.read(P1p);
		Reader<Texture1d> rp3("V3/segment/V3_Lwhite_sulci_Par_float.tex");
		rp3.read(P3p);
		Reader<Texture1d> rp4("V4/segment/V4_Lwhite_sulci_Par_float.tex");
		rp4.read(P4p);
		Reader<Texture1d> rp6("V6/segment/V6_Lwhite_sulci_Par_float.tex");
		rp6.read(P6p);
		Reader<Texture1d> rp7("V7/segment/V7_Lwhite_sulci_Par_float.tex");
		rp7.read(P7p);
		Reader<Texture1d> rp8("V8/segment/V8_Lwhite_sulci_Par_float.tex");
		rp8.read(P8p);
		Reader<Texture1d> rp10("V10/segment/V10_Lwhite_sulci_Par_float.tex");
		rp10.read(P10p);
		Reader<Texture1d> rp11("V11/segment/V11_Lwhite_sulci_Par_float.tex");
		rp11.read(P11p);
	
		
		
		size1=V1p[0].nItem();
		size3=V3p[0].nItem();
		size4=V4p[0].nItem();
		size6=V6p[0].nItem();
		size7=V7p[0].nItem();
		size8=V8p[0].nItem();
		size10=V10p[0].nItem();
		size11=V11p[0].nItem();
		
		result1=0;
		result5=0;
		result7=0;
		result21=0;
		result25=0;
		result27=0;
		
		result3=0;
		result17=0;
		result19=0;
		result23=0;

		cpt1=0;
		cpt5=0;
		cpt7=0;
		cpt21=0;
		cpt25=0;
		cpt27=0;
		
		cpt3=0;
		cpt17=0;
		cpt19=0;
		cpt23=0;
		int tab[10];
// 		std::cout<<"Tableau des valeurs = "<<std::endl;
// 		for(int i=0;i<10;i++)
// 			tab[i]=0;
// 		for(int i=0; i<size1;i++)
// 		{
// 			if(P1p[0].item(i)!=0)
// 			{
// 				for(int j=0;j<10;j++)
// 				{
// 					if( tab[j]==P1p[0].item(i) )
// 						j=10;
// 					else
// 						if( tab[j]==0 )
// 						{
// 							tab[j]=P1p[0].item(i);
// 							j=10;
// 						}
// 						else
// 						{}
// 				}
// 			}
// 		}
// 		for(int i=0;i<10;i++)
// 			std::cout<<"Valeur = "<<tab[i]<<std::endl;
			
	}
	
	
inline void InitValues::CountMeridian(  )
	{
	
		

		for(int i=0;i<size1;i++)
		{
			if( P1m[0].item(i)!=0 )
			{
				values_meridian[ P1m[0].item(i) ].push_back( V1m[0].item(i) );
				cpt_meridian[ P1m[0].item(i) ]++;
			}
			if( P1p[0].item(i)!=0 )
			{
				values_parallel[ P1p[0].item(i) ].push_back( V1p[0].item(i) );
				cpt_parallel[ P1p[0].item(i) ]++;
			}
		}
		
		for(int i=0;i<size3;i++)
		{
			if( P3m[0].item(i)!=0 )
			{
				values_meridian[ P3m[0].item(i) ].push_back( V3m[0].item(i) );
				cpt_meridian[ P3m[0].item(i) ]++;
			}
			if( P3p[0].item(i)!=0 )
			{
				values_parallel[ P3p[0].item(i) ].push_back( V3p[0].item(i) );
				cpt_parallel[ P3p[0].item(i) ]++;
			}
		}
		
		for(int i=0;i<size4;i++)
		{
			if( P4m[0].item(i)!=0 )
			{
				values_meridian[ P4m[0].item(i) ].push_back( V4m[0].item(i) );
				cpt_meridian[ P4m[0].item(i) ]++;
			}
			if( P4p[0].item(i)!=0 )
			{
				values_parallel[ P4p[0].item(i) ].push_back( V4p[0].item(i) );
				cpt_parallel[ P4p[0].item(i) ]++;
			}
		}
		
		for(int i=0;i<size6;i++)
		{
			if( P6m[0].item(i)!=0 )
			{
				values_meridian[ P6m[0].item(i) ].push_back( V6m[0].item(i) );
				cpt_meridian[ P6m[0].item(i) ]++;
			}
			if( P6p[0].item(i)!=0 )
			{
				values_parallel[ P6p[0].item(i) ].push_back( V6p[0].item(i) );
				cpt_parallel[ P6p[0].item(i) ]++;
			}
		}
		for(int i=0;i<size6;i++)
		{
			if( P6m[0].item(i)!=0 )
			{
				values_meridian[ P6m[0].item(i) ].push_back( V6m[0].item(i) );
				cpt_meridian[ P6m[0].item(i) ]++;
			}
			if( P6p[0].item(i)!=0 )
			{
				values_parallel[ P6p[0].item(i) ].push_back( V6p[0].item(i) );
				cpt_parallel[ P6p[0].item(i) ]++;
			}
		}
		
		for(int i=0;i<size7;i++)
		{
			if( P7m[0].item(i)!=0 )
			{
				values_meridian[ P7m[0].item(i) ].push_back( V7m[0].item(i) );
				cpt_meridian[ P7m[0].item(i) ]++;
			}
			if( P7p[0].item(i)!=0 )
			{
				values_parallel[ P7p[0].item(i) ].push_back( V7p[0].item(i) );
				cpt_parallel[ P7p[0].item(i) ]++;
			}
		}
		
		for(int i=0;i<size8;i++)
		{
			if( P8m[0].item(i)!=0 )
			{
				values_meridian[ P8m[0].item(i) ].push_back( V8m[0].item(i) );
				cpt_meridian[ P8m[0].item(i) ]++;
			}
			if( P8p[0].item(i)!=0 )
			{
				values_parallel[ P8p[0].item(i) ].push_back( V8p[0].item(i) );
				cpt_parallel[ P8p[0].item(i) ]++;
			}
		}
		
		for(int i=0;i<size10;i++)
		{
			if( P10m[0].item(i)!=0 )
			{
				values_meridian[ P10m[0].item(i) ].push_back( V10m[0].item(i) );
				cpt_meridian[ P10m[0].item(i) ]++;
			}
			if( P10p[0].item(i)!=0 )
			{
				values_parallel[ P10p[0].item(i) ].push_back( V10p[0].item(i) );
				cpt_parallel[ P10p[0].item(i) ]++;
			}
		}
		
		for(int i=0;i<size11;i++)
		{
			if( P11m[0].item(i)!=0 )
			{
				values_meridian[ P11m[0].item(i) ].push_back( V11m[0].item(i) );
				cpt_meridian[ P11m[0].item(i) ]++;
			}
			if( P11p[0].item(i)!=0 )
			{
				values_parallel[ P11p[0].item(i) ].push_back( V11p[0].item(i) );
				cpt_parallel[ P11p[0].item(i) ]++;
			}
		}
		
		
		
		
		
			
			/*			if( P1m[0].item(i)==1 )
			{
				cpt1++;
				result1 += V1m[0].item(i);
			}
			if( P1m[0].item(i)==5 )
			{
				cpt5++;
				result5 += V1m[0].item(i);
			}
			if( P1m[0].item(i)==7 )
			{
				cpt7++;
				result7 += V1m[0].item(i);
			}
			if( P1m[0].item(i)==21 )
			{
				cpt21++;
				result21 += V1m[0].item(i);
			}
			if( P1m[0].item(i)==25 )
			{
				cpt25++;
				result25 += V1m[0].item(i);
			}
			if( P1m[0].item(i)==27 )
			{
				cpt27++;
				result27 += V1m[0].item(i);
			}
		}
		
		for(int i=0;i<size3;i++)
		{
			if( P3m[0].item(i)==1 )
			{
				cpt1++;
				result1 += V3m[0].item(i);
			}
			if( P3m[0].item(i)==5 )
			{
				cpt5++;
				result5 += V3m[0].item(i);
			}
			if( P3m[0].item(i)==7 )
			{
				cpt7++;
				result7 += V3m[0].item(i);
			}
			if( P3m[0].item(i)==21 )
			{
				cpt21++;
				result21 += V3m[0].item(i);
			}
			if( P3m[0].item(i)==25 )
			{
				cpt25++;
				result25 += V3m[0].item(i);
			}
			if( P3m[0].item(i)==27 )
			{
				cpt27++;
				result27 += V3m[0].item(i);
			}
		}
		
		for(int i=0;i<size4;i++)
		{
			if( P4m[0].item(i)==1 )
			{
				cpt1++;
				result1 += V4m[0].item(i);
			}
			if( P4m[0].item(i)==5 )
			{
				cpt5++;
				result5 += V4m[0].item(i);
			}
			if( P4m[0].item(i)==7 )
			{
				cpt7++;
				result7 += V4m[0].item(i);
			}
			if( P4m[0].item(i)==21 )
			{
				cpt21++;
				result21 += V4m[0].item(i);
			}
			if( P4m[0].item(i)==25 )
			{
				cpt25++;
				result25 += V4m[0].item(i);
			}
			if( P4m[0].item(i)==27 )
			{
				cpt27++;
				result27 += V4m[0].item(i);
			}
		}
		
		for(int i=0;i<size6;i++)
		{
			if( P6m[0].item(i)==1 )
			{
				cpt1++;
				result1 += V6m[0].item(i);
			}
			if( P6m[0].item(i)==5 )
			{
				cpt5++;
				result5 += V6m[0].item(i);
			}
			if( P6m[0].item(i)==7 )
			{
				cpt7++;
				result7 += V6m[0].item(i);
			}
			if( P6m[0].item(i)==21 )
			{
				cpt21++;
				result21 += V6m[0].item(i);
			}
			if( P6m[0].item(i)==25 )
			{
				cpt25++;
				result25 += V6m[0].item(i);
			}
			if( P6m[0].item(i)==27 )
			{
				cpt27++;
				result27 += V6m[0].item(i);
			}
		}
		
		for(int i=0;i<size7;i++)
		{
			if( P7m[0].item(i)==1 )
			{
				cpt1++;
				result1 += V7m[0].item(i);
			}
			if( P7m[0].item(i)==5 )
			{
				cpt5++;
				result5 += V7m[0].item(i);
			}
			if( P7m[0].item(i)==7 )
			{
				cpt7++;
				result7 += V7m[0].item(i);
			}
			if( P7m[0].item(i)==21 )
			{
				cpt21++;
				result21 += V7m[0].item(i);
			}
			if( P7m[0].item(i)==25 )
			{
				cpt25++;
				result25 += V7m[0].item(i);
			}
			if( P7m[0].item(i)==27 )
			{
				cpt27++;
				result27 += V7m[0].item(i);
			}
		}
		
		
		for(int i=0;i<size8;i++)
		{
			if( P8m[0].item(i)==1 )
			{
				cpt1++;
				result1 += V8m[0].item(i);
			}
			if( P8m[0].item(i)==5 )
			{
				cpt5++;
				result5 += V8m[0].item(i);
			}
			if( P8m[0].item(i)==7 )
			{
				cpt7++;
				result7 += V8m[0].item(i);
			}
			if( P8m[0].item(i)==21 )
			{
				cpt21++;
				result21 += V8m[0].item(i);
			}
			if( P8m[0].item(i)==25 )
			{
				cpt25++;
				result25 += V8m[0].item(i);
			}
			if( P8m[0].item(i)==27 )
			{
				cpt27++;
				result27 += V8m[0].item(i);
			}
		}
		
		
		for(int i=0;i<size10;i++)
		{
			if( P10m[0].item(i)==1 )
			{
				cpt1++;
				result1 += V10m[0].item(i);
			}
			if( P10m[0].item(i)==5 )
			{
				cpt5++;
				result5 += V10m[0].item(i);
			}
			if( P10m[0].item(i)==7 )
			{
				cpt7++;
				result7 += V10m[0].item(i);
			}
			if( P10m[0].item(i)==21 )
			{
				cpt21++;
				result21 += V10m[0].item(i);
			}
			if( P10m[0].item(i)==25 )
			{
				cpt25++;
				result25 += V10m[0].item(i);
			}
			if( P10m[0].item(i)==27 )
			{
				cpt27++;
				result27 += V10m[0].item(i);
			}
		}
		
		
		for(int i=0;i<size11;i++)
		{
			if( P11m[0].item(i)==1 )
			{
				cpt1++;
				result1 += V11m[0].item(i);
			}
			if( P11m[0].item(i)==5 )
			{
				cpt5++;
				result5 += V11m[0].item(i);
			}
			if( P11m[0].item(i)==7 )
			{
				cpt7++;
				result7 += V11m[0].item(i);
			}
			if( P11m[0].item(i)==21 )
			{
				cpt21++;
				result21 += V11m[0].item(i);
			}
			if( P11m[0].item(i)==25 )
			{
				cpt25++;
				result25 += V11m[0].item(i);
			}
			if( P11m[0].item(i)==27 )
			{
				cpt27++;
				result27 += V11m[0].item(i);
			}
		}*/
		
		
		
		FILE *f_m;
		FILE *f_p;
		
		std::string adstr_m ="stats_meridian.st";
		std::string adstr_p ="stats_parallel.st";
		
		const char* ad_m;
		const char* ad_p;
		
		ad_m=adstr_m.c_str();
		ad_p=adstr_p.c_str();
		
		f_m=fopen( ad_m , "w" );
		f_p=fopen( ad_p , "w" );
		
		/****RESULTS****/
		
		std::map<float, std::vector<float> >::iterator itMap;
		itMap = values_meridian.begin();
		fprintf( f_m,"MERIDIANS!!!\n\n",(*itMap).first );
		for ( ; itMap != values_meridian.end(); ++itMap)
		{
			fprintf( f_m,"VALUE = %f\n",(*itMap).first );
			fprintf( f_m,"CPT = %f\n",cpt_meridian[(*itMap).first] );
			for(int i=0; i< (*itMap).second.size(); i++)
			{
				fprintf( f_m,"%f\n",(*itMap).second[i] );
			}
			
			/*
			std::cout<<"VALUE = "<<(*itMap).first<<std::endl;
			std::cout<<"CPT = "<<cpt_meridian[(*itMap).first]<<std::endl;
			for(int i=0; i< (*itMap).second.size(); i++)
			{
				std::cout<<(*itMap).second[i]<<std::endl;
			}*/
		}
		fprintf( f_m,"SSTTOOPPP!!!\n\n",(*itMap).first );
		
		
		itMap = values_parallel.begin();
		fprintf( f_p,"PARRALLELS!!!\n\n",(*itMap).first );
		for ( ; itMap != values_parallel.end(); ++itMap)
		{
			fprintf( f_p,"VALUE = %f\n",(*itMap).first );
			fprintf( f_p,"CPT = %f\n",cpt_parallel[(*itMap).first] );
			for(int i=0; i< (*itMap).second.size(); i++)
			{
				fprintf( f_p,"%f\n",(*itMap).second[i] );
			}
			/*
			std::cout<<"VALUE = "<<(*itMap).first<<std::endl;
			std::cout<<"CPT = "<<cpt_parallel[(*itMap).first]<<std::endl;
			for(int i=0; i< (*itMap).second.size(); i++)
			{
				std::cout<<(*itMap).second[i]<<std::endl;
			}*/
		}
		fprintf( f_p,"SSTTOOPPP!!!\n\n",(*itMap).first );
		
		fclose(f_m);
		fclose(f_p);
		
		/*
		result1=result1/cpt1;
		result5=result5/cpt5;
		result7=result7/cpt7;
		result21=result21/cpt21;
		result25=result25/cpt25;
		result27=result27/cpt27;
		
		std::cout<<"RESULTATS MERIDIENS!!!!"<<std::endl;
		std::cout<<"pour 1 : "<<result1<<std::endl;
		std::cout<<"pour 5 : "<<result5<<std::endl;
		std::cout<<"pour 7 : "<<result7<<std::endl;
		std::cout<<"pour 21 : "<<result21<<std::endl;
		std::cout<<"pour 25 : "<<result25<<std::endl;
		std::cout<<"pour 27 : "<<result27<<std::endl;*/
		
	}
	
	
inline void InitValues::CountParallels(  )
	{
	
		/****Lecture des surfaces******/
	
/*		AimsSurfaceTriangle S1;
		AimsSurfaceTriangle S3;
		AimsSurfaceTriangle S4;
		AimsSurfaceTriangle S6;
		AimsSurfaceTriangle S7;
		AimsSurfaceTriangle S8;
		AimsSurfaceTriangle S10;
		AimsSurfaceTriangle S11;
		Reader < AimsSurfaceTriangle > r1("V1/tri/V1_Lwhite.mesh");
		r1.read( S1 );
		Reader < AimsSurfaceTriangle > r3("V3/tri/V3_Lwhite.mesh");
		r3.read( S3 );
		Reader < AimsSurfaceTriangle > r4("V4/tri/V4_Lwhite.mesh");
		r4.read( S4 );
		Reader < AimsSurfaceTriangle > r6("V6/tri/V6_Lwhite.mesh");
		r6.read( S6 );
		Reader < AimsSurfaceTriangle > r7("V7/tri/V7_Lwhite.mesh");
		r7.read( S7 );
		Reader < AimsSurfaceTriangle > r8("V8/tri/V8_Lwhite.mesh");
		r8.read( S8 );
		Reader < AimsSurfaceTriangle > r10("V10/tri/V10_Lwhite.mesh");
		r10.read( S10 );
		Reader < AimsSurfaceTriangle > r11("V11/tri/V11_Lwhite.mesh");
		r11.read( S11 );
		*/
		
		//float et3,ec15,ec19,ec23;
		float min3,max3,min15,max15,min19,max19,min23,max23;
		min3=P1p[0].item(0);
		max3=P1p[0].item(0);
		min15=P1p[0].item(0);
		max15=P1p[0].item(0);
		min19=P1p[0].item(0);
		max19=P1p[0].item(0);
		min23=P1p[0].item(0);
		max23=P1p[0].item(0);
		
		std::vector<float> tab3;
		std::vector<float> tab15;
		std::vector<float> tab19;
		std::vector<float> tab23;
		
		for(int i=0;i<size1;i++)
		{
			if( P1p[0].item(i)==3 )
			{
				cpt3++;
				tab3.push_back(V1p[0].item(i));
				result3 += V1p[0].item(i);
				plotPar(V1m[0].item(i),V1p[0].item(i),3,size1);
				if(min3<V1p[0].item(i))
					min3=V1p[0].item(i);
				if(max3>V1p[0].item(i))
					max3=V1p[0].item(i);
			}
			if( P1p[0].item(i)==15 )
			{
				cpt17++;
				tab15.push_back(V1p[0].item(i));
				result17 += V1p[0].item(i);
				plotPar(V1m[0].item(i),V1p[0].item(i),17,size1);
				if(min15<V1p[0].item(i))
					min15=V1p[0].item(i);
				if(max15>V1p[0].item(i))
					max15=V1p[0].item(i);
			}
			if( P1p[0].item(i)==19 )
			{
				cpt19++;
				tab19.push_back(V1p[0].item(i));
				result19 += V1p[0].item(i);
				plotPar(V1m[0].item(i),V1p[0].item(i),19,size1);
				if(min19<V1p[0].item(i))
					min19=V1p[0].item(i);
				if(max19>V1p[0].item(i))
					max19=V1p[0].item(i);
			}
			if( P1p[0].item(i)==23 )
			{
				cpt23++;
				tab23.push_back(V1p[0].item(i));
				result23 += V1p[0].item(i);
				plotPar(V1m[0].item(i),V1p[0].item(i),23,size1);
				if(min23<V1p[0].item(i))
					min23=V1p[0].item(i);
				if(max23>V1p[0].item(i))
					max23=V1p[0].item(i);
			}
		}
		
		for(int i=0;i<size3;i++)
		{
			if( P3p[0].item(i)==3 )
			{
				cpt3++;
				tab3.push_back(V3p[0].item(i));
				result3 += V3p[0].item(i);
				plotPar(V3m[0].item(i),V3p[0].item(i),3,size3);
				if(min3<V3p[0].item(i))
					min3=V3p[0].item(i);
				if(max3>V3p[0].item(i))
					max3=V3p[0].item(i);
			}
			if( P3p[0].item(i)==15 )
			{
				cpt17++;
				tab15.push_back(V3p[0].item(i));
				result17 += V3p[0].item(i);
				plotPar(V3m[0].item(i),V3p[0].item(i),17,size3);
				if(min15<V3p[0].item(i))
					min15=V3p[0].item(i);
				if(max15>V3p[0].item(i))
					max15=V3p[0].item(i);
			}
			if( P3p[0].item(i)==19 )
			{
				cpt19++;
				tab19.push_back(V3p[0].item(i));
				result19 += V3p[0].item(i);
				plotPar(V3m[0].item(i),V3p[0].item(i),19,size3);
				if(min19<V3p[0].item(i))
					min19=V3p[0].item(i);
				if(max19>V3p[0].item(i))
					max19=V3p[0].item(i);
			}
			if( P3p[0].item(i)==23 )
			{
				cpt23++;
				tab23.push_back(V3p[0].item(i));
				result23 += V3p[0].item(i);
				plotPar(V3m[0].item(i),V3p[0].item(i),23,size3);
				if(min23<V3p[0].item(i))
					min23=V3p[0].item(i);
				if(max23>V3p[0].item(i))
					max23=V3p[0].item(i);
			}
		}
		
		
		
		for(int i=0;i<size4;i++)
		{
			if( P4p[0].item(i)==3 )
			{
				cpt3++;
				tab3.push_back(V4p[0].item(i));
				result3 += V4p[0].item(i);
				plotPar(V4m[0].item(i),V4p[0].item(i),3,size4);
				if(min3<V4p[0].item(i))
					min3=V4p[0].item(i);
				if(max3>V4p[0].item(i))
					max3=V4p[0].item(i);
			}
			if( P4p[0].item(i)==15 )
			{
				cpt17++;
				tab15.push_back(V4p[0].item(i));
				result17 += V4p[0].item(i);
				plotPar(V4m[0].item(i),V4p[0].item(i),17,size4);
				if(min15<V4p[0].item(i))
					min15=V4p[0].item(i);
				if(max15>V4p[0].item(i))
					max15=V4p[0].item(i);
			}
			if( P4p[0].item(i)==19 )
			{
				cpt19++;
				tab19.push_back(V4p[0].item(i));
				result19 += V4p[0].item(i);
				plotPar(V4m[0].item(i),V4p[0].item(i),19,size4);
				if(min19<V4p[0].item(i))
					min19=V4p[0].item(i);
				if(max19>V4p[0].item(i))
					max19=V4p[0].item(i);
			}
			if( P4p[0].item(i)==23 )
			{
				cpt23++;
				tab23.push_back(V4p[0].item(i));
				result23 += V4p[0].item(i);
				plotPar(V4m[0].item(i),V4p[0].item(i),23,size4);
				if(min23<V4p[0].item(i))
					min23=V4p[0].item(i);
				if(max23>V4p[0].item(i))
					max23=V4p[0].item(i);
			}
		}
		
		
		for(int i=0;i<size6;i++)
		{
			if( P6p[0].item(i)==3 )
			{
				cpt3++;
				tab3.push_back(V6p[0].item(i));
				result3 += V6p[0].item(i);
				plotPar(V6m[0].item(i),V6p[0].item(i),3,size6);
				if(min3<V6p[0].item(i))
					min3=V6p[0].item(i);
				if(max3>V6p[0].item(i))
					max3=V6p[0].item(i);
			}
			if( P6p[0].item(i)==15 )
			{
				cpt17++;
				tab15.push_back(V6p[0].item(i));
				result17 += V6p[0].item(i);
				plotPar(V6m[0].item(i),V6p[0].item(i),17,size6);
				if(min15<V6p[0].item(i))
					min15=V6p[0].item(i);
				if(max15>V6p[0].item(i))
					max15=V6p[0].item(i);
			}
			if( P6p[0].item(i)==19 )
			{
				cpt19++;
				tab19.push_back(V6p[0].item(i));
				result19 += V6p[0].item(i);
				plotPar(V6m[0].item(i),V6p[0].item(i),19,size6);
				if(min19<V6p[0].item(i))
					min19=V6p[0].item(i);
				if(max19>V6p[0].item(i))
					max19=V6p[0].item(i);
			}
			if( P6p[0].item(i)==23 )
			{
				cpt23++;
				tab23.push_back(V6p[0].item(i));
				result23 += V6p[0].item(i);
				plotPar(V6m[0].item(i),V6p[0].item(i),23,size6);
				if(min23<V6p[0].item(i))
					min23=V6p[0].item(i);
				if(max23>V6p[0].item(i))
					max23=V6p[0].item(i);
			}
		}
		
		
		for(int i=0;i<size7;i++)
		{
			if( P7p[0].item(i)==3 )
			{
				cpt3++;
				tab3.push_back(V7p[0].item(i));
				result3 += V7p[0].item(i);
				plotPar(V7m[0].item(i),V7p[0].item(i),3,size7);
				if(min3<V7p[0].item(i))
					min3=V7p[0].item(i);
				if(max3>V7p[0].item(i))
					max3=V7p[0].item(i);
			}
			if( P7p[0].item(i)==15 )
			{
				cpt17++;
				tab15.push_back(V7p[0].item(i));
				result17 += V7p[0].item(i);
				plotPar(V7m[0].item(i),V7p[0].item(i),17,size7);
				if(min15<V7p[0].item(i))
					min15=V7p[0].item(i);
				if(max15>V7p[0].item(i))
					max15=V7p[0].item(i);
			}
			if( P7p[0].item(i)==19 )
			{
				cpt19++;
				tab19.push_back(V7p[0].item(i));
				result19 += V7p[0].item(i);
				plotPar(V7m[0].item(i),V7p[0].item(i),19,size7);
				if(min19<V7p[0].item(i))
					min19=V7p[0].item(i);
				if(max19>V7p[0].item(i))
					max19=V7p[0].item(i);
			}
			if( P7p[0].item(i)==23 )
			{
				cpt23++;
				tab23.push_back(V7p[0].item(i));
				result23 += V7p[0].item(i);
				plotPar(V7m[0].item(i),V7p[0].item(i),23,size7);
				if(min23<V7p[0].item(i))
					min23=V7p[0].item(i);
				if(max23>V7p[0].item(i))
					max23=V7p[0].item(i);
			}
		}
		
		
		for(int i=0;i<size8;i++)
		{
			if( P8p[0].item(i)==3 )
			{
				cpt3++;
				tab3.push_back(V8p[0].item(i));
				result3 += V8p[0].item(i);
				plotPar(V8m[0].item(i),V8p[0].item(i),3,size8);
				if(min3<V8p[0].item(i))
					min3=V8p[0].item(i);
				if(max3>V8p[0].item(i))
					max3=V8p[0].item(i);
			}
			if( P8p[0].item(i)==15 )
			{
				cpt17++;
				tab15.push_back(V8p[0].item(i));
				result17 += V8p[0].item(i);
				plotPar(V8m[0].item(i),V8p[0].item(i),17,size8);
				if(min15<V8p[0].item(i))
					min15=V8p[0].item(i);
				if(max15>V8p[0].item(i))
					max15=V8p[0].item(i);
			}
			if( P8p[0].item(i)==19 )
			{
				cpt19++;
				tab19.push_back(V8p[0].item(i));
				result19 += V8p[0].item(i);
				plotPar(V8m[0].item(i),V8p[0].item(i),19,size8);
				if(min19<V8p[0].item(i))
					min19=V8p[0].item(i);
				if(max19>V8p[0].item(i))
					max19=V8p[0].item(i);
			}
			if( P8p[0].item(i)==23 )
			{
				cpt23++;
				tab19.push_back(V8p[0].item(i));
				result23 += V8p[0].item(i);
				plotPar(V8m[0].item(i),V8p[0].item(i),23,size8);
				if(min23<V8p[0].item(i))
					min23=V8p[0].item(i);
				if(max23>V8p[0].item(i))
					max23=V8p[0].item(i);
			}
		}
		
		
		for(int i=0;i<size10;i++)
		{
			if( P10p[0].item(i)==3 )
			{
				cpt3++;
				tab3.push_back(V10p[0].item(i));
				result3 += V10p[0].item(i);
				plotPar(V10m[0].item(i),V10p[0].item(i),3,size10);
				if(min3<V10p[0].item(i))
					min3=V10p[0].item(i);
				if(max3>V10p[0].item(i))
					max3=V10p[0].item(i);
			}
			if( P10p[0].item(i)==15 )
			{
				cpt17++;
				tab15.push_back(V10p[0].item(i));
				result17 += V10p[0].item(i);
				plotPar(V10m[0].item(i),V10p[0].item(i),17,size10);
				if(min15<V10p[0].item(i))
					min15=V10p[0].item(i);
				if(max15>V10p[0].item(i))
					max15=V10p[0].item(i);
			}
			if( P10p[0].item(i)==19 )
			{
				cpt19++;
				tab19.push_back(V10p[0].item(i));
				result19 += V10p[0].item(i);
				plotPar(V10m[0].item(i),V10p[0].item(i),19,size10);
				if(min19<V10p[0].item(i))
					min19=V10p[0].item(i);
				if(max19>V10p[0].item(i))
					max19=V10p[0].item(i);
			}
			if( P10p[0].item(i)==23 )
			{
				cpt23++;
				tab23.push_back(V10p[0].item(i));
				result23 += V10p[0].item(i);
				plotPar(V10m[0].item(i),V10p[0].item(i),23,size10);
				if(min23<V10p[0].item(i))
					min23=V10p[0].item(i);
				if(max23>V10p[0].item(i))
					max23=V10p[0].item(i);
			}
		}
		
		
		for(int i=0;i<size11;i++)
		{
			if( P11p[0].item(i)==3 )
			{
				cpt3++;
				tab3.push_back(V11p[0].item(i));
				result3 += V11p[0].item(i);
				plotPar(V11m[0].item(i),V11p[0].item(i),3,size11);
				if(min3<V11p[0].item(i))
					min3=V11p[0].item(i);
				if(max3>V11p[0].item(i))
					max3=V11p[0].item(i);
			}
			if( P11p[0].item(i)==15 )
			{
				cpt17++;
				tab15.push_back(V11p[0].item(i));
				result17 += V11p[0].item(i);
				plotPar(V11m[0].item(i),V11p[0].item(i),17,size11);
				if(min15<V11p[0].item(i))
					min15=V11p[0].item(i);
				if(max15>V11p[0].item(i))
					max15=V11p[0].item(i);
			}
			if( P11p[0].item(i)==19 )
			{
				cpt19++;
				tab19.push_back(V11p[0].item(i));
				result19 += V11p[0].item(i);
				plotPar(V11m[0].item(i),V11p[0].item(i),19,size11);
				if(min19<V11p[0].item(i))
					min19=V11p[0].item(i);
				if(max19>V11p[0].item(i))
					max19=V11p[0].item(i);
			}
			if( P11p[0].item(i)==23 )
			{
				cpt23++;
				tab23.push_back(V11p[0].item(i));
				result23 += V11p[0].item(i);
				plotPar(V11m[0].item(i),V11p[0].item(i),23,size11);
				if(min23<V11p[0].item(i))
					min23=V11p[0].item(i);
				if(max23>V11p[0].item(i))
					max23=V11p[0].item(i);
			}
		}
		
		
		
		/****RESULTS****/
		
		result3=result3/cpt3;
		result17=result17/cpt17;
		result19=result19/cpt19;
		result23=result23/cpt23;
		
		std::cout<<"RESULTATS PARALLELES!!!!"<<std::endl;
		std::cout<<"pour 3 : "<<result3<<"  Min="<<min3<<" et max="<<max3<<std::endl;
		std::cout<<"pour 17 : "<<result17<<"  Min="<<min15<<" et max="<<max15<<std::endl;
		std::cout<<"pour 19 : "<<result19<<"  Min="<<min19<<" et max="<<max19<<std::endl;
		std::cout<<"pour 23 : "<<result23<<"  Min="<<min23<<" et max="<<max23<<std::endl;
		
		Writer<TimeTexture<float> > latW( "atlasLatitude.tex");
		latW.write(texLati);
		
		//Calcul ecart type
		float som=0;
		for(int i=0;i<cpt3;i++)
		{
			std::cout<<tab3[i]<<std::endl;
			som+=pow( (tab3[i]-result3) , 2 );
		}
		som=sqrt(som/cpt3);
		std::cout<<"Ecart Type pour 3 = "<<som<<std::endl;
		
		
		som=0;
		for(int i=0;i<cpt17;i++)
		{
			som+=pow( (tab15[i]-result17) , 2 );
		}
		som=sqrt(som/cpt17);
		std::cout<<"Ecart Type pour 17 = "<<som<<std::endl;
		
		som=0;
		for(int i=0;i<cpt19;i++)
		{
			som+=pow( (tab19[i]-result19) , 2 );
		}
		som=sqrt(som/cpt19);
		std::cout<<"Ecart Type pour 19 = "<<som<<std::endl;
		
		som=0;
		for(int i=0;i<cpt23;i++)
		{
			som+=pow( (tab23[i]-result23) , 2 );
		}
		som=sqrt(som/cpt23);
		std::cout<<"Ecart Type pour 23 = "<<som<<std::endl;
	}
	
inline void InitValues::plotPar(float u, float v, float value, int n)
	{
		//std::cout<<"Plotting paralleles constraints on a template"<<std::endl;
		
//   std::vector<Point3df> vertices=atlas.vertex();
//   int nv=vertices.size();
	
		int plot=0;
		float s,t;
		int nn=n;
		double dmin, d;
		dmin=1000.0;
		
		//finding the closet point on the mesh
		for (int j=0; j<nv; j++)
		{
			s=longi.item(j);
			t=lati.item(j);
			d=((u-s)*(u-s) + (v-t)*(v-t));
			if (d<dmin)
			{
				dmin=d;
				plot=j;
			}
		}
		//generating little sphere
		
		texLati[0].item(plot)=value;
	
	}
	
	
inline void InitValues::findLimitsGyrus()
	{
		float xMin, xMax, yMin, yMax, xMin2,  label, label2, label3; //xMax2, yMin2, yMax2,
		TimeTexture<float> gyriMap;
		TimeTexture<float> gyriResult(1,size10);
		Reader <TimeTexture<float> > Rg("V10/segment/V10_Lwhite_gyri_float.tex");
		Rg.read(gyriMap);
// 		xMin=V3m[0].item(0);
// 		xMax=V3m[0].item(0);
// 		yMin=V3p[0].item(0);
// 		yMax=V3p[0].item(0);
// 		xMin2=V3m[0].item(0);
// 		xMax2=V3m[0].item(0);
// 		yMin2=V3p[0].item(0);
// 		yMax2=V3p[0].item(0);
		label=11;
		label2=12;
		label3=2;
		int c=0;
		int c2=0;
		for(int i=1;i<size10;i++)
		{
			float x=gyriMap[2].item(i);
			gyriResult[0].item(i)=0;
			if(x==label)
			{
				c++;
				gyriResult[0].item(i)=10;
// 				if(V3m[0].item(i)>xMax && V3m[0].item(i) < 358)
// 					xMax=V3m[0].item(i);
// 				if(V3m[0].item(i)<xMin && V3m[0].item(i) > 2)
// 					xMin=V3m[0].item(i);
// 				if(V3p[0].item(i)>yMax && V3p[0].item(i) < 178)
// 					yMax=V3p[0].item(i);
// 				if(V3p[0].item(i)<yMin && V3p[0].item(i) > 2)
// 					yMin=V3p[0].item(i);
			}
			if(x==label2)
			{
				c2++;
				gyriResult[0].item(i)=20;/*
				if(V3m[0].item(i)>xMax2 && V3m[0].item(i) < 358)
					xMax2=V3m[0].item(i);
				if(V3m[0].item(i)<xMin2 && V3m[0].item(i) > 2)
					xMin2=V3m[0].item(i);
				if(V3p[0].item(i)>yMax2 && V3p[0].item(i) < 178)
					yMax2=V3p[0].item(i);
				if(V3p[0].item(i)<yMin2 && V3p[0].item(i) > 2)
					yMin2=V3p[0].item(i);*/
			}
			if(x==label3)
				gyriResult[0].item(i)=30;
			if(x==3)
				gyriResult[0].item(i)=40;
			if(x==4)
				gyriResult[0].item(i)=50;
			if(x==5)
				gyriResult[0].item(i)=60;
			if(x==6)
				gyriResult[0].item(i)=70;
			if(x==7)
				gyriResult[0].item(i)=80;
			if(x==8)
				gyriResult[0].item(i)=90;
			if(x==9)
				gyriResult[0].item(i)=100;
			if(x==10)
				gyriResult[0].item(i)=110;
			if(x==13)
				gyriResult[0].item(i)=120;
			if(x==14)
				gyriResult[0].item(i)=130;
			if(x==15)
				gyriResult[0].item(i)=140;
			if(x==16)
				gyriResult[0].item(i)=150;
			if(x==17)
				gyriResult[0].item(i)=160;
			if(x==18)
				gyriResult[0].item(i)=170;
			if(x==0)
				gyriResult[0].item(i)=200;

		}
		std::cout<<"longitude min et max: "<<xMin<<" et "<<xMax<<std::endl;
		std::cout<<"latitude min et max : "<<yMin<<" et "<<yMax<<std::endl;
		std::cout<<"longitude min2 et max2: "<<xMin<<" et "<<xMax<<std::endl;
		std::cout<<"latitude min2 et max2 : "<<yMin<<" et "<<yMax<<std::endl;
		std::cout<<"compteur = "<<c<<std::endl;
		std::cout<<"compteur2 = "<<c2<<std::endl;
		Writer<TimeTexture<float> > latW( "V10/gyrusPreCentral.tex");
		latW.write(gyriResult);
	}

	
inline void InitValues::createGyrus()
	{
		TimeTexture<float> gyriMap;
		TimeTexture<float> texP;
		TimeTexture<float> texM;
		AimsSurfaceTriangle meshResult;
		AimsSurfaceTriangle mesh;
		SurfaceManip merging;
		
		std::cout<<"OK0"<<std::endl;
		
		Reader < AimsSurfaceTriangle > rm("V10/tri/V10_Lwhite_inflated_30.mesh");
		rm.read( mesh );
		
		texP=V10p;
		texM=V10m;
		
		int lon[]={2,10,30,100,280,300,340,358};
		int lat[]={2,27,65,80,100,178};
		
		for(int i=0;i<6;i++)
		{
			AimsSurfaceTriangle msh;
			IsoLine mt(mesh, texP, lat[i]);
			msh=mt.makeTubes();
		
			merging.meshMerge<3, Void>(meshResult, msh);
			//delete msh;
		}
		
		for(int i=0;i<8;i++)
		{
			AimsSurfaceTriangle msh;
			IsoLine mt(mesh, texM, lon[i]);
			msh=mt.makeTubes();
		
			merging.meshMerge<3, Void>(meshResult, msh);
			//delete msh;
		}
		
	
		Writer<AimsSurfaceTriangle> sphW("V10/meshGyri30.mesh");
		sphW.write(meshResult);
	
	}
	
inline void InitValues::drawGyri()
	{
	
		TimeTexture<float> longitude;
		TimeTexture<float> latitude;
		
		int size;
		
		std::string adLong=adress+"/longitude.tex";
		std::string adLat=adress+"/latitude.tex";
		
		//cout<<"Adress Longitude = "<<adLong<<endl;
		//cout<<"Adress Latitude = "<<adLat<<endl;
		
		Reader< TimeTexture<float> > longR( adLong );
		longR.read(longitude);
		Reader< TimeTexture<float> > latR( adLat );
		latR.read(latitude);
		
		size=longitude.nItem();
		
		TimeTexture<float> gyri(1,size);
		
		float u=0;
		float v=0;
		
		//float step=0;
		
		//int lon[]={0,16,43,58,100,287,303,340,360};
		//int lat[]={2,30,79,98,125,179};
		
		for(int i=0;i<size;i++)
		{		
			gyri[0].item(i)=0;
			
			u=longitude[0].item(i);
			v=latitude[0].item(i);
			
			//Pre-central
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
				gyri[0].item(i)=40;
		}
		
		Writer<Texture1d> gW(adress+"/gyriMap.tex");
		gW.write(gyri);
		
		FILE *f;
		
		std::string adstr = adress+"/stats.st";
		
		const char* ad;
		
		ad=adstr.c_str();
		
		f=fopen( ad , "w" );
		
		fprintf( f,"%f\n",areaTotal(gyri) );
		
		
				
		for(int i=1; i<=40; i++)
 			fprintf( f,"%f\n",calcAreas(gyri, i) );

		std::cout<<"fichier ecrit!"<<std::endl;

		fclose(f);
		
// 		String ad = adress+"stats.st";
// 		
// 		ifstream f(ad);
// 		if (!f)
// 		{
// 			cerr << "Impossible d'ouvrir le fichier!" << endl;
// 			return (-1);
// 		}


//		std::cout << "Aire totale de l'hemisphere: ";
// 		f<<areaTotal(gyri)<<std::endl;
		 		
		/*		
		for(int i=1; i<=40; i++)
 			f<<(int)calcAreas(gyri, i)<<std::endl;

		std::cout<<"fichier ecrit!"<<std::endl;

		f.close();*/
		
// //		std::cout<<"avant le pre-central: ";
// 		for(int i=6; i<=10; i++)
// 			std::cout<<(int)calcAreas(gyri, i)<<" ";
// 		std::cout<<std::endl;
// 			
// // 		std::cout<<"Aires des gyri pre-centraux: ";
// 		for(int i=1; i<=5; i++)
// 			std::cout<<(int)calcAreas(gyri, i)<<" ";
// 		std::cout<<std::endl;
// 		
// // 		std::cout<<"Aires des gyri post-centraux: ";
// 		for(int i=36; i<=40; i++)
// 			std::cout<<(int)calcAreas(gyri, i)<<" ";
// 		std::cout<<std::endl;
// 		
// 		
// // 		std::cout<<"apr� le post-central: ";
// 		for(int i=31; i<=35; i++)
// 			std::cout<<(int)calcAreas(gyri, i)<<" ";
// 		std::cout<<std::endl;
// 		
// 		
// 		std::cout<<std::endl;
		
	}	
	
inline void InitValues::drawGyriNew()
	{
	
	
		TimeTexture<float> longitude;
		TimeTexture<float> latitude;
		TimeTexture<short> gyr_arnaud_origin;
		
		int size;
		
		std::string adLong=adress+"/cortical_referential/longitude.tex";
		std::string adLat=adress+"/cortical_referential/latitude.tex";
		std::string adGyr=adress+"/segment/"+adress+"_Lwhite_gyri.tex";
		
		std::cout<<"Adress Longitude = "<<adLong<<std::endl;
		std::cout<<"Adress Latgitude = "<<adLat<<std::endl;
		
		Reader< TimeTexture<float> > longR( adLong );
		longR.read(longitude);
		Reader< TimeTexture<float> > latR( adLat );
		latR.read(latitude);
		Reader< TimeTexture<short> > gyrR( adGyr );
		gyrR.read(gyr_arnaud_origin);
		
		size=longitude.nItem();
		
		TimeTexture<float> gyr_arnaud(1,size);
		TimeTexture<float> gyri(1,size);
		
		float u=0;
		float v=0;
		
		for(int i=0;i<size;i++)
		{
			gyri[0].item(i)=0;
			
			gyr_arnaud[0].item(i)=(float)gyr_arnaud_origin[2].item(i);
			
			u=longitude[0].item(i);
			v=latitude[0].item(i);
			
			//Pre-central
			if( u<=16 && v>30 && v<= 179)
				gyri[0].item(i)=9;
				
			//Frontal inferior
			if( u>16 && u<=43 && v>30 && v<= 79)
				gyri[0].item(i)=2;
				
			//Frontal middle
			if( u>16 && u<=43 && v>79 && v<= 125)
				gyri[0].item(i)=11;
				
			//Frontal superior
			if( u>16 && u<=43 && v>125 && v<= 179)
				gyri[0].item(i)=5;
				
			//Orbital??
			if( u>43 && u<=58 && v>30 && v<= 179)
				gyri[0].item(i)=12;
				
			//?
			if( u>58 && u<=100 && v>30 && v<= 179)
				gyri[0].item(i)=7;
				
			//Parietal
			if( u>100 && u<=287 && v>30 && v<= 79)
				gyri[0].item(i)=8;
				
			//Temporal inf�ieur
			if( u>100 && u<=287 && v>79 && v<= 98)
				gyri[0].item(i)=13;
				
			//Temporal middle
			if( (u>100 && u<=287 && v>98 && v<= 125) || (u>287 && u<=303 && v>79 && v<= 125) )
				gyri[0].item(i)=14;
				
			//Temporal superior
			if( u>100 && u<=303 && v>125 && v<= 179)
				gyri[0].item(i)=6;
				
			//Cuneus
			if( (u>287 && u<=340 && v>30 && v<= 79) || (u>303 && u<=340 && v>79 && v<= 98) )
				gyri[0].item(i)=15;
				
			//?
			if( u>303 && u<=340 && v>98 && v<= 179)
				gyri[0].item(i)=4;
				
			//Post-central
			if( u>340 && u<=360 && v>30 && v<= 179)
				gyri[0].item(i)=10;
		}
		
		Writer<Texture1d> gW(adress+"/gyriMapNew.tex");
		gW.write(gyri);
		
		Writer<Texture1d> gaW(adress+"/gyriArnaud.tex");
		gaW.write(gyr_arnaud);
		
	}	
	
inline void InitValues::frontal()
	{
		TimeTexture<float> lat1;
		TimeTexture<float> lat2;
		TimeTexture<float> lat3;
		TimeTexture<float> lat4;
		
		TimeTexture<float> proj1;
		TimeTexture<float> proj2;
		TimeTexture<float> proj3;
		TimeTexture<float> proj4;
		
/*		int size1;
		int size2;
		int size3;
		int size4;*/
		
		Reader< TimeTexture<float> > l1R( "V1/cortical_referential/latitude.tex" );
		l1R.read(lat1);
		Reader< TimeTexture<float> > l2R( "V4/cortical_referential/latitude.tex" );
		l2R.read(lat2);
		Reader< TimeTexture<float> > l3R( "V7/cortical_referential/latitude.tex" );
		l3R.read(lat3);
		Reader< TimeTexture<float> > l4R( "V8/cortical_referential/latitude.tex" );
		l4R.read(lat4);
		
		Reader< TimeTexture<float> > l1p( "V1/segment/V1_Lwhite_sulci_float.tex" );
		l1p.read(lat1);
		Reader< TimeTexture<float> > l2p( "V4/segment/V4_Lwhite_sulci_float.tex" );
		l2p.read(lat2);
		Reader< TimeTexture<float> > l3p( "V7/segment/V7_Lwhite_sulci_float.tex" );
		l3p.read(lat3);
		Reader< TimeTexture<float> > l4p( "V8/segment/V8_Lwhite_sulci_float.tex" );
		l4p.read(lat4);
		
	}	
	
inline float InitValues::areaTotal(TimeTexture<float> gyrus)
	{
		std::vector< Point3df > vertex;
		std::vector< AimsVector< uint,3> > poly;
	
		vertex = mesh.vertex();
		poly = mesh.polygon();

		TimeTexture<float> tex;
		unsigned  i, p, t1, t2, t3;
		//AimsSurfaceTriangle::iterator itMesh;
		//unsigned		nedge = 0;
		
		float aire;
		
		p = poly.size();
	
		aire = 0;
		
		for( i=0; i<p; ++i )
		{
			t1 = poly[i][0];
			t2 = poly[i][1];
			t3 = poly[i][2];
			if( gyrus[0].item(t1)!=0 && gyrus[0].item(t2)!=0 && gyrus[0].item(t3)!=0 )
				aire += triangleAreas( vertex[t1], vertex[t2], vertex[t3] );
			
		}
	
		return aire;
	
	}
	
	
inline float InitValues::calcAreas(TimeTexture<float> gyrus, float value)
	{
		std::vector< Point3df > vertex;
		std::vector< AimsVector< uint,3> > poly;
	
		vertex = mesh.vertex();
		poly = mesh.polygon();

		TimeTexture<float> tex;
		unsigned  i, p, t1, t2, t3;
		//AimsSurfaceTriangle::iterator itMesh;
		//unsigned		nedge = 0;
		
		float aire;
		
		p = poly.size();
	
		aire = 0;
		
		for( i=0; i<p; ++i )
		{
			t1 = poly[i][0];
			t2 = poly[i][1];
			t3 = poly[i][2];
			
			if( gyrus[0].item(t1)==value && gyrus[0].item(t2)==value && gyrus[0].item(t3)==value )
			{
				aire += triangleAreas( vertex[t1], vertex[t2], vertex[t3] );
			}
			
		}
	
		//std::cout<<"Aire du gyrus "<<value<<" = "<<aire<<std::endl;
	
		return aire;
	
	}
	
inline float InitValues::triangleAreas(Point3df & A, Point3df & B, Point3df & C )
	{
	
		float aire = 0;
		float dista,distb,distc;
		
		distc = sqrt( ( pow( (B[0]-A[0]),2) + pow( (B[1]-A[1]),2) + pow( (B[2]-A[2]),2) ) );
		dista = sqrt( ( pow( (C[0]-B[0]),2) + pow( (C[1]-B[1]),2) + pow( (C[2]-B[2]),2) ) );
		distb = sqrt( ( pow( (A[0]-C[0]),2) + pow( (A[1]-C[1]),2) + pow( (A[2]-C[2]),2) ) );
	
		float inter1;
		
		inter1 = ( ( (distc*distc - distb*distb)/dista ) + dista ) / 2;
	
		aire = (dista/2) * sqrt ( distc*distc - inter1*inter1 );
	
		return aire;
	}
}

