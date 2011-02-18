#include <cstdlib>
#include <cortical_surface/surfacereferential/corticalReferential.h>

#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshvoronoi.h>
#include <aims/scalespace/meshDiffuse.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/distancemap/meshmorphomat.h>

using namespace aims::meshdistance;


namespace aims {
    /****************MAIN METHODS****************/
    
    //Process is the main method. It calls the parallels propagation
    //and the meridians propagation
    
    void CorticalReferential::process()
    {
    	std::cout << "Before preprocessing..." << std::endl;
    	constraintPreprocess();
    	
    	TimeTexture<float> roots ( 1, size );
    	Reader < TimeTexture<float> > rtp ( adr_par );
    	rtp.read ( roots );
    	
    	int tab [180];
    	
    	for ( int i = 0 ; i < 180 ; i++ )
    		tab[i] = 0;
    	
    	for ( int i = 0 ; i < size ; i++ )
    	{
    		if ( roots[0].item(i) != 0 )
    			tab[(int) roots[0].item(i) ]++;
    	}
    	
    	for ( int i = 0 ; i < 180 ; i++ )
    	{
    		if( tab[i] != 0 )
    			std::cout << "Val=" << i << " - size=" << tab[i] << std::endl;
    	}
    	std::cout << "After preprocessing..." << std::endl;
    	
    	switch ( choice_process )
    	{
    		case 1:
    		    std::cerr << "choice_process=" << choice_process << " - Latitude propagation only" << std::endl;
    			latitudePropagation();
    			break;
    		case 2:
    		    std::cerr << "choice_process=" << choice_process << " - Longitude propagation only" << std::endl;
    			longitudePropagation();
    			break;
    		case 3: 
    			break;
    		default:
    		    std::cerr << "choice_process=" << choice_process << " - Latitude and longitude propagations" << std::endl;
    			latitudePropagation();
    			longitudePropagation();
    			break;
    	}
    }
    
    //Preprocessing constraints (i.e. finding pole)
    
    void CorticalReferential::constraintPreprocess()
    {
    	
    	std::cout << " Constraint preprocess...reading cleaned lat/lon constraints textures" << std::endl;
    	Reader < TimeTexture<float> > rtp(adr_par);
    	rtp.read( constraint_lat_cleaned );
    	
    	Reader < TimeTexture<float> > rtm(adr_mer);
    	rtm.read( constraint_long_cleaned );
    	
    	std::cout << "  Find poles..." << std::endl;
    	poles_points = find_poles ( poles, mesh );
    	std::cout << "   ok..." << std::endl;
    	nord = poles_points.first;
    	sud = poles_points.second;
    	
    	std::cout << "done" << std::endl;
    }
    
    
    //Defining the new insula pole point
    void CorticalReferential::cingularPoint()
    {
    	TimeTexture<float> dist_orig ( 1, size );
    	TimeTexture<short> pole_ins_short ( 1, size );
    	TimeTexture<float> pole_insula ( 1, size );
    	TimeTexture<short> dist_erod ( 1, size );
    	TimeTexture<float> dist_cc ( 1,size );
    	for ( int i = 0 ; i < size ; i++ )
    	{
    		if ( poles[0].item(i)==1 )
    		{
    			dist_orig[0].item(i)=1;
    		}
    		else
    			dist_orig[0].item(i)=0;
    		if ( poles[0].item(i)==180 )
    		{
    			pole_ins_short[0].item(i)=(short)poles[0].item(i);
    		}
    		else
    		{
    			pole_ins_short[0].item(i)=0;
    		}
    	}
    	
    	std::cout << "meshdistance" << std::endl;
    	std::cout << "debut" << std::endl;
    	dist_cc[0]=meshdistance::MeshDistance( mesh[0], dist_orig[0],true);
    	dist_erod[0]=MeshErosion<short> ( mesh[0], pole_ins_short[0], short(0), -1, 3.0, true);
    	
    	for (int i=0;i<size; i++)
    	{
    		pole_insula[0].item(i)=0;
    		if(dist_erod[0].item(i)!=0)
    			pole_insula[0].item(i)=180;
    	}
    /*	Writer<Texture1d> wT1cr("point_poleinsula.tex");
    	wT1cr.write(pole_insula);
    	Writer<Texture1d> wT1c("point_distinsula.tex");
    	wT1c.write(dist_cc);*/
    	
    	float dist_min=100;
    	ind_min=0;
    	for (int i=0;i<size; i++)
    	{
    		if( (pole_insula[0].item(i)==180) && (dist_cc[0].item(i)<dist_min) )
    		{
    			dist_min=dist_cc[0].item(i);
    			ind_min=i;
    		}
    	}
    }
    
    //The whole latitude propagation process
    
    void CorticalReferential::latitudePropagation()
    {
    	
    	std::cout << std::endl;
    	std::cout << "LATITUDE" << std::endl;
    	
    	TimeTexture<float> diff_lat_result(1,size);
    	TimeTexture<float> diff_lat_dilate_float(1,size);
    	TimeTexture<short> diff_lat_dilate(1,size);
    	
    	//Nouveau pole insulaire
    	cingularPoint();
    	int sud_temp=0;
    	sud_temp=sud;
    	sud=ind_min;
    	
    	//pole_insula_single[0].item(ind_min)=180;
    	
    	//ebarbulage
    	std::cout << "Cleaning constraints.........";
    	fflush(stdout);
    	
    	TimeTexture<float> contraint_tmp ( 1, size );
    	
    	//A LA PLACE DE TOUT CE QUI SUIT!!
    	//ON LIT SIMPLEMENT LA TEXTURE CREEE PAR CONSTRAINT_CLEANER.H
    	
    	std::cout << "done" << std::endl;
    	
    	////////////////////////////////////////////////////////////
    	// CERCLES POLAIRES!!
    	////////////////////////////////////////////////////////////
    	TimeTexture<float> cercles(1,size);
    	
    	std::set<uint>::const_iterator itvoisin;
    	
    	for(int i=0; i<size; i++)
    	{
    		cercles[0].item(i)=0;
    		int cpt1=0, cpt2=0;
    		if(poles[0].item(i)!=0)
    		{
    			itvoisin=neigh[i].begin();
    			for(; itvoisin!=neigh[i].end(); itvoisin++)
    			{
    				if(poles[0].item(*itvoisin)!=0)
    					cpt1++;
    				if(poles[0].item(*itvoisin)==0)
    					cpt2++;
    			}
    			//Enleve les restes de contrainte a l'interieur des poles
    			constraint_lat_cleaned[0].item(i)=0;
    		}
    		if( cpt1!=0 && cpt2!=0 )
    			cercles[0].item(i)=poles[0].item(i);
    		
    		if( cercles[0].item(i)==1 )
    			cercle_polaire[0].push_back(31);
    		else
    		{
    			if( cercles[0].item(i)==180 )
    				cercle_polaire[0].push_back(151);
    			else
    				cercle_polaire[0].push_back(0);
    		}
    	}
    	
    	////////////////////////////////////////////////////////////
    	// FIN!!
    	////////////////////////////////////////////////////////////
    	
    	//Mise a jour de la texture d'origine pour la propagation
    	
    	for(int i=0; i<size; i++)
    	{
    		if( i==sud )
    			constraint_lat_cleaned[0].item(i)=181;
    		else if(i==nord )
    			constraint_lat_cleaned[0].item(i)=1;
    		
    		if( cercle_polaire[0].item(i)!=0 )
    			constraint_lat_cleaned[0].item(i)=cercle_polaire[0].item(i);
    	}
    	
    	//process de propagation des coordonnees
    	std::cout << "Diffusion....................";
    	fflush(stdout);
    	diff_lat_result = diffusionLatitudeRelax ( constraint_lat_cleaned );
    	
    	//Remise de 0 a 180
    	
    	for(int i=0; i<size; i++)
    		diff_lat_result[0].item(i)=diff_lat_result[0].item(i)-1;
    	
    	
    	diff_lat_result[0].item(nord)=0;
    	diff_lat_result[0].item(sud)=180;
    	
    	Writer<Texture1d> wT1c(adr_lat);
    	wT1c.write(diff_lat_result);
    	
    	
    	sud=sud_temp;
    	
    	std::cout << "done" << std::endl;
    
    }
    
    //The whole longitude propagation process
    void CorticalReferential::longitudePropagation()
    {
    
    	std::cout << std::endl;
    	std::cout << "LONGITUDE" << std::endl;
    
    	TimeTexture<float> te(1,size);
    	TimeTexture<float> constraint_long_side(1,size); //les 2 cotes du mer. d'origine
    	TimeTexture<float> diff_long_result(1,size);
    	TimeTexture<float> constraint_empty(1,size);
    	TimeTexture<float> new_constraint_pole(1,size);
    	
    	init_texture_single(te);
    
    	std::cout << "Origin Meridian.............." << std::endl;
    	constraint_long_side = origin_meridian ( constraint_long_cleaned, nord, sud, neigh, mesh, poles );
    	
    // 	Writer<Texture1d> sideB("sides_origine.tex");
    // 	sideB.write(constraint_long_side);
    
    	
    //***A VIRER
    /*	Writer<Texture1d> wT3mtt("meridien_long_side_avant.tex");
    	wT3mtt.write(constraint_long_side);*/
    // 	Writer<Texture1d> wT3mff("meridien_long_cleaned_avant.tex");
    // 	wT3mff.write(constraint_long_cleaned);
    
    	TimeTexture<float> constraint_long_side_long(1,size);  //sauvegarde des cotes y compris dans les poles, pour la diffusion dans le cingulaire
    	TimeTexture<float> constraint_long_side_temp(1,size);
    	TimeTexture<float> save_mer_origin(1,size);
    	std::cout << "done" << std::endl;
    	
    	
    	std::set<uint>::iterator itt;
    	
    	for(int i=0; i<size; i++)
    	{
    		if( poles[0].item(i)==0 )
    		{
    			float val=0;
    			int cptt=0;
    			itt=neigh[i].begin();
    			for(; itt!=neigh[i].end(); itt++)
    			{
    				if( poles[0].item(*itt)!=0 )
    					val=poles[0].item(*itt);
    				else
    					cptt++;
    			}
    			if(cptt==0)
    			{
    				//Permet d'eliminer les points solitaires dans les poles differents de la valeur de pole
    				poles[0].item(i)=val;
    				std::cout << "correction="<<i << std::endl;
    			}
    		}
    	}
    	
    	for(int i=0; i<size; i++)
    	{
    		save_mer_origin[0].item(i)=0;
    		if(constraint_long_cleaned[0].item(i)==360)
    			save_mer_origin[0].item(i)=360;
    	}
    	
    	for(int i=0; i<size; i++)
    	{
    		constraint_long_side_long[0].item(i)=constraint_long_side[0].item(i);
    		new_constraint_pole[0].item(i)=constraint_long_side[0].item(i);
    		if( poles[0].item(i)!=0 )
    		{
    			constraint_long_cleaned[0].item(i)=360;
    			constraint_long_side[0].item(i)=0;
    		}
    		constraint_long_side_temp[0].item(i)=constraint_long_side[0].item(i);
    		
    	}
    	
    // 	Writer<Texture1d> sideWss("sides.tex");
    // 	sideWss.write(constraint_long_side_temp);
    
    
    	constraint_long_side = defineSides( constraint_long_side_temp, constraint_long_cleaned, mesh, neigh);
    
    	//selon le cas:
    	if(context==1)
    		for(int i=0; i<size; i++)
    	{
    		if(constraint_long_side[0].item(i)==2)
    			constraint_long_side[0].item(i)=4;
    		else
    			if(constraint_long_side[0].item(i)==4)
    				constraint_long_side[0].item(i)=2;
    	}
    
    
    	te[0].item(nord)=50;
    	te[0].item(sud)=100;
    /*	Writer<Texture1d> sideWssp("points_lele.tex");
    	sideWssp.write(te);*/
    	
    	
    	for(int i=0; i<size; i++)
    	{
    		if( constraint_long_cleaned[0].item(i)==360 )
    			constraint_empty[0].item(i)=360;
    		else
    			constraint_empty[0].item(i)=0;
    	}
    	
    /*	Writer<Texture1d> wT3c("poles.tex");
    	wT3c.write(poles);*/
    	
    	int Save_Point_Insula=-1;
    	
    	for(int i=0; i<size; i++)
    	{
    		if( poles[0].item(i)==180)
    		{
    			
    			std::set<uint>::iterator iti;
    			iti=neigh[i].begin();
    			int yes1=0;
    			int yes2=0;
    			int yes3=0;
    			for(; iti!=neigh[i].end(); iti++)
    			{
    				if( poles[0].item(*iti)==0 && constraint_long_cleaned[0].item(*iti)==360 )
    					yes1++;
    				if( constraint_long_side[0].item(*iti)==2 )
    					yes2++;
    				if( constraint_long_side[0].item(*iti)==4 )
    					yes3++;
    			}
    			if( yes1!=0 && yes2!=0 && yes3!=0 )
    				Save_Point_Insula=i;
    		}
    		
    	}
    	
    	if(Save_Point_Insula==-1)
    	{
    		std::cout << "Point insula non trouvé...." << std::endl;
    		for(int i=0; i<size; i++)
    		{
    			if( poles[0].item(i)==180)
    			{
    				
    				std::set<uint>::iterator it;
    				it=neigh[i].begin();
    				for(; it!=neigh[i].end(); it++)
    				{
    					if( poles[0].item(*it)!=180 )
    					{
    						std::set<uint>::iterator iti;
    						iti=neigh[(*it)].begin();
    						int yes1=0;
    						int yes2=0;
    						int yes3=0;
    						for(; iti!=neigh[(*it)].end(); iti++)
    						{
    							if( poles[0].item(*iti)==0 && constraint_long_cleaned[0].item(*iti)==360 )
    								yes1++;
    							if( constraint_long_side[0].item(*iti)==2 )
    								yes2++;
    							if( constraint_long_side[0].item(*iti)==4 )
    								yes3++;
    						}
    						if( yes1!=0 && yes2!=0 && yes3!=0 )
    							Save_Point_Insula=i;
    					}
    				}
    			}
    		}
    	}
    	
    	std::cout << "Point insula="<<Save_Point_Insula << std::endl;
    	
    // 	Writer<Texture1d> wT3c("/home/olivier/constraint_long_cleaned.tex");
    // 	wT3c.write(constraint_long_cleaned);
    // 	
    // 	Writer<Texture1d> wT3cs("/home/olivier/onstraint_long_side.tex");
    // 	wT3cs.write(constraint_long_side);
    	
    	//Diffusion de la longitude
    	std::cout << "Diffusion....................";
    	fflush(stdout);
    	
    	//***********************************************************************
    	//Diffusion longitude****************************************
    	//***********************************************************************
    	
    	diff_long_result=diffusionLongitudeRelax(constraint_long_cleaned, constraint_long_side);
    	
    //	Writer<Texture1d> wT1c("diff_lon_sans_pole.tex");
    //	wT1c.write(diff_long_result);*/
    	
    /*	Writer<Texture1d> wT3mc("meridien_long_side.tex");
    	wT3mc.write(constraint_long_cleaned);*/
    /*	Writer<Texture1d> wT3m("meridien_long_cleaned.tex");
    	wT3m.write(constraint_long_side);*/
    	
    	
    	TimeTexture<float> temp(1,size);
    	init_texture_single(temp);
    	
    	for(int i=0; i<size; i++)
    	{
    		if( constraint_long_cleaned[0].item(i)==360 && poles[0].item(i)==0 )
    		{
    			diff_long_result[0].item(i)=360;
    			temp[0].item(i)=360;
    		}
    	}
    	
    	//A REMETTRE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	Writer<Texture1d> wT1kkc(adr_lon);
    	wT1kkc.write(diff_long_result);
    	//-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	
    /*	Writer<Texture1d> wT1eepo("poles.tex");
    	wT1eepo.write(poles);*/
    	//***********************************************************************
    	//Diffusion pole insula longitude****************************************
    	//***********************************************************************
    	
    	std::cout << "Poles parameterization..." << std::endl;
    	//****declarations et init textures
    	
    	TimeTexture<float>  poleNew(1,size);    //nouveaux poles
    	TimeTexture<float>  bord(1,size);    //bord du pole sud insula
    	TimeTexture<float>  bordValues(1,size);     //bords avec les valeurs obtenues lors de la diff
    	TimeTexture<float>  new_contraint(1,size);    //resultat diffusion sans les poles
    	TimeTexture<float>  bordExtreme(1,size);    //bords avec les extremites du meridien origine
    	TimeTexture<float>  bordSides(1,size);    //bords avec les valeurs sides
    	
    	TimeTexture<float>  new_contraint_insula(1,size);    //resultat diffusion sans les poles
    	
    	TimeTexture<float>  bordCing(1,size);    //bord du pole sud insula
    	TimeTexture<float>  bordValuesCing(1,size);     //bords avec les valeurs obtenues lors de la diff
    	TimeTexture<float>  bordExtremeCing(1,size);    //bords avec les extremites du meridien origine
    	TimeTexture<float>  bordSidesCing(1,size);    //bords avec les valeurs sides
    	
    	TimeTexture<float>  poleSave(1,size);    //sauvegarde des regions polaires
    	TimeTexture<float>  newOriginMeridian(1,size);    //new meridien d'origine
    	TimeTexture<float>  new_contraint_long_side(1,size);   //sides du nouveau meridien
    	TimeTexture<float>  new_long_insula(1,size);   //sides du nouveau meridien
    	
    	TimeTexture<float>  newOriginMeridianCing(1,size);    //new meridien d'origine
    	TimeTexture<float>  new_contraint_long_sideCing(1,size);   //sides du nouveau meridien
    	TimeTexture<float>  new_long_cing(1,size);   //sides du nouveau meridien
    	TimeTexture<float>  new_contraintCing(1,size);    //resultat diffusion sans les poles
    	
    	init_texture_single(bord);
    	init_texture_single(bordValues);
    	init_texture_single(bordExtreme);
    	init_texture_single(bordSides);
    	init_texture_single(poleNew);
    	init_texture_single(poleSave);
    	init_texture_single(newOriginMeridian);
    	init_texture_single(new_contraint_long_side);
    	init_texture_single(new_contraint);
    	init_texture_single(new_long_insula);
    	init_texture_single(bordCing);
    	init_texture_single(bordValuesCing);
    	init_texture_single(bordExtremeCing);
    	init_texture_single(bordSidesCing);
    	init_texture_single(newOriginMeridianCing);
    	init_texture_single(new_contraint_long_sideCing);
    	init_texture_single(new_long_cing);
    	init_texture_single(new_contraintCing);
    	
    	//OK!!
    	
    // 	Reader<Texture1d> wT1r("lon_old.tex");
    	Reader<Texture1d> wT1r(adr_lon);
    	wT1r.read(diff_long_result);
    	
    	//****mise a jour du nouveau point pole
    	cingularPoint();
    	int sud_temp=0;
    	sud_temp=sud;
    	sud=ind_min;
    	
    	TimeTexture<float> ttt(1,size);
    	init_texture_single(ttt);
    	int cpt=0,cpt1=0,cpt2=0,cpt3=0;
    	std::set<uint>::iterator it;
    	
    	
    	
    	//Sauvegarde des calottes polaires
    	for(int i=0; i<size; i++)
    		poleSave[0].item(i) = poles[0].item(i);
    	
    /*	Writer<Texture1d> wT3lpst("poles_origin.tex");
    	wT3lpst.write(poleSave);*/
    	
    	for(int i=0; i<size; i++)
    	{
    		if(poles[0].item(i)==0)
    		{
    			cpt++;
    			ttt[0].item(i)=0;
    			if(diff_long_result[0].item(i)==360)
    			{
    				//ttt[0].item(i)=100;
    				it=neigh[i].begin();
    				for(; it!=neigh[i].end(); it++)
    				{
    					if( diff_long_result[0].item(*it)!=360 )
    					{
    						diff_long_result[0].item(i)=diff_long_result[0].item(*it);
    					}
    				}
    			}
    		}
    		else
    		{
    			if(poles[0].item(i)==180)
    			{
    				cpt1++;
    				ttt[0].item(i)=30;
    			}
    			else
    			{
    				if(poles[0].item(i)==1)
    				{
    					cpt2++;
    					ttt[0].item(i)=10;
    				}
    				else
    				{
    					cpt3++;
    					ttt[0].item(i)=20;
    				}
    			}
    		}
    	}
    	//Second tour pour éliminer les trucs bizarres
    	for(int i=0; i<size; i++)
    	{
    		if(poles[0].item(i)==0)
    		{
    			cpt++;
    			ttt[0].item(i)=0;
    			if(diff_long_result[0].item(i)==360)
    			{
    				//ttt[0].item(i)=100;
    				it=neigh[i].begin();
    				for(; it!=neigh[i].end(); it++)
    				{
    					if( diff_long_result[0].item(*it)!=360 )
    					{
    						diff_long_result[0].item(i)=diff_long_result[0].item(*it);
    					}
    				}
    			}
    		}
    	}
    	
    	TimeTexture<float> ty(1,size);
    	init_texture_single(ty);
    	for(int i=0; i<size; i++)
    	{
    		//Pour l'insula
    		if( poles[0].item(i)==180 )
    		{
    			it=neigh[i].begin();
    			for(; it!=neigh[i].end(); it++)
    			{
    				if( poles[0].item(*it)!=180 )
    				{
    					bord[0].item(*it) = 1;
    					//bordValues[0].item(i) = new_contraint[0].item(i);
    					bordValues[0].item(*it) = diff_long_result[0].item(*it);
    				}
    			}
    		}
    		
    		//pour le cingulaire
    		
    		//if( poles[0].item(i)>0 && poles[0].item(i)<1.5 )
    		if( poles[0].item(i)==1 )
    		{
    			it=neigh[i].begin();
    			for(; it!=neigh[i].end(); it++)
    			{
    				//if(poles[0].item(*it)!=1 && poles[0].item(*it)<359)
    				if(poles[0].item(*it)==0)
    				{
    					bordCing[0].item(*it) = 1;
    					//bordValuesCing[0].item(i) = new_contraint[0].item(i);
    					bordValuesCing[0].item(*it) = diff_long_result[0].item(*it);
    					//ttt[0].item(*it)=180;
    				}
    			}
    		}
    		
    		if( bordValuesCing[0].item(i) != 0 || poleSave[0].item(i)==1 )
    		{
    			new_contraint_long_sideCing[0].item(i) = constraint_long_side_long[0].item(i);
    		}
    		else
    		{
    			new_contraintCing[0].item(i) = diff_long_result[0].item(i);
    			ty[0].item(i)=10;
    		}
    		
    		//if(  bordCing[0].item(i) == 1  )
    		if(  bordValuesCing[0].item(i) != 0  )
    			new_contraintCing[0].item(i) = bordValuesCing[0].item(i);
    	}
    	
    /*	Writer<Texture1d> wT1rbbvf("bordMer.tex");
    	wT1rbbvf.write(ty);*/
    	
    // 	std::cout << "cpt0="<<cpt<<" - cpt1="<<cpt2<<" - cpt180="<<cpt1<<"cpt autres="<<cpt3<<" - total="<<cpt+cpt1+cpt2+cpt3 << std::endl;
    	
    /*	Writer<Texture1d> wT1rbtbv("ttt.tex");
    	wT1rbtbv.write(ttt);*/
    	
    	
    	//meridien origine pour le cingulaire
    	for(int i=0; i<size; i++)
    	{
    		if( bordValuesCing[0].item(i) != 0 || poleSave[0].item(i)==1 )
    		{
    			int indTemp1=0,indTemp2=0;
    			it=neigh[i].begin();
    			for(; it!=neigh[i].end(); it++)
    			{
    				if(constraint_long_side_long[0].item(*it)==2)
    					indTemp1++;
    				if(constraint_long_side_long[0].item(*it)==4)
    					indTemp2++;
    			}
    			if(indTemp1!=0 && indTemp2!=0 && save_mer_origin[0].item(i)!=0)
    			{
    				new_contraintCing[0].item(i)=360;
    				temp[0].item(i)=360;
    			}
    		}
    	}
    	new_contraintCing[0].item(nord)=360;
    	
    	
    /*	Writer<Texture1d> w_bord("tmp_avant_avant.tex");
    	w_bord.write(temp);*/
    			
    	//Writer<Texture1d> wT1rbb("bord.tex");
    	//wT1rbb.write(bord);
    	
    /*	Writer<Texture1d> wT1rbbv("bordValuesCing.tex");
    	wT1rbbv.write(bordValuesCing);
    	
    	Writer<Texture1d> w_bord("bord.tex");
    	w_bord.write(bord);*/
    	
    	int ind_bord_insula=0;
    	int cpttt=0;
    	
    	for(int i=0; i<size; i++)
    	{
    		if(bord[0].item(i)!=0)
    		{
    			if(constraint_long_side[0].item(i)!=0 )
    			{
    				bordSides[0].item(i)=constraint_long_side[0].item(i);
    				cpttt++;
    			}
    		}
    	}
    	
    // 	std::cout << "cpttt="<<cpttt << std::endl;
    	
    	for(int i=0; i<size; i++)
    	{
    		if(bord[0].item(i)!=0)
    		{
    			int indTemp1=0,indTemp2=0;
    			it=neigh[i].begin();
    			for(; it!=neigh[i].end(); it++)
    			{
    // 				if(constraint_long_side[0].item(*it)==2)
    				if( (bord[0].item(*it)!=0) && (constraint_long_side[0].item(*it)==2) )
    					indTemp1++;
    // 				if(constraint_long_side[0].item(*it)==4)
    				if( (bord[0].item(*it)!=0) && (constraint_long_side[0].item(*it)==4) )
    					indTemp2++;
    			}
    			if(indTemp1!=0 && indTemp2!=0)
    			{
    				if( (indTemp1+indTemp2)>=cpttt)
    				{
    					bordExtreme[0].item(i)=180;
    					ind_bord_insula=i;
    				}
    			}
    		}
    	}
    	
    	std::cout << "ind="<<ind_bord_insula<<" - Save_Point_Insula="<<Save_Point_Insula << std::endl;
    /*	if(ind_bord_insula==0)
    		std::cout << "PROBLEM!!!!" << std::endl;*/
    	
    	//Test
    	if(Save_Point_Insula==-1)
    	{}
    	else
    	{
    		ind_bord_insula=Save_Point_Insula;
    	}
    	
    	std::cout << "ind new="<<ind_bord_insula << std::endl;
    	
    // 	Writer<Texture1d> w_constraint_long_side("constraint_long_side.tex");
    // 	w_constraint_long_side.write(constraint_long_side);
    // 	Writer<Texture1d> w_bord_values("bord_values.tex");
    // 	w_bord_values.write(bordValues);
    // 	Writer<Texture1d> w_bord("bord.tex");
    // 	w_bord.write(bord);
    // 	Writer<Texture1d> w_bord_side("bord_side.tex");
    // 	w_bord_side.write(bordSides);
    // 	Writer<Texture1d> w_bord_constraint_long_side("bord_constraint_long_side.tex");
    // 	w_bord_constraint_long_side.write(constraint_long_side);
    	
    	newOriginMeridian[0].item(sud)=360;
    	newOriginMeridian[0].item(ind_bord_insula)=360;
    	
    /*	Writer<Texture1d> w_bord_insula_point("bord_insula_point.tex");
    	w_bord_insula_point.write(newOriginMeridian);*/
    	
    	findNearNeigh(ind_bord_insula, sud, newOriginMeridian, 360, mesh, neigh );
    	
    // 	Writer<Texture1d> w_bord_insula_point_after("meridien_origine_insula.tex");
    // 	w_bord_insula_point_after.write(newOriginMeridian);
    	
    	new_contraint_long_side=originNeighbourgs(newOriginMeridian, ind_bord_insula, sud, mesh, neigh, poles);
    	
    	
    	new_contraint_long_side[0].item(sud)=0;
    	new_contraint_long_side[0].item(ind_bord_insula)=0;
    	
    // 	Writer<Texture1d> wT1ff("new_contraint_lon_side_AVANT.tex");
    // 	wT1ff.write(new_contraint_long_side);
    	
    	TimeTexture<float> bordtmp(1,size);
    	init_texture_single(bordtmp);
    	
    	//Rajoute des voisins au point du pole, pour décider des côtés
    	
    	for(int i=0; i<size; i++)
    	{
    		if( new_contraint_long_side[0].item(i)==0 && i!=ind_bord_insula )
    		{
    			float va=0;
    			
    			std::set<uint>::iterator its;
    			its=neigh[i].begin();
    			int ind1=0, ind2=0;
    			for(; its!=neigh[i].end(); its++)
    			{
    				if( new_contraint_long_side[0].item(*its)!=0 )
    				{
    					va=new_contraint_long_side[0].item(*its);
    					ind1++;
    				}
    				if( (*its)==(uint)ind_bord_insula )
    					ind2++;
    			}
    			if( ind1!=0 && ind2!=0 )
    				bordtmp[0].item(i)=va;
    		}
    	}
    	
    	
    	for(int i=0; i<size; i++)
    	{
    		if( bordtmp[0].item(i)!=0 )
    			new_contraint_long_side[0].item(i)=bordtmp[0].item(i);
    	}
    	
    	for(int i=0; i<size; i++)
    	{
    		if( new_contraint_long_side[0].item(i)==0 && i!=ind_bord_insula )
    		{
    			float va=0;
    			
    			std::set<uint>::iterator its;
    			its=neigh[i].begin();
    			int ind1=0, ind2=0;
    			for(; its!=neigh[i].end(); its++)
    			{
    				if( new_contraint_long_side[0].item(*its)!=0 )
    				{
    					va=new_contraint_long_side[0].item(*its);
    					ind1++;
    				}
    				if( (*its)==(uint)ind_bord_insula )
    					ind2++;
    			}
    			if( ind1!=0 && ind2!=0 )
    				bordtmp[0].item(i)=va;
    		}
    	}
    	
    /*	Writer<Texture1d> wT1btmp("bordtmp.tex");
    	wT1btmp.write(bordtmp);*/
    	for(int i=0; i<size; i++)
    	{
    		if( bordtmp[0].item(i)!=0 )
    			new_contraint_long_side[0].item(i)=bordtmp[0].item(i);
    	}
    	
    	
    	TimeTexture<float> bordtmp1(1,size);
    	init_texture_single(bordtmp1);
    	
    	for(int i=0; i<size; i++)
    	{
    		if( constraint_long_side[0].item(i)==0 && i!=ind_bord_insula )
    		{
    			float va=0;
    			
    			std::set<uint>::iterator its;
    			its=neigh[i].begin();
    			int ind1=0, ind2=0;
    			for(; its!=neigh[i].end(); its++)
    			{
    				if( constraint_long_side[0].item(*its)!=0 )
    				{
    					va=constraint_long_side[0].item(*its);
    					ind1++;
    				}
    				if( (*its)==(uint)ind_bord_insula )
    					ind2++;
    			}
    			if( ind1!=0 && ind2!=0 )
    				bordtmp1[0].item(i)=va;
    		}
    	}
    	for(int i=0; i<size; i++)
    	{
    		if( bordtmp1[0].item(i)!=0 )
    			constraint_long_side[0].item(i)=bordtmp1[0].item(i);
    	}
    	
    	for(int i=0; i<size; i++)
    	{
    		if( constraint_long_side[0].item(i)==0 && i!=ind_bord_insula )
    		{
    			float va=0;
    			
    			std::set<uint>::iterator its;
    			its=neigh[i].begin();
    			int ind1=0, ind2=0;
    			for(; its!=neigh[i].end(); its++)
    			{
    				if( constraint_long_side[0].item(*its)!=0 )
    				{
    					va=constraint_long_side[0].item(*its);
    					ind1++;
    				}
    				if( (*its)==(uint)ind_bord_insula )
    					ind2++;
    			}
    			if( ind1!=0 && ind2!=0 )
    				bordtmp1[0].item(i)=va;
    		}
    	}
    /*	Writer<Texture1d> wT1btmp1("bordtmp1.tex");
    	wT1btmp1.write(bordtmp1);*/
    	
    	for(int i=0; i<size; i++)
    	{
    		if( bordtmp1[0].item(i)!=0 )
    			constraint_long_side[0].item(i)=bordtmp1[0].item(i);
    	}
    	
    /*	Writer<Texture1d> wT1resinsntaa("contraint_long_side_AFTER.tex");
    	wT1resinsntaa.write(constraint_long_side);*/
    	
    	
    	int ind_change=0;
    	for(int i=0; i<size; i++)
    	{
    		if( bordSides[0].item(i)!=0 && new_contraint_long_side[0].item(i)!=0)
    		{
    			//std::cout << "IN 1" << std::endl;
    			if( (constraint_long_side[0].item(i)!= new_contraint_long_side[0].item(i)) && (new_contraint_long_side[0].item(i)!=0) )
    			{
    				ind_change++;
    				//std::cout << "CHANGE!!" << std::endl;
    			}
    		}
    	}
    	
    	TimeTexture<float> sides_temp(1,size);
    	
    	//Mise a jour de sides selon comment c'est dans la diff totale
    	if( ind_change != 0 )
    	{
    		for(int i=0; i<size; i++)
    		{
    			if( new_contraint_long_side[0].item(i)==2 )
    				new_contraint_long_side[0].item(i)=4;
    			else
    				if( new_contraint_long_side[0].item(i)== 4)
    					new_contraint_long_side[0].item(i)=2;
    				else{}
    		}
    	}
    	for(int i=0; i<size; i++)
    	{
    		sides_temp[0].item(i)=new_contraint_long_side[0].item(i);
    	}
    	
    /*	Writer<Texture1d> wSides("sides_avant.tex");
    	wSides.write(sides_temp);
    	
    	Writer<Texture1d> wSidesNew("bordValues.tex");
    	wSidesNew.write(bordValues);*/
    	
    	new_contraint_long_side = defineSidesPoles( sides_temp, bordValues, mesh, neigh );
    	
    	for(int i=0; i<size; i++)
    	{
    		new_contraint[0].item(i)=bordValues[0].item(i);
    		if( poleSave[0].item(i)==180)
    		{
    			poles[0].item(i)=0;
    			if( bord[0].item(i)==0)
    				new_contraint[0].item(i)=0;
    		}
    		if( poleSave[0].item(i)!=180)
    		{
    			if( bord[0].item(i)==0)
    			{
    				poles[0].item(i)=1;
    				new_contraint[0].item(i)=360;
    			}
    			if( bord[0].item(i)!=0)
    			{
    				poles[0].item(i)=0;
    			}
    			
    		}
    		if( newOriginMeridian[0].item(i)!=0 )
    		{
    			new_contraint[0].item(i)=360;
    			temp[0].item(i)=360;
    		}
    		//if(bordSides[0].item(i)!=0)
    		//	new_contraint_long_side[0].item(i)=bordSides[0].item(i);
    	}
    	
    	
    /*	Writer<Texture1d> wT1kdakc("tmp_avant.tex");
    	wT1kdakc.write(temp);*/
    	
    	new_contraint[0].item(sud)=360;
    	poles[0].item(sud)=180;
    	temp[0].item(sud)=360;
    	temp[0].item(nord)=360;
    	
    	
    /*	Writer<Texture1d> wT1kdkc("tmp.tex");
    	wT1kdkc.write(temp);*/
    	
    	
    /*	Writer<Texture1d> wT1resinsnt("new_contraint_lon_side.tex");
    	wT1resinsnt.write(new_contraint_long_side);*/
    /*	Writer<Texture1d> wT1resinsnn("new_contraint_lon.tex");
    	wT1resinsnn.write(new_contraint);*/
    /*	Writer<Texture1d> wT1resinsfr("new_contraintlongitude.tex");
    	wT1resinsfr.write(new_contraint);*/
    	//diff pole insula*************************************************
    	// A REMETTRE!!!!!!!!!!!!!!!!!!!!!!!!!
    	new_long_insula=diffusionLongitude(new_contraint, new_contraint_long_side, poleSave, 0);
    	
    /*	Writer<Texture1d> wT1resins("new_long_insula.tex");
    	wT1resins.write(new_long_insula);*/
    	
    	//*********************************************************************
    	
    	//poles pour le cingulaire
    	for(int i=0; i<size; i++)
    	{
    		poles[0].item(i)=0;
    		if( bordCing[0].item(i) == 1 || poleSave[0].item(i)==1  )
    		{}
    		else
    		{
    			poles[0].item(i)=180;
    		}
    	}
    	poles[0].item(nord)=1;
    // 	Writer<Texture1d> wT1eepeo("poles_cingulaire.tex");
    // 	wT1eepeo.write(poles);
    // 	
    // 	TimeTexture<float> t(1,size);
    	ind_change=0;
    	for(int i=0; i<size; i++)
    	{
    //		t[0].item(i)=0;
    		if( constraint_long_side_long[0].item(i)!=0 && poleSave[0].item(i)==0 )
    		{
    			if( constraint_long_side_long[0].item(i)!= constraint_long_side[0].item(i) )
    			{
    				ind_change++;
    // 				std::cout << constraint_long_side_long[0].item(i)<<" - "<<constraint_long_side[0].item(i) << std::endl;
    // 				t[0].item(i)=constraint_long_side_long[0].item(i);
    			}
    		}
    	}
    
    /*	Writer<Texture1d> wTt("t.tex");
    	wTt.write(t);*/
    	
    // 	std::cout << "nombre de differences: "<<ind_change << std::endl;
    	//Mise a jour de sides selon comment c'est dans la diff totale
    // 	if( ind_change != 0 )
    // 	{
    // 		//std::cout << "changement sides" << std::endl;
    // 		for(int i=0; i<size; i++)
    // 		{
    // 			if( new_contraint_long_sideCing[0].item(i)==2 )
    // 				new_contraint_long_sideCing[0].item(i)=4;
    // 			else
    // 				if( new_contraint_long_sideCing[0].item(i)== 4)
    // 					new_contraint_long_sideCing[0].item(i)=2;
    // 		}
    // 	}
    	//Mise a jour de sides selon comment c'est dans la diff totale
    	if( ind_change != 0 )
    	{
    		//std::cout << "changement sides" << std::endl;
    		for(int i=0; i<size; i++)
    		{
    			if( constraint_long_side_long[0].item(i)==2 )
    				constraint_long_side_long[0].item(i)=4;
    			else
    				if( constraint_long_side_long[0].item(i)== 4)
    					constraint_long_side_long[0].item(i)=2;
    		}
    	}
    
    /*	Writer<Texture1d> wT1rescingor("constraint_long_side.tex");
    	wT1rescingor.write(constraint_long_side);*/
    // 	Writer<Texture1d> wT1rescingc("new_contraintCing.tex");
    // 	wT1rescingc.write(new_contraintCing);
    /*	Writer<Texture1d> wT1rescingff("new_contraint_long_side_long.tex");
    	wT1rescingff.write(constraint_long_side_long);*/
    // 	Writer<Texture1d> wT1rescingfq("poles_diff_cing.tex");
    // 	wT1rescingfq.write(poles);
    	//Diff pole cing*************************************************************
    //	new_long_cing=diffusionLongitude(new_contraintCing, new_contraint_long_sideCing, poleSave, 1);
    	new_long_cing=diffusionLongitude(new_contraintCing, constraint_long_side_long, poleSave, 1);
    
    /*	Writer<Texture1d> wT1rescing("new_long_cing.tex");
    	wT1rescing.write(new_long_cing);*/
    	
    	//******************************************************************************
    
    /*	it=neigh[sud].begin();
    	for(; it!=neigh[sud].end(); it++)
    	{
    	poles_save[0].item(*it)=180;
    }
    */
    
    	for(int i=0; i<size; i++)
    	{
    		poles[0].item(i) = poleSave[0].item(i);
    	}
    	
    	
    	TimeTexture<float> diff_total(1,size);
    	init_texture_single(diff_total);
    	
    	for(int i=0;i<size; i++)
    	{
    		if( poles[0].item(i) == 180 )
    			diff_total[0].item(i) = new_long_insula[0].item(i);
    		else
    			if( poles[0].item(i) == 1 )
    				diff_total[0].item(i) = new_long_cing[0].item(i);
    		else
    		{
    			diff_total[0].item(i) = diff_long_result[0].item(i);
    				//std::cout << "YEAH";
    		}
    		if(temp[0].item(i)>359)
    			diff_total[0].item(i) = 360;
    	}
    	diff_total[0].item(nord)=360;
    	diff_total[0].item(sud)=360;
    	
    	Writer<Texture1d> wcipfp(adr_lon);
    	wcipfp.write(diff_total);
    	
    	sud=sud_temp;
    	
    	std::cout << "done" << std::endl;
    	
    }
    
    
    
    
    /****************DIFFUSION METHODS****************/
    
    TimeTexture<float> CorticalReferential::diffusionLatitudeRelax( TimeTexture<float> & tex ) 
    {
    
    	TimeTexture<float>  result(1,size);
    	Texture<float>     smooth(size), lapl;
    	float			s;
        int maxIter=90000;
    
    	//float Beta=0.5;
    
    	std::vector<unsigned>::iterator itlist;
    
    
    	TimeTexture<float> lat_forbid(1,size);
    	TimeTexture<float> contraint(1,size);
    	TimeTexture<float> init1(1,size);
    	TimeTexture<float> init2(1,size);
    	TimeTexture<float> init(1,size);
    	TimeTexture<float> temp(1,size);
    	TimeTexture<float> cInit(1,size);
    	TimeTexture<float> betaInit(1,size);
    	TimeTexture<float> betaMap(1,size);
    	TimeTexture<float> cMap(1,size);
    	TimeTexture<short> pole1(1,size);
    	TimeTexture<short> pole2(1,size);
    	init_texture_single(lat_forbid);
    	init_texture_single(contraint);
    	init_texture_single(init1);
    	init_texture_single(init2);
    	init_texture_single(init);
    	init_texture_single(betaMap);
    	init_texture_single(betaInit);
    	init_texture_single(cMap);
    	init_texture_single(cInit);
    // 	init_texture_single(pole1);
    // 	init_texture_single(pole2);
    	
    	
    	//place la source dans une texture
         for(int i=0;i<size; i++)
         {
              pole1[0].item(i)=0;
              pole2[0].item(i)=0;
         }
         
    	pole1[0].item(nord)=100;
    	pole2[0].item(nord)=-1;
    	
    	pole1[0].item(sud)=-1;
    	pole2[0].item(sud)=100;
    	
    	/*A REMETTRE!!*/
    	std::cout << "meshdistance1" << std::endl;
    	init1[0]=meshdistance::MeshDistance( mesh[0], pole1[0], false );
    	init2[0]=meshdistance::MeshDistance( mesh[0], pole2[0], false );
    	std::cout << "meshdistance1 OK" << std::endl;
    	float max_init1=0;
    	float max_init2=0;
    	for(int i=0;i<size; i++)
    	{
              if ((i!=nord) && (i!=sud))
              {
                   if(init1[0].item(i)>max_init1)
                        max_init1=init1[0].item(i);
                   if(init2[0].item(i)>max_init2)
                        max_init2=init2[0].item(i);
              }
    	}
     	
    	for(int i=0;i<size; i++)
    	{
              if ((i!=nord) && (i!=sud))
              {
                   init1[0].item(i) = (init1[0].item(i)*180)/max_init1  ;
                   init2[0].item(i) = (init2[0].item(i)*180)/max_init2  ;
                   init2[0].item(i) = sqrt( pow( (init2[0].item(i)-180), 2));
    // 		init[0].item(i) = 0  ;
              }
    	}
    
    	
    	for(int i=0;i<size; i++)
    	{
              if ((i!=nord) && (i!=sud))
              {
                   init[0].item(i)=(init1[0].item(i)+init2[0].item(i))/2;
    		/*if(init1[0].item(i)>180)
    		init[0].item(i)=init2[0].item(i);*/
              }
    	}
    	
    	float max_init=0;
    	float min_init=50;
    	for(int i=0;i<size; i++)
    	{
              if ((i!=nord) && (i!=sud))
              {
                   if(init[0].item(i)>max_init)
                        max_init=init[0].item(i);
                   if(init[0].item(i)<min_init)
                        min_init=init[0].item(i);
              }
    	}
    	for(int i=0;i<size; i++)
    	{
              if ((i!=nord) && (i!=sud))
                   init[0].item(i) =( ( (init[0].item(i)-min_init)/(max_init-min_init) ) * 180 ) +1 ;
    	}
    	
    	init[0].item(nord) = 1;
    	init[0].item(sud) = 181;
    	
    /*	Writer< Texture1d > wTiricpf("initLat.tex");
    	wTiricpf.write(init);*/
    	
    	std::set<uint>::const_iterator it;
    
    	for(int i=0;i<size; i++)
    	{
    		//smooth = diffusion & tex = contraintes
    		smooth.item(i) = init[0].item(i);
    
    		if(tex[0].item(i) != 0)
    		{
    			contraint[0].item(i)=tex[0].item(i);
    			smooth.item(i)=tex[0].item(i);
    		}
    		
    		if( cercle_polaire[0].item(i)!=0 )
    			lat_forbid[0].item(i)=1;
    
    		//cree la texture sans les bords pour la condition d'arret
    		
    // 		it=neigh[i].begin();
    // 		for(; it!=neigh[i].end(); it++)
    // 			if( poles[0].item(*it)!=0 )
    // 				lat_forbid[0].item(*it)=1;
    		
    	}
    	
    	//std::cout << "ok1" << std::endl;
    	for(int i=0;i<size; i++)
    	{
    		if( (tex[0].item(i)!=0) && (lat_forbid[0].item(i)==0) )
    		{
    			cInit[0].item(i)=tex[0].item(i);
    			betaInit[0].item(i)=_Beta;
    		}
    		if(cInit[0].item(i)!=0)
    		{
    			smooth.item(i)=cInit[0].item(i);
    			contraint[0].item(i)=cInit[0].item(i);
    		}
    	}
    	
    	//std::cout << "ok2" << std::endl;
    	cInit[0].item(nord)=0;
    	betaInit[0].item(nord)=0;
    	cInit[0].item(sud)=0;
    	betaInit[0].item(sud)=0;
    	
    	
    	if(typeBeta==0)
    	{
    		betaMap=dilate_texture(betaInit,0.5, neigh, mesh);
    		//std::cout << "ok3" << std::endl;
    		cMap=dilate_texture(cInit,1, neigh, mesh);
    		//std::cout << "ok4" << std::endl;
    	}
    	else
    	{
    		for(int i=0;i<size;i++)
    		{
    			betaMap[0].item(i)=betaInit[0].item(i);
    			cMap[0].item(i)=cInit[0].item(i);
    		}
    	}
    	
    	//A VIRER!!
    	smooth.item(nord) = 1;
    	smooth.item(sud) = 181;
    	//STOP!!
    	lat_forbid[0].item(nord)=1;
    	lat_forbid[0].item(sud)=1;
    
    // 	Writer<Texture1d> wbetainit("betaMap.tex");
    // 	wbetainit.write(betaMap);
    // 	Writer<Texture1d> wcinit("cMap.tex");
    // 	wcinit.write(cMap);
    // 	Writer< Texture1d > wTy1c("lat_forbid_lat.tex");
    // 	wTy1c.write(lat_forbid);
    
    	TimeTexture<float> diff_laplacian(1,size);
    // 	TimeTexture<float> diff_temp(1,size);
    // 	TimeTexture<float> diff_multi(1000,size);
    	//init_texture_single(
    	//int cp=0;
    				
    	float max=1;
    	//float moy=0;
    	float moy_previous=0;
    	std::cout << "Processing, avec _dt=" << _dt << std::endl;
    	std::cout << "max=" << max << std::endl;
    
    	//DIFFUSION PROCESS
    
    	int iter=0;
    // 	int cpt_multi=0;
    	do
    	{
    		lapl =  AimsMeshLaplacian(smooth, weightLapl);
    /*		if( iter<100)
    		{
    			for ( int i=0; i<size; i++)
    			{
    				diff_multi[cpt_multi].item(i)=smooth.item(i);
    			}
    			cpt_multi++;
    		}*/
    		if ( (iter%100)==0  && (iter>1) )
    		{
    			//cp++;
    			float moy=0;
    			
    			//Calcul de la condition d'arret
    			int cpt=0;
    			for (int i=0; i<size; i++)
    			{
    				if ( lat_forbid[0].item(i)==0 )
    				{
    					moy+=fabs( lapl.item(i) );
    					cpt++;
    				}
    			}
    			moy=moy/cpt;
    			max=fabs(moy-moy_previous);
    			moy_previous=moy;
    			
    			std::cout << "iter " << iter << " avec max = " << max  << std::endl;
    /*			if(cpt_multi<1000)
    			{
    				for ( int i=0; i<size; i++)
    				{
    					diff_multi[cpt_multi].item(i)=smooth.item(i);
    				}
    				Writer< Texture1d > wT1m("diff_lat_multi.tex");
    				wT1m.write(diff_multi);
    				cpt_multi++;
    			}*/
    /*			for ( int i=0; i<size; i++)
    			{
    				diff_temp[0].item(i)=smooth.item(i);
    			}
    			Writer< Texture1d > wT1c(adr_lat);
    			wT1c.write(diff_temp);*/
    			
    		}
    
    		
    		//Mise a jour des textures
    		float epsilon=0.02;
    		for (int i=0; i<size; ++i)
    		{
    			if( fabs(cMap[0].item(i))>epsilon ) 
    			{
    				temp[0].item(i) = smooth.item(i) - cMap[0].item(i);
    			}
    			else
    				temp[0].item(i) = 0;
    			
    			if ( lat_forbid.item(i)==0 )
    			{
    				s = smooth.item(i) + _dt * ( lapl.item(i) - betaMap[0].item(i)*temp[0].item(i) );
    				smooth.item(i) = s;
    			}
    		}
    		iter++;
    	}
    	while((max > criterium)&&(iter<maxIter));
    // 	while(iter<300000);
    	
    	std::cout  << std::endl;
    	//Le resultat est stocke dans une TimeTexture
    	for ( int i=0; i<size; ++i)
    	{
    		result[0].item(i)=smooth.item(i);
    	}
    //	diff_laplacian[0]=lapl;
    	//std::cout << "Max : " << max << std::endl;
    
    // 	std::string address_lapl="./diff_laplacian.tex";
    // 	Writer<Texture1d> w_lapl(address_lapl);
    // 	if( !w_lapl.write( diff_laplacian ) ) {}
    
    	return(result);
    }
    
    
    
    TimeTexture<float> CorticalReferential::diffusionLatitude( TimeTexture<float> & tex ) 
    {
    
    	TimeTexture<float>  result(1,size);
    	Texture<float>     smooth(size), lapl;
    	float			s;
    
    	std::vector<unsigned>::iterator itlist;
    
    
    	
    	TimeTexture<float> lat_forbid(1,size);
    	init_texture_single(lat_forbid);
    	
    	
    	//place la source dans une texture
    
    	std::set<uint>::const_iterator it;
    	
    	for(int i=0;i<size; i++)
    	{
    		//smooth = diffusion & tex = contraintes
    		smooth.item(i) = 90;
    		if(tex[0].item(i) != 0)
    			smooth.item(i) = tex[0].item(i);
    		
    		//cree la texture sans les bords pour la condition d'arret
    		if( poles[0].item(i)!=0 || tex[0].item(i)!=0 )
    			lat_forbid[0].item(i)=1;
    		it=neigh[i].begin();
    		for(; it!=neigh[i].end(); it++)
    			if( poles[0].item(*it)!=0 || tex[0].item(*it)!=0 )
    				lat_forbid[0].item(*it)=1;
    	}
    
    
    	TimeTexture<float> diff_laplacian(1,size);
    	TimeTexture<float> diff_temp(1,size);
    // 	TimeTexture<float> diff_temp_big(1000,size);
    	//init_texture_single(
    	//int cp=0;
    				
    	float max=0;
    	//float moy=0;
    	float moy_previous=0;
    	std::cout << "Processing, avec _dt=" << _dt << std::endl;
    
    	//DIFFUSION PROCESS
    //	for (int iter=0; iter< t; ++iter)
    
    	int iter=0;
    // 	int cpt=0;
    	do
    	{
    		
    		lapl =  AimsMeshLaplacian(smooth, weightLapl);
    /*		if( iter<100)
    		{
    			for ( int i=0; i<size; i++)
    			{
    				diff_temp_big[cpt].item(i)=smooth.item(i);
    			}
    			cpt++;
    		}*/
    		if ( (iter%100)==0 )
    		{
    			//cp++;
    			float moy=0;
    // 			std::cout << " moy avant process="<< moy  << std::endl;
    
    			//Calcul de la condition d'arret
    			int cpt=0;
    			for (int i=0; i<size; i++)
    			{
    				if ( lat_forbid[0].item(i)==0 )
    				{
    					moy+=fabs( lapl.item(i) );
    					cpt++;
    				}
    			}
    // 			std::cout << " cpt="<< cpt  << std::endl;
    // 			std::cout << " moy_previous="<< moy_previous  << std::endl;
    			moy=moy/cpt;
    // 			std::cout << " moy="<< moy  << std::endl;
    			max=fabs(moy-moy_previous);
    			moy_previous=moy;
    			
    			std::cout << "iter " << iter << " avec max = " << max  << std::endl;
    /*			for ( int i=0; i<size; i++)
    			{
    			diff_temp[0].item(i)=smooth.item(i);
    		}
    			Writer< Texture1d > wT1c("diff_latitude_processing.tex");
    			wT1c.write(diff_temp);*/
    			
    /*			for ( int i=0; i<size; i++)
    			{
    				diff_temp_big[cpt].item(i)=smooth.item(i);
    			}*/
    // 			std::cout << "ECRIT!!!!" << std::endl;
    // 			Writer< Texture1d > wT1f("diff_latitude_big.tex");
    // 			wT1f.write(diff_temp_big);
    			cpt++;
    		}
    
    		
    		//Mise a jour des textures
    		for (int i=0; i<size; i++)
    		{
    			if ( tex.item(i)==0 )
    			{
    				s = smooth.item(i) + _dt * lapl.item(i);
    				smooth.item(i) = s;
    			}
    		}
    		
    		iter++;
    	}
    	while(max > criterium);
    	//while(iter<300000);
    	
    	std::cout  << std::endl;
    	//Le resultat est stocke dans une TimeTexture
    	for ( int i=0; i<size; ++i)
    	{
    		result[0].item(i)=smooth.item(i);
    	}
    //	diff_laplacian[0]=lapl;
    	//std::cout << "Max : " << max << std::endl;
    
    // 	std::string address_lapl="./diff_laplacian.tex";
    // 	Writer<Texture1d> w_lapl(address_lapl);
    // 	if( !w_lapl.write( diff_laplacian ) ) {}
    
    /*	Writer< Texture1d > wT1f("diff_latitude_big.tex");
    	wT1f.write(diff_temp_big);*/
    	
    	return(result);
    }
    
    TimeTexture<float> CorticalReferential::diffusionLongitudeRelax( TimeTexture<float> & tex, TimeTexture<float> & side) 
    {
    
    	TimeTexture<float>  result(1,size);
    	Texture<float>     smooth(size), lapl;
    	float			s;
         float epsi=0.0001;
         int maxIter=90000;
    
    	//Remove weights from poles
    	unsigned  uneigh;
    	float weight;	
    	TimeTexture<float> contraint(1,size);
    	TimeTexture<float> long_forbid(1,size);
    	TimeTexture<float> long_forbid_comp(1,size);
    	TimeTexture<float> temp(1,size);
    	TimeTexture<float> cInit(1,size);
    	TimeTexture<float> betaInit(1,size);
    	TimeTexture<float> betaMap(1,size);
    	TimeTexture<float> cMap(1,size);
    	TimeTexture<float> init1(1,size);
    	TimeTexture<float> init2(1,size);
    	TimeTexture<float> init(1,size);
    	TimeTexture<short> pole1(1,size);
    	TimeTexture<short> pole2(1,size);
    	init_texture_single(contraint);
    	init_texture_single(long_forbid);
    	init_texture_single(long_forbid_comp);
    	init_texture_single(temp);
    	init_texture_single(betaMap);
    	init_texture_single(betaInit);
    	init_texture_single(cMap);
    	init_texture_single(cInit);
    	init_texture_single(init1);
    	init_texture_single(init2);
    	init_texture_single(init);
    /*	init_texture_single(pole1);
    	init_texture_single(pole2);*/
    	
    	
    	//place la source dans une texture
    
    	for(int i=0;i<size; i++)
    	{
    		if(tex[0].item(i) == 360)
    		{
    			pole1[0].item(i)=(short) -1;
    			pole2[0].item(i)=(short)-1;
    		}
    		else
    		{
    			if (side[0].item(i) == 4)
    				pole1[0].item(i)=(short) 100;
    			else
    				pole1[0].item(i)=(short) 0;
    			
    			if (side[0].item(i) == 2)
    				pole2[0].item(i)=(short) 100;
    			else
    				pole2[0].item(i)=(short) 0;
    		}
    	}
    
    //      Writer< TimeTexture<float> > wT0("/home/olivier/t0_before.tex");
    //      wT0.write(tex);
    // 	Writer< TimeTexture<short> > wTp1("/home/olivier/pole1_lon.tex");
    // 	wTp1.write(pole1);
    // 	Writer< TimeTexture<short> > wTp2("/home/olivier/pole2_lon.tex");
    // 	wTp2.write(pole2);
    	
    	/*A REMETTRE!!*/
    	
    	init1[0]=meshdistance::MeshDistance( mesh[0], pole1[0], false );
    	init2[0]=meshdistance::MeshDistance( mesh[0], pole2[0], false );
    
    //      Writer< Texture1d > wTiric1("/home/olivier/init1_lon.tex");
    //      wTiric1.write(init1);
    //      Writer< Texture1d > wTiric2("/home/olivier/init2_lon.tex");
    //      wTiric2.write(init2);
         
    	float max_init1=0;
    	float min_init1=360;
    	float max_init2=0;
    	float min_init2=360;
    	for(int i=0;i<size; i++)
    	{
              if (pole1[0].item(i) > -1)
              {
                   if(init1[0].item(i)>max_init1)
                        max_init1=init1[0].item(i);
                   if(init1[0].item(i)<min_init1)
                        min_init1=init1[0].item(i);
    		
                   if(init2[0].item(i)>max_init2)
                        max_init2=init2[0].item(i);
                   if(init2[0].item(i)<min_init2)
                        min_init2=init2[0].item(i);
              }
    	}
    //      std::cerr << "distanceMap1 : min=" << min_init1 << ", max=" << max_init1 << std::endl;
    //      std::cerr << "distanceMap2 : min=" << min_init1 << ", max=" << max_init1 << std::endl;
    
    	
    	for(int i=0;i<size; i++)
    	{
              if (pole1[0].item(i) > -1)
              {
                   init1[0].item(i) = (init1[0].item(i)*360)/max_init1  ;
                   init2[0].item(i) = (init2[0].item(i)*360)/max_init2  ;
                   init2[0].item(i) = (float) sqrt( pow( (init2[0].item(i)-360), 2));
              }
         }
    	
    //      Writer< Texture1d > wTiricp1("/home/olivier/init1.1_lon.tex");
    //      wTiricp1.write(init1);
    //      Writer< Texture1d > wTiricp2("/home/olivier/init2.1_lon.tex");
    //      wTiricp2.write(init2);
    	
    	for(int i=0;i<size; i++)
    	{
    		init[0].item(i)=(init1[0].item(i)+init2[0].item(i))/2;
    		/*if(init1[0].item(i)>180)
    			init[0].item(i)=init2[0].item(i);*/
    	}
    
     /*    Writer< Texture1d > wTiricpp("/home/olivier/init1+2_lon.tex");
         wTiricpp.write(init)*/;
         
    	float max_init=0;
    	float min_init=50;
    	for(int i=0;i<size; i++)
    	{
    		if(init[0].item(i)>max_init)
    			max_init=init[0].item(i);
    		if(init[0].item(i)<min_init)
    			min_init=init[0].item(i);
    	}
    	for(int i=0;i<size; i++)
    	{
    		//init[0].item(i) = (init[0].item(i)*360)/max_init  ;
    		init[0].item(i) = ( (init[0].item(i)-min_init)/(max_init-min_init) ) * 360  ;
    	}
    
    
    // EDIT Olivier : la texture init contient l'initialisation de la carte avant diffusion.
    // Elle est construite à partir de cartes de distances
    // Ici j'introduis une modif qui la rescale par morceau pour déjà correspondre le mieux
    // possible aux contraintes. Les contraintes, et leurs valeurs, sont dans tex[0]
    
    
    //      Writer< Texture1d > wtex0("/home/olivier/init.tex");
    //      wtex0.write(init);
    //      std::cerr << "init written" << std::endl;
    
         std::map<float, std::vector<float> > mapCon;
         std::map<float, float> mapRescale;
         for (int i=0; i<size; i++)
         {
              if ((tex[0].item(i) >0) && (poles[0].item(i) < epsi) && (tex[0].item(i) < 360)) // a bit tricky here because we are handling floats (why ?)
              {
                   if (mapCon.find(tex[0].item(i)) == mapCon.end())
                   {
                        mapCon[tex[0].item(i)]=std::vector<float>();
                        mapCon[tex[0].item(i)].push_back(init[0].item(i));
                   }
                   else
                   {
                        mapCon[tex[0].item(i)].push_back(init[0].item(i));
                   }
              }
         }
    
    
         std::map<float, std::vector<float> >::iterator mapConIt;
         std::vector<float> vecCon;
         std::vector<float>::iterator vecConIt;
         float con, mean;
    
         for (mapConIt=mapCon.begin(); mapConIt!=mapCon.end(); ++mapConIt)
         {
              vecCon=(*mapConIt).second;
              con=(*mapConIt).first;
              mean=0;
              int count=0;
    
              for ( vecConIt=vecCon.begin(); vecConIt!=vecCon.end(); ++vecConIt)
              {
                   mean += (*vecConIt);
                   count ++;
              }
              mapRescale[con]=(float) (mean/float(count));
              std::cerr << "\t" << con << "->" << mapRescale[con] << std::endl;
    
         }
         mapRescale[360]=float(360);
    
         float bound1, bound2, con1, con2;
    
    
         for (int i=0; i<size; i++)
         {
              float val=init[0].item(i);
              std::map<float, float>::iterator rescaleIt=mapRescale.begin();
              bound1=0.0; con1=0.0;
              if ((val!=0) && (val !=360))
              {
                   for ( ; rescaleIt!=mapRescale.end(); ++rescaleIt)
                   {
                        con2=(*rescaleIt).first; bound2=(*rescaleIt).second;
                        if ((val>=bound1) && (val<=bound2))
                        {
                             init[0].item(i)=(((con2-con1)/float(bound2-bound1))*(val-bound1)) + con1;
                             con1=con2; bound1=bound2;
                        }
                        else
                        {
                             con1=con2; bound1=bound2;
                        }
                   }
              }
         }
    
    
    // EDIT Olivier : fin de la modification de l'initialisation
    
    
    //      Writer< Texture1d > wTiricp("/home/olivier/initFinal_lon.tex");
    //      wTiricp.write(init);
    
    
    
    	
    // 	TimeTexture<float> diff_multi(1000,size);
    	std::map<unsigned, std::set< std::pair<unsigned,float> > >  weightLapl1;
    	
    	weightLapl1=weightLapl;
    	
    	std::set< std::pair<unsigned,float> >::iterator isp,esp;
    	std::map<unsigned, std::set< std::pair<unsigned,float> > >::iterator  il, el;
    	for (il = weightLapl1.begin(), el = weightLapl1.end(); il != el; ++il)
    	{
    		for ( isp = (il->second).begin(), esp = (il->second).end(); isp != esp; ++isp)
    		{
    			uneigh = isp->first;
    			weight = isp->second;
    			//itlist=liste.begin();
    			if( poles[0].item(uneigh)!=0 )
    			{
    				il->second.erase(std::pair<unsigned,float>(uneigh,weight) );
    				il->second.insert(std::pair<unsigned,float>(uneigh,0) );
    			}
    			else {}
    		}
    	}
    
    	std::set<uint>::const_iterator it;
    
    	for(int i=0; i<size; i++)
    	{
    		//smooth = diffusion & tex = contraintes
    		smooth.item(i) = init[0].item(i);
    		if(tex[0].item(i) != 0)
    		{
    			contraint[0].item(i)=tex[0].item(i);
    			//if( poles[0].item(i)!=0 )
    			//smooth.item(i) = tex[0].item(i);
    		}
    
    		//cree la texture sans les bords pour la condition d'arret
    		if( poles[0].item(i)!=0 || tex[0].item(i)==360 )
    		{
    			long_forbid[0].item(i)=1;
    			smooth.item(i) = 360;
    		}
    		it=neigh[i].begin();
    
    		for(; it!=neigh[i].end(); it++)
    			if( poles[0].item(*it)!=0 || tex[0].item(*it)==360 )
    				long_forbid[0].item(*it)=1;
    
    		//Cree la texture compl�entaire
    		if(long_forbid[0].item(i)==1)
    			long_forbid_comp[0].item(i)=0;
    		else
    			long_forbid_comp[0].item(i)=1;
    	}
    
    // 	Writer< Texture1d > wTiic("tex.tex");
    // 	wTiic.write(tex);
    
    	for(int i=0;i<size; i++)
    	{
    		if( (tex[0].item(i)!=0) && (long_forbid[0].item(i)==0) )
    		{
    			cInit[0].item(i)=tex[0].item(i);
    			betaInit[0].item(i)=_Beta;
    		}
    		if(cInit[0].item(i)!=0)
    		{
    			smooth.item(i)=cInit[0].item(i);
    			contraint[0].item(i)=cInit[0].item(i);
    		}
    	}
    
    //      std::cerr << "DEBUG : Writing betaaaaaaa" << std::endl; 
    //      Writer< Texture1d > wbetInit("/home/olivier/initBeta.tex");
    //      wbetInit.write(betaInit);
    
    
    	if(typeBeta==0)
    	{
    		betaMap=dilate_texture(betaInit,0.5, neigh, mesh);
    		//std::cout << "ok3" << std::endl;
    		cMap=dilate_texture(cInit,1, neigh, mesh);
    		//std::cout << "ok4" << std::endl;
    	}
    	else
    	{
    		for(int i=0;i<size;i++)
    		{
    			betaMap[0].item(i)=betaInit[0].item(i);
    			cMap[0].item(i)=cInit[0].item(i);
    		}
    	}
    
    // 	Writer<Texture1d> wbetainit("/home/olivier/betaMap_lon.tex");
    // 	wbetainit.write(betaMap);
    // 	Writer<Texture1d> wcinit("/home/olivier/cMap_lon.tex");
    // 	wcinit.write(cMap);
    //      Writer<Texture1d> wtex0("/home/olivier/texZero_lon.tex");
    //      wtex0.write(tex);
    
         TimeTexture<float> t_smooth(1,size);
         t_smooth[0]=smooth;
    //      Writer<TimeTexture<float> > wsmooth("/home/olivier/smooth_lon.tex");
    //      wsmooth.write(t_smooth);
    //      std::cerr << "Wrote debug textures..." << std::endl;
    	
    // 	TimeTexture<float> diff_temp(1,size);
    
    	TimeTexture<float> diff_laplacian(1,size);
    	float moy_previous=0;
    	std::cout << "Processing" << std::endl;
    
    /*    TimeTexture<float> tt;*/
    	float max=10;
    	int iter=0;
    	int cpt_multi=0;
    	//for (int iter=0; iter< t; ++iter)
    	do
    	{
    		lapl =  AimsMeshLaplacian_meridian(smooth, tex, weightLapl1, side);
    		
    /*		if( iter<100 )
    		{
    			for ( int i=0; i<size; i++)
    			{
    				diff_multi[cpt_multi].item(i)=smooth.item(i);
    			}
    			cpt_multi++;
    		}*/
    		
    		if ((iter%100)==0 && iter>1 )
    		{
    			//cp++;
    			float moy=0;
    			std::cout << "iter " << iter << " avec max = " << max << std::endl;
    			
    /*			if(cpt_multi<1000)
    			{
    				for ( int i=0; i<size; i++)
    				{
    					diff_multi[cpt_multi].item(i)=smooth.item(i);
    				}
    				Writer< Texture1d > wT1m("diff_lon_multi.tex");
    				wT1m.write(diff_multi);
    				cpt_multi++;
    			}*/
    			
    /*			for ( int i=0; i<size; i++)
    			{
    				diff_temp[0].item(i)=smooth.item(i);
    			}
    			Writer< Texture1d > wT1c(adr_lon);
    			wT1c.write(diff_temp);*/
    			
    			//Calcul de la condition d'arret
    			int cpt=0;
    			for (int i=0; i<size; i++)
    			{
    				if ( long_forbid[0].item(i)==0 )
    				{
    					moy+=fabs( lapl.item(i) );    //* long_forbid[0].item(i) );
    					cpt++;                        // = cpt + cpt * long_forbid[0].item(i) ;
    				}
    			}
    			//std::cout << "moy " << moy << std::endl;
    			moy=moy/cpt;
    			max=fabs(moy-moy_previous);
    			//std::cout << "moy " << moy << " avec cpt = " << cpt << std::endl;
    			moy_previous=moy;
    /*			for ( int i=0; i<size; i++)
    			{
    				diff_temp_big[cpt].item(i)=smooth.item(i);
    			}
    			std::cout << "ECRIT!!!!" << std::endl;
    			Writer< Texture1d > wT1f("diff_latitude_big.tex");
    			wT1f.write(diff_temp_big);
    			cpt++;*/
    			
    		}
    		
    		float epsilon=0.02;
    		for ( int i=0; i<size; i++)
    		{
    			if( fabs(cMap[0].item(i))>epsilon ) 
    			{
    				temp[0].item(i) = smooth.item(i) - cMap[0].item(i);
    			}
    			else
    				temp[0].item(i) = 0;
    			
    			if ( long_forbid.item(i)==0 )
    			{
    				s = smooth.item(i) + _dt * ( lapl.item(i) - betaMap[0].item(i)*temp[0].item(i) );
    				if(s>360)
    					s=360;
    				smooth.item(i) = s;
    			}
    		}
    //     if( iter % 10 == 0 )
    //     tt[iter/10] = smooth;
    //     if( iter == 1000 )
    //     {
    //     Writer<TimeTexture<float> > wtemp("/home/olivier/source.tex");
    //     wtemp.write(tt);
    //     std::cin >> iter;
    //     }
              iter++;
        }
         while((max > criterium)&&(iter<maxIter));
    	//while(iter<300000);
    	
    	std::cout  << std::endl;
    	
    	//Le resultat est stocke dans une TimeTexture
    	for ( int i=0; i<size; i++)
    	{
    		result[0].item(i)=smooth.item(i);
    	}
    	diff_laplacian[0]=lapl;
    
    /*	std::string address_lapl="./diff_laplacian.tex";
    	Writer<Texture1d> w_lapl(address_lapl);
    	if( !w_lapl.write( diff_laplacian ) ) {}*/
    		//return( false );
    
    // 	std::cout << "Max : " << max << std::endl;
    	return(result);
    
    }
    
    TimeTexture<float> CorticalReferential::diffusionLongitude( TimeTexture<float> & tex, TimeTexture<float> & side, TimeTexture<float> & poleSave, int ind)
    {
    
    	//ind=0 : pole insulaire; ind=1 : pole cingulaire
    	
    	TimeTexture<float>  result(1,size);
    	Texture<float>     smooth(size), lapl;
    	float			s;
    
    	//Remove weights from poles
    	unsigned  uneigh;
    	float weight;	
    	TimeTexture<float> long_forbid(1,size);
    	init_texture_single(long_forbid);
    	
    	TimeTexture<float> long_forbid_comp(1,size);
    	init_texture_single(long_forbid_comp);
    	
    	TimeTexture<float> noLapl(1,size);
    	init_texture_single(noLapl);
    	
    	
    	
    	TimeTexture<float> interdit(1,size);
    	init_texture_single(interdit);
    	std::set<uint>::const_iterator itr;
    	
    	int poleT=sud, cpt=0;
    	
    	itr=neigh[nord].begin();
    	for(; itr!=neigh[nord].end(); itr++)
    		if(poles[0].item(*itr)!=0)
    			cpt++;
    	if(cpt==0)
    		poleT=nord;
    	
    	
    	itr=neigh[poleT].begin();
    	for(; itr!=neigh[poleT].end(); itr++)
    		interdit[0].item(*itr)=1;
    	
    	
    	std::map<unsigned, std::set< std::pair<unsigned,float> > >  weightLapl2;
    	
    	weightLapl2=weightLapl;
    	
    	std::set< std::pair<unsigned,float> >::iterator isp,esp;
    	std::map<unsigned, std::set< std::pair<unsigned,float> > >::iterator  il, el;
    	for (il = weightLapl2.begin(), el = weightLapl2.end(); il != el; ++il)
    	{
    		for ( isp = (il->second).begin(), esp = (il->second).end(); isp != esp; ++isp)
    		{
    			uneigh = isp->first;
    			weight = isp->second;
    			//itlist=liste.begin();
    			if( interdit[0].item(uneigh)!=0 )
    			{
    // 				std::cout << uneigh<<" à 0" << std::endl;
    				il->second.erase(std::pair<unsigned,float>(uneigh,weight) );
    				il->second.insert(std::pair<unsigned,float>(uneigh,0) );
    				noLapl[0].item(uneigh)=1;
    			}
    			else {}
    			if( uneigh==(uint)poleT )
    			{
    // 				std::cout << "Pole à 0" << std::endl;
    				il->second.erase(std::pair<unsigned,float>(uneigh,weight) );
    				il->second.insert(std::pair<unsigned,float>(uneigh,0) );
    				noLapl[0].item(uneigh)=1;
    			}
    		}
    	}
    	
    	
    /*	Writer<Texture1d> wT3lttn("test_no_lapl.tex");
    	wT3lttn.write(noLapl);
    
    	Writer<Texture1d> wT3ltt("interdit.tex");
    	wT3ltt.write(interdit);
    	
    	Writer<Texture1d> wT3lptt("long_pole.tex");
    	wT3lptt.write(poles);*/
    	
    	std::set<uint>::const_iterator it;
    	
    	for(int i=0; i<size; i++)
    	{
    		//smooth = diffusion & tex = contraintes
    		smooth.item(i) = 180;
    		if( tex[0].item(i) != 0 )
    			smooth.item(i) = tex[0].item(i);
    			
    		//cree la texture sans les bords pour la condition d'arret
    		if( poles[0].item(i)!=0 )
    			long_forbid[0].item(i)=1;
    		
    		//pour l'insula
    		if( ind==0 && poleSave[0].item(i)!=180 )// && tex[0].item(i)!=0 )
    			long_forbid[0].item(i)=1;
    		//pour le cingulaire
    		if( ind==1 && poleSave[0].item(i)!=1 )// && tex[0].item(i)!=0 )
    			long_forbid[0].item(i)=1;
    		
    		if( noLapl[0].item(i)!=0 )
    			long_forbid[0].item(i)=1;
    // 		it=neigh[i].begin();
    // 		for(; it!=neigh[i].end(); it++)
    // 			if( (poles[0].item(*it)!=0 || tex[0].item(*it)!=0) && (side[0].item(i)==0) )
    // 				long_forbid[0].item(*it)=1;
    				
    		//Cree la texture compl�entaire
    		if(long_forbid[0].item(i)==1)
    			long_forbid_comp[0].item(i)=0;
    		else
    			long_forbid_comp[0].item(i)=1;
    	}
    
    // 	TimeTexture<float> diff_temp(1,size);
    	//int cp=0;
    	
    // 	Writer<Texture1d> wT3ls("test_long_side.tex");
    // 	wT3ls.write(side);
    
    /*	Writer<Texture1d> wT3fd("test_long_tex_proc.tex");
    	wT3fd.write(tex);
    	for ( int i=0; i<size; ++i)
    	{
    	diff_temp[0].item(i)=smooth.item(i);
    }
    // 			Writer<Texture1d> wT1c("diff_long_smooth.tex");
    // 			wT1c.write(diff_temp);
    */
    /*	Writer<Texture1d> wT3lf("test_long_forbid.tex");
    	wT3lf.write(long_forbid);
    	
    	Writer<Texture1d> wT3ltst("long_pole_side.tex");
    	wT3ltst.write(side);
    	
    	Writer<Texture1d> wT3ltstss("tex.tex");
    	wT3ltstss.write(tex);*/
    	
    /*	for ( int i=0; i<size; i++)
    	{
    		diff_temp[0].item(i)=smooth.item(i);
    	}
    	Writer<Texture1d> wT3ltsttt("smooth.tex");
    	wT3ltsttt.write(diff_temp);*/
    	
    	TimeTexture<float> diff_laplacian(1,size);
    	
    	//TimeTexture<float> diff_temp(1,size);
    	float moy_previous=0;
    	std::cout << "Processing" << std::endl;
    
    	float max=10;
    	int iter=0;
    	//for (int iter=0; iter< t; ++iter)
    	do
    	{
    		lapl =  AimsMeshLaplacian_meridian(smooth, tex, weightLapl2, side);
    		
    		if ((iter%100)==0 && iter!=0 )
    		{
    			//cp++;
    			float moy=0;
    			std::cout << "iter " << iter << " avec max = " << max << std::endl;
    			/* for ( int i=0; i<size; i++)
    			{
    				diff_temp[0].item(i)=smooth.item(i);
    			}
    			Writer< Texture1d > wT1c(adr_lon);
    			wT1c.write(diff_temp);*/
    			
    			//Calcul de la condition d'arret
    			int cpt=0;
    			for (int i=0; i<size; i++)
    			{
    				if ( long_forbid[0].item(i)==0 )
    				{
    					moy+=fabs( lapl.item(i) );    //* long_forbid[0].item(i) );
    					cpt++;                        // = cpt + cpt * long_forbid[0].item(i) ;
    				}
    			}
    			//std::cout << "moy " << moy << std::endl;
    			moy=moy/cpt;
    			max=fabs(moy-moy_previous);
    			//std::cout << "moy " << moy << " avec cpt = " << cpt << std::endl;
    			moy_previous=moy;
    			
    		}
    		for ( int i=0; i<size; i++)
    		{
    			if ( tex.item(i)==0 )
    			{
    				s = smooth.item(i) + _dt * lapl.item(i);
    				if(s>360)
    					s=360;
    				smooth.item(i) = s;
    			}
    		}
    		iter++;
    
    	}
    	while(max > criterium);
    	// 	while(iter<300000);
    	
    	std::cout  << std::endl;
    	
    	//Le resultat est stocke dans une TimeTexture
    	for ( int i=0; i<size; i++)
    	{
    		result[0].item(i)=smooth.item(i);
    	}
    	diff_laplacian[0]=lapl;
    
        /*	std::string address_lapl="./diff_laplacian.tex";
        	Writer<Texture1d> w_lapl(address_lapl);
        	if( !w_lapl.write( diff_laplacian ) ) {}*/
        		//return( false );
        
        // 	std::cout << "Max : " << max << std::endl;
    	return(result);
    
    }
    
    
    Texture<float> CorticalReferential::AimsMeshLaplacian_meridian( const Texture<float> &smooth, TimeTexture<float> &source, const std::map<unsigned, std::set< std::pair<unsigned,float> > > &lapl, TimeTexture<float> & side)
    {
    	unsigned				neighb,node;
    	Texture<float>				tex(size);
    	std::map<unsigned, std::set< std::pair<unsigned,float> > >::const_iterator il,el;
    	std::set< std::pair<unsigned,float> >::iterator      ip,ep;
    	float					L, weight;
    	Point3df				temp, temp_current;
    	int 					cpt = 0;
    	ASSERT ( lapl.size() == (unsigned)size);
    
    
    	for (il = lapl.begin(), el = lapl.end(); il != el; ++il, cpt++)
    	{
    		node = il->first;
    		L = 0;
    		temp=vert[cpt];
    
    		//Pondered sum on the neighbour of the node
    		for ( ip = (il->second).begin(), ep = (il->second).end(); ip != ep; ++ip    )
    		{
    			temp_current=vert[(ip->first)];
    			neighb = ip->first;
    			weight = ip->second;
    
    			//We're on the meridian
    			//Side chnge with left or right hemisphere
    			if( (side[0].item(node)==2) && (source[0].item(neighb)==360) && ( !(poles[0].item(neighb)!=0) ) )		//2 for left hemisphere, 4 for right hemisphere
    				L += weight * (- smooth.item(node) + 360 );
    			else
    			//if we are NOT going in the trigonometric direction
    				if( (side[0].item(node)==4) && (source[0].item(neighb)==360) && ( !(poles[0].item(neighb)!=0) ) )		//4 for left hemisphere, 2 for right hemisphere
    					L += weight * (-smooth.item(node));
    					
    				//for all other nodes
    			else
    				L += weight * (- smooth.item(node) + smooth.item(neighb));
    		}
    		tex.item(node) = L;
    	}
    	return(tex);
    }
    
}
    
