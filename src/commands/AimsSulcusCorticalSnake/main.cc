/*
 *  Copyright (C) 1997-2005 CEA
 *
 *  This software and supporting documentation were developed by
 *
 *      Laboratoire LSIS, �uipe LXAO,
 *      Marseille, France
 *
 *   CEA/DSV/SHFJ
 *   4 place du General Leclerc
 *   91401 Orsay cedex
 *   France
 */
#include <iostream>
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/io/process.h>
#include <aims/io/finder.h>
#include <string.h>
#include <aims/io/aimsGraphW.h>
#include <aims/def/path.h>
#include <cartobase/stream/fileutil.h>

#include <aims/surfacereferential/sulcusCorticalSnake.h>

using namespace aims;
using namespace carto;
using namespace std;


int main( int argc, const char** argv )
{
	try
	{
		TimeTexture<short> constraint;
		TimeTexture<float> constraint_float;
		//std::string constraint;
		float value;
		float alpha1, alpha2, alpha3;
		AimsSurfaceTriangle mesh;
		std::string adr_constraint;
		std::string adr_mesh;
		std::string adr_graph;
		std::string label_p;
		std::string adr_output;
		TimeTexture<float> result;

		AimsApplication    app( argc, argv, "Project an entire sulcus using multi-scale snake" );
		app.addOption( adr_constraint, "-i", "input constraint texture" );
		app.alias( "--input", "-i" );
		app.addOption( adr_mesh, "-m", "hemisphere mesh" );
		app.alias( "--mesh", "-m" );
		app.addOption( value, "-v", "value of the constraint" );
		app.alias( "--value", "-v" );
		app.addOption( alpha1, "-a", "Constraint distance parameter" );
		app.alias( "--alpha1", "-a" );
		app.addOption( alpha2, "-b", "curvature parameter" );
		app.alias( "--alpha2", "-b" );
		app.addOption( alpha3, "-c", "elasticity parameter" );
		app.alias( "--alpha3", "-c" );
		app.addOption( adr_output, "-o", "output texture" );
		app.alias( "--output", "-o" );
		app.initialize();
		
		
		Reader < TimeTexture<short> > r1(adr_constraint);
		r1.read( constraint );
		
			//Opening brain mesh
		Reader < AimsSurfaceTriangle > r2(adr_mesh);
		r2.read( mesh );
		
// 		std::cout<<"Copie!"<<std::endl;
		for( uint i=0; i<constraint[0].nItem(); i++)
		{
			constraint_float[0].push_back( (float)(constraint[0].item(i)) );
		//	std::cout<<(float)constraint[0].item(i)<<" ";
		}
		
		SulcusCorticalSnake k( constraint_float, value, alpha1, alpha2, alpha3, mesh);
		result = k.compute_snake();
// 		k.create_bottom_volume();
		
		Writer<Texture1d> wr1(adr_output);
		wr1.write(result);
		
		return( 0 );
	}

	/*  catch( user_interruption & )
	{
	}*/
	catch( exception & e )
	{
		cerr << e.what() << endl;
	}
	return 1;
}


