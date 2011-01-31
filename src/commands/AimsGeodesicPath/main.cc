/*
 *  Copyright (C) 1997-2005 CEA
 *
 *  This software and supporting documentation were developed by
 *
 *      Laboratoire LSIS, equipe LXAO,
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
//#include <aims/io/writer.h>
#include <aims/io/process.h>
//#include <aims/io/finder.h>
#include <string.h>

//#include <aims/io/aimsGraphW.h>
//#include <aims/def/path.h>
//#include <cartobase/stream/fileutil.h>

#include <aims/geodesicpath/geodesic_mesh.h>
#include <aims/geodesicpath/geodesic_algorithm_dijkstra.h>
#include <aims/geodesicpath/geodesic_algorithm_subdivision.h>
#include <aims/geodesicpath/geodesic_algorithm_exact.h>
//#include <cortical_surface/geodesicpath/geodesic_algorithm_graph_base.h>


#include <fstream>

#include <aims/mesh/curv.h>
//#include <aims/data/data.h>
#include <aims/mesh/surfaceOperation.h>
//#include <aims/mesh/geometric.h>
//#include <aims/data/data_g.h>
//#include <aims/io/io_g.h>
//#include <aims/math/math_g.h>
//#include <aims/vector/vector.h>
//#include <aims/mesh/texture.h>
//#include <aims/distancemap/meshdistance_d.h>
//#include <aims/distancemap/distancemap_g.h>
//#include <aims/morphology/morphology_g.h>

#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace geodesic;



int main( int argc, const char** argv )
{
	try
	{
	  string meshFileIn;
	  string FileOut= "./out";
	  int method = 0;
	  int constraintType = 0;
	  int strain = 5;
	  unsigned source,target;

		AimsApplication    app( argc, argv, "Compute the shortest path between two vertex" );

		app.addOption( meshFileIn, "-i", "mesh" );
		app.alias( "--input", "-i" );

		app.addOption( source, "-s", "index of source vertex" );
		app.alias( "--source", "-s" );

		app.addOption( target, "-t", "index of target vertex" );
		app.alias( "--target", "-t" );

		app.addOption( FileOut, "-o", "output file without extension file specified (.tex or .mesh)" );
		app.alias( "--output", "-o" );

		app.addOption( constraintType, "-c", "constraintType:\n\"0\" -> no constraint\n"
        "\"1\" -> constrained sulci\n\"2\" -> constrained gyri\n\"3\" -> Exact geodesic path",true);
		app.alias( "--constraint", "-c" );

	  app.addOption( strain, "-st", "strain parameter (5 by default)",true );
	  app.alias( "--strain", "-st" );

		app.initialize();

    // read triangulation
    cout << "reading triangulation   : " << flush;
    AimsSurfaceTriangle surface;
    Reader<AimsSurfaceTriangle> triR( meshFileIn );
    triR >> surface;
    cout << "done" << endl;

    // compute and copy curvature
    TimeTexture<float> texCurv;
    cout << "compute texture curvature : ";
    texCurv = TimeTexture<float>(1, surface.vertex().size());
    texCurv = AimsMeshCurvature(surface[0]);
    cout << "done" << endl;

    float *f = (float*) malloc (texCurv[0].nItem() * sizeof(float));
    for( uint i = 0; i < texCurv[0].nItem(); i++)
    {
    f[i] = (float)(texCurv[0].item(i));
    }

    // copy vertex and faces vector
    std::vector<double> pointsSP;
    std::vector<unsigned> facesSP;
    vector<Point3df> & vert = surface.vertex();
    vector<AimsVector<uint, 3> > & tri = surface.polygon();
    pointsSP.resize(3*vert.size());
    facesSP.resize(3*tri.size());

    for (uint j = 0; j < (int) vert.size(); j++)
    {
      pointsSP[3*j] = vert[j][0];
      pointsSP[3*j+1] = vert[j][1];
      pointsSP[3*j+2] = vert[j][2];
    }
    for (uint j = 0; j < (int) tri.size(); j++)
    {
      facesSP[3*j] = tri[j][0];
      facesSP[3*j+1] = tri[j][1];
      facesSP[3*j+2] = tri[j][2];
    }

    // compute adjacence graph
    geodesic::Mesh meshSP;
    cout << "compute adjacence graph : ";
    if (constraintType != 3)
      meshSP.initialize_mesh_data(pointsSP,facesSP, f,constraintType,strain);
    else
      meshSP.initialize_mesh_data(pointsSP,facesSP, NULL ,constraintType,0);

    cout << "done" << endl;

    // compute shortest path
    cout << "compute shortest path : ";

    std::vector<geodesic::SurfacePoint> sources;
    sources.push_back(geodesic::SurfacePoint(&meshSP.vertices()[source]));

    std::vector<geodesic::SurfacePoint> targets;
    targets.push_back(geodesic::SurfacePoint(&meshSP.vertices()[target]));

    printf("indice source = %d target = %d \n",source, target);

    // clear path
    std::vector<geodesic::SurfacePoint> SPath;
    SPath.clear();

    if (constraintType != 3)
    {
      //writing path in the output texture
      TimeTexture<float> texOut(1, surface.vertex().size() );

      // dijkstra method
      geodesic::GeodesicAlgorithmDijkstra *dijkstra_algorithm;
      dijkstra_algorithm = new geodesic::GeodesicAlgorithmDijkstra(&meshSP);

      std::vector<int> listIndexVertexPathSP;
      listIndexVertexPathSP.clear();

      geodesic::SurfacePoint short_sources(&meshSP.vertices()[source]);
      geodesic::SurfacePoint short_targets(&meshSP.vertices()[target]);

      dijkstra_algorithm->geodesic(short_sources,short_targets, SPath, listIndexVertexPathSP);

      //std::vector<int>::iterator ite;
      reverse(listIndexVertexPathSP.begin(),listIndexVertexPathSP.end());
      listIndexVertexPathSP.push_back((int)target);

      cout << "shortest path (index vertex) = ";
      for (unsigned i = 0; i < listIndexVertexPathSP.size(); i++)
        cout << listIndexVertexPathSP[i] << " " ;

      cout << endl;

      for (unsigned i = 0; i < listIndexVertexPathSP.size(); i++)
        texOut[0].item(listIndexVertexPathSP[i]) = 1;

      FileOut = FileOut + ".tex";
      Writer<TimeTexture<float> > texW(FileOut);
      texW << texOut;
    }
    else
    {
      geodesic::GeodesicAlgorithmExact *exact_algorithm;
      exact_algorithm = new geodesic::GeodesicAlgorithmExact(&meshSP);

      //geodesic::GeodesicAlgorithmSubdivision *subdivision_algorithm;
      //subdivision_algorithm = new geodesic::GeodesicAlgorithmSubdivision(&meshSP,2);

      exact_algorithm->propagate(sources);    //cover the whole mesh
      exact_algorithm->print_statistics();
      exact_algorithm->trace_back(targets[0], SPath);

      //geodesic::print_info_about_path(SPath);

      AimsSurfaceTriangle meshOut, *tmpMeshOut;

      tmpMeshOut = new AimsSurfaceTriangle;

      std::vector< Point3df > vertexList;
      Point3df newVertex;

      int i;

      for (i = 0; i < SPath.size(); ++i)
      {
        newVertex[0] = SPath[i].x();
        newVertex[1] = SPath[i].y();
        newVertex[2] = SPath[i].z();
        cout << "(" << SPath[i].x() << ',' << SPath[i].y() << ',' <<  SPath[i].z() << ")\n";
        vertexList.push_back(newVertex);
      }

      for (i = 0; i < vertexList.size() - 1; ++i)
      {
        tmpMeshOut =  SurfaceGenerator::sphere(vertexList[i], 0.25 ,20 );
        SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
        delete tmpMeshOut;

        tmpMeshOut = SurfaceGenerator::cylinder( vertexList[i],vertexList[i+1], 0.2, 0.2, 12, false, true );
        SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
        delete tmpMeshOut;
      }

      tmpMeshOut =  SurfaceGenerator::sphere(vertexList[i], 0.2 ,10 );
      SurfaceManip::meshMerge( meshOut, *tmpMeshOut );
      delete tmpMeshOut;

      FileOut = FileOut + ".mesh";
      Writer<AimsSurfaceTriangle> wm(FileOut);
      wm.write(meshOut);
    }

    cout << "writing " << FileOut << " done\n";

		return( 0 );
	}

	catch( user_interruption & )
	{
	}
	catch( exception & e )
	{
		cerr << e.what() << endl;
	}
	return 1;
}


