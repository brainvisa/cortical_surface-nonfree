//AimsMeshPaint -im /home/arnaud/Bureau/vrac/AYo1_t2_nobias_cortex_hemil.mesh -it /home/arnaud/Bureau/vrac/AYo1_blobs.tex -ot /home/arnaud/Bureau/vrac/texout.gii

//AimsMeshPaint -im /home/arnaud/Bureau/vrac/subject01_Lwhite_inflated.mesh -it /home/arnaud/Bureau/vrac/subject01_Lwhite_curv.tex -ic blue_red_bis.rgb -ot /home/arnaud/Bureau/vrac/texout.gii
//anatomist --userLevel 4  /home/arnaud/Bureau/vrac/subject01_Lwhite_inflated.mesh /home/arnaud/Bureau/vrac/subject01_Lwhite_curv.tex
#include "meshpaint.h"

int main(int argc, const char **argv)
{

//DECLARATIONS
std::string adressTexIn="./";
std::string adressMeshIn="./";
std::string adressTexOut="./";
std::string colorMap="./";

AimsSurfaceTriangle meshIn;
TimeTexture<float> texIn;

AimsApplication     app( argc, argv, "MeshPaint : draw a texture on mesh");

try
{
app.addOption( adressMeshIn, "-im", "input mesh");
app.alias( "--inputMesh", "-im" );
app.addOption( adressTexIn, "-it", "input texture");
app.alias( "--inputTex", "-it" );
app.addOption( colorMap, "-ic", "input colormap");
app.alias( "--inputColorMap", "-ic" );
app.addOption( adressTexOut, "-ot", "output texture");
app.alias( "--outputTex", "-ot" );
app.initialize();

std::cout << "Reading mesh and texture" << endl;

Reader < AimsSurfaceTriangle > rmeshIn(adressMeshIn);
rmeshIn.read( meshIn );

Reader < TimeTexture<float> > rtexIn(adressTexIn);
rtexIn.read( texIn );

//for (uint i=0; i<20; i++)
//{
//cout << (float) texIn[0].item(i) << endl;
//}
//cout << "computing neighbours  " << endl;
//vector<set<uint> >  neigh = SurfaceManip::surfaceNeighbours(meshIn);
//uint sizeV= meshIn.vertex().size();

////texture copy
//TimeTexture<short> texOut(1, sizeV);
//for (uint i=0; i<sizeV; i++)
//{
//texOut[0].item(i)= (short) floor(texIn[0].item(i));
//}
//
//std::cout << "Writing new texture" << endl;
//Writer< TimeTexture<short> > wtexOut(adressTexOut);
//
//Object da_attrTex = Object::value( IntDictionary() );
//Object da_sub_attrTex = Object::value( Dictionary() );
//
//da_sub_attrTex->setProperty( "intent", "NIFTI_INTENT_LABEL" );
//da_sub_attrTex->setProperty( "data_type", "U32" );
//da_sub_attrTex->setProperty( "encoding",3);
//
//if( da_sub_attrTex->size() != 0 )
//   da_attrTex->setArrayItem( 0, da_sub_attrTex );
//
//if( da_attrTex->size() != 0 )
// texOut.header().setProperty( "GIFTI_dataarrays_info", da_attrTex );
//
//wtexOut.write( texOut );
//
////test writer mesh
//AimsSurfaceTriangle meshOut;
//Writer < AimsSurfaceTriangle > wmeshOut(adressMeshOut);
//
//Object da_attrMesh = Object::value( IntDictionary() );
//Object da_sub_attrMesh = Object::value( Dictionary() );
//
//da_sub_attrMesh->setProperty( "encoding",3);
//
//if( da_sub_attrMesh->size() != 0 )
// da_attrMesh->setArrayItem( 0, da_sub_attrMesh );
//
//if( da_attrMesh->size() != 0 )
// meshOut.header().setProperty( "GIFTI_dataarrays_info", da_attrMesh );
//
//for (uint i=0; i<sizeV; i++)
//{
//meshOut.vertex() = meshIn.vertex();
//}
//
//meshOut.header().setProperty( "vertex_number", sizeV );
//
//uint sizeP= meshIn.polygon().size();
//
//for (uint i=0; i<sizeP; i++)
//{
//meshOut.polygon() = meshIn.polygon();
//}
//
//meshOut.header().setProperty( "polygon_number", sizeP );
//
//wmeshOut.write( meshOut );

}
catch( user_interruption & )
   {
      return EXIT_FAILURE;
   }
   catch( exception & e )
   {
      cerr << e.what() << endl;
      return EXIT_FAILURE;
   }

QApplication a(argc, (char **)argv);
MeshPaint w(meshIn,texIn,colorMap);
w.show();
return a.exec();
}

