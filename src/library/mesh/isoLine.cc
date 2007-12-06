
#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/isoLine.h>

using namespace aims;


AimsSurfaceTriangle IsoLine::makeTubes()
{
	AimsSurfaceTriangle meshResult;
	TimeTexture<short> tex;
	unsigned  i, t=0, p, t1, t2, t3;
	//AimsSurfaceTriangle::iterator itMesh;
	//unsigned		nedge = 0;
	AimsSurfaceTriangle	*msh, *tmpMesh;

	tmpMesh=new AimsSurfaceTriangle;
	
	tex=setVertices();
	
	p = poly.size();

	for( i=0; i<p; ++i )
	{
		Point3df v1,v2;
		t1 = poly[i][0];
		t2 = poly[i][1];
		t3 = poly[i][2];
		
		short res=tex[0].item(t1)+tex[0].item(t2)+tex[0].item(t3);
		
		if( (res==22) || (res==12) ) //if 2 vertices are = value
		{
			if( (tex[0].item(t1)==10) || (tex[0].item(t1)==20) )
				msh = SurfaceGenerator::cylinder( vertex[t2], vertex[t3], radius1, radius2, facets, closed, smooth );
			else
			if( (tex[0].item(t2)==10) || (tex[0].item(t2)==20) )
				msh = SurfaceGenerator::cylinder( vertex[t1], vertex[t3], radius1, radius2, facets, closed, smooth );
			else
			if( (tex[0].item(t3)==10) || (tex[0].item(t3)==20) )
				msh = SurfaceGenerator::cylinder( vertex[t1], vertex[t2], radius1, radius2, facets, closed, smooth );
			
			SurfaceManip::meshMerge( *tmpMesh, *msh );
			delete msh;
		}
		if(res==40) //if 2 vertices are = value
		{
			if (tex[0].item(t1)==20)
			{
				v1=createNewVertex(vertex[t1],vertex[t2],t1,t2);
				v2=createNewVertex(vertex[t1],vertex[t3],t1,t3);
				msh = SurfaceGenerator::cylinder( v1, v2, radius1, radius2, facets, closed, smooth );
			}
			if (tex[0].item(t2)==20)
			{
				v1=createNewVertex(vertex[t2],vertex[t1],t2,t1);
				v2=createNewVertex(vertex[t2],vertex[t3],t2,t3);
				msh = SurfaceGenerator::cylinder( v1, v2, radius1, radius2, facets, closed, smooth );
			}
			if (tex[0].item(t3)==20)
			{
				v1=createNewVertex(vertex[t3],vertex[t1],t3,t1);
				v2=createNewVertex(vertex[t3],vertex[t2],t3,t2);
				msh = SurfaceGenerator::cylinder( v1, v2, radius1, radius2, facets, closed, smooth );
			}
			SurfaceManip::meshMerge( *tmpMesh, *msh );
			delete msh;
		}
			
		if(res==50)
		{
			if (tex[0].item(t1)==10)
			{
				v1=createNewVertex(vertex[t1],vertex[t2],t1,t2);
				v2=createNewVertex(vertex[t1],vertex[t3],t1,t3);
				msh = SurfaceGenerator::cylinder( v1, v2, radius1, radius2, facets, closed, smooth );
			}
			if (tex[0].item(t2)==10)
			{
				v1=createNewVertex(vertex[t2],vertex[t1],t2,t1);
				v2=createNewVertex(vertex[t2],vertex[t3],t2,t3);
				msh = SurfaceGenerator::cylinder( v1, v2, radius1, radius2, facets, closed, smooth );
			}
			if (tex[0].item(t3)==10)
			{
				v1=createNewVertex(vertex[t3],vertex[t1],t3,t1);
				v2=createNewVertex(vertex[t3],vertex[t2],t3,t2);
				msh = SurfaceGenerator::cylinder( v1, v2, radius1, radius2, facets, closed, smooth );
			}
			SurfaceManip::meshMerge( *tmpMesh, *msh );
			delete msh;
		}
		if(res==31)
		{
			if (tex[0].item(t1)==1)
			{
				v1=createNewVertex(vertex[t2],vertex[t3],t2,t3);
				msh = SurfaceGenerator::cylinder( vertex[t1], v1, radius1, radius2, facets, closed, smooth );
			}
			if (tex[0].item(t2)==1)
			{
				v1=createNewVertex(vertex[t1],vertex[t3],t1,t3);
				msh = SurfaceGenerator::cylinder( vertex[t2], v1, radius1, radius2, facets, closed, smooth );
			}
			if (tex[0].item(t3)==1)
			{
				v1=createNewVertex(vertex[t1],vertex[t2],t1,t2);
				msh = SurfaceGenerator::cylinder( vertex[t3], v1, radius1, radius2, facets, closed, smooth );
			}
			SurfaceManip::meshMerge( *tmpMesh, *msh );
			delete msh;
			
		}
	}
	meshResult[t]=(*tmpMesh)[0];
	delete tmpMesh;
	return meshResult;
}

TimeTexture<short> IsoLine::setVertices()
{
	int size = texOriginal[0].nItem();
	TimeTexture<short> marked_texture(1,size);
	for(int i=0;i<size;i++)
	{
		if(int(rint(texOriginal[0].item(i)==value)))
		{
			marked_texture[0].item(i)=1;
		}
		if(int(rint(texOriginal[0].item(i)<value)))
		{
			marked_texture[0].item(i)=20;
		}
		if(int(rint(texOriginal[0].item(i)>value)))
		{
			marked_texture[0].item(i)=10;
		}
	}
	return marked_texture;
}


Point3df IsoLine::createNewVertex(Point3df & v1, Point3df & v2, int t1, int t2)
{
	Point3df newV;
	float p1,p2;
	p1= (float)value - texOriginal[0].item(t1);
	p2= (float)value - texOriginal[0].item(t2);
	
	if( p1< 0 )
		p1=-p1;
	if( p2< 0 )
		p2=-p2;
	
	newV[0] = ( (v1[0]*p2) + (v2[0]*p1) ) / (p1+p2);
	newV[1] = ( (v1[1]*p2) + (v2[1]*p1) ) / (p1+p2);
	newV[2] = ( (v1[2]*p2) + (v2[2]*p1) ) / (p1+p2);
	return newV;
}


