
#include <cstdlib>
#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation.h>
#include <aims/mesh/surfacegen.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <cortical_surface/mesh/isoLine.h>

using namespace aims;


AimsSurfaceTriangle IsoLine::makeTubes(int val)
{
	AimsSurfaceTriangle meshResult;
	TimeTexture<short> tex;
	unsigned  i, t=0, p, t1, t2, t3;
	//AimsSurfaceTriangle::iterator itMesh;
	//unsigned		nedge = 0;
	AimsSurfaceTriangle *msh, *tmpMesh;
        
        value=val;

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
		if(fabs(texOriginal[0].item(i)-(float)value)<0.0001)
		{
			marked_texture[0].item(i)=1;
		}
		if(texOriginal[0].item(i)<(float)value)
		{
			marked_texture[0].item(i)=20;
		}
		if(texOriginal[0].item(i)>(float)value)
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

AimsSegments IsoLine::makeLine(int val)
{
        AimsSegments line;
	TimeTexture<short> tex;
	unsigned  i, p, t1, t2, t3;
	
        value=val;
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
				addSegment(vertex[t2], vertex[t3], &line);
			else
			if( (tex[0].item(t2)==10) || (tex[0].item(t2)==20) )
				addSegment( vertex[t1], vertex[t3], &line);
			else
			if( (tex[0].item(t3)==10) || (tex[0].item(t3)==20) )
				addSegment( vertex[t1], vertex[t2], &line);
		}
		if(res==40) //if 2 vertices are = value
		{
			if (tex[0].item(t1)==20)
			{
				v1=createNewVertex(vertex[t1],vertex[t2],t1,t2);
				v2=createNewVertex(vertex[t1],vertex[t3],t1,t3);
				addSegment(v1, v2, &line);
			}
			if (tex[0].item(t2)==20)
			{
				v1=createNewVertex(vertex[t2],vertex[t1],t2,t1);
				v2=createNewVertex(vertex[t2],vertex[t3],t2,t3);
				addSegment(v1, v2, &line);
			}
			if (tex[0].item(t3)==20)
			{
				v1=createNewVertex(vertex[t3],vertex[t1],t3,t1);
				v2=createNewVertex(vertex[t3],vertex[t2],t3,t2);
				addSegment(v1, v2, &line);
			}
		}
			
		if(res==50)
		{
			if (tex[0].item(t1)==10)
			{
				v1=createNewVertex(vertex[t1],vertex[t2],t1,t2);
				v2=createNewVertex(vertex[t1],vertex[t3],t1,t3);
				addSegment(v1, v2, &line);
			}
			if (tex[0].item(t2)==10)
			{
				v1=createNewVertex(vertex[t2],vertex[t1],t2,t1);
				v2=createNewVertex(vertex[t2],vertex[t3],t2,t3);
				addSegment(v1, v2, &line);
			}
			if (tex[0].item(t3)==10)
			{
				v1=createNewVertex(vertex[t3],vertex[t1],t3,t1);
				v2=createNewVertex(vertex[t3],vertex[t2],t3,t2);
				addSegment(v1, v2, &line);
			}
		}
		if(res==31)
		{
			if (tex[0].item(t1)==1)
			{
				v1=createNewVertex(vertex[t2],vertex[t3],t2,t3);
				addSegment(vertex[t1], v1, &line);
			}
			if (tex[0].item(t2)==1)
			{
				v1=createNewVertex(vertex[t1],vertex[t3],t1,t3);
				addSegment(vertex[t2], v1, &line);
			}
			if (tex[0].item(t3)==1)
			{
				v1=createNewVertex(vertex[t1],vertex[t2],t1,t2);
				addSegment(vertex[t3], v1, &line);
			}			
		}
	}
	return line;
}

void IsoLine::addSegment(Point3df v1, Point3df v2, AimsSegments *line)
{
     std::vector<Point3df>  & vert=line->vertex();
     std::vector<AimsVector<uint,2> > & poly=line->polygon();
     uint i;
     uint ind_v1=1000000, ind_v2=1000000, ind;
     double diff1, diff2;
     
     for (i=0; (i<vert.size()) && (ind_v1==1000000) && (ind_v2==1000000); i++)
     {
          diff1=dnorm(vert[i]-v1);
          diff2=dnorm(vert[i]-v2);
          if (diff1<0.001) ind_v1=i;
          if (diff2<0.001) ind_v2=i;
      }
      ind=vert.size()-1;
      if (ind_v1==1000000)
      {
          vert.push_back(v1);
          ind=ind+1;
          ind_v1=ind;
      }
      if (ind_v2==1000000)
      {
          vert.push_back(v2);
          ind=ind+1;
          ind_v2=ind;
      }
      poly.push_back(AimsVector<uint,2>(ind_v1, ind_v2));
}





