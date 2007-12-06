/*Automatic constraints projection*/


#ifndef AIMS_AUTO_CONSTRAINTS_H
#define AIMS_AUTO_CONSTRAINTS_H


#include <iostream>
#include <stdlib.h>
#include <aims/mesh/texture.h>
#include <aims/mesh/curv.h>
#include <aims/mesh/surfaceOperation_d.h>
#include <aims/distancemap/meshdistance_d.h>
#include <aims/distancemap/meshvoronoi_d.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/scalespace/meshDiffuse.h>
#include <graph/graph/graph.h>
#include <graph/graph/greader.h>
#include <graph/tree/tree.h>
#include <graph/tree/treader.h>
#include <cartobase/object/sreader.h>
#include <cartobase/object/syntax.h>
#include <aims/io/aimsGraphR.h>
#include <aims/io/aimsGraphW.h>
#include <aims/graph/graphmanip.h>
#include <set>
#include <map>
#include <vector>
#include <iomanip>
#include <assert.h>
#include <aims/def/path.h>
#include <aims/roi/hie.h>
#include <stdio.h>
#include <cartobase/object/syntax.h>

using namespace aims;
using namespace carto;
using namespace std;

class autoConstraints
{

public:
	const std::string & graphfile;
	const std::string & hiefile;
	const std::string & trlfile;
	const std::string & gyrfile;
	const std::string & syntfile;
	const std::string & type;
	SyntaxSet	*_syntax;
	Graph *gr;
	//Hierarchy *hi;
	Tree *tr;


	autoConstraints(const std::string & _graphfile, const std::string & _hiefile, 
				const std::string & _trlfile, const std::string & _gyrfile, 
				const std::string & _syntfile, std::string & _type):
				graphfile(_graphfile),
				hiefile(_hiefile),
				trlfile(_trlfile),
				gyrfile(_gyrfile),
				syntfile(_syntfile),
				type(_type)
	{

	
		cout<<"adress graph = "<<graphfile<<endl;
		cout<<"adress hiefile = "<<hiefile<<endl;
		cout<<"adress trlfile = "<<trlfile<<endl;
		cout<<"adress gyrfile = "<<gyrfile<<endl;
		cout<<"adress syntfile = "<<syntfile<<endl;
		
		Reader<Graph> grd( graphfile.c_str() );
		gr = new Graph;
		gr = grd.read();

		tr = new Tree;
		_syntax = new SyntaxSet;
		string sname = syntfile.c_str();//Path::singleton().syntax() + "/hierarchy.stx";

		SyntaxReader	sr( sname );
		sr.read( *_syntax );

		try
		{
			TreeReader	trd( hiefile.c_str(), *_syntax );
			trd >> *tr;
		}
		catch( exception & e )
		{
			cerr << e.what() << endl;
			delete tr;
			cout<<"ERROR!!"<<endl;
		}

	}

	void runCommands();
	void graph2trl(const std::string & filename);
	void hie2gyr(const std::string & filename);
	TimeTexture<float> argValues2TimeTexture();
	char* removeSide(std::string original);
	bool returnEmptyTexture(set<TimeTexture<float>*> textures);
	TimeTexture<float> argValues2TimeTexture(int number);
};

			
inline
void autoConstraints::runCommands()
{
	graph2trl(trlfile);
	hie2gyr(gyrfile);

}

inline
char * autoConstraints::removeSide(std::string _original)
{
	char *result;
	const char *original;
	int i=0;
	original = _original.c_str();
	do
	{
		if(original[i]=='_')
		{
			result[i]='\0';
		}
		else result[i]=original[i];
		i++;
	}
 	while(original[i]!='\0');
	_original[i]='\0';
	return (result);
}


inline
void autoConstraints::graph2trl(const std::string & filename)
{
	std::string name;
	char *tmp;
	std::ofstream *_stream;
	_stream = new ofstream( filename.c_str(), ios::out | ios::trunc ) ;
//	if((*_stream).bad())
//		return(1); 
	(*_stream) << "#Translation file, taking sulcus names from the current graph (no template)" << '\n';    
	
	set<Vertex*>	sv = gr->vertices();
	set<Vertex*>::iterator	im, em = sv.end();
	for( im=sv.begin(); im!=em; ++im )
	{
		(*im)->getProperty("label", name);
		if(name=="unknown")
			(*_stream) << "%unknown";
		else
		{
			tmp=removeSide( name );
			int i;
			(*_stream) << "%" << tmp;
		}
		(*_stream) << "\n";
		name=' ';
// 		cout<<"wrote %"<<*tmp<<" / ";
	}
	
	cout<<endl;
	(*_stream).close();
}

inline
void autoConstraints::hie2gyr(const std::string & filename)
{
	std::string name;
	char * tmp;
	std::ofstream *_stream;
	_stream = new ofstream( filename.c_str(), ios::out | ios::trunc );
	(*_stream) << "/ "<<type<<" / ";
	Tree::const_iterator	in, fn=tr->end();
	//	éléments
	for( in=tr->begin(); in!=fn; ++in )
	{
		((Tree *)(*in))->getProperty("label", name);
		tmp=removeSide( name );
		(*_stream) << tmp <<" ";
	}
	(*_stream) << "/"<<endl;
	(*_stream).close();
}
/*
inline
bool autoConstraints::returnEmptyTexture(set<TimeTexture<float>*> textures, int value)
{
	set<Vertex*>::iterator	im, em = textures.end();
	//Pour chaque texture, on regarde si la valeur "value" est utilisée
	for( im=textures.begin(); im!=em; ++im )
	{
		for(int j=0;j<
	}
}
*/

inline
TimeTexture<float> autoConstraints::argValues2TimeTexture(int number)
{
			/*
	set<TimeTexture<float>*> textures;
	TimeTexture<float> tex1;
	TimeTexture<float> tex2;
	TimeTexture<float> tex3;
	TimeTexture<float> tex4;
	TimeTexture<float> tex5;
	TimeTexture<float> result;
	
	//Pour chaque sulci, on crée une nouvelle dimension dans result
	for(int i=0;i<number;i++)
	{
	
	}
	
*/

}

// inline
// TimeTexture<float> autoConstraints::argValues2TimeTexture()
// {
// 	TimeTexture<float> result;
//
// 	std::ifstream *_stream;
//
//
// 	return result;
//
// }
// 


#endif

