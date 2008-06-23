#include <aims/getopt/getopt2.h>

using namespace aims;
using namespace carto;
using namespace std;

void ConcatenerGraphes( const vector<Graph*> in, Graph & out, const string & subjectatt ){
  vector<Graph*>::const_iterator  ig, eg = in.end();
  Graph::iterator     iv, ev;
  string        subject;
  for( ig=in.begin(); ig!=eg; ++ig )
  {
    (*ig)->getProperty( subjectatt, subject );
    for( iv=(*ig)->begin(), ev=(*ig)->end(); iv!=ev; ++iv ){
      (*iv)->setProperty( "subject", subject );
    }
    (*ig)->extract( out, (*ig)->begin(), (*ig)->end() );
  }
  out.setProperty( "subject", subjectatt );
}


vector<string> splitGraphFile(string graphFile){
  vector<string>  graphFiles;
  string::size_type pos = graphFile.find( '|' );
  string::size_type pos2 = 0;

  while( pos != string::npos ) {
    graphFiles.push_back( graphFile.substr( pos2, pos - pos2 ) );
    pos2 = pos + 1;
    pos = graphFile.find( '|', pos2 );
  }
  graphFiles.push_back( graphFile.substr( pos2, graphFile.length() - pos2 ) );
  return graphFiles;
}

void LireGraphes(string graphFile, string output, Graph &primal ){
  vector<string> graphFiles = splitGraphFile(graphFile);
  if( output.empty() )
    output = graphFile;
  vector<string> outputs = splitGraphFile(output);

  if( graphFiles.size() != outputs.size() )
    throw runtime_error( "Input and output graphs numbers differ" );
  if( graphFiles.size() == 1 ){
    Reader<Graph> fr( graphFile );
    try{
      cout << "Lecture du graphe sillons " << graphFile << endl;
      fr.read(primal);
      Graph::iterator iv;
      string        subject;
      primal.getProperty("sujet", subject);
      for( iv=primal.begin(); iv!=primal.end(); ++iv )
        (*iv)->setProperty( "subject", subject );
      cout << "Lecture FGraph OK." << endl;
    }
    catch( parse_error & e )   {
      cerr << e.what() << " : " << e.filename() << ", line " << e.line() << endl;
      throw;
    }
    catch( exception & e )  {
      cerr << graphFile << ": " << e.what() << endl;
      throw;
    }
  }
  else  {
    vector<string>::iterator  i, e = graphFiles.end();
    vector<Graph *> vg(1);
    for( i=graphFiles.begin(); i!=e; ++i )   {
      Reader<Graph> fr( *i );
      try  {
        Graph test;
        cout << "Lecture du graphe " << *i << endl;
        fr.read(test);
        cout << "Lecture FGraph OK.\n";
        vg[0] = &test;
        ConcatenerGraphes( vg, primal, "sujet" ); // -----*****
      }
      catch( parse_error & e ) {
        cerr << e.what() << " : " << e.filename() << ", line "  << e.line() << endl;
        throw;
      }
      catch( exception & e ){
        cerr << *i << ": " << e.what() << endl;
        throw;
      }
    }
  }
}

void RecupererGraphesIndividuels( const std::vector<Graph*> subjects, Graph & multi, const std::string & vertexatt ){
  Graph::const_iterator       imv, emv=multi.end();
  vector<Graph*>::const_iterator    ig, eg = subjects.end();
  map<string, map<int, Vertex *> >    smap;
  map<string, map<int, Vertex *> >::iterator  is, es = smap.end();
  string          subjatt, s;
  int           index;
  Graph::const_iterator       isv, esv;
  map<int, Vertex *>::iterator      iv;
  string          l;

  multi.getProperty( "subject", subjatt );

  for( ig=subjects.begin(); ig!=eg; ++ig )
  {
    (*ig)->getProperty( subjatt, s );
    map<int, Vertex *>  & m = smap[ s ];
    for( isv=(*ig)->begin(), esv=(*ig)->end(); isv!=esv; ++isv )
    {
      (*isv)->getProperty( vertexatt, index );
      m[index] = *isv;
    }
  }

  for( imv=multi.begin(); imv!=emv; ++imv )
  {
    (*imv)->getProperty( "subject", s );
    is = smap.find( s );
    if( is != es )
    {
      (*imv)->getProperty( vertexatt, index );
      iv = is->second.find( index );
      if( iv == is->second.end() )
        cerr << "warning: unfortunately, index " << index
            << " in subject " << s << " doesn't exist\n";
      else
      {
        (*imv)->getProperty( "label", l );
        iv->second->setProperty( "label", l );
        iv->second->setProperty( "name", l );
      }
    }
  }
}

void SauvegarderGraphes(Graph &primal, string graphFile, string output){
  vector<string> graphFiles = splitGraphFile(graphFile);
  if( output.empty() )
    output = graphFile;
  vector<string> outputs = splitGraphFile(output);
  if( outputs.size() == 1 ){
    Writer<Graph> fr( output );
    fr.write(primal);
  }
  else  {
    unsigned  i, n = graphFiles.size();
    vector<Graph *> vg( 1 );
    for( i=0; i<n; ++i ){
      try  {
        Reader<Graph> fr( graphFiles[i] );
        Graph tmpfg;
        fr.read(tmpfg);

        vg[0] = &tmpfg;
        RecupererGraphesIndividuels( vg, primal, "index" ); // warning...
        cout << "Writing :" << outputs[i] << endl;
        Writer<Graph> fw( outputs[i] );
        fw.write(tmpfg);
      }
      catch( parse_error & e ) {
        cerr << e.what() << " : " << e.filename() << ", line "
            << e.line() << endl;
      }
      catch( exception & e ) {
        cerr << e.what() << endl;
      }
    }
  }
}



