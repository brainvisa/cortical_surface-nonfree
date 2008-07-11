#    Hopital Pitie Salpetriere
#    91 boulevard de l'Hopital
#    75634 Paris cedex 13
#    France
#    --
#    CEA/DSV/SHFJ
#    4 place du General Leclerc
#    91401 Orsay cedex
#    France
#    --
#    CNRS UPR640-LENA
#    Hopital Pitie Salpetriere
#    47 boulevard de l'Hopital
#    75651 Paris cedex 13
#    France
#
#  $Id$
#

from neuroProcesses import *
import shfjGlobals
from soma import aims
import sys
from qt import *
import os
import anatomist.direct.api as anatomist

name = 'Average Sulci'
userLevel = 2


signature = Signature(
    'sulci', ListOf(ReadDiskItem('Mesh', 'MESH mesh' )),
    'mean', WriteDiskItem('Mesh', 'MESH mesh' ),
    'dev', WriteDiskItem('Texture', 'Texture' )
)

def initialization( self ):
     pass

def execution ( self, context ):

     reader = aims.Reader()
     ns=len(self.sulci)
     sulcus=[]
     v=[]
     av=[]
     print 'Found ', ns, ' sulci'
     for i in xrange(ns) :
          print 'Reading ', self.sulci[i].fullPath()
          sulcus.append(reader.read(self.sulci[i].fullPath()))
          v.append(sulcus[i].vertex(0))
     nv=v[0].size()
     for i in xrange( nv ) :
          p=aims.vector_POINT3DF()
          p=[0.0, 0.0, 0.0]
          av.append(p)
          for j in xrange( ns ) :
               av[i]+=v[j][i]
          av[i]=av[i]/ns

     m = aims.AimsTimeSurface_3()
     m.vertex(0).assign( av )
    # m.normal(0).assign( sulcus[0].normal(0) )
     m.polygon(0).assign( sulcus[0].polygon(0) )
     m.updateNormals()
     a = anatomist.Anatomist()
     am = a.toAObject( m )
     aw=a.createWindow( '3D' )
     a.addObjects( [ am ], [ aw ] )
     W = aims.Writer()
     W.write(m, self.mean.fullPath)



