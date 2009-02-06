# Copyright CEA and IFR 49 (2000-2005)
#
#  This software and supporting documentation were developed by
#      CEA/DSV/SHFJ and IFR 49
#      4 place du General Leclerc
#      91401 Orsay cedex
#      France
#
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license version 2 and that you accept its terms.

from neuroProcesses import *
import shfjGlobals

name = 'Primal Sketch from Texture'

userLevel = 2

signature = Signature(
    'texture', ReadDiskItem( 'Texture','Texture'),
    'white', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh'),
    'primalsketch', WriteDiskItem( 'Primalsketch graph', 'Graph' ),
    'tMin', Float(),
    'tMax', Float(),
    'whiteAux', ReadDiskItem( 'Hemisphere White Mesh', 'MESH mesh'),
    'subject', String(),
    'filterout', Float(),
    'intersectioncriterium', Integer(),
    'latitude', ReadDiskItem( 'Texture','Texture'),
    'longitude', ReadDiskItem('Texture','Texture')
    )

def initialization( self ):
     self.setOptional('whiteAux','subject','filterout','intersectioncriterium','latitude','longitude')
     self.linkParameters( 'primalsketch', 'texture' )
     self.linkParameters( 'white', 'texture' )
     self.tMin = 1.0
     self.tMax = 64.0
     self.filterout = 2.0
     self.intersectioncriterium = 10
     

def execution( self, context ):
     scales=context.temporary( 'GIS image' )
     blobs=context.temporary( 'GIS image' )

     call_list = [ 'AimsTexturePrimalSketch',
                   '-t', self.texture.fullName(),
                   '-o', self.primalsketch.fullName(),
		   '-m', self.white.fullName(),
                   '-os', scales.fullName(),
                   '-ob', blobs.fullName(),
                   '-t1', self.tMin,
                   '-t2', self.tMax,
		   ]
     if (self.whiteAux is not None):
	     call_list += ['-mX', self.whiteAux.fullName()]
     if (self.subject is not ""):
	     call_list += ['-sj', self.subject]
     if (self.filterout is not ""):
	     call_list += ['-f', self.filterout]
     if (self.latitude is not None):
	     call_list += ['-l', self.latitude.fullName()]
     if (self.longitude is not None):
	     call_list += ['-L', self.longitude.fullName()]
     if (self.intersectioncriterium is not ""):
	     call_list += ['-iP', self.intersectioncriterium]
     
 
     context.write('Starting primal sketch computation')
     apply( context.system, call_list )
     context.write('Finished')
