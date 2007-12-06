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

name = 'Sulcus Parameterization'

userLevel = 2

signature = Signature(
    'graph', ReadDiskItem( 'Cortical folds graph', 'Graph' ),
    'mri', ReadDiskItem( 'T1 MRI Bias Corrected', 'GIS image' ),
    'label_attributes', Choice( 'label', 'name' ),
    'label_values', String(),
    'orientation', Choice( 'Top->Bottom', 'Front->Back' ),
    'sulcus_mesh', WriteDiskItem( 'Mesh', 'MESH mesh' ),
    'texture_param1', WriteDiskItem( 'Coordinate Texture', 'Texture' ),
    'texture_param2', WriteDiskItem( 'Coordinate Texture', 'Texture' ),
    'coordinates_grid', WriteDiskItem( 'Mesh', 'MESH mesh'),
    'dilation', Float()
)

def initialization( self ):
     self.linkParameters( 'mri', 'graph' )
     self.label_attributes = 'name'
     self.dilation = 1.0


def execution( self, context ):
     sulcusIm=context.temporary( 'GIS image' )
     hullIm=context.temporary(  'GIS image' )
     bottomIm=context.temporary(  'GIS image' )
     simplesurf=context.temporary( 'GIS image' )
     dilatedIm=context.temporary(  'GIS image' )
     isoL=context.temporary( 'MESH mesh')
     transform=''

     if (self.orientation=='Top->Bottom'):
          orient=0
     else :
          orient=1

     context.write('Extracting sulcus and buckets from graph')

     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         compress= 'No',
                         bucket= 'Sulci',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )
     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         bottom = bottomIm.fullPath(),
                         compress= 'No',
                         bucket= 'Bottoms',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )
     context.runProcess( 'Create Sulcus Label Volume',
                         graph = self.graph,
                         mri = self.mri,
                         sulci =  sulcusIm.fullPath(),
                         simple_surface = simplesurf.fullPath(),
                         hull_junction = hullIm.fullPath(),
                         compress= 'No',
                         bucket= 'Junctions with brain hull',
                         label_attributes = self.label_attributes,
                         label_values = self.label_values,
                         node_edge_types='All' )
     dilating = [ 'AimsDilation',
                 '-i', sulcusIm.fullPath(),
                 '-o', dilatedIm.fullPath(),
                 '-e', self.dilation ]

     context.write('Remeshing sulcus')

     meshing = [ 'AimsMesh',
                 '-i', dilatedIm.fullPath(),
                 '-o', self.sulcus_mesh.fullPath(),
                 '-l', '32767',
                 '--smooth',
                 '--smoothIt', '20' ]

     apply( context.system, dilating )
     apply( context.system, meshing )

     test=self.sulcus_mesh.fullName()
     sulcusMname=test + '_32767_0.mesh'
     moving = [ 'mv', sulcusMname, self.sulcus_mesh.fullPath()]
     movingminf = [ 'mv', sulcusMname + '.minf', self.sulcus_mesh.fullPath() + '.minf']
     apply( context.system, moving )
     apply( context.system, movingminf )

     context.write('Starting parameterisation')

     parameterising = [ 'AimsParameterizeSulcus',
                        '-i', self.sulcus_mesh.fullPath(),
                        '-b', bottomIm.fullPath(),
                        '-t', hullIm.fullPath(),
                        '-o', orient,
                        '-ox', self.texture_param1.fullPath(),
                        '-oy', self.texture_param2.fullPath() ]

     apply(context.system, parameterising )

     context.write('Computing coordinate grid')

     i=0;
     iso = [ 'AimsMeshIsoLine',
              '-i', self.sulcus_mesh.fullPath(),
              '-t', self.texture_param2.fullPath(),
              '-o', isoL.fullPath(),
              '-v', i ]
     conc = [ 'AimsZCat',
              '-i', isoL.fullPath(),
              '-o', self.coordinates_grid.fullPath() ]
     apply(context.system, iso)
     apply(context.system, conc)
     iso = [ 'AimsMeshIsoLine',
              '-i', self.sulcus_mesh.fullPath(),
              '-t', self.texture_param1.fullPath(),
              '-o', isoL.fullPath(),
              '-v', i ]
     conc = [ 'AimsZCat',
              '-i', isoL.fullPath(),
              '-o', self.coordinates_grid.fullPath() ]
     apply(context.system, iso)
     apply(context.system, conc)
     i=i+10
     while (i<=100):
          iso = [ 'AimsMeshIsoLine',
                   '-i', self.sulcus_mesh.fullPath(),
                   '-t', self.texture_param2.fullPath(),
                   '-o', isoL.fullPath(),
                   '-v', i ]
          conc = [ 'AimsZCat',
                   '-i', isoL.fullPath(), self.coordinates_grid.fullPath(),
                   '-o', self.coordinates_grid.fullPath() ]
          apply(context.system, iso)
          apply(context.system, conc)
          i=i+10
     i=10
     while (i<=200):
          iso = [ 'AimsMeshIsoLine',
                   '-i', self.sulcus_mesh.fullPath(),
                   '-t', self.texture_param1.fullPath(),
                   '-o', isoL.fullPath(),
                   '-v', i ]
          conc = [ 'AimsZCat',
                   '-i', isoL.fullPath(), self.coordinates_grid.fullPath(),
                   '-o', self.coordinates_grid.fullPath() ]
          apply(context.system, iso)
          apply(context.system, conc)
          i=i+10

     context.write('Finished')

