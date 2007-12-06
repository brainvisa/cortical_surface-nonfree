

from neuroProcesses import *
import shfjGlobals     

name = 'Cortical Surface Parameterization Pipeline'
userLevel = 2

signature = Signature( 
  'Lgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'left' } ),
  'Rgraph', ReadDiskItem( 'Cortical folds graph', 'Graph',requiredAttributes={ 'side': 'right' } )
  )


def initialization( self ):
    self.linkParameters( 'Lgraph','Rgraph')
    self.linkParameters( 'Rgraph','Lgraph')
    eNode = SerialExecutionNode( self.name, parameterized=self )

    eNode.addChild( 'LeftHemisphere_Process',
                    ProcessExecutionNode( 'LeftHemisphereProcess', optional = 1 ) )
    eNode.addChild( 'RightHemisphere_Process',
                    ProcessExecutionNode( 'RightHemisphereProcess', optional = 1 ) )
    

    eNode.addLink( 'LeftHemisphere_Process.Lgraph', 'Lgraph' )
    eNode.addLink( 'RightHemisphere_Process.Rgraph', 'Rgraph' )
    
    self.setExecutionNode( eNode )
