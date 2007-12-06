
include( 'base' )

insert( 'nomenclature', 'surfaceanalysis',
  SetContent(
    'constraint_correspondance', SetType( 'Constraint coordinates values' ),
  )
)

insert( 'nomenclature', 'surfaceanalysis',
  SetContent(
    'surfaceReferential', SetType( 'Surface Label Translation' ),
  )
)

insert( 'nomenclature', 'surfaceanalysis',
  SetContent(
    'surfaceRefModel_par', SetType( 'Latitude Constraint Gyri Model' ),
  )
)

insert( 'nomenclature', 'surfaceanalysis',
  SetContent(
    'surfaceRefModel_mer', SetType( 'Longitude Constraint Gyri Model' ),
  )
)
