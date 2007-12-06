include( 'base' )

insert( '{protocol}/{subject}',
  'surface', SetWeakAttr( 'category', 'surface coordinates'), SetContent(
    '<subject>_Lhippo', SetType( 'Left hippocampus pole texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rhippo', SetType( 'Right hippocampus pole texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_poles', SetType( 'Left poles texture' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_poles', SetType( 'Right poles texture' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lat_cleaned', SetType( 'Left hemisphere latitude cleaned constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lat_cleaned', SetType( 'Right hemisphere latitude cleaned constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lon_cleaned', SetType( 'Left hemisphere longitude cleaned constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lon_cleaned', SetType( 'Right hemisphere longitude cleaned constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lat_roots', SetType( 'Left hemisphere latitude constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lat_roots', SetType( 'Right hemisphere latitude constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lon_roots', SetType( 'Left hemisphere longitude constraints texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lon_roots', SetType( 'Right hemisphere longitude constraints texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lat', SetType( 'Left hemisphere latitude texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lat', SetType( 'Right hemisphere latitude texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_lon', SetType( 'Left hemisphere longitude texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_lon', SetType( 'Right hemisphere longitude texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_rootsValues', SetType( 'Constraint coordinates values'),
    '<subject>_L_grid', SetType( 'Left hemisphere coordinate grid'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_grid', SetType( 'Right hemisphere coordinate grid'), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_gyri', SetType( 'Left hemisphere gyri parcellation texture'), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyri', SetType( 'Right hemisphere gyri parcellation texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lhippo_Volume', SetType( 'Left Cingular Pole Template Subject' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_Rhippo_Volume', SetType( 'Right Cingular Pole Template Subject' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Talairach_To_Subject_Transformation', SetType( 'Talairach To Subject Transformation'),
    '<subject>_Subject_To_Template_Transformation', SetType( 'Subject To Template Transformation' ),
    '<subject>_L_gyriGraph', SetType( 'Left Gyri Graph' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyriGraph', SetType( 'Right Gyri Graph' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_L_gyriVolume', SetType( 'Left Gyri Volume' ), SetWeakAttr( 'side', 'left' ),
    '<subject>_R_gyriVolume', SetType( 'Right Gyri Volume' ), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_KERNEL', SetType( 'Projection convolution kernels'), SetWeakAttr( 'side', 'left' ) ,
    '<subject>_Rwhite_KERNEL', SetType( 'Projection convolution kernels'), SetWeakAttr( 'side', 'right' ) ,
    '{volume}_<subject>_Lwhite_projection', SetType( 'Functional texture'), SetWeakAttr( 'side', 'left' ) ,
    '{volume}_<subject>_Rwhite_projection', SetType( 'Functional texture'), SetWeakAttr( 'side', 'right' ),
    '<subject>_Lwhite_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'left', 'mesh', 'white'),
    '<subject>_Rwhite_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'right', 'mesh', 'white'),
    '<subject>_Lhemi_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'left', 'mesh', 'hemi'),
    '<subject>_Rhemi_thickness', SetType('Cortical thickness'), SetWeakAttr('side', 'right', 'mesh', 'hemi')
  ),
)


