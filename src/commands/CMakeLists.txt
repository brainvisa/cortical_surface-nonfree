#------------------------------------------------------------------------------
# Build targets for commands by looking into directories listed in connectomist_commands.pro
#------------------------------------------------------------------------------
BRAINVISA_GET_FILE_LIST_FROM_PRO(  ${CMAKE_CURRENT_SOURCE_DIR}/cortical_surface_commands.pro "PSUBDIRS" _commands )
foreach( _command ${_commands} )
  file( GLOB _pro "${_command}/*.pro" )
  BRAINVISA_GET_FILE_LIST_FROM_PRO( "${_pro}"
                              TARGET _target
                              SOURCES _proSources )
  set( _sources )
  foreach( _f ${_proSources} )
    set( _sources ${_sources} ${_command}/${_f} )
  endforeach()
  brainvisa_add_executable( ${_target} ${_sources} )
  target_link_libraries( ${_target} cortical_surface ${AIMS-FREE_LIBRARIES} )
  BRAINVISA_INSTALL( TARGETS ${_target}
                     DESTINATION bin
                     COMPONENT ${PROJECT_NAME} )
endforeach()
