BRAINVISA_GET_FILE_LIST_FROM_PRO(  ${CMAKE_CURRENT_SOURCE_DIR}/cortical_surface_library.pro "HEADERS" _headers "SOURCES" _sources )
BRAINVISA_COPY_AND_INSTALL_HEADERS( _headers cortical_surface headersTarget )
add_library( cortical_surface SHARED ${_sources} )
add_dependencies( cortical_surface ${headersTarget} )
set_property( TARGET cortical_surface PROPERTY VERSION ${${PROJECT_NAME}_VERSION} )
target_link_libraries( cortical_surface ${AIMS-FREE_LIBRARIES} )
BRAINVISA_INSTALL( TARGETS cortical_surface
                   DESTINATION lib
                   COMPONENT ${PROJECT_NAME} )

