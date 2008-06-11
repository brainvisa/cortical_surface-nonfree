
getFromPro(cortical_surface_library.pro HEADERS SOURCES)

SHFJ_COPY_HEADERS(include/cortical_surface ${HEADERS})


ADD_LIBRARY(cortical_surface SHARED
${SOURCES}
)

TARGET_LINK_LIBRARIES(cortical_surface
aimsalgo
)

INSTALL(TARGETS cortical_surface DESTINATION lib/)
