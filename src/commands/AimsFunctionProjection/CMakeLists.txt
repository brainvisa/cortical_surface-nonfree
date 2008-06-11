
ADD_EXECUTABLE(AimsFunctionProjection
createKernels.cc
projectOntoMesh.cc
main.cc
)

TARGET_LINK_LIBRARIES(AimsFunctionProjection
cortical_surface
)

INSTALL(TARGETS AimsFunctionProjection DESTINATION bin/)
