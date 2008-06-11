
ADD_EXECUTABLE(AimsParameterizeGyri
main.cc
vector_operations.cc
basictex_operations.cc
model_operations.cc
vertices_operations.cc
constraints_operations.cc
gyri_operations.cc
intersec_operations.cc
compconn_operations.cc
mesh_operations.cc
verif_operations.cc
misctex_operations.cc
gyri_parameterization.cc
)

TARGET_LINK_LIBRARIES(AimsParameterizeGyri
cortical_surface
)

INSTALL(TARGETS AimsParameterizeGyri DESTINATION bin/)
