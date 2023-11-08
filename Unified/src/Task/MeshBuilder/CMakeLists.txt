
set(MeshBuilderSrc
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshBuilder.cpp
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshDistributor.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitter.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitterCommon.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitter2.inl
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitter3.inl
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitterCommon.inl
)

set(SpecificMeshBuildersSrc
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/SpecificMeshBuilders/BasicMeshBuilder.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/SpecificMeshBuilders/SalomeMeshBuilder.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/SpecificMeshBuilders/SpecificMeshBuilders.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/SpecificMeshBuilders/UnstructuredMeshBuilder.h
)

set(ElasticSystemSrc 
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystem.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/IniStates.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/VectorFunctors.h
)

set(TriGeomSrc
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/DistributedMeshIO.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIO.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMesh.h
)

add_executable(meshbuilder
  ${MeshBuilderSrc}
  ${SpecificMeshBuildersSrc}
  ${VtkSrc}
  ${MathsSrc}
  ${TinyXmlSrc}
  ${ElasticSystemSrc}
  ${TriGeomSrc}
  ${TriangleSrc}
  ${MetisGKlibSrc}
  ${MetisIncludeSrc}
  ${MetisLibSrc}
  ${MetisProgramSrc}
)

target_link_libraries(meshbuilder
    PUBLIC ZLIB::ZLIB 
    PUBLIC OpenMP::OpenMP_CXX
    PUBLIC metis
    PUBLIC tinyxml
    PUBLIC Eigen3::Eigen
)

if (COPY_BINARIES_TO_BIN_FOLDER)
    add_custom_command(TARGET meshbuilder POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/.bin
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:meshbuilder> ${CMAKE_SOURCE_DIR}/.bin/$<TARGET_FILE_BASE_NAME:meshbuilder>-$<CONFIG>$<TARGET_FILE_SUFFIX:meshbuilder>)
endif()
