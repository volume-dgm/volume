set(TestSrc
  ${CMAKE_SOURCE_DIR}/src/Task/Testing/TestMain.cpp
  ${CMAKE_SOURCE_DIR}/src/Task/Task/Task.h
  ${CMAKE_SOURCE_DIR}/src/Task/Testing/TestRun.inl
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/MeshSettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/ParserUtil.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/ScheduleSettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/SettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/SnapshotSettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/SolverSettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/TaskSettingsParser.h
)

set(ElasticSystemSrc 
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystem.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/IniStates.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/VectorFunctors.h

  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/SourceFunctors.h

  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/BoundaryInfos.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ContactInfos.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystemCommon.h
)

set(ElasticVolumeMeshSrc
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/ElasticVolumeMesh.h

  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/DistributedElasticVolumeMesh.h

  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/DistributedElasticVolumeMeshCommon.h
)

set(GeomMeshSrc
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/AABBTree.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/BucketStorage.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/DisjointSetBuilder.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/TimeHierarchyLevelsManager.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/TopologyReconstructor.h

  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/DistributedMeshIO.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/TransitionInfo.h

  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIO.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIOBase.h

  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMeshCommon.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMeshBase.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMesh.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/AdditionalCellInfo.h
)

set(VolumeMethodSrc 
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/Spaces/LagrangeSpace.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/Spaces/PolynomialPrecomputer.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/Spaces/PolynomialSpace.h

  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/FunctionGetters/SourceFunctionGetter.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/FunctionGetters/GhostFunctionGetter.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/FunctionGetters/FunctionGetters.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/FunctionGetters/BoundaryFunctionGetter.h

  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/VolumeMesh/VolumeMesh.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/Cell.h
)

set(BINDIR ${CMAKE_SOURCE_DIR}/bin)

MAKE_DIRECTORY(${BINDIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BINDIR})

IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    SET(BinaryFileName testd)
ENDIF()

IF (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    SET(BinaryFileName test)
ENDIF()

enable_testing()

add_executable(${BinaryFileName}
  ${TestSrc}
  ${NetworkSrc}
  ${MathsSrc}
  ${TinyXmlSrc}
  ${QuadratureIntegrationSrc}
  ${DifferentialSolversSrc}
  ${VtkSrc}
  ${ElasticSystemSrc}
  ${ElasticVolumeMeshSrc}
  ${GeomMeshSrc}
  ${VolumeMethodSrc}
)

target_link_libraries(test
  PRIVATE ZLIB::ZLIB
  PRIVATE MPI::MPI_CXX
  PRIVATE OpenMP::OpenMP_CXX
  PRIVATE tinyxml
  PRIVATE quadrature_integration
  PRIVATE Eigen3::Eigen
  PRIVATE gtest_main
)

include(GoogleTest)
gtest_discover_tests(
  ${BinaryFileName}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/Task/Testing/
)
