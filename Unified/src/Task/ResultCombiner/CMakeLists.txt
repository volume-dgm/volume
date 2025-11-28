set(ResultCombinerSrc
  ${CMAKE_SOURCE_DIR}/src/Task/ResultCombiner/ResultCombiner.cpp
  ${CMAKE_SOURCE_DIR}/src/Task/ResultCombiner/SampleIO.h
)

set(ElasticSystemSrc 
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystem.h
)

set(BINDIR ${CMAKE_SOURCE_DIR}/bin)

MAKE_DIRECTORY(${BINDIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BINDIR})

IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    SET(BinaryFileName resultcombinerd)
ENDIF()

IF (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    SET(BinaryFileName resultcombiner)
ENDIF()

add_executable(${BinaryFileName}
${ResultCombinerSrc}
${VtkSrc}
${MathsSrc}
${TinyXmlSrc}
${ElasticSystemSrc}
${TriGeomSrc}
)

target_link_libraries(resultcombiner
  PUBLIC ZLIB::ZLIB 
  PUBLIC OpenMP::OpenMP_CXX
  PUBLIC tinyxml
  PUBLIC Eigen3::Eigen
)		
