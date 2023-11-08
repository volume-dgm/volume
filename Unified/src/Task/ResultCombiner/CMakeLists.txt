
set(ResultCombinerSrc
  ${CMAKE_SOURCE_DIR}/src/Task/ResultCombiner/ResultCombiner.cpp
  ${CMAKE_SOURCE_DIR}/src/Task/ResultCombiner/SampleIO.h
)

set(ElasticSystemSrc 
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystem.h
)

add_executable(resultcombiner
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

if (COPY_BINARIES_TO_BIN_FOLDER)
    add_custom_command(TARGET resultcombiner POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/.bin
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:resultcombiner> ${CMAKE_SOURCE_DIR}/.bin/$<TARGET_FILE_BASE_NAME:resultcombiner>-$<CONFIG>$<TARGET_FILE_SUFFIX:resultcombiner>)
endif()
