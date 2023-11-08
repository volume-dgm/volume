cmake_minimum_required(VERSION 3.15)
project(Volume)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}     -mtune=native -funroll-loops -ftree-vectorize -mfpmath=sse -ffast-math -ftree-vectorizer-verbose=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native -funroll-loops -ftree-vectorize -mfpmath=sse -ffast-math -ftree-vectorizer-verbose=0")
if (MINGW)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message("Generating makefiles for: " ${CMAKE_SYSTEM})
message("Config: " ${CMAKE_BUILD_TYPE})

option(COPY_BINARIES_TO_BIN_FOLDER "After compiling copy binaries to .bin folder." YES)


find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_compile_definitions(Eigen3::Eigen
    INTERFACE EIGEN_DEFAULT_TO_ROW_MAJOR 
    INTERFACE EIGEN_DONT_PARALLELIZE 
    INTERFACE EIGEN_NO_DEBUG 
    INTERFACE EIGEN_NO_AUTOMATIC_RESIZING 
    INTERFACE EIGEN_FAST_MATH EIGEN_NO_MALLOC
)

set(QuadratureIntegrationSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/tetrahedron_arbq_rule.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/triangle_fekete_rule.cpp
)

set(TinyXmlSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinystr.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxml.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxmlerror.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxmlparser.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinystr.h
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxml.h
)

set(MetisGKlibSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/b64.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/blas.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/csr.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/error.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/evaluate.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/fkvkselect.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/fs.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/getopt.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_arch.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_externs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_getopt.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/GKlib.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_macros.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkblas.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkmemory.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkpqueue2.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkpqueue.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkrandom.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mksort.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkutils.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_proto.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gkregex.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gkregex.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_struct.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_types.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_graph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/htable.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/io.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/itemsets.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/mcore.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/memory.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/ms_inttypes.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/ms_stat.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/ms_stdint.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/omp.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/pdb.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/pqueue.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/random.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/rw.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/seq.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/sort.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/string.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/timers.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/tokenizer.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_util.c
)

set(MetisIncludeSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/include/metis.h
)

set(MetisLibSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/auxapi.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/balance.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/bucketsort.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/checkgraph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/coarsen.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/compress.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/contig.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/debug.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/fm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/fortran.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/frename.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/gklib.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/gklib_defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/gklib_rename.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/graph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/initpart.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/kmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/kwayfm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/kwayrefine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/macros.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/mcutil.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/mesh.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/meshpart.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/metislib.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/minconn.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/mincover.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/mmd.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/ometis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/options.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/parmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/pmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/proto.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/refine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/rename.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/separator.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/sfm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/srefine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/stat.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/stdheaders.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/struct.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/timing.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/util.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/wspace.c
)

set(MetisProgramSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/programs/metisbin.h
)

add_library(quadrature_integration STATIC ${QuadratureIntegrationSrc})
target_include_directories(quadrature_integration PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/)

add_library(tinyxml STATIC ${TinyXmlSrc})
target_compile_definitions(tinyxml PUBLIC TIXML_USE_STL)


add_library(metis STATIC ${MetisGKlibSrc} ${MetisIncludeSrc} ${MetisLibSrc})
target_include_directories(metis 
    PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib 
    PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis 
    INTERFACE ${CMAKE_SOURCE_DIR}/3rdparty/metis/programs
)
target_compile_definitions(metis PUBLIC USE_GKREGEX)


set(VtkSrc
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/BasicVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/ContactVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/ContactVtkWriter.inl
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/MeshVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/SnapshotVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/VtkReader.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/CellInfoVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/CellInfoVtkWriter.inl
  ${CMAKE_SOURCE_DIR}/src/IO/Segy/SegySeismo.h
  ${CMAKE_SOURCE_DIR}/src/IO/Segy/SegySeismo.inl
)

set(UtilsSrc
  ${CMAKE_SOURCE_DIR}/src/Utils/Base64.h
  ${CMAKE_SOURCE_DIR}/src/Utils/Utils.h
)

set(MathsSrc
  ${CMAKE_SOURCE_DIR}/src/Maths/AABB2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/AABB3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Coordinates2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Coordinates3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Matrix2x2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Matrix3x3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/MatrixInv.h
  ${CMAKE_SOURCE_DIR}/src/Maths/MatrixMaths.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Random.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Spaces.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Util.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Vector2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Vector3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Tensor2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Tensor3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/VectorMaths.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Collisions.h
)

set(DifferentialSolversSrc
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/DifferentialSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/DifferentialSystem.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/AdaptiveRKSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/DormandPrinceSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/EulerSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/RK2Solver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/RungeSolver4.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/TrapezoidSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/RungeSolver8_10.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/RungeSolver12_14.h
)

set(NonlinearSolversSrc
  ${CMAKE_SOURCE_DIR}/src/NonlinearSolvers/NonlinearSystemSolver.h
  ${CMAKE_SOURCE_DIR}/src/NonlinearSolvers/NonlinearFunctionSystem.h
  ${CMAKE_SOURCE_DIR}/src/NonlinearSolvers/Solvers/SimpleIterations.h
  ${CMAKE_SOURCE_DIR}/src/NonlinearSolvers/Solvers/KrylovProjector.h
)

set(NetworkSrc
  ${CMAKE_SOURCE_DIR}/src/Network/NetworkInterface.h
  ${CMAKE_SOURCE_DIR}/src/Network/NetworkInterface.cpp
)

add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/Task)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/ResultCombiner)