cmake_minimum_required(VERSION 3.30)
project(Volume)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_COMPILER /usr/bin/gcc)
set(CMAKE_CXX_COMPILER /usr/bin/g++)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -funroll-loops -Wno-ignored-attributes -Wno-switch")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native -funroll-loops -Wno-ignored-attributes -Wno-switch")
if(DEBUG_SYMBOLS)
  message(STATUS "Compiling with -g")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()
if(MINGW)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif()

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

message("Generating makefiles for: " ${CMAKE_SYSTEM} "\n")
message("Build type: " ${CMAKE_BUILD_TYPE} "\n")

FIND_PACKAGE(MPI REQUIRED)
if (MPI_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)

FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif (OPENMP_FOUND)

FIND_PACKAGE(ZLIB REQUIRED)
#IF (ZLIB_FOUND)
#    MESSAGE(STATUS "Setting up ZLIB")    
#    include_directories( ${ZLIB_INCLUDE_DIRS} )    
#ENDIF (ZLIB_FOUND)

#FIND_PACKAGE(Eigen REQUIRED NO_MODULE HINTS "${CMAKE_SOURCE_DIR}/3rdparty/eigen-5.0.0/cmake")
#IF (EIGEN_FOUND)
#   MESSAGE(STATUS "Setting up eigen")    
#ENDIF (EIGEN_FOUND)


set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/bin)

MAKE_DIRECTORY(${CMAKE_SOURCE_DIR}/bin)

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

set(TinyXmlSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinystr.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxml.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxmlerror.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxmlparser.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinystr.h
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxml.h
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

set(QuadratureIntegrationSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/tetrahedron_arbq_rule.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/triangle_fekete_rule.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/legendre_rule.h
)

set(GKlibSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/b64.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/blas.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/cache.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/csr.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/error.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/evaluate.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/fkvkselect.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/fs.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/getopt.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/graph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/htable.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/io.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/itemsets.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/mcore.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/memory.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/pqueue.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/random.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/rw.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/seq.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/sort.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/string.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/timers.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/tokenizer.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/gk_util.c
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src/gkregex.c
)

set(GKlibInclude
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_macros.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/GKlib.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_mkblas.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_mkmemory.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_mkpqueue2.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_mkpqueue.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_mkrandom.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_mksort.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_mkutils.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_ms_inttypes.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_ms_stat.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_ms_stdint.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_proto.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gkregex.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_struct.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_types.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_arch.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_externs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include/gk_getopt.h
)

set(MetisLibSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/auxapi.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/balance.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/bucketsort.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/checkgraph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/coarsen.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/compress.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/contig.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/debug.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/fm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/fortran.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/frename.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/gklib.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/gklib_defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/gklib_rename.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/graph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/initpart.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/kmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/kwayfm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/kwayrefine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/macros.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/mcutil.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/mesh.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/meshpart.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/metislib.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/minconn.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/mincover.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/mmd.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/ometis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/options.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/parmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/pmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/proto.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/refine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/rename.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/separator.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/sfm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/srefine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/stat.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/stdheaders.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/struct.h
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/timing.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/util.c
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis/wspace.c
)

set(MetisProgram
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/programs/metisbin.h
)

set(MetisInclude
  ${CMAKE_SOURCE_DIR}/3rdparty/METIS/include/metis.h
)

set(EigenSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Dense}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Householder}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/LU}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/PaStiXSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/SparseCholesky}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/SPQRSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/SuperLUSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Version}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Cholesky}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Eigen}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/IterativeLinearSolvers}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/MetisSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Core}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/QR}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/SparseCore}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/StdDeque}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/SVD}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/CholmodSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Eigenvalues}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Jacobi}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/OrderingMethods}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/QtAlignedMalloc}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/SparseLU}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/ThreadPool}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/StdList}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/AccelerateSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/UmfPackSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/StdVector}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/KLUSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/PardisoSupport}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/SparseQR}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Sparse}
  ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/Eigen/Geometry}
)

add_library(quadrature_integration STATIC ${QuadratureIntegrationSrc})
target_include_directories(quadrature_integration PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/)

add_library(tinyxml STATIC ${TinyXmlSrc})
target_include_directories(tinyxml PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/)
#target_compile_definitions(tinyxml PUBLIC TIXML_USE_STL)

add_library(metis STATIC ${GKlibSrc} ${GKlibInclude} ${MetisLibSrc} ${MetisInclude})
target_include_directories(metis
  PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/METIS/libmetis
  INTERFACE ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/include
  PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/GKlib/src
  INTERFACE ${CMAKE_SOURCE_DIR}/3rdparty/METIS/programs
  INTERFACE ${CMAKE_SOURCE_DIR}/3rdparty/METIS/include
)
target_compile_definitions(metis
    PUBLIC IDXTYPEWIDTH=32
    PUBLIC REALTYPEWIDTH=32
)

add_library(eigen INTERFACE)
add_library (Eigen3::Eigen ALIAS eigen)
target_include_directories(eigen INTERFACE ${CMAKE_SOURCE_DIR}/3rdparty/eigen5/)
target_compile_definitions(eigen
    INTERFACE EIGEN_DEFAULT_TO_ROW_MAJOR 
    INTERFACE EIGEN_DONT_PARALLELIZE 
    INTERFACE EIGEN_NO_DEBUG 
    INTERFACE EIGEN_NO_AUTOMATIC_RESIZING 
    INTERFACE EIGEN_NO_MALLOC
)

add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/Task)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/ResultCombiner)
if(COMPILE_TESTS)
  message(STATUS "Compiling with googletest testcases")
  add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/googletest)
  add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/Testing)
endif()
