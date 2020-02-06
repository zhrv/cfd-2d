cmake_minimum_required(VERSION 2.6)
project(cfd_2d)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")


set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

#
#  https://github.com/hypre-space/hypre
#
set(HYPRE_INCLUDE "${CONTRIB_DIR}/hypre/src/hypre/include")
set(HYPRE_LIB "${CONTRIB_DIR}/hypre/src/hypre/lib")


#
#  http://glaros.dtc.umn.edu/gkhome/metis/metis/download
#
set(METIS_INCLUDE "${CONTRIB_DIR}/metis-5.1.0/local/include")
set(METIS_LIB "${CONTRIB_DIR}/metis-5.1.0/local/lib")


include_directories(
		${HYPRE_INCLUDE}
		${METIS_INCLUDE}
		src
		src/methods
		src/methods/solvers
		src/mesh
		src/tinyxml)
link_directories(
		${MPI_LIB}
		${HYPRE_LIB}
		${METIS_LIB}
)

set(METHODS_SRC
		src/methods/decomp.cpp
		src/methods/solver.cpp
		src/methods/fvm_tvd.cpp
		src/methods/fvm_tvd_implicit.cpp
		src/methods/fem_rkdg.cpp
		src/methods/fem_dg_implicit.cpp
		src/methods/LimiterDG.cpp
		src/methods/LimiterDGCockburn.cpp
)

set(SOLVERS_SRC
		src/methods/solvers/CSR.cpp
		src/methods/solvers/MatrixSolver.cpp
		src/methods/solvers/SolverHypre.cpp
		src/methods/solvers/SolverHypreBoomerAmg.cpp
		src/methods/solvers/SolverHypreCustomSeidel.cpp
		src/methods/solvers/SolverHypreFlexGmres.cpp
		src/methods/solvers/SolverHypreFlexGmresPrecAMG.cpp
		src/methods/solvers/SolverHypreGmres.cpp
		src/methods/solvers/SolverHyprePcg.cpp
		src/methods/solvers/SolverZeidel.cpp
)

set(MESH_SRC
		src/mesh/grid.cpp
		src/mesh/MeshReader.cpp
		src/mesh/MeshReaderBerkleyTriangle.cpp
		src/mesh/MeshReaderSalomeUnv.cpp
		src/mesh/MeshReaderMsh2.cpp
		src/mesh/MeshReaderGmsh.cpp
)

set(TINYXML_SRC
		src/tinyxml/tinystr.cpp
		src/tinyxml/tinyxml.cpp
		src/tinyxml/tinyxmlerror.cpp
		src/tinyxml/tinyxmlparser.cpp
)

add_executable(	cfd-2d
		src/main.cpp
		src/global.cpp
		src/bnd_cond.cpp
		${SOLVERS_SRC}
		${METHODS_SRC}
		${MESH_SRC}
		${TINYXML_SRC}
)

target_link_libraries(cfd-2d HYPRE metis)
install(TARGETS cfd-2d RUNTIME DESTINATION bin)
