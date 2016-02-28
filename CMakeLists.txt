cmake_minimum_required(VERSION 2.6)
project(cfd_2d)

# необходимо прописать путь к библиотеке MPI
set(MPI_INCLUDE "/usr/include/mpich")
set(MPI_LIB "${CONTRIB_DIR}/boxlib/lib")

set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

include_directories(
		${MPI_INCLUDE}
		src src/methods
		src/methods/solvers/
		src/mesh
		src/tinyxml)
link_directories(${MPI_LIB})

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
)

set(TINYXML_SRC
		/src/tinyxml/tinystr.cpp
		/src/tinyxml/tinyxml.cpp
)

add_executable(	cfd-2d 
		src/main.cpp 
		src/global.cpp
		src/bnd_cond.cpp
		${SOLVERS_SRC}
		${METHODS_SRC}
		${MESH_SRC}
)

target_link_libraries(cfd-2d  mpi)
install(TARGETS cfd-2d RUNTIME DESTINATION bin)
