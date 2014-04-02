cmake_minimum_required(VERSION 2.8)
project(PEANOCLAW)

#########################################
# Peano Settings
#########################################

set(PEANO_DIMENSIONS 2 CACHE STRING "number of dimensions")

find_package(HDF5 COMPONENTS HL C)
find_package(VTK)

# for mekka scenario
find_package(PkgConfig)
pkg_search_module(NETCDF netcdf>=4)

#include_directories(${HDF5_INCLUDE_DIRS})

#add_definitions(-DDim${PEANO_DIMENSIONS} -DParallel -DSharedTBB -DSharedMemoryParallelisation -DAsserts -DDebug -DTrace)
add_definitions(-DDim${PEANO_DIMENSIONS})
add_definitions(-std=c++11)

option(PEANO_PARALLEL "shall we build a MPI enabled Peano")
if (PEANO_PARALLEL)
    find_package(MPI)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    add_definitions(-DParallel)
endif (PEANO_PARALLEL)

option(PEANO_TBB "shall we build a TBB enabled Peano")
if (PEANO_TBB)
    set(TBB_ROOT "/work/wittmanr/mytools/tbb41_20121003oss")
    set(TBB_INCLUDE_DIRS "${TBB_ROOT}/include")
    set(TBB_LIBRARY_DIRS "${TBB_ROOT}/build/linux_intel64_gcc_cc4.6_libc2.15_kernel3.7.1_release")
    set(TBB_LINK_FLAGS "-ltbb")

    include_directories(${TBB_INCLUDE_DIRS})
    link_directories(${TBB_LIBRARY_DIRS})
    add_definitions(-DSharedTBB)
endif (PEANO_TBB)

option(PEANO_OMP "shall we build an OMP enabled Peano")
if (PEANO_OMP)
    add_definitions(-DSharedOMP)
endif (PEANO_OMP)

option(PEANO_ASSERTS "shall we build with Assertions enabled")
if (PEANO_ASSERTS)
    add_definitions(-DAsserts)
endif (PEANO_ASSERTS)

option(PEANO_TRACE "shall we build with trace enabled")
if (PEANO_TRACE)
     add_definitions(-DTrace)
endif (PEANO_TRACE)

option(PEANO_DEBUG "shall we build with Debug enabled")
if (PEANO_DEBUG)
    add_definitions(-DDebug)
endif (PEANO_DEBUG)

option(PEANO_NORECURSIONUNROLLING "disable recursion unrolling")
if (PEANO_NORECURSIONUNROLLING)
    add_definitions(-DnoUseRecursionUnrollingOnRegularPatches)
endif (PEANO_NORECURSIONUNROLLING)

option(PEANO_NOCACHEACTIONSETS "do not cache action sets")
if (PEANO_NOCACHEACTIONSETS)
    add_definitions(-DnoCacheActionSets)
endif (PEANO_NOCACHEACTIONSETS)

option(PEANO_NODISTINGUISHSUBTREES "do not distinguish stationary and instationary subtrees")
if (PEANO_NODISTINGUISHSUBTREES)
    add_definitions(-DnoDistinguishStationaryAndInstationarySubtrees)
endif (PEANO_NODISTINGUISHSUBTREES)

option(PEANO_LOGOFF "disable noise")
if (PEANO_LOGOFF)
    add_definitions(-DLogOff)
endif (PEANO_LOGOFF)

option (SCALASCA "instrument with scalasca?")
if (SCALASCA)
    set(CMAKE_C_COMPILER skin_mpicc)
    set(CMAKE_CXX_COMPILER skin_mpicxx)
    add_definitions(-DParallel)
endif (SCALASCA)

if (HDF5_FOUND)
    add_definitions(-DWITH_HDF5)
endif(HDF5_FOUND)

#########################################

include_directories(peanokernel)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(peanokernel)

include_directories(peanotoolboxes/ControlLoopLoadBalancer)
add_subdirectory(peanotoolboxes/ControlLoopLoadBalancer)

#########################################

set(Python_ADDITIONAL_VERSIONS "2.7")
find_package(PythonLibs)

option (PEANOCLAW_PYCLAW "enable PyClaw support")
option (PEANOCLAW_SWE "enable SWE support")
option (PEANOCLAW_FULLSWOF2D "enable FullSWOF2D support")

if (PEANOCLAW_PYCLAW)
    add_definitions(-DPEANOCLAW_PYCLAW)

    set (PEANOCLAW_PYCLAW_SOURCES
        src/peanoclaw/pyclaw/PyClaw.cpp
        src/peanoclaw/pyclaw/PyClawState.cpp
        src/peanoclaw/pyclaw/pyclawBindings.cpp
        src/peanoclaw/pyclaw/FluxCorrectionCallbackWrapper.cpp
        src/peanoclaw/pyclaw/InterpolationCallbackWrapper.cpp
        src/peanoclaw/pyclaw/RestrictionCallbackWrapper.cpp
    )
endif (PEANOCLAW_PYCLAW)

if (PEANOCLAW_SWE)
    add_definitions(-DPEANOCLAW_SWE)

    set(PEANOCLAW_SWE_SOURCES
        src/peanoclaw/native/SWEKernel.cpp   
        src/peanoclaw/native/SWE_WavePropagationBlock_patch.cpp

        # scenarios
        src/peanoclaw/native/BreakingDam.cpp
        src/peanoclaw/native/MekkaFlood.cpp
        src/peanoclaw/native/BathymetryHelper.cpp

        # swe part
        src/peanoclaw/native/SWE/src/blocks/SWE_Block.cpp
        src/peanoclaw/native/SWE/src/blocks/SWE_WavePropagationBlock.cpp
    )
endif(PEANOCLAW_SWE)

if (PEANOCLAW_FULLSWOF2D)
    add_definitions(-DPEANOCLAW_FULLSWOF2D)

    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libboundaryconditions)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libflux)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libfrictions)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libinitializations)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/liblimitations)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libparameters)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libparser)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/librain_infiltration)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libreconstructions)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libsave)
    include_directories(src/peanoclaw/native/FullSWOF_2D/Headers/libschemes)

    set(PEANOCLAW_FULLSWOF2D_SOURCES
        src/peanoclaw/native/FullSWOF_2D/Sources/libboundaryconditions/bc_imp_discharge.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libboundaryconditions/bc_imp_height.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libboundaryconditions/bc_neumann.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libboundaryconditions/bc_periodic.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libboundaryconditions/bc_wall.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libboundaryconditions/boundary_condition.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libboundaryconditions/choice_condition.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libflux/choice_flux.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libflux/f_hll2.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libflux/f_hll.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libflux/flux.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libflux/f_rusanov.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libfrictions/choice_friction.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libfrictions/fr_darcy_weisbach.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libfrictions/friction.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libfrictions/fr_manning.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libfrictions/no_friction.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/choice_init_huv.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/choice_init_topo.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/huv_generated.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/huv_generated_radial_dam_dry.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/huv_generated_radial_dam_wet.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/huv_generated_thacker.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/huv_read.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/initialization_huv.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/initialization_topo.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/topo_generated_flat.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/topo_generated_thacker.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libinitializations/topo_read.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/liblimitations/choice_limiter.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/liblimitations/limiter.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/liblimitations/minmod.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/liblimitations/vanalbada.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/liblimitations/vanleer.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libparameters/parameters.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libparser/parser.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/choice_infiltration.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/choice_rain.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/greenampt.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/infiltration.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/no_infiltration.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/no_rain.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/rain.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/rain_generated.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/librain_infiltration/rain_read.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libreconstructions/choice_reconstruction.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libreconstructions/eno.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libreconstructions/eno_mod.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libreconstructions/hydrostatic.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libreconstructions/muscl.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libreconstructions/reconstruction.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libsave/choice_output.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libsave/gnuplot.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libsave/no_evolution_file.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libsave/output.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libsave/vtk_out.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libschemes/choice_scheme.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libschemes/order1.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libschemes/order2.cpp
        src/peanoclaw/native/FullSWOF_2D/Sources/libschemes/scheme.cpp

        src/peanoclaw/native/FullSWOF2D.cpp

        src/peanoclaw/native/dem.cpp

        # "standalone" of FullSWOF2D combination
        src/peanoclaw/native/MekkaFlood_solver.cpp
    )
endif(PEANOCLAW_FULLSWOF2D)

set(PEANOCLAW_SOURCES
    src/peanoclaw/Area.cpp
    src/peanoclaw/Cell.cpp
    src/peanoclaw/State.cpp
    src/peanoclaw/Patch.cpp
    src/peanoclaw/PatchPlotter.cpp
    src/peanoclaw/Vertex.cpp
    src/peanoclaw/Numerics.cpp
    src/peanoclaw/NumericsFactory.cpp
    src/peanoclaw/ParallelSubgrid.cpp

    src/peanoclaw/interSubgridCommunication/GhostLayerCompositor.cpp
    src/peanoclaw/interSubgridCommunication/GhostlayerCompositorFunctors.cpp
    src/peanoclaw/interSubgridCommunication/DefaultRestriction.cpp
    src/peanoclaw/interSubgridCommunication/DefaultInterpolation.cpp
    src/peanoclaw/interSubgridCommunication/DefaultFluxCorrection.cpp
    src/peanoclaw/interSubgridCommunication/Extrapolation.cpp
    src/peanoclaw/interSubgridCommunication/GridLevelTransfer.cpp
    src/peanoclaw/interSubgridCommunication/aspects/AdjacentSubgrids.cpp
    src/peanoclaw/interSubgridCommunication/aspects/AdjacentVertices.cpp

    src/peanoclaw/adapters/Cleanup.cpp
    src/peanoclaw/adapters/GatherCurrentSolution.cpp
    src/peanoclaw/adapters/GatherCurrentSolutionAndValidateGrid.cpp
    src/peanoclaw/adapters/InitialiseGrid.cpp
    src/peanoclaw/adapters/InitialiseAndValidateGrid.cpp
    src/peanoclaw/adapters/Plot.cpp
    src/peanoclaw/adapters/PlotAndValidateGrid.cpp
    src/peanoclaw/adapters/Remesh.cpp
    src/peanoclaw/adapters/SolveTimestep.cpp
    src/peanoclaw/adapters/SolveTimestepAndValidateGrid.cpp
    src/peanoclaw/adapters/SolveTimestepAndPlot.cpp
    src/peanoclaw/adapters/SolveTimestepAndPlotAndValidateGrid.cpp
    src/peanoclaw/configurations/PeanoClawConfigurationForSpacetreeGrid.cpp

    src/peanoclaw/grid/TimeIntervals.cpp

    src/peanoclaw/mappings/Cleanup.cpp
    src/peanoclaw/mappings/GatherCurrentSolution.cpp
    src/peanoclaw/mappings/InitialiseGrid.cpp
    src/peanoclaw/mappings/Plot.cpp
    src/peanoclaw/mappings/Remesh.cpp
    src/peanoclaw/mappings/SolveTimestep.cpp
    src/peanoclaw/mappings/ValidateGrid.cpp

    src/peanoclaw/records/Cell.cpp
    src/peanoclaw/records/CellDescription.cpp
    src/peanoclaw/records/Data.cpp
    src/peanoclaw/records/RepositoryState.cpp
    src/peanoclaw/records/State.cpp
    src/peanoclaw/records/Vertex.cpp
    src/peanoclaw/records/VertexDescription.cpp
    src/peanoclaw/records/PatchDescription.cpp

    src/peanoclaw/repositories/RepositoryArrayStack.cpp
    src/peanoclaw/repositories/RepositoryFactory.cpp
    src/peanoclaw/repositories/RepositorySTDStack.cpp
 
    src/peanoclaw/runners/PeanoClawLibraryRunner.cpp

    #src/peanoclaw/tests/TestCase.cpp

    src/peanoclaw/statistics/LevelInformation.cpp
    src/peanoclaw/statistics/LevelStatistics.cpp
    src/peanoclaw/statistics/Probe.cpp
    src/peanoclaw/statistics/PatchDescriptionDatabase.cpp
    src/peanoclaw/statistics/ParallelGridValidator.cpp
    src/peanoclaw/statistics/SubgridStatistics.cpp
    src/peanoclaw/statistics/ParallelStatistics.cpp
    src/peanoclaw/statistics/TimeIntervalStatistics.cpp

    src/peanoclaw/parallel/MasterWorkerAndForkJoinCommunicator.cpp
    src/peanoclaw/parallel/NeighbourCommunicator.cpp
    src/peanoclaw/parallel/SubgridCommunicator.cpp
    src/peanoclaw/parallel/PeanoClawControlLoopLoadBalancer.cpp

    src/peanoclaw/native/NativeKernel.cpp   
)

execute_process ( COMMAND python2 -c "from numpy import get_include; print get_include()" OUTPUT_VARIABLE NUMPY_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/peanoclaw/native/SWE/src
        ${PYTHON_INCLUDE_DIRS}
        ${NUMPY_INCLUDE_DIR}
)


if (PEANOCLAW_PYCLAW)
    set (PEANOCLAW_SOURCES ${PEANOCLAW_SOURCES} ${PEANOCLAW_PYCLAW_SOURCES})
    set (PEANOCLAW_LIBRARIES ${PEANOCLAW_LIBRARIES} ${PYTHON_LIBRARIES})
endif (PEANOCLAW_PYCLAW)

if (PEANOCLAW_SWE)
    add_definitions(-DSWE -DWAVE_PROPAGATION_SOLVER=1) # TODO: 4 was slow, why?
    set (PEANOCLAW_SOURCES ${PEANOCLAW_SOURCES} ${PEANOCLAW_SWE_SOURCES})
    set (PEANOCLAW_LIBRARIES ${PEANOCLAW_LIBRARIES})
endif (PEANOCLAW_SWE) 

if (PEANOCLAW_FULLSWOF2D)
    set (PEANOCLAW_SOURCES ${PEANOCLAW_SOURCES} ${PEANOCLAW_FULLSWOF2D_SOURCES})
    set (PEANOCLAW_LIBRARIES ${PEANOCLAW_LIBRARIES})
endif (PEANOCLAW_FULLSWOF2D) 

if (HDF5_FOUND)
    set (PEANOCLAW_LIBRARIES ${PEANOCLAW_LIBRARIES} -lhdf5)
endif (HDF5_FOUND)

if (PEANO_TBB)
    set (PEANO_LIBRARIES ${PEANOCLAW_LIBRARIES} ${TBB_LINK_FLAGS} -lrt)
endif (PEANO_TBB)

set (TOOLBOXES_LIBRARIES
    peano_controllooploadbalancer_shared
)

add_library (peano-claw-${PEANO_DIMENSIONS}d SHARED ${PEANOCLAW_SOURCES})
add_library (peanoclaw_static STATIC ${PEANOCLAW_SOURCES})
add_library (peanoclaw_static_norle STATIC ${PEANOCLAW_SOURCES})
SET_TARGET_PROPERTIES(peanoclaw_static_norle PROPERTIES COMPILE_FLAGS "-DnoPackedEmptyHeapMessages" )

target_link_libraries (peano-claw-${PEANO_DIMENSIONS}d peano_shared
    peano_tarch_shared ${PEANOCLAW_LIBRARIES} ${TOOLBOXES_LIBRARIES}
    ${PEANO_LIBRARIES} ${NETCDF_LIBRARIES} -lpng)
target_link_libraries (peanoclaw_static peano_static peano_tarch_static ${PEANOCLAW_LIBRARIES} ${TOOLBOXES_LIBRARIES} ${PEANO_LIBRARIES})
target_link_libraries (peanoclaw_static_norle peano_static peano_tarch_static ${PEANOCLAW_LIBRARIES} ${TOOLBOXES_LIBRARIES} ${PEANO_LIBRARIES})

if (PEANO_OMP)
    SET_TARGET_PROPERTIES(peano-claw-${PEANO_DIMENSIONS}d PROPERTIES COMPILE_FLAGS "-fopenmp")
endif (PEANO_OMP)
SET_TARGET_PROPERTIES(peano-claw-${PEANO_DIMENSIONS}d PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

if (PEANOCLAW_PYCLAW)
 
    # placeholder

endif (PEANOCLAW_PYCLAW)

if (PEANOCLAW_SWE OR PEANOCLAW_FULLSWOF2D)
    add_executable (peanoclaw_run src/peanoclaw/native/main.cpp)
    add_executable (peanoclaw_run_norle src/peanoclaw/native/main.cpp)

    target_link_libraries (peanoclaw_run peanoclaw_static ${NETCDF_LIBRARIES} -lpng)
    target_link_libraries (peanoclaw_run_norle peanoclaw_static_norle ${NETCDF_LIBRARIES} -lpng)
endif (PEANOCLAW_SWE OR PEANOCLAW_FULLSWOF2D)



## tools

if (PEANOCLAW_FULLSWOF2D)
    add_executable (FullSWOF2D_Kernel_Benchmark tools/FullSWOF2D_Kernel_Benchmark.cpp)
    target_link_libraries (FullSWOF2D_Kernel_Benchmark peanoclaw_static
        peano_static peano_tarch_static ${PEANOCLAW_LIBRARIES}  ${TOOLBOXES_LIBRARIES} ${PEANO_LIBRARIES} ${NETCDF_LIBRARIES})
endif (PEANOCLAW_FULLSWOF2D)

