find_package(PythonLibs)

option (PEANOCLAW_PYCLAW "enable PyClaw support")
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

    src/peanoclaw/parallel/MasterWorkerAndForkJoinCommunicator.cpp
    src/peanoclaw/parallel/NeighbourCommunicator.cpp

    src/peanoclaw/native/NativeKernel.cpp   
    src/peanoclaw/native/SWEKernel.cpp   
    src/peanoclaw/native/SWE_WavePropagationBlock_patch.cpp

    # swe part
    src/peanoclaw/native/SWE/src/blocks/SWE_Block.cpp
    src/peanoclaw/native/SWE/src/blocks/SWE_WavePropagationBlock.cpp

    # load balancing part
    src/peanoclaw/PeanoClawControlLoopLoadBalancer.cpp
)

execute_process ( COMMAND python -c "from numpy import get_include; print get_include()" OUTPUT_VARIABLE NUMPY_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/peanoclaw/native/SWE/src
        ${PYTHON_INCLUDE_DIRS}
        ${NUMPY_INCLUDE_DIR}
)

add_definitions(-DSWE -DWAVE_PROPAGATION_SOLVER=1) # TODO: 4 was slow, why?
#add_definitions(-std=c++0x)

if (PEANOCLAW_PYCLAW)
    set (PEANOCLAW_SOURCES ${PEANOCLAW_SOURCES} ${PEANOCLAW_PYCLAW_SOURCES})
    set (PEANOCLAW_LIBRARIES ${PEANOCLAW_LIBRARIES} ${PYTHON_LIBRARIES})
endif (PEANOCLAW_PYCLAW) 

if (HDF5_FOUND)
    set (PEANOCLAW_LIBRARIES ${PEANOCLAW_LIBRARIES} -lhdf5)
endif (HDF5_FOUND)

if (PEANO_TBB)
    set (PEANO_LIBRARIES ${PEANOCLAW_LIBRARIES} ${TBB_LINK_FLAGS} -lrt)
endif (PEANO_TBB)

add_library (peanoclaw SHARED ${PEANOCLAW_SOURCES})
add_library (peanoclaw_static STATIC ${PEANOCLAW_SOURCES})
add_library (peanoclaw_static_norle STATIC ${PEANOCLAW_SOURCES})
SET_TARGET_PROPERTIES(peanoclaw_static_norle PROPERTIES COMPILE_FLAGS "-DnoPackedEmptyHeapMessages" )

target_link_libraries (peanoclaw peano_shared peano_tarch_shared peano_toolboxes_shared ${PEANOCLAW_LIBRARIES} ${PEANO_LIBRARIES})
target_link_libraries (peanoclaw_static peano_static peano_tarch_static peano_toolboxes_static ${PEANOCLAW_LIBRARIES} ${PEANO_LIBRARIES})
target_link_libraries (peanoclaw_static_norle peano_static peano_tarch_static peano_toolboxes_static ${PEANOCLAW_LIBRARIES} ${PEANO_LIBRARIES})


if (PEANO_OMP)
    SET_TARGET_PROPERTIES(peanoclaw PROPERTIES COMPILE_FLAGS "-fopenmp")
endif (PEANO_OMP)
SET_TARGET_PROPERTIES(peanoclaw PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

add_executable (peanoclaw_run src/peanoclaw/native/main.cpp)
add_executable (peanoclaw_run_norle src/peanoclaw/native/main.cpp)

target_link_libraries (peanoclaw_run peanoclaw_static)
target_link_libraries (peanoclaw_run_norle peanoclaw_static_norle)
