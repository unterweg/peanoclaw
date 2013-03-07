find_package(PythonLibs)

set(PEANOCLAW_SOURCES
    src/peanoclaw/Area.cpp
    src/peanoclaw/Cell.cpp
    src/peanoclaw/State.cpp
    src/peanoclaw/Patch.cpp
    src/peanoclaw/PatchPlotter.cpp
    src/peanoclaw/Vertex.cpp
 
    src/peanoclaw/pyclaw/PyClaw.cpp
    src/peanoclaw/pyclaw/pyclawBindings.cpp

    src/peanoclaw/interSubgridCommunication/GhostLayerCompositor.cpp
    src/peanoclaw/interSubgridCommunication/Restriction.cpp
    src/peanoclaw/interSubgridCommunication/Interpolation.cpp
    src/peanoclaw/interSubgridCommunication/GridLevelTransfer.cpp

    src/peanoclaw/adapters/Cleanup.cpp
    src/peanoclaw/adapters/GatherCurrentSolution.cpp
    src/peanoclaw/adapters/InitialiseGrid.cpp
    src/peanoclaw/adapters/Plot.cpp
    src/peanoclaw/adapters/Remesh.cpp
    src/peanoclaw/adapters/SolveTimestep.cpp
    src/peanoclaw/adapters/SolveTimestepAndPlot.cpp
 
    src/peanoclaw/configurations/PeanoClawConfigurationForSpacetreeGrid.cpp

    src/peanoclaw/mappings/Cleanup.cpp
    src/peanoclaw/mappings/GatherCurrentSolution.cpp
    src/peanoclaw/mappings/InitialiseGrid.cpp
    src/peanoclaw/mappings/Plot.cpp
    src/peanoclaw/mappings/Remesh.cpp
    src/peanoclaw/mappings/SolveTimestep.cpp

    src/peanoclaw/records/Cell.cpp
    src/peanoclaw/records/CellDescription.cpp
    src/peanoclaw/records/Data.cpp
    src/peanoclaw/records/RepositoryState.cpp
    src/peanoclaw/records/State.cpp
    src/peanoclaw/records/Vertex.cpp
    src/peanoclaw/records/VertexDescription.cpp

    src/peanoclaw/repositories/RepositoryArrayStack.cpp
    src/peanoclaw/repositories/RepositoryFactory.cpp
    src/peanoclaw/repositories/RepositorySTDStack.cpp
 
    src/peanoclaw/runners/PeanoClawLibraryRunner.cpp

    #src/peanoclaw/tests/TestCase.cpp

    src/peanoclaw/statistics/LevelInformation.cpp
    src/peanoclaw/statistics/Probe.cpp
    )



include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${PYTHON_INCLUDE_DIRS}
        /usr/lib64/python2.7/site-packages/numpy/core/include # numpy.get_include()
)

add_library (peanoclaw SHARED ${PEANOCLAW_SOURCES})
target_link_libraries (peanoclaw peano_shared peano_tarch_shared ${PYTHON_LIBRARIES} ${TBB_LINK_FLAGS} -lrt)

SET_TARGET_PROPERTIES(peanoclaw PROPERTIES LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
