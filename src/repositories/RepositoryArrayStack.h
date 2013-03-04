// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _PEANOCLAW_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 
#define _PEANOCLAW_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 


#include "repositories/Repository.h"
#include "records/RepositoryState.h"

#include "State.h"
#include "Vertex.h"
#include "Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/VertexArrayStack.h"


 #include "adapters/InitialiseGrid.h" 
 #include "adapters/Plot.h" 
 #include "adapters/Remesh.h" 
 #include "adapters/SolveTimestep.h" 
 #include "adapters/SolveTimestepAndPlot.h" 
 #include "adapters/GatherCurrentSolution.h" 
 #include "adapters/Cleanup.h" 



namespace peanoclaw {
      namespace repositories {
        class RepositoryArrayStack;  
      }
}


class peanoclaw::repositories::RepositoryArrayStack: public peanoclaw::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellArrayStack<peanoclaw::Cell>       CellStack;
    typedef peano::stacks::VertexArrayStack<peanoclaw::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    peanoclaw::State          _solverState;
    peano::grid::RegularGridContainer<peanoclaw::Vertex,peanoclaw::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<peanoclaw::Vertex,peanoclaw::Cell,peanoclaw::State,VertexStack,CellStack,peanoclaw::adapters::InitialiseGrid> _gridWithInitialiseGrid;
    peano::grid::Grid<peanoclaw::Vertex,peanoclaw::Cell,peanoclaw::State,VertexStack,CellStack,peanoclaw::adapters::Plot> _gridWithPlot;
    peano::grid::Grid<peanoclaw::Vertex,peanoclaw::Cell,peanoclaw::State,VertexStack,CellStack,peanoclaw::adapters::Remesh> _gridWithRemesh;
    peano::grid::Grid<peanoclaw::Vertex,peanoclaw::Cell,peanoclaw::State,VertexStack,CellStack,peanoclaw::adapters::SolveTimestep> _gridWithSolveTimestep;
    peano::grid::Grid<peanoclaw::Vertex,peanoclaw::Cell,peanoclaw::State,VertexStack,CellStack,peanoclaw::adapters::SolveTimestepAndPlot> _gridWithSolveTimestepAndPlot;
    peano::grid::Grid<peanoclaw::Vertex,peanoclaw::Cell,peanoclaw::State,VertexStack,CellStack,peanoclaw::adapters::GatherCurrentSolution> _gridWithGatherCurrentSolution;
    peano::grid::Grid<peanoclaw::Vertex,peanoclaw::Cell,peanoclaw::State,VertexStack,CellStack,peanoclaw::adapters::Cleanup> _gridWithCleanup;

  
   peanoclaw::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureInitialiseGridCPUTime;
    tarch::timing::Measurement _measurePlotCPUTime;
    tarch::timing::Measurement _measureRemeshCPUTime;
    tarch::timing::Measurement _measureSolveTimestepCPUTime;
    tarch::timing::Measurement _measureSolveTimestepAndPlotCPUTime;
    tarch::timing::Measurement _measureGatherCurrentSolutionCPUTime;
    tarch::timing::Measurement _measureCleanupCPUTime;

    tarch::timing::Measurement _measureInitialiseGridCalendarTime;
    tarch::timing::Measurement _measurePlotCalendarTime;
    tarch::timing::Measurement _measureRemeshCalendarTime;
    tarch::timing::Measurement _measureSolveTimestepCalendarTime;
    tarch::timing::Measurement _measureSolveTimestepAndPlotCalendarTime;
    tarch::timing::Measurement _measureGatherCurrentSolutionCalendarTime;
    tarch::timing::Measurement _measureCleanupCalendarTime;


  public:
    RepositoryArrayStack(
      peano::geometry::Geometry&                   geometry,
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset,
      int                                          maximumSizeOfCellInOutStack,
      int                                          maximumSizeOfVertexInOutStack,
      int                                          maximumSizeOfVertexTemporaryStack
    );
    
    /**
     * Parallel Constructor
     *
     * Used in parallel mode only where the size of the domain is not known 
     * when the type of repository is determined.  
     */
    RepositoryArrayStack(
      peano::geometry::Geometry&                   geometry,
      int                                          maximumSizeOfCellInOutStack,
      int                                          maximumSizeOfVertexInOutStack,
      int                                          maximumSizeOfVertexTemporaryStack
    );
    
    virtual ~RepositoryArrayStack();

    virtual void restart(
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
      int                                          domainLevel
    );
         
    virtual void terminate();
        
    virtual peanoclaw::State& getState();

    virtual void iterate( bool reduceState = true );
    
    virtual void writeCheckpoint(peano::grid::Checkpoint<peanoclaw::Vertex, peanoclaw::Cell> * const checkpoint); 
    virtual void readCheckpoint( peano::grid::Checkpoint<peanoclaw::Vertex, peanoclaw::Cell> const * const checkpoint );
    virtual peano::grid::Checkpoint<peanoclaw::Vertex, peanoclaw::Cell>* createEmptyCheckpoint(); 

    virtual void switchToInitialiseGrid();    
    virtual void switchToPlot();    
    virtual void switchToRemesh();    
    virtual void switchToSolveTimestep();    
    virtual void switchToSolveTimestepAndPlot();    
    virtual void switchToGatherCurrentSolution();    
    virtual void switchToCleanup();    

    virtual bool isActiveAdapterInitialiseGrid() const;
    virtual bool isActiveAdapterPlot() const;
    virtual bool isActiveAdapterRemesh() const;
    virtual bool isActiveAdapterSolveTimestep() const;
    virtual bool isActiveAdapterSolveTimestepAndPlot() const;
    virtual bool isActiveAdapterGatherCurrentSolution() const;
    virtual bool isActiveAdapterCleanup() const;

     
    #ifdef Parallel
    virtual bool continueToIterate();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics() const;
    virtual void clearIterationStatistics();
};


#endif