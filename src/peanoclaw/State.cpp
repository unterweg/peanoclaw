#include "State.h"
#include "Cell.h"
#include "Vertex.h"

#include "peanoclaw/statistics/SubgridStatistics.h"

#include "peano/grid/Checkpoint.h"

tarch::logging::Log peanoclaw::State::_log("peanoclaw::State");

peanoclaw::State::State():
  Base(),
  _numerics(0),
  _subgridStatistics(new peanoclaw::statistics::SubgridStatistics)
{
  resetGlobalTimeIntervals();
  resetMinimalTimestep();
  resetTotalNumberOfCellUpdates();

  //Flags that might not be initialized by Peano
  _stateData.setHasChangedVertexOrCellState(false);
  _stateData.setHasRefined(false);
  _stateData.setHasErased(false);
  _stateData.setHasTriggeredRefinementForNextIteration(false);
  _stateData.setHasTriggeredEraseForNextIteration(false);
  _stateData.setAllPatchesEvolvedToGlobalTimestep(false);
  #ifdef Parallel
  _stateData.setCouldNotEraseDueToDecompositionFlag(false);
  _stateData.setSubWorkerIsInvolvedInJoinOrFork(false);
  #endif
  assertion(_subgridStatistics != 0);
}


peanoclaw::State::State(const Base::PersistentState& argument):
  Base(argument),
  _numerics(0),
  _subgridStatistics(new peanoclaw::statistics::SubgridStatistics)
{
  // @todo Insert your code here
  assertion(_subgridStatistics != 0);
  _subgridStatistics->setGlobalTimestepEndTime(argument.getGlobalTimestepEndTime());
}


void peanoclaw::State::writeToCheckpoint( peano::grid::Checkpoint<peanoclaw::Vertex,peanoclaw::Cell>& checkpoint ) const {
  // @todo Insert your code here
}

    
void peanoclaw::State::readFromCheckpoint( const peano::grid::Checkpoint<peanoclaw::Vertex,peanoclaw::Cell>& checkpoint ) {
  // @todo Insert your code here
}

int peanoclaw::State::getPlotNumber() const {
  return _stateData.getPlotNumber();
}

void peanoclaw::State::setPlotNumber(int plotNumber) {
  _stateData.setPlotNumber(plotNumber);
}

std::string peanoclaw::State::getPlotName() const {
  std::string plotName;
  int i = 0;
  while(i < 32 && _stateData.getPlotName(i) > 0) {
    plotName.push_back((char)_stateData.getPlotName(i));
    i++;
  }
  return plotName;
}

void peanoclaw::State::setPlotName(const std::string& plotName) {
  for(int i = 0; i < std::min(32, (int)plotName.length()); i++) {
    _stateData.setPlotName(i, (int)plotName[i]);
  }
  if(plotName.length() < 32) {
    _stateData.setPlotName(std::min(32, (int)plotName.length()), 0);
  }
}

void peanoclaw::State::setUnknownsPerSubcell(int unknownsPerSubcell) {
  _stateData.setUnknownsPerSubcell(unknownsPerSubcell);
}

int peanoclaw::State::getUnknownsPerSubcell() const {
  return _stateData.getUnknownsPerSubcell();
}

void peanoclaw::State::setNumberOfParametersWithoutGhostlayerPerSubcell(int numberOfParametersWithoutGhostlayerPerSubcell) {
  _stateData.setNumberOfParametersWithoutGhostlayerPerSubcell(numberOfParametersWithoutGhostlayerPerSubcell);
}

int peanoclaw::State::getNumberOfParametersWithoutGhostlayerPerSubcell() const {
  return _stateData.getNumberOfParametersWithoutGhostlayerPerSubcell();
}

void peanoclaw::State::setNumberOfParametersWithGhostlayerPerSubcell(int numberOfParametersWithGhostlayerPerSubcell) {
  _stateData.setNumberOfParametersWithGhostlayerPerSubcell(numberOfParametersWithGhostlayerPerSubcell);
}

int peanoclaw::State::getNumberOfParametersWithGhostlayerPerSubcell() const {
  return _stateData.getNumberOfParametersWithGhostlayerPerSubcell();
}

void peanoclaw::State::setDefaultSubdivisionFactor(const tarch::la::Vector<DIMENSIONS, int>& defaultSubdivisionFactor) {
  _stateData.setDefaultSubdivisionFactor(defaultSubdivisionFactor);
}

tarch::la::Vector<DIMENSIONS, int> peanoclaw::State::getDefaultSubdivisionFactor() const {
  return _stateData.getDefaultSubdivisionFactor();
}

void peanoclaw::State::setDefaultGhostLayerWidth(int defaultGhostLayerWidth) {
  _stateData.setDefaultGhostWidthLayer(defaultGhostLayerWidth);
}

int peanoclaw::State::getDefaultGhostLayerWidth() const {
  return _stateData.getDefaultGhostWidthLayer();
}

tarch::la::Vector<DIMENSIONS, double> peanoclaw::State::getInitialMaximalSubgridSize() const {
  return _stateData.getInitialMaximalSubgridSize();
}

void peanoclaw::State::setInitialMaximalSubgridSize(const tarch::la::Vector<DIMENSIONS, double>& h) {
  _stateData.setInitialMaximalSubgridSize(h);
}

void peanoclaw::State::setInitialTimestepSize(double initialTimestepSize) {
  _stateData.setInitialTimestepSize(initialTimestepSize);
}

double peanoclaw::State::getInitialTimestepSize() const {
  return _stateData.getInitialTimestepSize();
}

void peanoclaw::State::setNumerics(peanoclaw::Numerics& numerics) {
  _numerics = &numerics;
}

peanoclaw::Numerics* peanoclaw::State::getNumerics() const {
  return _numerics;
}

void peanoclaw::State::setProbeList(std::vector<peanoclaw::statistics::Probe> probeList) {
  _probeList = probeList;
}

std::vector<peanoclaw::statistics::Probe>& peanoclaw::State::getProbeList() {
  return _probeList;
}

void peanoclaw::State::setGlobalTimestepEndTime(double globalTimestepEndTime) {
  _stateData.setGlobalTimestepEndTime(globalTimestepEndTime);
  assertion(_subgridStatistics != 0);
  _subgridStatistics->setGlobalTimestepEndTime(globalTimestepEndTime);
}

double peanoclaw::State::getGlobalTimestepEndTime() const {
  return _stateData.getGlobalTimestepEndTime();
}

void peanoclaw::State::setAllPatchesEvolvedToGlobalTimestep(bool value) {
  _stateData.setAllPatchesEvolvedToGlobalTimestep(value);
}

bool peanoclaw::State::getAllPatchesEvolvedToGlobalTimestep() const {
  return _stateData.getAllPatchesEvolvedToGlobalTimestep();
}

void peanoclaw::State::setDomain(const tarch::la::Vector<DIMENSIONS, double>& offset, const tarch::la::Vector<DIMENSIONS, double>& size) {
  _stateData.setDomainOffset(offset);
  _stateData.setDomainSize(size);
}

tarch::la::Vector<DIMENSIONS, double> peanoclaw::State::getDomainOffset() {
  return _stateData.getDomainOffset();
}

tarch::la::Vector<DIMENSIONS, double> peanoclaw::State::getDomainSize() {
  return _stateData.getDomainSize();
}

void peanoclaw::State::updateGlobalTimeIntervals(
  double startMaximumLocalTimeInterval,
  double endMaximumLocalTimeInterval,
  double startMinimumLocalTimeInterval,
  double endMinimumLocalTimeInterval
) {
  _stateData.setStartMaximumGlobalTimeInterval(
    std::min(startMaximumLocalTimeInterval, _stateData.getStartMaximumGlobalTimeInterval())
  );
  _stateData.setEndMaximumGlobalTimeInterval(
    std::max(endMaximumLocalTimeInterval, _stateData.getEndMaximumGlobalTimeInterval())
  );
  _stateData.setStartMinimumGlobalTimeInterval(
    std::max(startMinimumLocalTimeInterval, _stateData.getStartMinimumGlobalTimeInterval())
  );
  _stateData.setEndMinimumGlobalTimeInterval(
    std::min(endMinimumLocalTimeInterval, _stateData.getEndMinimumGlobalTimeInterval())
  );
}

void peanoclaw::State::resetGlobalTimeIntervals() {
  _stateData.setStartMaximumGlobalTimeInterval(std::numeric_limits<double>::max());
  _stateData.setEndMaximumGlobalTimeInterval(-std::numeric_limits<double>::max());
  _stateData.setStartMinimumGlobalTimeInterval(-std::numeric_limits<double>::max());
  _stateData.setEndMinimumGlobalTimeInterval(std::numeric_limits<double>::max());
}

double peanoclaw::State::getStartMaximumGlobalTimeInterval() const {
  return _stateData.getStartMaximumGlobalTimeInterval();
}

double peanoclaw::State::getEndMaximumGlobalTimeInterval() const {
  return _stateData.getEndMaximumGlobalTimeInterval();
}

double peanoclaw::State::getStartMinimumGlobalTimeInterval() const {
  return _stateData.getStartMinimumGlobalTimeInterval();
}

double peanoclaw::State::getEndMinimumGlobalTimeInterval() const {
  return _stateData.getEndMinimumGlobalTimeInterval();
}

void peanoclaw::State::resetTotalNumberOfCellUpdates() {
  _stateData.setTotalNumberOfCellUpdates(0.0);
}

double peanoclaw::State::getTotalNumberOfCellUpdates() const {
  return _stateData.getTotalNumberOfCellUpdates();
}

void peanoclaw::State::resetMinimalTimestep() {
  _stateData.setMinimalTimestep(std::numeric_limits<double>::max());
//  _stateData.setMinimalEstimatedNextTimestepSize(std::numeric_limits<double>::max());
}

void peanoclaw::State::updateMinimalTimestep(double timestep) {
  _stateData.setMinimalTimestep(std::min(timestep, _stateData.getMinimalTimestep()));
}

double peanoclaw::State::getMinimalTimestep() const {
  return _stateData.getMinimalTimestep();
}

void peanoclaw::State::prepareGridIteration () {
  delete _subgridStatistics;
  _subgridStatistics = new peanoclaw::statistics::SubgridStatistics;
  _subgridStatistics->setGlobalTimestepEndTime(_stateData.getGlobalTimestepEndTime());
}

void peanoclaw::State::finalizeGridIteration () {
  if(_subgridStatistics != 0)
  {
    _subgridStatistics->finalizeIteration(*this);
    _subgridStatistics->setGlobalTimestepEndTime(_stateData.getGlobalTimestepEndTime());
    if(tarch::parallel::Node::getInstance().isGlobalMaster()) {
      _subgridStatisticsHistory.push_back(*_subgridStatistics);
    }
  }
}

peanoclaw::statistics::SubgridStatistics* peanoclaw::State::getSubgridStatistics() const {
  return _subgridStatistics;
}

std::list<peanoclaw::statistics::SubgridStatistics> peanoclaw::State::getSubgridStatisticsHistory() const {
  return _subgridStatisticsHistory;
}

void peanoclaw::State::plotStatisticsForLastGridIteration() const {
  assertion(_subgridStatisticsHistory.size() > 0);
  _subgridStatisticsHistory.back().logLevelStatistics("Statistics for last grid iteration");
  _subgridStatisticsHistory.back().logProcessStatistics("Statistics for last grid iteration");
}

void peanoclaw::State::plotTotalStatistics() const {
  peanoclaw::statistics::SubgridStatistics totalStatistics;
  for(std::list<peanoclaw::statistics::SubgridStatistics>::const_iterator i = _subgridStatisticsHistory.begin();
      i != _subgridStatisticsHistory.end(); i++) {
    totalStatistics.merge(*i);
  }
  totalStatistics.averageTotalSimulationValues(_subgridStatisticsHistory.size());
  totalStatistics.logLevelStatistics("Total Statistics");
  totalStatistics.logProcessStatistics("Total Statistics");
}

void peanoclaw::State::setIsInitializing(bool isInitializing) {
  _stateData.setIsInitializing(isInitializing);
}

bool peanoclaw::State::getIsInitializing() const {
  return _stateData.getIsInitializing();
}

void peanoclaw::State::enableRefinementCriterion(bool enabled) {
  _stateData.setIsRefinementCriterionEnabled(enabled);
}

bool peanoclaw::State::isRefinementCriterionEnabled() const {
  return _stateData.getIsRefinementCriterionEnabled();
}


void peanoclaw::State::setUseDimensionalSplittingExtrapolation(bool useDimensionalSplittingExtrapolation) {
  _stateData.setUseDimensionalSplittingExtrapolation(useDimensionalSplittingExtrapolation);
}

bool peanoclaw::State::useDimensionalSplittingExtrapolation() const {
  return _stateData.getUseDimensionalSplittingExtrapolation();
}

//void peanoclaw::State::resetLocalHeightOfWorkerTree() {
//  #ifdef Parallel
//  _stateData.setGlobalHeightOfWorkerTreeDuringLastIteration(
//    _stateData.getLocalHeightOfWorkerTree()
//  );
//  _stateData.setLocalHeightOfWorkerTree(0);
//  #endif
//}
//
//void peanoclaw::State::increaseLocalHeightOfWorkerTree() {
//  #ifdef Parallel
//  _stateData.setLocalHeightOfWorkerTree(_stateData.getLocalHeightOfWorkerTree() + 1);
//  #endif
//}
//
//void peanoclaw::State::updateLocalHeightOfWorkerTree(int localHeightOfWorkerTree) {
//  #ifdef Parallel
//  _stateData.setLocalHeightOfWorkerTree(
//    std::max(_stateData.getLocalHeightOfWorkerTree(), localHeightOfWorkerTree)
//  );
//  #endif
//}
//
//int peanoclaw::State::getLocalHeightOfWorkerTree() const {
//  #ifdef Parallel
//  return _stateData.getLocalHeightOfWorkerTree();
//  #else
//  return 0;
//  #endif
//}
//
//int peanoclaw::State::getGlobalHeightOfWorkerTreeDuringLastIteration() const {
//  #ifdef Parallel
//  return _stateData.getGlobalHeightOfWorkerTreeDuringLastIteration();
//  #else
//  return 0;
//  #endif
//}


void peanoclaw::State::setReduceReductions(bool reduceReductions) {
  #ifdef Parallel
  _stateData.setReduceReductions(reduceReductions);
  #endif
}

bool peanoclaw::State::shouldReduceReductions() const {
  #ifdef Parallel
  return _stateData.getReduceReductions();
  #else
  return true;
  #endif
}

void peanoclaw::State::setRestrictStatistics(bool restrictStatistics) {
  #ifdef Parallel
  return _stateData.setRestrictStatistics(restrictStatistics);
  #endif
}

bool peanoclaw::State::shouldRestrictStatistics() const {
  #ifdef Parallel
  return _stateData.getRestrictStatistics();
  #else
  return false;
  #endif
}

//void peanoclaw::State::updateMinimalEstimatedNextTimestepSize(double estimatedNextTimestepSize) {
//  _stateData.setMinimalEstimatedNextTimestepSize(
//    std::min(_stateData.getMinimalEstimatedNextTimestepSize(), estimatedNextTimestepSize)
//  );
//}
//
//double peanoclaw::State::getMinimalEstimatedNextTimestepSize() const {
//  return _stateData.getMinimalEstimatedNextTimestepSize();
//}
//
//void peanoclaw::State::setEnforceGlobalTimestepping(bool enforceGlobalTimestepping) {
//  _stateData.setEnforceGlobalTimestepping(enforceGlobalTimestepping);
//}
//
//bool peanoclaw::State::shouldEnforceGlobalTimestepping() const {
//  return _stateData.getEnforceGlobalTimestepping();
//}

void peanoclaw::State::enableFluxCorrection(bool enableFluxCorrection) {
  _stateData.setEnableFluxCorrection(enableFluxCorrection);
}

bool peanoclaw::State::isFluxCorrectionEnabled() const {
  return _stateData.getEnableFluxCorrection();
}

#ifdef Parallel
void peanoclaw::State::setEstimatedNumberOfIterationsToGlobalTimestep(int estimatedNumber) {
  _stateData.setEstimatedNumberOfIterationsToGlobalTimestep(estimatedNumber);
}

int peanoclaw::State::getEstimatedNumberOfIterationsToGlobalTimestep() const {
  return _stateData.getEstimatedNumberOfIterationsToGlobalTimestep();
}
#endif

void peanoclaw::State::shouldEstimateNeighborInducedMaximumTimestep(bool shouldEstimate) {
  _stateData.setEstimateNeighborInducedMaximumTimestep(shouldEstimate);
}

bool peanoclaw::State::estimateNeighborInducedMaximumTimestep() const {
  return _stateData.getEstimateNeighborInducedMaximumTimestep();
}
