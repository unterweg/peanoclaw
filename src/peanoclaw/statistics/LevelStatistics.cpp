#include "peanoclaw/statistics/LevelStatistics.h"

peanoclaw::statistics::LevelStatistics::PersistentRecords::PersistentRecords() {
   
}


peanoclaw::statistics::LevelStatistics::PersistentRecords::PersistentRecords(const double& area, const int& level, const double& numberOfPatches, const double& numberOfCells, const double& numberOfCellUpdates, const double& createdPatches, const double& destroyedPatches, const double& patchesBlockedDueToNeighbors, const double& patchesBlockedDueToGlobalTimestep, const double& patchesSkippingIteration, const double& patchesCoarsening, const int& estimatedNumberOfRemainingIterationsToGlobalTimestep):
_area(area),
_level(level),
_numberOfPatches(numberOfPatches),
_numberOfCells(numberOfCells),
_numberOfCellUpdates(numberOfCellUpdates),
_createdPatches(createdPatches),
_destroyedPatches(destroyedPatches),
_patchesBlockedDueToNeighbors(patchesBlockedDueToNeighbors),
_patchesBlockedDueToGlobalTimestep(patchesBlockedDueToGlobalTimestep),
_patchesSkippingIteration(patchesSkippingIteration),
_patchesCoarsening(patchesCoarsening),
_estimatedNumberOfRemainingIterationsToGlobalTimestep(estimatedNumberOfRemainingIterationsToGlobalTimestep) {
   
}

peanoclaw::statistics::LevelStatistics::LevelStatistics() {
   
}


peanoclaw::statistics::LevelStatistics::LevelStatistics(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._area, persistentRecords._level, persistentRecords._numberOfPatches, persistentRecords._numberOfCells, persistentRecords._numberOfCellUpdates, persistentRecords._createdPatches, persistentRecords._destroyedPatches, persistentRecords._patchesBlockedDueToNeighbors, persistentRecords._patchesBlockedDueToGlobalTimestep, persistentRecords._patchesSkippingIteration, persistentRecords._patchesCoarsening, persistentRecords._estimatedNumberOfRemainingIterationsToGlobalTimestep) {
   
}


peanoclaw::statistics::LevelStatistics::LevelStatistics(const double& area, const int& level, const double& numberOfPatches, const double& numberOfCells, const double& numberOfCellUpdates, const double& createdPatches, const double& destroyedPatches, const double& patchesBlockedDueToNeighbors, const double& patchesBlockedDueToGlobalTimestep, const double& patchesSkippingIteration, const double& patchesCoarsening, const int& estimatedNumberOfRemainingIterationsToGlobalTimestep):
_persistentRecords(area, level, numberOfPatches, numberOfCells, numberOfCellUpdates, createdPatches, destroyedPatches, patchesBlockedDueToNeighbors, patchesBlockedDueToGlobalTimestep, patchesSkippingIteration, patchesCoarsening, estimatedNumberOfRemainingIterationsToGlobalTimestep) {
   
}


peanoclaw::statistics::LevelStatistics::~LevelStatistics() { }



std::string peanoclaw::statistics::LevelStatistics::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void peanoclaw::statistics::LevelStatistics::toString (std::ostream& out) const {
   out << "("; 
   out << "area:" << getArea();
   out << ",";
   out << "level:" << getLevel();
   out << ",";
   out << "numberOfPatches:" << getNumberOfPatches();
   out << ",";
   out << "numberOfCells:" << getNumberOfCells();
   out << ",";
   out << "numberOfCellUpdates:" << getNumberOfCellUpdates();
   out << ",";
   out << "createdPatches:" << getCreatedPatches();
   out << ",";
   out << "destroyedPatches:" << getDestroyedPatches();
   out << ",";
   out << "patchesBlockedDueToNeighbors:" << getPatchesBlockedDueToNeighbors();
   out << ",";
   out << "patchesBlockedDueToGlobalTimestep:" << getPatchesBlockedDueToGlobalTimestep();
   out << ",";
   out << "patchesSkippingIteration:" << getPatchesSkippingIteration();
   out << ",";
   out << "patchesCoarsening:" << getPatchesCoarsening();
   out << ",";
   out << "estimatedNumberOfRemainingIterationsToGlobalTimestep:" << getEstimatedNumberOfRemainingIterationsToGlobalTimestep();
   out <<  ")";
}


peanoclaw::statistics::LevelStatistics::PersistentRecords peanoclaw::statistics::LevelStatistics::getPersistentRecords() const {
   return _persistentRecords;
}

peanoclaw::statistics::LevelStatisticsPacked peanoclaw::statistics::LevelStatistics::convert() const{
   return LevelStatisticsPacked(
      getArea(),
      getLevel(),
      getNumberOfPatches(),
      getNumberOfCells(),
      getNumberOfCellUpdates(),
      getCreatedPatches(),
      getDestroyedPatches(),
      getPatchesBlockedDueToNeighbors(),
      getPatchesBlockedDueToGlobalTimestep(),
      getPatchesSkippingIteration(),
      getPatchesCoarsening(),
      getEstimatedNumberOfRemainingIterationsToGlobalTimestep()
   );
}

#ifdef Parallel
   tarch::logging::Log peanoclaw::statistics::LevelStatistics::_log( "peanoclaw::statistics::LevelStatistics" );
   
   MPI_Datatype peanoclaw::statistics::LevelStatistics::Datatype = 0;
   MPI_Datatype peanoclaw::statistics::LevelStatistics::FullDatatype = 0;
   
   
   void peanoclaw::statistics::LevelStatistics::initDatatype() {
      {
         LevelStatistics dummyLevelStatistics[2];
         
         const int Attributes = 13;
         MPI_Datatype subtypes[Attributes] = {
            MPI_DOUBLE,		 //area
            MPI_INT,		 //level
            MPI_DOUBLE,		 //numberOfPatches
            MPI_DOUBLE,		 //numberOfCells
            MPI_DOUBLE,		 //numberOfCellUpdates
            MPI_DOUBLE,		 //createdPatches
            MPI_DOUBLE,		 //destroyedPatches
            MPI_DOUBLE,		 //patchesBlockedDueToNeighbors
            MPI_DOUBLE,		 //patchesBlockedDueToGlobalTimestep
            MPI_DOUBLE,		 //patchesSkippingIteration
            MPI_DOUBLE,		 //patchesCoarsening
            MPI_INT,		 //estimatedNumberOfRemainingIterationsToGlobalTimestep
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //area
            1,		 //level
            1,		 //numberOfPatches
            1,		 //numberOfCells
            1,		 //numberOfCellUpdates
            1,		 //createdPatches
            1,		 //destroyedPatches
            1,		 //patchesBlockedDueToNeighbors
            1,		 //patchesBlockedDueToGlobalTimestep
            1,		 //patchesSkippingIteration
            1,		 //patchesCoarsening
            1,		 //estimatedNumberOfRemainingIterationsToGlobalTimestep
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._area))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._level))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._numberOfPatches))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._numberOfCells))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._numberOfCellUpdates))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._createdPatches))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._destroyedPatches))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._patchesBlockedDueToNeighbors))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._patchesBlockedDueToGlobalTimestep))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._patchesSkippingIteration))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._patchesCoarsening))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._estimatedNumberOfRemainingIterationsToGlobalTimestep))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[1]._persistentRecords._area))), 		&disp[12] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &LevelStatistics::Datatype );
         MPI_Type_commit( &LevelStatistics::Datatype );
         
      }
      {
         LevelStatistics dummyLevelStatistics[2];
         
         const int Attributes = 13;
         MPI_Datatype subtypes[Attributes] = {
            MPI_DOUBLE,		 //area
            MPI_INT,		 //level
            MPI_DOUBLE,		 //numberOfPatches
            MPI_DOUBLE,		 //numberOfCells
            MPI_DOUBLE,		 //numberOfCellUpdates
            MPI_DOUBLE,		 //createdPatches
            MPI_DOUBLE,		 //destroyedPatches
            MPI_DOUBLE,		 //patchesBlockedDueToNeighbors
            MPI_DOUBLE,		 //patchesBlockedDueToGlobalTimestep
            MPI_DOUBLE,		 //patchesSkippingIteration
            MPI_DOUBLE,		 //patchesCoarsening
            MPI_INT,		 //estimatedNumberOfRemainingIterationsToGlobalTimestep
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //area
            1,		 //level
            1,		 //numberOfPatches
            1,		 //numberOfCells
            1,		 //numberOfCellUpdates
            1,		 //createdPatches
            1,		 //destroyedPatches
            1,		 //patchesBlockedDueToNeighbors
            1,		 //patchesBlockedDueToGlobalTimestep
            1,		 //patchesSkippingIteration
            1,		 //patchesCoarsening
            1,		 //estimatedNumberOfRemainingIterationsToGlobalTimestep
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._area))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._level))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._numberOfPatches))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._numberOfCells))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._numberOfCellUpdates))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._createdPatches))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._destroyedPatches))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._patchesBlockedDueToNeighbors))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._patchesBlockedDueToGlobalTimestep))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._patchesSkippingIteration))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._patchesCoarsening))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[0]._persistentRecords._estimatedNumberOfRemainingIterationsToGlobalTimestep))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatistics[1]._persistentRecords._area))), 		&disp[12] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &LevelStatistics::FullDatatype );
         MPI_Type_commit( &LevelStatistics::FullDatatype );
         
      }
      
   }
   
   
   void peanoclaw::statistics::LevelStatistics::shutdownDatatype() {
      MPI_Type_free( &LevelStatistics::Datatype );
      MPI_Type_free( &LevelStatistics::FullDatatype );
      
   }
   
   void peanoclaw::statistics::LevelStatistics::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message peanoclaw::statistics::LevelStatistics "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         
      }
      else {
      
         MPI_Request* sendRequestHandle = new MPI_Request();
         MPI_Status   status;
         int          flag = 0;
         int          result;
         
         clock_t      timeOutWarning   = -1;
         clock_t      timeOutShutdown  = -1;
         bool         triggeredTimeoutWarning = false;
         
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            result = MPI_Isend(
               this, 1, Datatype, destination,
               tag, tarch::parallel::Node::getInstance().getCommunicator(),
               sendRequestHandle
            );
            
         }
         else {
            result = MPI_Isend(
               this, 1, FullDatatype, destination,
               tag, tarch::parallel::Node::getInstance().getCommunicator(),
               sendRequestHandle
            );
            
         }
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message peanoclaw::statistics::LevelStatistics "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         result = MPI_Test( sendRequestHandle, &flag, &status );
         while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            if (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "testing for finished send task for peanoclaw::statistics::LevelStatistics "
               << toString()
               << " sent to node " << destination
               << " failed: " << tarch::parallel::MPIReturnValueToString(result);
               _log.error("send(int)", msg.str() );
            }
            
            // deadlock aspect
            if (
               tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
               (clock()>timeOutWarning) &&
               (!triggeredTimeoutWarning)
            ) {
               tarch::parallel::Node::getInstance().writeTimeOutWarning(
               "peanoclaw::statistics::LevelStatistics",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "peanoclaw::statistics::LevelStatistics",
               "send(int)", destination,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         #ifdef Debug
         _log.debug("send(int,int)", "sent " + toString() );
         #endif
         
      }
      
   }
   
   
   
   void peanoclaw::statistics::LevelStatistics::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive peanoclaw::statistics::LevelStatistics from node "
            << source << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "receive(int)", msg.str() );
         }
         
      }
      else {
      
         MPI_Request* sendRequestHandle = new MPI_Request();
         MPI_Status   status;
         int          flag = 0;
         int          result;
         
         clock_t      timeOutWarning   = -1;
         clock_t      timeOutShutdown  = -1;
         bool         triggeredTimeoutWarning = false;
         
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            result = MPI_Irecv(
               this, 1, Datatype, source, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
            );
            
         }
         else {
            result = MPI_Irecv(
               this, 1, FullDatatype, source, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
            );
            
         }
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive peanoclaw::statistics::LevelStatistics from node "
            << source << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "receive(int)", msg.str() );
         }
         
         result = MPI_Test( sendRequestHandle, &flag, &status );
         while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            if (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "testing for finished receive task for peanoclaw::statistics::LevelStatistics failed: "
               << tarch::parallel::MPIReturnValueToString(result);
               _log.error("receive(int)", msg.str() );
            }
            
            // deadlock aspect
            if (
               tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
               (clock()>timeOutWarning) &&
               (!triggeredTimeoutWarning)
            ) {
               tarch::parallel::Node::getInstance().writeTimeOutWarning(
               "peanoclaw::statistics::LevelStatistics",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "peanoclaw::statistics::LevelStatistics",
               "receive(int)", source,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         
         #ifdef Debug
         _log.debug("receive(int,int)", "received " + toString() ); 
         #endif
         
      }
      
   }
   
   
   
   bool peanoclaw::statistics::LevelStatistics::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Status status;
      int  flag        = 0;
      MPI_Iprobe(
         MPI_ANY_SOURCE, tag,
         tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
      );
      if (flag) {
         int  messageCounter;
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Get_count(&status, Datatype, &messageCounter);
         }
         else {
            MPI_Get_count(&status, FullDatatype, &messageCounter);
         }
         return messageCounter > 0;
      }
      else return false;
      
   }
   
   
#endif


peanoclaw::statistics::LevelStatisticsPacked::PersistentRecords::PersistentRecords() {
   
}


peanoclaw::statistics::LevelStatisticsPacked::PersistentRecords::PersistentRecords(const double& area, const int& level, const double& numberOfPatches, const double& numberOfCells, const double& numberOfCellUpdates, const double& createdPatches, const double& destroyedPatches, const double& patchesBlockedDueToNeighbors, const double& patchesBlockedDueToGlobalTimestep, const double& patchesSkippingIteration, const double& patchesCoarsening, const int& estimatedNumberOfRemainingIterationsToGlobalTimestep):
_area(area),
_level(level),
_numberOfPatches(numberOfPatches),
_numberOfCells(numberOfCells),
_numberOfCellUpdates(numberOfCellUpdates),
_createdPatches(createdPatches),
_destroyedPatches(destroyedPatches),
_patchesBlockedDueToNeighbors(patchesBlockedDueToNeighbors),
_patchesBlockedDueToGlobalTimestep(patchesBlockedDueToGlobalTimestep),
_patchesSkippingIteration(patchesSkippingIteration),
_patchesCoarsening(patchesCoarsening),
_estimatedNumberOfRemainingIterationsToGlobalTimestep(estimatedNumberOfRemainingIterationsToGlobalTimestep) {
   
}

peanoclaw::statistics::LevelStatisticsPacked::LevelStatisticsPacked() {
   
}


peanoclaw::statistics::LevelStatisticsPacked::LevelStatisticsPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._area, persistentRecords._level, persistentRecords._numberOfPatches, persistentRecords._numberOfCells, persistentRecords._numberOfCellUpdates, persistentRecords._createdPatches, persistentRecords._destroyedPatches, persistentRecords._patchesBlockedDueToNeighbors, persistentRecords._patchesBlockedDueToGlobalTimestep, persistentRecords._patchesSkippingIteration, persistentRecords._patchesCoarsening, persistentRecords._estimatedNumberOfRemainingIterationsToGlobalTimestep) {
   
}


peanoclaw::statistics::LevelStatisticsPacked::LevelStatisticsPacked(const double& area, const int& level, const double& numberOfPatches, const double& numberOfCells, const double& numberOfCellUpdates, const double& createdPatches, const double& destroyedPatches, const double& patchesBlockedDueToNeighbors, const double& patchesBlockedDueToGlobalTimestep, const double& patchesSkippingIteration, const double& patchesCoarsening, const int& estimatedNumberOfRemainingIterationsToGlobalTimestep):
_persistentRecords(area, level, numberOfPatches, numberOfCells, numberOfCellUpdates, createdPatches, destroyedPatches, patchesBlockedDueToNeighbors, patchesBlockedDueToGlobalTimestep, patchesSkippingIteration, patchesCoarsening, estimatedNumberOfRemainingIterationsToGlobalTimestep) {
   
}


peanoclaw::statistics::LevelStatisticsPacked::~LevelStatisticsPacked() { }



std::string peanoclaw::statistics::LevelStatisticsPacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void peanoclaw::statistics::LevelStatisticsPacked::toString (std::ostream& out) const {
   out << "("; 
   out << "area:" << getArea();
   out << ",";
   out << "level:" << getLevel();
   out << ",";
   out << "numberOfPatches:" << getNumberOfPatches();
   out << ",";
   out << "numberOfCells:" << getNumberOfCells();
   out << ",";
   out << "numberOfCellUpdates:" << getNumberOfCellUpdates();
   out << ",";
   out << "createdPatches:" << getCreatedPatches();
   out << ",";
   out << "destroyedPatches:" << getDestroyedPatches();
   out << ",";
   out << "patchesBlockedDueToNeighbors:" << getPatchesBlockedDueToNeighbors();
   out << ",";
   out << "patchesBlockedDueToGlobalTimestep:" << getPatchesBlockedDueToGlobalTimestep();
   out << ",";
   out << "patchesSkippingIteration:" << getPatchesSkippingIteration();
   out << ",";
   out << "patchesCoarsening:" << getPatchesCoarsening();
   out << ",";
   out << "estimatedNumberOfRemainingIterationsToGlobalTimestep:" << getEstimatedNumberOfRemainingIterationsToGlobalTimestep();
   out <<  ")";
}


peanoclaw::statistics::LevelStatisticsPacked::PersistentRecords peanoclaw::statistics::LevelStatisticsPacked::getPersistentRecords() const {
   return _persistentRecords;
}

peanoclaw::statistics::LevelStatistics peanoclaw::statistics::LevelStatisticsPacked::convert() const{
   return LevelStatistics(
      getArea(),
      getLevel(),
      getNumberOfPatches(),
      getNumberOfCells(),
      getNumberOfCellUpdates(),
      getCreatedPatches(),
      getDestroyedPatches(),
      getPatchesBlockedDueToNeighbors(),
      getPatchesBlockedDueToGlobalTimestep(),
      getPatchesSkippingIteration(),
      getPatchesCoarsening(),
      getEstimatedNumberOfRemainingIterationsToGlobalTimestep()
   );
}

#ifdef Parallel
   tarch::logging::Log peanoclaw::statistics::LevelStatisticsPacked::_log( "peanoclaw::statistics::LevelStatisticsPacked" );
   
   MPI_Datatype peanoclaw::statistics::LevelStatisticsPacked::Datatype = 0;
   MPI_Datatype peanoclaw::statistics::LevelStatisticsPacked::FullDatatype = 0;
   
   
   void peanoclaw::statistics::LevelStatisticsPacked::initDatatype() {
      {
         LevelStatisticsPacked dummyLevelStatisticsPacked[2];
         
         const int Attributes = 13;
         MPI_Datatype subtypes[Attributes] = {
            MPI_DOUBLE,		 //area
            MPI_INT,		 //level
            MPI_DOUBLE,		 //numberOfPatches
            MPI_DOUBLE,		 //numberOfCells
            MPI_DOUBLE,		 //numberOfCellUpdates
            MPI_DOUBLE,		 //createdPatches
            MPI_DOUBLE,		 //destroyedPatches
            MPI_DOUBLE,		 //patchesBlockedDueToNeighbors
            MPI_DOUBLE,		 //patchesBlockedDueToGlobalTimestep
            MPI_DOUBLE,		 //patchesSkippingIteration
            MPI_DOUBLE,		 //patchesCoarsening
            MPI_INT,		 //estimatedNumberOfRemainingIterationsToGlobalTimestep
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //area
            1,		 //level
            1,		 //numberOfPatches
            1,		 //numberOfCells
            1,		 //numberOfCellUpdates
            1,		 //createdPatches
            1,		 //destroyedPatches
            1,		 //patchesBlockedDueToNeighbors
            1,		 //patchesBlockedDueToGlobalTimestep
            1,		 //patchesSkippingIteration
            1,		 //patchesCoarsening
            1,		 //estimatedNumberOfRemainingIterationsToGlobalTimestep
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._area))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._level))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._numberOfPatches))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._numberOfCells))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._numberOfCellUpdates))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._createdPatches))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._destroyedPatches))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._patchesBlockedDueToNeighbors))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._patchesBlockedDueToGlobalTimestep))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._patchesSkippingIteration))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._patchesCoarsening))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._estimatedNumberOfRemainingIterationsToGlobalTimestep))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[1]._persistentRecords._area))), 		&disp[12] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &LevelStatisticsPacked::Datatype );
         MPI_Type_commit( &LevelStatisticsPacked::Datatype );
         
      }
      {
         LevelStatisticsPacked dummyLevelStatisticsPacked[2];
         
         const int Attributes = 13;
         MPI_Datatype subtypes[Attributes] = {
            MPI_DOUBLE,		 //area
            MPI_INT,		 //level
            MPI_DOUBLE,		 //numberOfPatches
            MPI_DOUBLE,		 //numberOfCells
            MPI_DOUBLE,		 //numberOfCellUpdates
            MPI_DOUBLE,		 //createdPatches
            MPI_DOUBLE,		 //destroyedPatches
            MPI_DOUBLE,		 //patchesBlockedDueToNeighbors
            MPI_DOUBLE,		 //patchesBlockedDueToGlobalTimestep
            MPI_DOUBLE,		 //patchesSkippingIteration
            MPI_DOUBLE,		 //patchesCoarsening
            MPI_INT,		 //estimatedNumberOfRemainingIterationsToGlobalTimestep
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //area
            1,		 //level
            1,		 //numberOfPatches
            1,		 //numberOfCells
            1,		 //numberOfCellUpdates
            1,		 //createdPatches
            1,		 //destroyedPatches
            1,		 //patchesBlockedDueToNeighbors
            1,		 //patchesBlockedDueToGlobalTimestep
            1,		 //patchesSkippingIteration
            1,		 //patchesCoarsening
            1,		 //estimatedNumberOfRemainingIterationsToGlobalTimestep
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._area))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._level))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._numberOfPatches))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._numberOfCells))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._numberOfCellUpdates))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._createdPatches))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._destroyedPatches))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._patchesBlockedDueToNeighbors))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._patchesBlockedDueToGlobalTimestep))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._patchesSkippingIteration))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._patchesCoarsening))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[0]._persistentRecords._estimatedNumberOfRemainingIterationsToGlobalTimestep))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyLevelStatisticsPacked[1]._persistentRecords._area))), 		&disp[12] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &LevelStatisticsPacked::FullDatatype );
         MPI_Type_commit( &LevelStatisticsPacked::FullDatatype );
         
      }
      
   }
   
   
   void peanoclaw::statistics::LevelStatisticsPacked::shutdownDatatype() {
      MPI_Type_free( &LevelStatisticsPacked::Datatype );
      MPI_Type_free( &LevelStatisticsPacked::FullDatatype );
      
   }
   
   void peanoclaw::statistics::LevelStatisticsPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message peanoclaw::statistics::LevelStatisticsPacked "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         
      }
      else {
      
         MPI_Request* sendRequestHandle = new MPI_Request();
         MPI_Status   status;
         int          flag = 0;
         int          result;
         
         clock_t      timeOutWarning   = -1;
         clock_t      timeOutShutdown  = -1;
         bool         triggeredTimeoutWarning = false;
         
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            result = MPI_Isend(
               this, 1, Datatype, destination,
               tag, tarch::parallel::Node::getInstance().getCommunicator(),
               sendRequestHandle
            );
            
         }
         else {
            result = MPI_Isend(
               this, 1, FullDatatype, destination,
               tag, tarch::parallel::Node::getInstance().getCommunicator(),
               sendRequestHandle
            );
            
         }
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message peanoclaw::statistics::LevelStatisticsPacked "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         result = MPI_Test( sendRequestHandle, &flag, &status );
         while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            if (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "testing for finished send task for peanoclaw::statistics::LevelStatisticsPacked "
               << toString()
               << " sent to node " << destination
               << " failed: " << tarch::parallel::MPIReturnValueToString(result);
               _log.error("send(int)", msg.str() );
            }
            
            // deadlock aspect
            if (
               tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
               (clock()>timeOutWarning) &&
               (!triggeredTimeoutWarning)
            ) {
               tarch::parallel::Node::getInstance().writeTimeOutWarning(
               "peanoclaw::statistics::LevelStatisticsPacked",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "peanoclaw::statistics::LevelStatisticsPacked",
               "send(int)", destination,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         #ifdef Debug
         _log.debug("send(int,int)", "sent " + toString() );
         #endif
         
      }
      
   }
   
   
   
   void peanoclaw::statistics::LevelStatisticsPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive peanoclaw::statistics::LevelStatisticsPacked from node "
            << source << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "receive(int)", msg.str() );
         }
         
      }
      else {
      
         MPI_Request* sendRequestHandle = new MPI_Request();
         MPI_Status   status;
         int          flag = 0;
         int          result;
         
         clock_t      timeOutWarning   = -1;
         clock_t      timeOutShutdown  = -1;
         bool         triggeredTimeoutWarning = false;
         
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            result = MPI_Irecv(
               this, 1, Datatype, source, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
            );
            
         }
         else {
            result = MPI_Irecv(
               this, 1, FullDatatype, source, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
            );
            
         }
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive peanoclaw::statistics::LevelStatisticsPacked from node "
            << source << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "receive(int)", msg.str() );
         }
         
         result = MPI_Test( sendRequestHandle, &flag, &status );
         while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            if (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "testing for finished receive task for peanoclaw::statistics::LevelStatisticsPacked failed: "
               << tarch::parallel::MPIReturnValueToString(result);
               _log.error("receive(int)", msg.str() );
            }
            
            // deadlock aspect
            if (
               tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
               (clock()>timeOutWarning) &&
               (!triggeredTimeoutWarning)
            ) {
               tarch::parallel::Node::getInstance().writeTimeOutWarning(
               "peanoclaw::statistics::LevelStatisticsPacked",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "peanoclaw::statistics::LevelStatisticsPacked",
               "receive(int)", source,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         
         #ifdef Debug
         _log.debug("receive(int,int)", "received " + toString() ); 
         #endif
         
      }
      
   }
   
   
   
   bool peanoclaw::statistics::LevelStatisticsPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Status status;
      int  flag        = 0;
      MPI_Iprobe(
         MPI_ANY_SOURCE, tag,
         tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
      );
      if (flag) {
         int  messageCounter;
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Get_count(&status, Datatype, &messageCounter);
         }
         else {
            MPI_Get_count(&status, FullDatatype, &messageCounter);
         }
         return messageCounter > 0;
      }
      else return false;
      
   }
   
   
#endif



