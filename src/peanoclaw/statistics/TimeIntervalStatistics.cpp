#include "peanoclaw/statistics/TimeIntervalStatistics.h"

peanoclaw::statistics::TimeIntervalStatistics::PersistentRecords::PersistentRecords() {
   
}


peanoclaw::statistics::TimeIntervalStatistics::PersistentRecords::PersistentRecords(const int& minimalPatchIndex, const int& minimalPatchParentIndex, const double& minimalPatchTime, const double& startMaximumLocalTimeInterval, const double& endMaximumLocalTimeInterval, const double& startMinimumLocalTimeInterval, const double& endMinimumLocalTimeInterval, const double& minimalTimestep, const bool& allPatchesEvolvedToGlobalTimestep, const double& averageGlobalTimeInterval, const double& globalTimestepEndTime, const bool& minimalPatchBlockedDueToCoarsening, const bool& minimalPatchBlockedDueToGlobalTimestep):
_minimalPatchIndex(minimalPatchIndex),
_minimalPatchParentIndex(minimalPatchParentIndex),
_minimalPatchTime(minimalPatchTime),
_startMaximumLocalTimeInterval(startMaximumLocalTimeInterval),
_endMaximumLocalTimeInterval(endMaximumLocalTimeInterval),
_startMinimumLocalTimeInterval(startMinimumLocalTimeInterval),
_endMinimumLocalTimeInterval(endMinimumLocalTimeInterval),
_minimalTimestep(minimalTimestep),
_allPatchesEvolvedToGlobalTimestep(allPatchesEvolvedToGlobalTimestep),
_averageGlobalTimeInterval(averageGlobalTimeInterval),
_globalTimestepEndTime(globalTimestepEndTime),
_minimalPatchBlockedDueToCoarsening(minimalPatchBlockedDueToCoarsening),
_minimalPatchBlockedDueToGlobalTimestep(minimalPatchBlockedDueToGlobalTimestep) {
   
}

peanoclaw::statistics::TimeIntervalStatistics::TimeIntervalStatistics() {
   
}


peanoclaw::statistics::TimeIntervalStatistics::TimeIntervalStatistics(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._minimalPatchIndex, persistentRecords._minimalPatchParentIndex, persistentRecords._minimalPatchTime, persistentRecords._startMaximumLocalTimeInterval, persistentRecords._endMaximumLocalTimeInterval, persistentRecords._startMinimumLocalTimeInterval, persistentRecords._endMinimumLocalTimeInterval, persistentRecords._minimalTimestep, persistentRecords._allPatchesEvolvedToGlobalTimestep, persistentRecords._averageGlobalTimeInterval, persistentRecords._globalTimestepEndTime, persistentRecords._minimalPatchBlockedDueToCoarsening, persistentRecords._minimalPatchBlockedDueToGlobalTimestep) {
   
}


peanoclaw::statistics::TimeIntervalStatistics::TimeIntervalStatistics(const int& minimalPatchIndex, const int& minimalPatchParentIndex, const double& minimalPatchTime, const double& startMaximumLocalTimeInterval, const double& endMaximumLocalTimeInterval, const double& startMinimumLocalTimeInterval, const double& endMinimumLocalTimeInterval, const double& minimalTimestep, const bool& allPatchesEvolvedToGlobalTimestep, const double& averageGlobalTimeInterval, const double& globalTimestepEndTime, const bool& minimalPatchBlockedDueToCoarsening, const bool& minimalPatchBlockedDueToGlobalTimestep):
_persistentRecords(minimalPatchIndex, minimalPatchParentIndex, minimalPatchTime, startMaximumLocalTimeInterval, endMaximumLocalTimeInterval, startMinimumLocalTimeInterval, endMinimumLocalTimeInterval, minimalTimestep, allPatchesEvolvedToGlobalTimestep, averageGlobalTimeInterval, globalTimestepEndTime, minimalPatchBlockedDueToCoarsening, minimalPatchBlockedDueToGlobalTimestep) {
   
}


peanoclaw::statistics::TimeIntervalStatistics::~TimeIntervalStatistics() { }



std::string peanoclaw::statistics::TimeIntervalStatistics::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void peanoclaw::statistics::TimeIntervalStatistics::toString (std::ostream& out) const {
   out << "("; 
   out << "minimalPatchIndex:" << getMinimalPatchIndex();
   out << ",";
   out << "minimalPatchParentIndex:" << getMinimalPatchParentIndex();
   out << ",";
   out << "minimalPatchTime:" << getMinimalPatchTime();
   out << ",";
   out << "startMaximumLocalTimeInterval:" << getStartMaximumLocalTimeInterval();
   out << ",";
   out << "endMaximumLocalTimeInterval:" << getEndMaximumLocalTimeInterval();
   out << ",";
   out << "startMinimumLocalTimeInterval:" << getStartMinimumLocalTimeInterval();
   out << ",";
   out << "endMinimumLocalTimeInterval:" << getEndMinimumLocalTimeInterval();
   out << ",";
   out << "minimalTimestep:" << getMinimalTimestep();
   out << ",";
   out << "allPatchesEvolvedToGlobalTimestep:" << getAllPatchesEvolvedToGlobalTimestep();
   out << ",";
   out << "averageGlobalTimeInterval:" << getAverageGlobalTimeInterval();
   out << ",";
   out << "globalTimestepEndTime:" << getGlobalTimestepEndTime();
   out << ",";
   out << "minimalPatchBlockedDueToCoarsening:" << getMinimalPatchBlockedDueToCoarsening();
   out << ",";
   out << "minimalPatchBlockedDueToGlobalTimestep:" << getMinimalPatchBlockedDueToGlobalTimestep();
   out <<  ")";
}


peanoclaw::statistics::TimeIntervalStatistics::PersistentRecords peanoclaw::statistics::TimeIntervalStatistics::getPersistentRecords() const {
   return _persistentRecords;
}

peanoclaw::statistics::TimeIntervalStatisticsPacked peanoclaw::statistics::TimeIntervalStatistics::convert() const{
   return TimeIntervalStatisticsPacked(
      getMinimalPatchIndex(),
      getMinimalPatchParentIndex(),
      getMinimalPatchTime(),
      getStartMaximumLocalTimeInterval(),
      getEndMaximumLocalTimeInterval(),
      getStartMinimumLocalTimeInterval(),
      getEndMinimumLocalTimeInterval(),
      getMinimalTimestep(),
      getAllPatchesEvolvedToGlobalTimestep(),
      getAverageGlobalTimeInterval(),
      getGlobalTimestepEndTime(),
      getMinimalPatchBlockedDueToCoarsening(),
      getMinimalPatchBlockedDueToGlobalTimestep()
   );
}

#ifdef Parallel
   tarch::logging::Log peanoclaw::statistics::TimeIntervalStatistics::_log( "peanoclaw::statistics::TimeIntervalStatistics" );
   
   MPI_Datatype peanoclaw::statistics::TimeIntervalStatistics::Datatype = 0;
   MPI_Datatype peanoclaw::statistics::TimeIntervalStatistics::FullDatatype = 0;
   
   
   void peanoclaw::statistics::TimeIntervalStatistics::initDatatype() {
      {
         TimeIntervalStatistics dummyTimeIntervalStatistics[2];
         
         const int Attributes = 14;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //minimalPatchIndex
            MPI_INT,		 //minimalPatchParentIndex
            MPI_DOUBLE,		 //minimalPatchTime
            MPI_DOUBLE,		 //startMaximumLocalTimeInterval
            MPI_DOUBLE,		 //endMaximumLocalTimeInterval
            MPI_DOUBLE,		 //startMinimumLocalTimeInterval
            MPI_DOUBLE,		 //endMinimumLocalTimeInterval
            MPI_DOUBLE,		 //minimalTimestep
            MPI_CHAR,		 //allPatchesEvolvedToGlobalTimestep
            MPI_DOUBLE,		 //averageGlobalTimeInterval
            MPI_DOUBLE,		 //globalTimestepEndTime
            MPI_CHAR,		 //minimalPatchBlockedDueToCoarsening
            MPI_CHAR,		 //minimalPatchBlockedDueToGlobalTimestep
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //minimalPatchIndex
            1,		 //minimalPatchParentIndex
            1,		 //minimalPatchTime
            1,		 //startMaximumLocalTimeInterval
            1,		 //endMaximumLocalTimeInterval
            1,		 //startMinimumLocalTimeInterval
            1,		 //endMinimumLocalTimeInterval
            1,		 //minimalTimestep
            1,		 //allPatchesEvolvedToGlobalTimestep
            1,		 //averageGlobalTimeInterval
            1,		 //globalTimestepEndTime
            1,		 //minimalPatchBlockedDueToCoarsening
            1,		 //minimalPatchBlockedDueToGlobalTimestep
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchIndex))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchParentIndex))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchTime))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._startMaximumLocalTimeInterval))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._endMaximumLocalTimeInterval))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._startMinimumLocalTimeInterval))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._endMinimumLocalTimeInterval))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalTimestep))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._allPatchesEvolvedToGlobalTimestep))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._averageGlobalTimeInterval))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._globalTimestepEndTime))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchBlockedDueToCoarsening))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchBlockedDueToGlobalTimestep))), 		&disp[12] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[1]._persistentRecords._minimalPatchIndex))), 		&disp[13] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &TimeIntervalStatistics::Datatype );
         MPI_Type_commit( &TimeIntervalStatistics::Datatype );
         
      }
      {
         TimeIntervalStatistics dummyTimeIntervalStatistics[2];
         
         const int Attributes = 14;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //minimalPatchIndex
            MPI_INT,		 //minimalPatchParentIndex
            MPI_DOUBLE,		 //minimalPatchTime
            MPI_DOUBLE,		 //startMaximumLocalTimeInterval
            MPI_DOUBLE,		 //endMaximumLocalTimeInterval
            MPI_DOUBLE,		 //startMinimumLocalTimeInterval
            MPI_DOUBLE,		 //endMinimumLocalTimeInterval
            MPI_DOUBLE,		 //minimalTimestep
            MPI_CHAR,		 //allPatchesEvolvedToGlobalTimestep
            MPI_DOUBLE,		 //averageGlobalTimeInterval
            MPI_DOUBLE,		 //globalTimestepEndTime
            MPI_CHAR,		 //minimalPatchBlockedDueToCoarsening
            MPI_CHAR,		 //minimalPatchBlockedDueToGlobalTimestep
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //minimalPatchIndex
            1,		 //minimalPatchParentIndex
            1,		 //minimalPatchTime
            1,		 //startMaximumLocalTimeInterval
            1,		 //endMaximumLocalTimeInterval
            1,		 //startMinimumLocalTimeInterval
            1,		 //endMinimumLocalTimeInterval
            1,		 //minimalTimestep
            1,		 //allPatchesEvolvedToGlobalTimestep
            1,		 //averageGlobalTimeInterval
            1,		 //globalTimestepEndTime
            1,		 //minimalPatchBlockedDueToCoarsening
            1,		 //minimalPatchBlockedDueToGlobalTimestep
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchIndex))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchParentIndex))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchTime))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._startMaximumLocalTimeInterval))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._endMaximumLocalTimeInterval))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._startMinimumLocalTimeInterval))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._endMinimumLocalTimeInterval))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalTimestep))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._allPatchesEvolvedToGlobalTimestep))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._averageGlobalTimeInterval))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._globalTimestepEndTime))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchBlockedDueToCoarsening))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[0]._persistentRecords._minimalPatchBlockedDueToGlobalTimestep))), 		&disp[12] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatistics[1]._persistentRecords._minimalPatchIndex))), 		&disp[13] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &TimeIntervalStatistics::FullDatatype );
         MPI_Type_commit( &TimeIntervalStatistics::FullDatatype );
         
      }
      
   }
   
   
   void peanoclaw::statistics::TimeIntervalStatistics::shutdownDatatype() {
      MPI_Type_free( &TimeIntervalStatistics::Datatype );
      MPI_Type_free( &TimeIntervalStatistics::FullDatatype );
      
   }
   
   void peanoclaw::statistics::TimeIntervalStatistics::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message peanoclaw::statistics::TimeIntervalStatistics "
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
            msg << "was not able to send message peanoclaw::statistics::TimeIntervalStatistics "
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
               msg << "testing for finished send task for peanoclaw::statistics::TimeIntervalStatistics "
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
               "peanoclaw::statistics::TimeIntervalStatistics",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "peanoclaw::statistics::TimeIntervalStatistics",
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
   
   
   
   void peanoclaw::statistics::TimeIntervalStatistics::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive peanoclaw::statistics::TimeIntervalStatistics from node "
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
            msg << "failed to start to receive peanoclaw::statistics::TimeIntervalStatistics from node "
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
               msg << "testing for finished receive task for peanoclaw::statistics::TimeIntervalStatistics failed: "
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
               "peanoclaw::statistics::TimeIntervalStatistics",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "peanoclaw::statistics::TimeIntervalStatistics",
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
   
   
   
   bool peanoclaw::statistics::TimeIntervalStatistics::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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


peanoclaw::statistics::TimeIntervalStatisticsPacked::PersistentRecords::PersistentRecords() {
   assertion((3 < (8 * sizeof(short int))));
   
}


peanoclaw::statistics::TimeIntervalStatisticsPacked::PersistentRecords::PersistentRecords(const int& minimalPatchIndex, const int& minimalPatchParentIndex, const double& minimalPatchTime, const double& startMaximumLocalTimeInterval, const double& endMaximumLocalTimeInterval, const double& startMinimumLocalTimeInterval, const double& endMinimumLocalTimeInterval, const double& minimalTimestep, const bool& allPatchesEvolvedToGlobalTimestep, const double& averageGlobalTimeInterval, const double& globalTimestepEndTime, const bool& minimalPatchBlockedDueToCoarsening, const bool& minimalPatchBlockedDueToGlobalTimestep):
_minimalPatchIndex(minimalPatchIndex),
_minimalPatchParentIndex(minimalPatchParentIndex),
_minimalPatchTime(minimalPatchTime),
_startMaximumLocalTimeInterval(startMaximumLocalTimeInterval),
_endMaximumLocalTimeInterval(endMaximumLocalTimeInterval),
_startMinimumLocalTimeInterval(startMinimumLocalTimeInterval),
_endMinimumLocalTimeInterval(endMinimumLocalTimeInterval),
_minimalTimestep(minimalTimestep),
_averageGlobalTimeInterval(averageGlobalTimeInterval),
_globalTimestepEndTime(globalTimestepEndTime) {
   setAllPatchesEvolvedToGlobalTimestep(allPatchesEvolvedToGlobalTimestep);
   setMinimalPatchBlockedDueToCoarsening(minimalPatchBlockedDueToCoarsening);
   setMinimalPatchBlockedDueToGlobalTimestep(minimalPatchBlockedDueToGlobalTimestep);
   assertion((3 < (8 * sizeof(short int))));
   
}

peanoclaw::statistics::TimeIntervalStatisticsPacked::TimeIntervalStatisticsPacked() {
   assertion((3 < (8 * sizeof(short int))));
   
}


peanoclaw::statistics::TimeIntervalStatisticsPacked::TimeIntervalStatisticsPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._minimalPatchIndex, persistentRecords._minimalPatchParentIndex, persistentRecords._minimalPatchTime, persistentRecords._startMaximumLocalTimeInterval, persistentRecords._endMaximumLocalTimeInterval, persistentRecords._startMinimumLocalTimeInterval, persistentRecords._endMinimumLocalTimeInterval, persistentRecords._minimalTimestep, persistentRecords.getAllPatchesEvolvedToGlobalTimestep(), persistentRecords._averageGlobalTimeInterval, persistentRecords._globalTimestepEndTime, persistentRecords.getMinimalPatchBlockedDueToCoarsening(), persistentRecords.getMinimalPatchBlockedDueToGlobalTimestep()) {
   assertion((3 < (8 * sizeof(short int))));
   
}


peanoclaw::statistics::TimeIntervalStatisticsPacked::TimeIntervalStatisticsPacked(const int& minimalPatchIndex, const int& minimalPatchParentIndex, const double& minimalPatchTime, const double& startMaximumLocalTimeInterval, const double& endMaximumLocalTimeInterval, const double& startMinimumLocalTimeInterval, const double& endMinimumLocalTimeInterval, const double& minimalTimestep, const bool& allPatchesEvolvedToGlobalTimestep, const double& averageGlobalTimeInterval, const double& globalTimestepEndTime, const bool& minimalPatchBlockedDueToCoarsening, const bool& minimalPatchBlockedDueToGlobalTimestep):
_persistentRecords(minimalPatchIndex, minimalPatchParentIndex, minimalPatchTime, startMaximumLocalTimeInterval, endMaximumLocalTimeInterval, startMinimumLocalTimeInterval, endMinimumLocalTimeInterval, minimalTimestep, allPatchesEvolvedToGlobalTimestep, averageGlobalTimeInterval, globalTimestepEndTime, minimalPatchBlockedDueToCoarsening, minimalPatchBlockedDueToGlobalTimestep) {
   assertion((3 < (8 * sizeof(short int))));
   
}


peanoclaw::statistics::TimeIntervalStatisticsPacked::~TimeIntervalStatisticsPacked() { }



std::string peanoclaw::statistics::TimeIntervalStatisticsPacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void peanoclaw::statistics::TimeIntervalStatisticsPacked::toString (std::ostream& out) const {
   out << "("; 
   out << "minimalPatchIndex:" << getMinimalPatchIndex();
   out << ",";
   out << "minimalPatchParentIndex:" << getMinimalPatchParentIndex();
   out << ",";
   out << "minimalPatchTime:" << getMinimalPatchTime();
   out << ",";
   out << "startMaximumLocalTimeInterval:" << getStartMaximumLocalTimeInterval();
   out << ",";
   out << "endMaximumLocalTimeInterval:" << getEndMaximumLocalTimeInterval();
   out << ",";
   out << "startMinimumLocalTimeInterval:" << getStartMinimumLocalTimeInterval();
   out << ",";
   out << "endMinimumLocalTimeInterval:" << getEndMinimumLocalTimeInterval();
   out << ",";
   out << "minimalTimestep:" << getMinimalTimestep();
   out << ",";
   out << "allPatchesEvolvedToGlobalTimestep:" << getAllPatchesEvolvedToGlobalTimestep();
   out << ",";
   out << "averageGlobalTimeInterval:" << getAverageGlobalTimeInterval();
   out << ",";
   out << "globalTimestepEndTime:" << getGlobalTimestepEndTime();
   out << ",";
   out << "minimalPatchBlockedDueToCoarsening:" << getMinimalPatchBlockedDueToCoarsening();
   out << ",";
   out << "minimalPatchBlockedDueToGlobalTimestep:" << getMinimalPatchBlockedDueToGlobalTimestep();
   out <<  ")";
}


peanoclaw::statistics::TimeIntervalStatisticsPacked::PersistentRecords peanoclaw::statistics::TimeIntervalStatisticsPacked::getPersistentRecords() const {
   return _persistentRecords;
}

peanoclaw::statistics::TimeIntervalStatistics peanoclaw::statistics::TimeIntervalStatisticsPacked::convert() const{
   return TimeIntervalStatistics(
      getMinimalPatchIndex(),
      getMinimalPatchParentIndex(),
      getMinimalPatchTime(),
      getStartMaximumLocalTimeInterval(),
      getEndMaximumLocalTimeInterval(),
      getStartMinimumLocalTimeInterval(),
      getEndMinimumLocalTimeInterval(),
      getMinimalTimestep(),
      getAllPatchesEvolvedToGlobalTimestep(),
      getAverageGlobalTimeInterval(),
      getGlobalTimestepEndTime(),
      getMinimalPatchBlockedDueToCoarsening(),
      getMinimalPatchBlockedDueToGlobalTimestep()
   );
}

#ifdef Parallel
   tarch::logging::Log peanoclaw::statistics::TimeIntervalStatisticsPacked::_log( "peanoclaw::statistics::TimeIntervalStatisticsPacked" );
   
   MPI_Datatype peanoclaw::statistics::TimeIntervalStatisticsPacked::Datatype = 0;
   MPI_Datatype peanoclaw::statistics::TimeIntervalStatisticsPacked::FullDatatype = 0;
   
   
   void peanoclaw::statistics::TimeIntervalStatisticsPacked::initDatatype() {
      {
         TimeIntervalStatisticsPacked dummyTimeIntervalStatisticsPacked[2];
         
         const int Attributes = 12;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //minimalPatchIndex
            MPI_INT,		 //minimalPatchParentIndex
            MPI_DOUBLE,		 //minimalPatchTime
            MPI_DOUBLE,		 //startMaximumLocalTimeInterval
            MPI_DOUBLE,		 //endMaximumLocalTimeInterval
            MPI_DOUBLE,		 //startMinimumLocalTimeInterval
            MPI_DOUBLE,		 //endMinimumLocalTimeInterval
            MPI_DOUBLE,		 //minimalTimestep
            MPI_DOUBLE,		 //averageGlobalTimeInterval
            MPI_DOUBLE,		 //globalTimestepEndTime
            MPI_SHORT,		 //_packedRecords0
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //minimalPatchIndex
            1,		 //minimalPatchParentIndex
            1,		 //minimalPatchTime
            1,		 //startMaximumLocalTimeInterval
            1,		 //endMaximumLocalTimeInterval
            1,		 //startMinimumLocalTimeInterval
            1,		 //endMinimumLocalTimeInterval
            1,		 //minimalTimestep
            1,		 //averageGlobalTimeInterval
            1,		 //globalTimestepEndTime
            1,		 //_packedRecords0
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._minimalPatchIndex))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._minimalPatchParentIndex))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._minimalPatchTime))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._startMaximumLocalTimeInterval))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._endMaximumLocalTimeInterval))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._startMinimumLocalTimeInterval))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._endMinimumLocalTimeInterval))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._minimalTimestep))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._averageGlobalTimeInterval))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._globalTimestepEndTime))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._packedRecords0))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[1]._persistentRecords._minimalPatchIndex))), 		&disp[11] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &TimeIntervalStatisticsPacked::Datatype );
         MPI_Type_commit( &TimeIntervalStatisticsPacked::Datatype );
         
      }
      {
         TimeIntervalStatisticsPacked dummyTimeIntervalStatisticsPacked[2];
         
         const int Attributes = 12;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //minimalPatchIndex
            MPI_INT,		 //minimalPatchParentIndex
            MPI_DOUBLE,		 //minimalPatchTime
            MPI_DOUBLE,		 //startMaximumLocalTimeInterval
            MPI_DOUBLE,		 //endMaximumLocalTimeInterval
            MPI_DOUBLE,		 //startMinimumLocalTimeInterval
            MPI_DOUBLE,		 //endMinimumLocalTimeInterval
            MPI_DOUBLE,		 //minimalTimestep
            MPI_DOUBLE,		 //averageGlobalTimeInterval
            MPI_DOUBLE,		 //globalTimestepEndTime
            MPI_SHORT,		 //_packedRecords0
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //minimalPatchIndex
            1,		 //minimalPatchParentIndex
            1,		 //minimalPatchTime
            1,		 //startMaximumLocalTimeInterval
            1,		 //endMaximumLocalTimeInterval
            1,		 //startMinimumLocalTimeInterval
            1,		 //endMinimumLocalTimeInterval
            1,		 //minimalTimestep
            1,		 //averageGlobalTimeInterval
            1,		 //globalTimestepEndTime
            1,		 //_packedRecords0
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._minimalPatchIndex))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._minimalPatchParentIndex))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._minimalPatchTime))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._startMaximumLocalTimeInterval))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._endMaximumLocalTimeInterval))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._startMinimumLocalTimeInterval))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._endMinimumLocalTimeInterval))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._minimalTimestep))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._averageGlobalTimeInterval))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._globalTimestepEndTime))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[0]._persistentRecords._packedRecords0))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyTimeIntervalStatisticsPacked[1]._persistentRecords._minimalPatchIndex))), 		&disp[11] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &TimeIntervalStatisticsPacked::FullDatatype );
         MPI_Type_commit( &TimeIntervalStatisticsPacked::FullDatatype );
         
      }
      
   }
   
   
   void peanoclaw::statistics::TimeIntervalStatisticsPacked::shutdownDatatype() {
      MPI_Type_free( &TimeIntervalStatisticsPacked::Datatype );
      MPI_Type_free( &TimeIntervalStatisticsPacked::FullDatatype );
      
   }
   
   void peanoclaw::statistics::TimeIntervalStatisticsPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message peanoclaw::statistics::TimeIntervalStatisticsPacked "
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
            msg << "was not able to send message peanoclaw::statistics::TimeIntervalStatisticsPacked "
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
               msg << "testing for finished send task for peanoclaw::statistics::TimeIntervalStatisticsPacked "
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
               "peanoclaw::statistics::TimeIntervalStatisticsPacked",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "peanoclaw::statistics::TimeIntervalStatisticsPacked",
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
   
   
   
   void peanoclaw::statistics::TimeIntervalStatisticsPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive peanoclaw::statistics::TimeIntervalStatisticsPacked from node "
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
            msg << "failed to start to receive peanoclaw::statistics::TimeIntervalStatisticsPacked from node "
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
               msg << "testing for finished receive task for peanoclaw::statistics::TimeIntervalStatisticsPacked failed: "
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
               "peanoclaw::statistics::TimeIntervalStatisticsPacked",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "peanoclaw::statistics::TimeIntervalStatisticsPacked",
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
   
   
   
   bool peanoclaw::statistics::TimeIntervalStatisticsPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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



