#ifndef _PEANOCLAW_STATISTICS_PROCESSSTATISTICSENTRY_H
#define _PEANOCLAW_STATISTICS_PROCESSSTATISTICSENTRY_H

#include "peano/utils/Globals.h"
#include "tarch/compiler/CompilerSpecificSettings.h"
#include "peano/utils/PeanoOptimisations.h"
#ifdef Parallel
	#include "tarch/parallel/Node.h"
#endif
#ifdef Parallel
	#include <mpi.h>
#endif
#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"
#include <bitset>
#include <complex>
#include <string>
#include <iostream>

namespace peanoclaw {
   namespace statistics {
      class ProcessStatisticsEntry;
      class ProcessStatisticsEntryPacked;
   }
}

/**
 * @author This class is generated by DaStGen
 * 		   DataStructureGenerator (DaStGen)
 * 		   2007-2009 Wolfgang Eckhardt
 * 		   2012      Tobias Weinzierl
 *
 * 		   build date: 09-02-2014 14:40
 *
 * @date   02/06/2014 13:19
 */
class peanoclaw::statistics::ProcessStatisticsEntry { 
   
   public:
      
      typedef peanoclaw::statistics::ProcessStatisticsEntryPacked Packed;
      
      struct PersistentRecords {
         int _rank;
         int _numberOfCellUpdates;
         /**
          * Generated
          */
         PersistentRecords();
         
         /**
          * Generated
          */
         PersistentRecords(const int& rank, const int& numberOfCellUpdates);
         
         
         inline int getRank() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _rank;
         }
         
         
         
         inline void setRank(const int& rank) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _rank = rank;
         }
         
         
         
         inline int getNumberOfCellUpdates() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _numberOfCellUpdates;
         }
         
         
         
         inline void setNumberOfCellUpdates(const int& numberOfCellUpdates) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _numberOfCellUpdates = numberOfCellUpdates;
         }
         
         
         
      };
      
   private: 
      PersistentRecords _persistentRecords;
      
   public:
      /**
       * Generated
       */
      ProcessStatisticsEntry();
      
      /**
       * Generated
       */
      ProcessStatisticsEntry(const PersistentRecords& persistentRecords);
      
      /**
       * Generated
       */
      ProcessStatisticsEntry(const int& rank, const int& numberOfCellUpdates);
      
      /**
       * Generated
       */
      ~ProcessStatisticsEntry();
      
      
      inline int getRank() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         return _persistentRecords._rank;
      }
      
      
      
      inline void setRank(const int& rank) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         _persistentRecords._rank = rank;
      }
      
      
      
      inline int getNumberOfCellUpdates() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         return _persistentRecords._numberOfCellUpdates;
      }
      
      
      
      inline void setNumberOfCellUpdates(const int& numberOfCellUpdates) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         _persistentRecords._numberOfCellUpdates = numberOfCellUpdates;
      }
      
      
      /**
       * Generated
       */
      std::string toString() const;
      
      /**
       * Generated
       */
      void toString(std::ostream& out) const;
      
      
      PersistentRecords getPersistentRecords() const;
      /**
       * Generated
       */
      ProcessStatisticsEntryPacked convert() const;
      
      
   #ifdef Parallel
      protected:
         static tarch::logging::Log _log;
         
      public:
         
         /**
          * Global that represents the mpi datatype.
          * There are two variants: Datatype identifies only those attributes marked with
          * parallelise. FullDatatype instead identifies the whole record with all fields.
          */
         static MPI_Datatype Datatype;
         static MPI_Datatype FullDatatype;
         
         /**
          * Initializes the data type for the mpi operations. Has to be called
          * before the very first send or receive operation is called.
          */
         static void initDatatype();
         
         static void shutdownDatatype();
         
         void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
         
         void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
         
         static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
         
         #endif
            
         };
         
         /**
          * @author This class is generated by DaStGen
          * 		   DataStructureGenerator (DaStGen)
          * 		   2007-2009 Wolfgang Eckhardt
          * 		   2012      Tobias Weinzierl
          *
          * 		   build date: 09-02-2014 14:40
          *
          * @date   02/06/2014 13:19
          */
         class peanoclaw::statistics::ProcessStatisticsEntryPacked { 
            
            public:
               
               struct PersistentRecords {
                  int _rank;
                  int _numberOfCellUpdates;
                  /**
                   * Generated
                   */
                  PersistentRecords();
                  
                  /**
                   * Generated
                   */
                  PersistentRecords(const int& rank, const int& numberOfCellUpdates);
                  
                  
                  inline int getRank() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     return _rank;
                  }
                  
                  
                  
                  inline void setRank(const int& rank) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     _rank = rank;
                  }
                  
                  
                  
                  inline int getNumberOfCellUpdates() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     return _numberOfCellUpdates;
                  }
                  
                  
                  
                  inline void setNumberOfCellUpdates(const int& numberOfCellUpdates) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     _numberOfCellUpdates = numberOfCellUpdates;
                  }
                  
                  
                  
               };
               
            private: 
               PersistentRecords _persistentRecords;
               
            public:
               /**
                * Generated
                */
               ProcessStatisticsEntryPacked();
               
               /**
                * Generated
                */
               ProcessStatisticsEntryPacked(const PersistentRecords& persistentRecords);
               
               /**
                * Generated
                */
               ProcessStatisticsEntryPacked(const int& rank, const int& numberOfCellUpdates);
               
               /**
                * Generated
                */
               ~ProcessStatisticsEntryPacked();
               
               
               inline int getRank() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  return _persistentRecords._rank;
               }
               
               
               
               inline void setRank(const int& rank) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  _persistentRecords._rank = rank;
               }
               
               
               
               inline int getNumberOfCellUpdates() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  return _persistentRecords._numberOfCellUpdates;
               }
               
               
               
               inline void setNumberOfCellUpdates(const int& numberOfCellUpdates) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  _persistentRecords._numberOfCellUpdates = numberOfCellUpdates;
               }
               
               
               /**
                * Generated
                */
               std::string toString() const;
               
               /**
                * Generated
                */
               void toString(std::ostream& out) const;
               
               
               PersistentRecords getPersistentRecords() const;
               /**
                * Generated
                */
               ProcessStatisticsEntry convert() const;
               
               
            #ifdef Parallel
               protected:
                  static tarch::logging::Log _log;
                  
               public:
                  
                  /**
                   * Global that represents the mpi datatype.
                   * There are two variants: Datatype identifies only those attributes marked with
                   * parallelise. FullDatatype instead identifies the whole record with all fields.
                   */
                  static MPI_Datatype Datatype;
                  static MPI_Datatype FullDatatype;
                  
                  /**
                   * Initializes the data type for the mpi operations. Has to be called
                   * before the very first send or receive operation is called.
                   */
                  static void initDatatype();
                  
                  static void shutdownDatatype();
                  
                  void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
                  
                  void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
                  
                  static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
                  
                  #endif
                     
                  };
                  
                  #endif
                  
