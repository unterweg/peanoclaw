#ifndef _PEANOCLAW_RECORDS_DOUBLEDATA_H
#define _PEANOCLAW_RECORDS_DOUBLEDATA_H

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
   namespace records {
      class DoubleData;
      class DoubleDataPacked;
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
 * @date   12/05/2015 20:06
 */
class peanoclaw::records::DoubleData { 
   
   public:
      
      typedef peanoclaw::records::DoubleDataPacked Packed;
      
      struct PersistentRecords {
         double _u;
         /**
          * Generated
          */
         PersistentRecords();
         
         /**
          * Generated
          */
         PersistentRecords(const double& u);
         
         
         inline double getU() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _u;
         }
         
         
         
         inline void setU(const double& u) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _u = u;
         }
         
         
         
      };
      
   private: 
      public:

      PersistentRecords _persistentRecords;
      private:

      
   public:
      /**
       * Generated
       */
      DoubleData();
      
      /**
       * Generated
       */
      DoubleData(const PersistentRecords& persistentRecords);
      
      /**
       * Generated
       */
      DoubleData(const double& u);
      
      /**
       * Generated
       */
      ~DoubleData();
      
      
      inline double getU() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         return _persistentRecords._u;
      }
      
      
      
      inline void setU(const double& u) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
         _persistentRecords._u = u;
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
      DoubleDataPacked convert() const;
      
      
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
          * @date   12/05/2015 20:06
          */
         class peanoclaw::records::DoubleDataPacked { 
            
            public:
               
               struct PersistentRecords {
                  double _u;
                  /**
                   * Generated
                   */
                  PersistentRecords();
                  
                  /**
                   * Generated
                   */
                  PersistentRecords(const double& u);
                  
                  
                  inline double getU() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     return _u;
                  }
                  
                  
                  
                  inline void setU(const double& u) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                     _u = u;
                  }
                  
                  
                  
               };
               
            private: 
               PersistentRecords _persistentRecords;
               
            public:
               /**
                * Generated
                */
               DoubleDataPacked();
               
               /**
                * Generated
                */
               DoubleDataPacked(const PersistentRecords& persistentRecords);
               
               /**
                * Generated
                */
               DoubleDataPacked(const double& u);
               
               /**
                * Generated
                */
               ~DoubleDataPacked();
               
               
               inline double getU() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  return _persistentRecords._u;
               }
               
               
               
               inline void setU(const double& u) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
                  _persistentRecords._u = u;
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
               DoubleData convert() const;
               
               
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
                  
