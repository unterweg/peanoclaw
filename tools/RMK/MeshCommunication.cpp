#include <unistd.h>

#include <tarch/parallel/NodePool.h>

#include <peano/parallel/MeshCommunication.h>

#ifdef Parallel
tarch::logging::Log MeshCommunication::_log("MeshCommunication");

void MeshCommunication::sendAsynchronousAndBlocking(void *data, int size, const MPI_Datatype& datatype, int destination, int tag, bool exchangeDataBlocking) {
    if (exchangeDataBlocking) {
      
        const int result = MPI_Send(data, size, datatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
        assertionEquals1(result, MPI_SUCCESS, tarch::parallel::MPIReturnValueToString(result));
     }
     else {
        MPI_Request* sendRequestHandle = new MPI_Request();
        MPI_Status   status;
        int          flag = 0;
        int          result;
      
        clock_t      timeOutWarning   = -1;
        clock_t      timeOutShutdown  = -1;
        bool         triggeredTimeoutWarning = false;
      
        result = MPI_Isend(
            data, size, datatype, destination,
            tag, tarch::parallel::Node::getInstance().getCommunicator(),
            sendRequestHandle
        );
        assertionEquals1(result, MPI_SUCCESS, tarch::parallel::MPIReturnValueToString(result));
        
        result = MPI_Test( sendRequestHandle, &flag, &status );
        assertionEquals1(result, MPI_SUCCESS, tarch::parallel::MPIReturnValueToString(result));
        while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            assertionEquals1(result, MPI_SUCCESS, tarch::parallel::MPIReturnValueToString(result));
         
            // deadlock aspect
            if (
                tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                (clock()>timeOutWarning) &&
                (!triggeredTimeoutWarning)
            ) {
                tarch::parallel::Node::getInstance().writeTimeOutWarning(
                    "poissonwithjacobi_wittmanr::records::VertexPacked",
                    "send(int)", destination,tag,1
                );
                triggeredTimeoutWarning = true;
            }

            if (
                tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                (clock()>timeOutShutdown)
            ) {
                tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                    "poissonwithjacobi_wittmanr::records::VertexPacked",
                    "send(int)", destination,tag,1
                );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
        }
      
        delete sendRequestHandle;
    }
}

void MeshCommunication::receiveAsynchronousAndBlocking(void *data, int size, const MPI_Datatype& datatype, int source, int tag, bool exchangeDataBlocking) {
    if (exchangeDataBlocking) {
      
        MPI_Status  status;
        const int   result = MPI_Recv(data, size, datatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
        assertionEquals1(result, MPI_SUCCESS, tarch::parallel::MPIReturnValueToString(result));
     }
     else {
        MPI_Request* sendRequestHandle = new MPI_Request();
        MPI_Status   status;
        int          flag = 0;
        int          result;
          
        clock_t      timeOutWarning   = -1;
        clock_t      timeOutShutdown  = -1;
        bool         triggeredTimeoutWarning = false;
          
        result = MPI_Irecv(
          data, size, datatype, source, tag,
          tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
        );
        assertionEquals1(result, MPI_SUCCESS, tarch::parallel::MPIReturnValueToString(result));

        int sleepDelay = 10;
        result = MPI_Test( sendRequestHandle, &flag, &status );
        while (!flag) {
           if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
           if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
           result = MPI_Test( sendRequestHandle, &flag, &status );
           assertionEquals1(result, MPI_SUCCESS, tarch::parallel::MPIReturnValueToString(result));
             
           // deadlock aspect
           if (
              tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
              (clock()>timeOutWarning) &&
              (!triggeredTimeoutWarning)
           ) {
              tarch::parallel::Node::getInstance().writeTimeOutWarning(
              "poissonwithjacobi_wittmanr::records::VertexPacked",
              "receive(int)", source,tag,1
              );
              triggeredTimeoutWarning = true;
           }
           if (
              tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
              (clock()>timeOutShutdown)
           ) {
              tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
              "poissonwithjacobi_wittmanr::records::VertexPacked",
              "receive(int)", source,tag,1
              );
           }
           tarch::parallel::Node::getInstance().receiveDanglingMessages();

           usleep(sleepDelay);
           sleepDelay = std::min(3000, sleepDelay + 10);
        }
          
        delete sendRequestHandle;
     }
}
#endif
