/**
 * @file parallel_compuation.h
 * @author Konrad Merkel
 * @brief Routines to perform parallel compuations using mpi.
 *
 */

#include "mpi_tags.h"
#include "potential.h"
#include "solver.h"
#include "filehandler.h"
#include "scheduler.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <chrono>
#include <filesystem>
#include <mpi.h>

/**
 * @brief master can send a termination signal to all workers after work is done.
 *
 */
void sendTerminationSignal() {
    int N = -1;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank != 0) return;  // should only executed by master!

    cout << "Terminate all workers.\n";

    // send termination signal to workers
    for (int worker=1; worker < size; worker++)
        MPI_Send(&N, 1, MPI_INT, worker, TAG_NUM_TASKS, MPI_COMM_WORLD);
}


/**
 * @brief Sends a list of tasks to a single worker. The number of
 * tasks might be different compared to other workers.
 *
 */
int sendTasks(int worker, vector<Integral> &tasks)
{
    // convert to array (we only want one block of memory)
    int N = tasks.size();
    int M = 13;

    if (N<=0) {
        throw runtime_error("Cannot send empty tasklist. Do you want to terminate the worker instead?");
    }

    int l = 0;
    int* data = (int*) malloc(sizeof(int)*N*M);
    for (int n=0; n<N; n++){
        for (int m=0; m<M; m++) {
            data[l] = tasks[n].indexes[m];
            l++;
        }
    }

    // cout << "send " << N << " tasks to worker " << worker << endl;

    // first send number of tasks
    MPI_Send(&N, 1, MPI_INT, worker, TAG_NUM_TASKS, MPI_COMM_WORLD);

    return MPI_Send(data, N*M, MPI_INT, worker, TAG_DATA, MPI_COMM_WORLD);
}


/**
 * @brief Receives a list of tasks for the worker. The number of tasks can be flexible and
 * will also be received.
 *
 */
vector<Integral> recvTasks()
{
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // cout << rank << ": wait for new tasks\n";

    int N_tasks = 0;
    MPI_Recv(&N_tasks, 1, MPI_INT, 0, TAG_NUM_TASKS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // cout << rank << ": recv N_tasks=" << N_tasks << endl;

    if (N_tasks < 0) {  // termination signal
        return vector<Integral>(0);
    }

    // cout << "try to recv " << N_tasks << " tasks\n";
    int* data=(int*) malloc(sizeof(int)*N_tasks*13);
    MPI_Recv(data, N_tasks*13, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    // if (rank == 1) {
    //     cout << rank << ": RECV tasks:" << endl;
    //     for (int i=0; i<N_tasks*13; i++) {
    //         cout << data[i] << " ";
    //     }
    //     cout << endl;
    // }


    // convert to vector<vector<int>>
    int l = 0;
    vector<Integral>tasks(N_tasks);
    for (int n=0; n<N_tasks; n++) {
        //tasks[n].indexes = vector<int>(13);
        for (int i=0; i<13; i++) {
            tasks[n].indexes[i] = data[l];
            l++;
        }
    }
    return tasks;
}


/**
 * @brief Sends all results from the worker back to the master. It is not necessary
 * to communicate the length of the data array (number of integral) because the master
 * already knows.
 *
 */
int sendResults(vector<Integral>& tasks)
{
    // create array of results
    int N = tasks.size();
    double* data = (double*) malloc(sizeof(double)*N);

    for (int i=0; i<N; i++) {
        data[i] = tasks[i].value;
    }
    int ret = MPI_Send(data, N, MPI_DOUBLE, 0, TAG_DATA, MPI_COMM_WORLD);

    free(data);
    return ret;
}

/**
 * @brief Receives result from worker and adds the information to the tasks vector (for this specific worker)
 *
 */
void recvResults(int worker, vector<Integral>& tasks)
{
    int N = tasks.size();
    double* data = (double*) malloc(sizeof(double)*N);
    MPI_Recv(data, N, MPI_DOUBLE, worker, TAG_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i=0; i<N; i++){
        tasks[i].value = data[i];
    }
    free(data);
}

/**
 * @brief Calculates integrals in parallel.
 *
 * Integrals are obtained from a Scheduler at the master process and will be solved by all
 * MPI processes using a Solver. Results are send back to the master process and are stored
 * inside the Scheduler.
 *
 * @param impl
 * @param myScheduler
 * @param OMP_threads
 * @param outfile
 * @param restartfile
 */
void run_parallel_calculations(Solver& impl, Scheduler* myScheduler, const int OMP_threads, string const& outfile, string const& restartfile)
{
    int rank, num_worker;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);
    vector<Integral> myTasks;  // working packet for each worker (and master)


    /**
     *  Solve integrals
     **/
    if (rank==0) {  // master
        vector<vector<Integral>> tasks_worker;

        cout << "\nSave data in " << outfile << endl;
        cout << "Save progress in restart file " << restartfile << endl;
        cout << "Start evaluation of integrals...\n";

        if (myScheduler==nullptr) {
            throw runtime_error("Scheduler is null pointer.");
        }

        int counter = 0;
        while (! myScheduler->isFinished()) {

            // generate list of tasks for each worker
            // cout << "# Assign tasks to workers...";
            // chrono::steady_clock::time_point time_start = chrono::steady_clock::now();
            tasks_worker = myScheduler->getNextTasks();

            // print sizes (for debugging)
            // for (int worker=0; worker<num_worker; worker++){
            //     cout << "worker " << worker << ": size=" << tasks_worker[worker].size() << endl;
            // }
            // for (int worker=0; worker<num_worker; worker++){
            //     for (Integral const& i : tasks_worker[worker]) {
            //         cout << "worker " << worker << ": " << i.toString() << endl;
            //     }
            // }

            // distribute to workers
            myTasks = tasks_worker[0];
            for (int worker=1; worker<num_worker; worker++){
                if (tasks_worker[worker].size() > 0)
                    sendTasks(worker, tasks_worker[worker]);  // send tasks
            }

            // chrono::steady_clock::time_point time_stop = chrono::steady_clock::now();
            // cout << " done. Time spent (ms): " << chrono::duration_cast<chrono::microseconds>(time_stop - time_start).count()/1000. << "\n";
            // cout << "# Calculate batch of integrals...\n";

            // master performs its tasks
            impl.calculate(myTasks, false,1, OMP_threads);

            // collect all data
            // cout << "# Wait for workers to finish and collect results...\n";
            tasks_worker[0] = myTasks;
            for (int worker=1; worker<num_worker; worker++){
                if (tasks_worker[worker].size()) {  // only recive results if tasks where assigned
                    recvResults(worker, tasks_worker[worker]);
                }
            }

            // add integrals to list of known integrals and get maxValue of this batch
            myScheduler->update(tasks_worker, true);

            // save intermediate results to file
            counter++;
            if (counter >= 1) {
                counter = 0;

                // cout << "#\n# Save data in " << outfile << endl;
                OutputService::writeResults(outfile, myScheduler->getKnownIntegrals(), "intermediate result", myScheduler->getMaxDist());
                // cout << "#\n# Save state of scheduler" << endl;
                myScheduler->saveState(restartfile);
            }
        }
        // terminate workers
        cout << "\nCalculation finished.\n";
        sendTerminationSignal();

        // cout << "\nWrite data in " << outfile << endl;
        OutputService::writeResults(outfile, myScheduler->getKnownIntegrals(), "", myScheduler->getMaxDist());
        // cout << "#\n# Save state of scheduler" << endl;
        myScheduler->saveState(restartfile);

    } else {  // worker processes

        while (true) {
            myTasks = recvTasks();
            if (myTasks.size() == 0) {
                // cout << rank << " recv. termination signal!\n";
                break;
            }

            // all workers perform their tasks
            impl.calculate(myTasks, false,1, OMP_threads);

            // send back results
            sendResults(myTasks);
        }
    }
}




void gather_results(list<OpticalDipole>& dipoles) {

    int rank, num_worker;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);

    if (rank!=0) {  // worker

        //serialize data to send via MPI
        int N = dipoles.size();
        if (N==0) {  // worker has no data
            MPI_Send(&N, 1, MPI_INT, 0, TAG_DIPOLES_NUM, MPI_COMM_WORLD);
            return;
        }

        double* data = (double*)malloc(N*3*sizeof(double));
        int* indexes = (int*)malloc(N*5*sizeof(int));

        int i = 0;
        for (const auto& dipole : dipoles) {

            for (int j=0; j<5; j++) {
                indexes[5*i+j] = dipole.indexes[j];
            }

            vector<double> d = dipole.getDipole();
            data[3*i+0] = d[0];
            data[3*i+1] = d[1];
            data[3*i+2] = d[2];

            i++;
        }

        MPI_Send(&N, 1, MPI_INT, 0, TAG_DIPOLES_NUM, MPI_COMM_WORLD);
        MPI_Send(indexes, N*5, MPI_INT, 0, TAG_DIPOLES_INDEXES, MPI_COMM_WORLD);
        MPI_Send(data, N*3, MPI_DOUBLE, 0, TAG_DIPOLES_DATA, MPI_COMM_WORLD);

        free(data);
        free(indexes);
        return;

    } else {  // master process
        int* indexes = nullptr;
        double* data = nullptr;

        for (int worker = 1; worker < num_worker; worker++) {
            //cout << "worker " << worker << " of " << num_worker << endl;
            int N = 0;
            MPI_Recv(&N, 1, MPI_INT, worker, TAG_DIPOLES_NUM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //cout << "rank=" << rank <<" worker " << worker << " N=" << N << endl;

            if (N==0) continue; // worker has no data

            indexes=(int*) malloc(sizeof(int)*N*5);
            MPI_Recv(indexes, N*5, MPI_INT, worker, TAG_DIPOLES_INDEXES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            data=(double*) malloc(sizeof(double)*N*3);
            MPI_Recv(data, N*3, MPI_DOUBLE, worker, TAG_DIPOLES_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            // add recived data to map
            for (int i=0; i<N; i++) {

                vector<int> ids(5);
                for (int j=0; j<5; j++) ids[j] = indexes[5*i+j];

                vector<double> value(3);
                for (int j=0; j<3; j++) value[j] = data[3*i+j];

                dipoles.push_back( OpticalDipole(ids, value) );
            }

            // free memory
            free(data);
            free(indexes);
            data = nullptr;
            indexes = nullptr;
        }
        return;
    }
}
