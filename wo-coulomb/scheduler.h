#ifndef SCHEDULER_H
#define SCHEDULER_H

/**
 * @file scheduler.h
 * @author Konrad Merkel
 * @brief Managing known and unkown integrals
 *
 */

#include "coulombIntegral.h"
#include "filehandler.h"
#include "density.h"

#include <vector>
#include <list>
#include <map>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <filesystem>

using namespace std;

/**
 * Helper function for creating shells (especially for sorting each shell)
 */
bool compare_R_nosign(const vector<int>& first, const vector<int>& second) {

    int sum_first=0;
    int sum_second=0;

    for (auto& n : first) sum_first += abs(n);
    for (auto& n : second) sum_second += abs(n);

    return sum_first<sum_second;
}

double getMaxRadiusInsideCell(vector< vector<double> > const& lattice, vector<bool> const& CRYSTAL_PERIODIC,
    vector<vector<double>> const& pos_c, vector<vector<double>> const& pos_v, vector<double> supercell_d,
    vector<double> const& origin, vector<vector<double>> const& unitcell)
{
    //
    // 1st estimate the largest shere that fits inside the supercell
    //
    double volume = det3x3(lattice);

    // normal vectors of planes
    vector<double> axb = crossProduct(lattice[0], lattice[1]);
    vector<double> axc = crossProduct(lattice[0], lattice[2]);
    vector<double> bxc = crossProduct(lattice[1], lattice[2]);

    // distances of the middle point M=(a+b+c)/2 to each plane
    vector<double> distances{};
    if (CRYSTAL_PERIODIC[2]) distances.push_back(volume/(2.*L2Norm(axb)));
    if (CRYSTAL_PERIODIC[1]) distances.push_back(volume/(2.*L2Norm(axc)));
    if (CRYSTAL_PERIODIC[0]) distances.push_back(volume/(2.*L2Norm(bxc)));

    if (distances.size() == 0) throw runtime_error("Cannot dermine MAX_ALLOWED_DISTANCE since no crystal-direction is marked as periodic!");

    double max_allowed_distance = *min_element(distances.begin(), distances.end());


    //
    // 2nd check if actual distances with Wannier functions outside the supercell are smaller
    //

    assert(supercell_d.size() == 3);
    // bounds of the inner supercell
    vector<int> minR{
        int(round(origin[0])),
        int(round(origin[1])),
        int(round(origin[2])),
    };

    vector<int> maxR{
        int(round(supercell_d[0] + origin[0])),
        int(round(supercell_d[1] + origin[1])),
        int(round(supercell_d[2] + origin[2])),
    };

    vector<int> supercell_g = vector<int>{maxR[0], maxR[1], maxR[2]};  // upper limit
    vector<int> supercell_l = vector<int>{minR[0]-1, minR[1]-1, minR[2]-1}; // lower limit

    // make sure periodic boundary conditions are considered
    for (size_t i=0; i<3; i++) {
        if (CRYSTAL_PERIODIC[i] == 0) {
            supercell_g[i] = 0;
            supercell_l[i] = 0;
        }
    }

    // cout << "minR = " << minR[0] << " " << minR[1] << " " << minR[2] << endl;
    // cout << "maxR = " << maxR[0] << " " << maxR[1] << " " << maxR[2] << endl;
    // cout << "supercell = " << supercell_d[0] << " " << supercell_d[1] << " " << supercell_d[2] << endl;
    // cout << "supercell_g = " << supercell_g[0] << " " << supercell_g[1] << " " << supercell_g[2] << endl;
    // cout << "supercell_l = " << supercell_l[0] << " " << supercell_l[1] << " " << supercell_l[2] << endl;

    // go through all Wannier functions outside the supercell and check the distance
    for (int Sx=supercell_l[0]; Sx<=supercell_g[0]; Sx++) {
        for (int Sy=supercell_l[1]; Sy<=supercell_g[1]; Sy++) {
            for (int Sz=supercell_l[2]; Sz<=supercell_g[2]; Sz++) {

                if (Sx >= minR[0] && Sx < maxR[0] &&
                    Sy >= minR[1] && Sy < maxR[1] &&
                    Sz >= minR[2] && Sz < maxR[2]) continue;  // skip if inside of supercell

                vector<double> vec_shift = {double(Sx), double(Sy), double(Sz)};
                vec_shift = matVecMul3x3(unitcell, vec_shift);

                for (size_t v=0; v<pos_v.size(); v++) {
                    for (size_t c=0; c<pos_c.size(); c++) {
                        double dist = L2Norm(vector<double>{
                            pos_c[c][0]-pos_v[v][0]-vec_shift[0],
                            pos_c[c][1]-pos_v[v][1]-vec_shift[1],
                            pos_c[c][2]-pos_v[v][2]-vec_shift[2]
                        });

                        if (dist < max_allowed_distance) {
                            max_allowed_distance = dist - 1e-5;  // -1e-5 to avoid numerical problems
                            // cout << "New max_allowed_distance: " << max_allowed_distance << " Angström" << endl;
                            // cout << "S = " << Sx << " " << Sy << " " << Sz << endl;
                            // cout << "c = " << c << " v = " << v << endl;
                        }
                    }
                }
            }
        }
    }

    return max_allowed_distance;
}

/**
 * @brief Create a shells of real space lattice vectors.
 *
 * Real space lattice vectors like Rc,Rv or RD are iterated in shells starting
 * at the origin unit cell (0,0,0) and its neighbors. Every shell contains a set
 * of vectors (unit cells) that have a specific distance to the origin. Shells
 * do not overlap, i.e. they form disjunct sets of unit cells.
 *
 * @param origin    origin of the meshgrid in units of lattice vectors (c.f. WannierFunction::getOriginInUnitcellBasis())
 * @param supercell supercell of in units of lattice vectors (c.f. WannierFunction::getLatticeInUnitcellBasis())
 * @param unitcell  in Angstrom
 * @param r_v       position of one valence Wannier function in Angstrom (Wannier center)
 * @param r_c       position of one conduction Wannier function in Angstrom (Wannier center)
 * @param factor    scaling factor
 * @return list<list<vector<int>>>
 */
vector<vector<vector<int>>> createShells(vector<double> const& origin, vector<double> const& supercell, vector<vector<double>> const& unitcell,
                    vector<double>const& r_v, vector<double>const& r_c, vector<bool> const& CRYSTAL_PERIODIC=vector<bool>{true, true, true}, double factor=1.0)
{
    vector<int> minR{
        int(round(origin[0])),
        int(round(origin[1])),
        int(round(origin[2])),
    };

    vector<int> maxR{
        int(round(supercell[0] + origin[0])),
        int(round(supercell[1] + origin[1])),
        int(round(supercell[2] + origin[2])),
    };

    // turn off directions where we dont have periodic boundary conditions
    int N_unitcells = 1;
    for (size_t i=0; i<3; i++) {
        if (! CRYSTAL_PERIODIC[i]) {
            maxR[i] = 1;
            minR[i] = 0;
        }else{
            N_unitcells = N_unitcells * round(supercell[i]);
        }
    }

    vector<vector<double>> unitcell_T = transpose3x3(unitcell);
    vector<vector<vector<int>>> shells;
    double dx= sqrt( pow(unitcell[0][0],2) + pow(unitcell[0][1],2) + pow(unitcell[0][2],2) )/factor;
    double minDist=-0.1, maxDist=dx;

    int n = 0;
    int sum = 0;
    while (sum < N_unitcells)
    {
        // create new shell
        vector<vector<int>> currentShell;

        // go over entire supercell and filter for distances between minDist and maxDist
        for (int x=minR[0]; x<maxR[0]; x++) {
            for (int y=minR[1]; y<maxR[1]; y++) {
                for (int z=minR[2]; z<maxR[2]; z++) {
                    vector<double> r_vec = matVecMul3x3(unitcell_T, vector<double>{double(x),double(y),double(z)});
                    double dist = sqrt(pow(r_c[0]- r_v[0] - r_vec[0],2)+pow( r_c[1]- r_v[1] - r_vec[1],2)+pow( r_c[2]- r_v[2] - r_vec[2],2));

                    if ((dist <= maxDist) && (dist > minDist)) {
                        currentShell.push_back(vector<int>{x,y,z});
                        sum++;
                    }

                }
            }
        }
        if (currentShell.size() > 0) {
            sort(currentShell.begin(), currentShell.end(), compare_R_nosign);
            shells.push_back(currentShell);
        }

        //cout << "Shell " << n << " size: " << currentShell.size() << " minDist: " << minDist << " maxDist: " << maxDist << endl;

        // update distance range for new shell
        n++;
        minDist = maxDist;
        maxDist += dx/sqrt(n);

    }
    //cout << "sum: " << sum << endl;
    if (sum != N_unitcells) {
        throw runtime_error("Creation of shells: Number of unitcells does not match!");
    }

    return shells;

}


/**
 * @brief General interface of any Scheduler object
 *
 * A scheduler has two central duties:
 *   - It figures out which integrals are relevant to calculate, i.e. it
 *     uses estimates and relations between integrals to generate a list
 *     of integrals that might be worth calculating (see getNextTasks()).
 *   - It stores calculated Coulomb/LFE integrals and is therefore the
 *     main data hub (see update() and getKnownIntegrals()).
 *
 * There should only be one scheduler. In case of MPI parallelization the
 * scheduler should only be on the master node (rank=0). Worker nodes do not
 * have a scheduler.
 *
 * A typical workflow would be done iteratively:
 *   - In every iteration the scheduler proposes a list of batch_size
 *     integrals for each worker that need to be calculated using the
 *     getNextTasks() method.
 *   - After all integrals are calculated the scheduler is updated using
 *     the update() method, i.e. all known integrals are saved in the
 *     scheduler and can be used for further estimations.
 *   - Using the isFinished() method the scheduler can determine if all
 *     relevant integrals have been found.
 *
 * The distribution of each list to the workers and calculation of integrals
 * is NOT part of the scheduler and need to be done somewhere else.
 *
 */
class Scheduler
{
protected:
    int batch_size;     //!< number of generated tasks for each worker
    int num_worker;     //!< number of workers (usually number of MPI processes)
    map<vector<int>, Integral> solved_ints; //!< all integrals that are known (maybe also from previous runs or different schedulers)
                                            //!< key = Integral::indexes

public:
    Scheduler(int batch_size, int num_worker) : batch_size(batch_size), num_worker(num_worker) {}
    virtual ~Scheduler() {}

    /**
     * @brief Generates a list of integrals for each worker that need to be solved in the next iteration.
     *
     * Which integrals are considered as relevant depends on the type of scheduler and its parameters.
     *
     * @return vector<vector<Integral>> unknown integrals, dimensions=(num_workers, batch_size or smaller)
     */
    virtual vector<vector<Integral>> getNextTasks() =0;  // returns a list of Integrals for each worker

    /**
     * @brief Updates / adds solved integrals to internal data structure (solved_ints).
     *
     * @param solved_integrals  dimensions=(num_workers, batch_size)
     * @param verbose           prints integrals and additional information
     */
    virtual void update( vector<vector<Integral>> const& solved_integrals, bool verbose=true ) =0;

    /**
     * @brief Determines if calculations are finished and all relevant integrals are calculated
     *
     * Which integrals are relevant depends on the type of scheduler and its parameters.
     *
     * @return true  all relevant integrals are calculated
     * @return false some relevant integrals might be unknown
     */
    virtual bool isFinished() const=0;

    /**
     * @brief Get a list of all known Integrals.
     */
    virtual const list<Integral> getKnownIntegrals() const
    {
        // convert to list
        list<Integral> knownInts;
        for (auto itr : solved_ints) {
            knownInts.push_back(itr.second);
        }
        return knownInts;
    }

    /**
     * @brief Get a map of all known Integrals.
     */
    virtual const map<vector<int>, Integral> getKnownIntegrals_asmap() const {return solved_ints;}

    /**
     * @brief Set the known Integrals manually
     */
    virtual void setKnownIntegrals(const vector<Integral>& knownInts) {
        this->solved_ints = map<vector<int>, Integral>();
        for (const Integral& i : knownInts) {
            this->solved_ints.insert({i.indexes, i});
        }
    }
    virtual void setKnownIntegrals(const map<vector<int>, Integral>& knownInts) {
        this->solved_ints = knownInts;  // TODO
    }

    /**
     * @brief Get the maximal electron-hole distance
     *
     * This is only important for Coulomb integrals (not for LFE).
     * The returned value can be used for the head of the COULOMB file.
     *
     */
    virtual double getMaxDist() const {return -1;}


    int getBatchSize() const {return this->batch_size;}

    /**
     * @brief Provides a short summary of what the scheduler tries to do
     *
     * It should contain all used parameters and a detailed description
     * of the purpose. This should help the user to understand what is
     * happending.
     *
     */
    virtual void printStartupSummary() const =0;

    /**
     * @brief Tries to restart the scheduler from a previous run (RESTART file)
     *
     * This is only a wrapper for the loadState() method. It gives more user output
     * and a propper error handling. The program exits when the restart_file exists
     * but cannot be loaded (e.g. wrong scheduler).
     *
     * @param restart_file  file that contains the state
     * @param coulomb_file  file that contains all known integrals
     * @return true         success
     * @return false        restart file does not exist
     */
    virtual bool tryRestartScheduler(string restart_file, string coulomb_file) {
        // load from previous calculation (plan file is not used)
        if (filesystem::exists(restart_file)) {
            cout << "Found restart file. Try to continue from previous calculation using \n";
            cout << "\trestart file: " << restart_file << endl;
            cout << "\tcoulomb file: " << coulomb_file << endl;
            cout << "Plan file will be ignored. Use information from restart file only." << endl;
            if (this->loadState(restart_file, coulomb_file)) {
                cout << "Successfully loaded previous state of the scheduler.\n";
                return true;
            }else{
                cout << "Could not load previous state of the scheduler. --> Quit this sick job!!!\n\n";
                cout << "You can start from scratch and use the plan file if you remove the restart file first.\n";
                exit(1);
            }
        }else {
            return false;
        }
    }

    /**
     * @brief Writes the internal state of the scheduler to a file.
     *
     * This enables the program to stop and resume later. The restart
     * file does not contain the known integrals. These are stored in
     * the COULOMB file.
     * The restart_file can be very different for each scheduler.
     *
     * @param restart_file  the file where the state is saved
     * @return true         success
     * @return false        failed
     */
    virtual bool saveState(string restart_file) const =0;

    /**
     * @brief Loads the internal state from a restart file.
     *
     * The restart_file can be very different for each scheduler.
     * It is usually not possible to load from restart_files from
     * other schedulers.
     *
     * @param restart_file  file that contains the state
     * @param coulomb_file  file that contains all known integrals
     * @return true         success
     * @return false        failed
     */
    virtual bool loadState(string restart_file, string coulomb_file)=0;
};


/**
 * @brief Scheduler that processes a fixed list of integrals.
 *
 * This scheduler uses a user defined, fixed list of integrals and splits it
 * into batches for each worker.
 * All Integrals of the list are calculated. Integrals that are not in the list
 * are not calculated. The scheduler is finished when all integrals of the list
 * are calculated.
 *
 * The list can contain Coulomb or LFE integrals but not both at the
 * same time. This scheduler does not use any estimates / relation between
 * integrals or heuristics.
 *
 */
class FixedListScheduler : public Scheduler
{
protected:
    vector<Integral> plan;   //!< fixed list of relevant integrals
    int curser;              //!< current position in the list (next unknown integral)
    double maxDistance;      //!< maximal electron-hole distance

    /**
     * @brief Helper function that reshapes a vector of tasks into multiple vectors for each worker
     *
     * If the tasks.size() is not dividable by num_workers then empty integrals are added such that
     * every worker has the same amount of integrals (padding).
     *
     * @param tasks         list of tasks in this iteration
     * @param num_worker    number of workers
     * @return vector<vector<Integral>> dimensions: (num_workers, N / num_workers)
     */
    vector<vector<Integral>> splitTasksToWorkers(vector<Integral> tasks)
    {
        vector<vector<Integral>> worker_lists(num_worker);
        int N = tasks.size();
        int N_per_worker = N/num_worker;

        vector<int> sizes(num_worker);
        for (int i =0; i<num_worker; i++) {
            sizes[i] = N_per_worker;
        }
        // rest if num_worker is not a divisor of N
        for (int i=0; i< N-N_per_worker*num_worker; i++) {
            sizes[i]++;
        }

        // print sizes
        // cout << "sizes of distributed tasks: ";
        // for (int i =0; i<num_worker; i++) {
        //     cout << sizes[i] << " ";
        // }
        // cout << endl;

        int l = 0;
        for (int i=0; i<num_worker; i++) {

            // create a vector of empty integrals for each worker
            worker_lists[i] = vector<Integral>(sizes[i]);
            // worker_lists[i] = vector<Integral>(this->batch_size);

            // assign tasks (integrals that need to be solved)
            for (int j=0; j<sizes[i]; j++) {
                worker_lists[i][j] = tasks[l];
                l++;
            }

            // create padding such that we always have the same dimensions of the output vectors
            // for (int j=sizes[i]; j<this->batch_size; j++) {
            //     worker_lists[i][j] = Integral::CreateEmptyIntegral();
            // }
        }
        return worker_lists;
    }

public:

    /**
     * @brief Construct a new FixedListScheduler from a list of unsolved integrals.
     *
     * @param batch_size
     * @param num_worker
     * @param plan_         fixed list of integrals that need to be solved
     */
    FixedListScheduler(int batch_size, int num_worker, vector<Integral> const& plan_=vector<Integral>{})
    : Scheduler(batch_size, num_worker), plan(plan_), curser(0), maxDistance(-1)
    {}

    /**
     * @brief Set the vector of integrals that needs to be calculated.
     */
    void setIntegrals(vector<Integral> const& new_plan) {
        this->plan = new_plan;
        curser = 0;  // reset curser for new list of integrals
    }


    vector<vector<Integral>> getNextTasks() override {
        int N = min(batch_size*num_worker, int(plan.size())-curser);
        vector<Integral> tasks(N);
        tasks.assign(plan.begin()+curser, plan.begin()+curser+N);
        curser = curser + N;

        return splitTasksToWorkers(tasks);
    }

    void update( vector<vector<Integral>> const& solved_integrals, bool verbose=true ) override {
        double maxValue = 0.0;
        for (size_t i=0; i<solved_integrals.size(); i++){
            for (size_t j=0; j<solved_integrals[i].size(); j++){

                if (solved_integrals[i][j].isEmpty()) continue;  // skip No-Integrals
                solved_ints.insert({solved_integrals[i][j].indexes, solved_integrals[i][j]});

                maxValue = max(maxValue, abs(solved_integrals[i][j].value));
                // if (verbose) solved_integrals[i][j].print();
            }
        }
        if (verbose) {
            // cout << "# End of mini-batch.\n";
            // cout << "# Max value of this batch: " << maxValue << endl;
            // cout << "# Remaining size of plan = " << plan.size()-curser << endl;
            cout << "\r"
             << "Solved integrals: " << solved_ints.size() << " of " << plan.size()
             << " (" << solved_ints.size()*100.0/ plan.size()<< "%);  "
             << "max value of batch: " << maxValue << "eV"
             << flush;
        }
    }

    bool isFinished() const override {
        return curser >= int(plan.size());
    }

    bool saveState(string state_file) const override {
        ofstream file(state_file);
        if (!file.is_open()) {
            cout << "Cannot open file!" << endl;
            return false;
        }

        file << "FixedListScheduler\n";
        file << maxDistance << endl;
        file << curser << "\t" << plan.size() << endl;
        file << "Plan:\n";
        for (const Integral& e: plan) {
            file << e.toString() << endl;
        }
        file << "Solved integrals:" << endl;
        file << this->solved_ints.size() << endl;

        file.close();
        return true;
    }

    bool loadState(string state_file, string coulomb_file) override {


        // open Coulomb file if exits
        if (!filesystem::exists(coulomb_file)) {
            cerr << "Cannot open file: " << coulomb_file << " because it does not exist!"<< endl;
            return false;
        }
        cout << "Read previous calculated integrals from " << coulomb_file << endl;
        auto f = CoulombFileReader(coulomb_file);
        if (! f.wasSuccessful()) {
            cerr << "Cannot read integrals from file.\n";
            return false;
        }
        const vector<Integral>& integrals = f.getElements();

        if (!filesystem::exists(state_file)) {
            cerr << "Cannot open file: " << state_file << " because it does not exist!"<< endl;
            return false;
        }

        ifstream file(state_file);
        if (!file.is_open()) {
            cerr << "Cannot open file: " << state_file << endl;
            file.close();
            return false;
        }

        string line;
        getline(file, line);
        if (line != "FixedListScheduler") {
            cerr << "Cannot load state of the scheduler of type " << line << " (should be FixedListScheduler)\n";
            file.close();
            return false;
        }

        getline(file, line); // maxDistance
        double new_maxDistance = stod(line);

        getline(file, line);  // curser and plan.size()
        vector<string> v = splitBySpacesAndTabs(line);//split(line, string("\t"));
        if (v.size() != 2) {
            cerr << "Cannot load state of the scheduler: 2 arguments required but "<< v.size() << " are given.\n";
            file.close();
            return false;
        }

        int new_curser = stoi(v[0]);
        int Num_plan = stoi(v[1]);

        getline(file, line);
        if (line != "Plan:") {
            cerr << "Error while parsing plan in restart file (cannot find keyword)\n";
            file.close();
            return false;
        }

        vector<Integral> new_plan(Num_plan);
        for (int i=0; i<Num_plan; i++) {
            getline(file, line);

            vector<string> v = splitBySpacesAndTabs(line); //split(line, string("\t"));
            if (v.size() < 13) {
                cerr << "Error while parsing plan in restart file\n";
                file.close();
                return false;
            }
            new_plan[i] = Integral(stoi(v[0]), stoi(v[1]), stoi(v[2]), stoi(v[3]), stoi(v[4]), stoi(v[5]),
                            stoi(v[6]), stoi(v[7]), stoi(v[8]), stoi(v[9]), stoi(v[10]), stoi(v[11]), stoi(v[12]));
        }


        // read solved ints
        getline(file, line);
        if (line != "Solved integrals:") {
            cerr << "Error while parsing solved integrals (cannot find keyword)\n";
            file.close();
            return false;
        }
        getline(file, line);
        int Num_ints = stoi(line);
        cout << "Number of solved ints: " << Num_ints << endl;
        file.close();

        if (Num_ints != int(integrals.size())) {
            cerr << "Cannot load state of the scheduler: Expected "<< Num_ints << " calculated integrals but found " << integrals.size() << ".\n";
            return false;
        }

        // actually set state
        this->setKnownIntegrals(integrals);
        this->curser = new_curser;
        this->maxDistance = new_maxDistance;
        this->plan = new_plan;
        cout << "set curser to " << this->curser << endl;

        return true;
    }

    void printStartupSummary() const override {
        cout << "You calculate integrals from a fixed list:\n";
        printPlan(this->plan);
        cout << "All threshold parameters (DIFF_MONO_THRESHOLD, DISTANCE_THRESHOLD,\n";
        cout << "ENERGY_THRESHOLD, ABSCHARGE_THRESHOLD) are ignored!\n";
    }

    double getMaxDist() const override {return this->maxDistance;}
};


/**
 * @brief Creates a list of all Density-Density Coulomb integrals up to a given electron-hole distance (maxDistance).
 *
 * Like the FixedListScheduler, this scheduler does not use any estimates, relation between integrals or
 * heuristics. It can only be used for density-density Coulomb integrals.
 */
class DensityFixedRScheduler : public FixedListScheduler
{
public:
    /**
     * @brief Creates a (fixed) list of all density-density Coulomb integrals up to an electron-hole distance of maxDist.
     *
     * Already known integrals (provided by solved_ints_) are skiped.
     *
     * @param batch_size
     * @param num_worker
     * @param maxDist   maximal electron-hole distance
     * @param id_c      list of all conduction Wannier function ids
     * @param id_v      list of all valence Wannier function ids
     * @param shells    shells of RD (e.g. obtained from createShells())
     * @param pos_c     conduction Wannier center
     * @param pos_v     valence Wannier center
     * @param unitcell  primitive unit cell of the material
     * @param solved_ints_ already known integrals that can be skiped.
     */
    DensityFixedRScheduler(int batch_size, int num_worker, double maxDist, vector<int> id_c, vector<int> id_v,
        vector<vector<vector<int>>> const& shells, vector<vector<double>> const& pos_c, vector<vector<double>> const& pos_v,
        vector<vector<double>> const& unitcell, map<vector<int>, Integral> const& solved_ints_)
        : FixedListScheduler(batch_size, num_worker)
    {
        // set all attributes
        this->setKnownIntegrals(solved_ints_);
        this->maxDistance = maxDist;

        if ((id_c.size() != pos_c.size()) || (id_v.size() != pos_v.size())) {
            throw runtime_error("ids and positons arrays need to have the same length!");
        }

        //cout << "solved_ints.size() = " << this->solved_ints.size() << endl;

        // create plan
        auto unitcell_T = transpose3x3(unitcell);
        list<Integral> plan_tmp;

        // go over all density-density integrals and create a plan
        for (size_t c=0; c<id_c.size(); c++) {
            for(size_t v=0; v<id_v.size(); v++) {
                for(auto shell: shells) {
                    for(auto RD : shell) {
                        // It is intended to go over each possible RD!!!
                        // Even if previous shells don't have contributions there can
                        // be situations where we would miss some integrals. To go over
                        // each RD without any shortcut is necessary for make it save.

                        vector<double> vec_shift = matVecMul3x3(unitcell_T, RD);
                        double dist = L2Norm(vector<double>{
                                pos_c[c][0]-pos_v[v][0]-vec_shift[0],
                                pos_c[c][1]-pos_v[v][1]-vec_shift[1],
                                pos_c[c][2]-pos_v[v][2]-vec_shift[2]
                            });

                        // cout << "dist: " << dist << " : " << c << "  " << v << "  " << RD[0] << "  " << RD[1] << "  " << RD[2] << endl;
                        // TODO: check if maxDistance goes beyond supercell

                        if (dist<=maxDistance + 1e-5) {  // +1e-5 to avoid numerical issues in wo-optics.x
                            auto index = vector<int>{id_c[c], id_c[c], id_v[v], id_v[v], RD[0],RD[1],RD[2], 0,0,0, 0,0,0};

                            // check if integral is known
                            auto itr = solved_ints.find(index);
                            if(itr == solved_ints.end())
                                plan_tmp.push_back(Integral(index));  // add integral to plan
                        }
                    }
                }
            }
        }

        // convert to vector
        this->plan = vector<Integral>{ begin(plan_tmp), end(plan_tmp) };

        // print plan
        // cout << "maxDist = " << maxDistance << endl;
        // cout << "Plan:\n";
        // printPlan(this->plan,false);
    }
};



/**
 * @brief Base class for all Coulomb integral schedulers.
 *
 * This class provides shared functionality for any Coulomb integral scheduler.
 * It is not intended for use with LFE integrals.
 *
 * The class stores all combinations of valence and conduction densities in a vector named tuples.
 * Each entry in tuples contains the Wannier indices and the shifts Rc (conduction) and Rv (valence) of the densities.
 * The tuples vector has dimensions (N, 2), where N is the number of combinations.
 * Child classes generate this list of combinations depending on the type of integral.
 *
 * A Coulomb integral is defined by a tuple (a pair of conduction and valence densities) combined with an RD vector,
 * which represents the shift between conduction and valence densities in terms of lattice vectors.
 * RD vectors are processed in shells. The cursors object stores the current shell and element within that shell
 * for each tuple.
 *
 * Each worker (MPI process) is assigned a specific tuple to work on, as recorded in worker2Tuple.
 * Different workers may work on the same tuple.
 *
 * This class provides functionalities to generate and update integrals but does not enforce thresholds,
 * estimates, or heuristics. It therefore processes all possible integrals based on the tuples.
 *
 * To apply estimates, thresholds, relationships between integrals, or other heuristics, implement the
 * check_relevance() and shouldTupleMarkedAsFinished() methods in the derived classes.
 */
class CoulombScheduler : public Scheduler
{
private:
    vector<vector<vector<int>>> shells;     //!< all possible RD vectors organized in shells
    vector<vector<int>> cursers;            //!< progress for each tuple:
                                            //!< tuble-id --> (shell, element of this shell)
                                            //!< this only stores which integrals have been assigned to workers. It does not contain
                                            //!< information about which results are already known. This is stored in the status array.

    vector<int> worker2Tuple;               //!< maps each worker to a tuple that he is working on.
                                            //!< This makes sure that workers can reuse the same density
                                            //!< over and over again.

    int nextTuple;                          //!< Candidate for the next tuple that might be picked if a worker finished
                                            //!< his current tuple. However, this tuple might already be finished! Please use
                                            //!< getNextTuple() to get a proper next tuple.

    vector<vector<vector<int>>> status;     //!< stores information about which integrals (results) are known.
                                            //!< (tuple_id, shell_id) --> (at least one integral is relevant in this shell, number of known integrals for this shell)
                                            //!< It stores if isRelevant() criterion is satisfied by any the integrals of the
                                            //!< shell or not and how many values of integrals for this shell are already known.

    vector<vector<vector<int>>> lastTasks;

    /**
     * @brief Returns next integral for a given tuple_id and increments the curser
     *
     * If no integral can be found (tuple is finished) then the curser is set to (-1,0)
     * and a NO-integral is returned.
     *
     * @param tuple_id     index of the tuples vector
     * @return Integral
     */
    Integral getNextIntegral(size_t tuple_id, int& shell, int& element_in_shell)
    {
        // check if tuple is finished
        if (isFinished(tuple_id)) {
            shell = -1;
            element_in_shell = -1;
            return Integral::CreateEmptyIntegral();
        }

        shell = cursers[tuple_id][0];
        element_in_shell = cursers[tuple_id][1];
        vector<int> RD = shells[shell][element_in_shell];

        Integral i(
            tuples[tuple_id][0],  // cDensity
            tuples[tuple_id][1],  // vDensity
            RD
        );

        // update curser
        cursers[tuple_id][1]++;
        if (cursers[tuple_id][1] >= int(shells[cursers[tuple_id][0]].size())) { // start new shell
            cursers[tuple_id][0]++;
            cursers[tuple_id][1] = 0;
            if (cursers[tuple_id][0] >= int(shells.size())) { setFinished(tuple_id); }  // finished all shells
        }

        return i;
    }

    /**
     * @brief Returns next tuple that can be assigned to a worker.
     *
     * It uses the private attribute nextTuple but also checks if the tuple is already
     * finished. If all tuples are finished it returns 0.
     * The value of nextTuple is updated (however it might point to a finished tuple).
     *
     * @return int   index of tuples[] vector
     */
    size_t getNextTuple()
    {
        size_t j=nextTuple;
        for (size_t i=0; i<tuples.size(); i++) {
            if (! this->isFinished(j)) {
                nextTuple = j+1;
                if (nextTuple>=int(tuples.size())) nextTuple=0;
                return j;
            }
            j++;
            if (j>=tuples.size()) j=0;
        }
        return 0;
    }


    void updateStatus(size_t tuple_id, size_t shell_id, bool relevant) {
        status[tuple_id][shell_id][0] = bool(status[tuple_id][shell_id][0]) || relevant;  // is at least one integral relevant in shell
        status[tuple_id][shell_id][1]++;  // number of known integrals

        // check if tuple can be marked as finished
        if (! status[tuple_id][shell_id][0]) {  // no known integral is relevant
            if (status[tuple_id][shell_id][1] == int(shells[shell_id].size())) {  // all integrals are known
                setFinished(tuple_id);
            }
        }
    }

    /**
     * @brief Marks a tuple as finished
     *
     * @param tuple_id
     */
    void setFinished(size_t tuple_id) {
        cursers[tuple_id][0] = -1;
        cursers[tuple_id][1] = 0;
    }


protected:
    vector<vector<Density_descr>> tuples;   //!< all combinations of valence and conduction densities that are of interest regardless of RD
                                            //!< dimensions: (N,2)

    /**
     * @brief Checks if a tuple is marked as finished
     *
     * It does not check if the tuple is finished due to some criterion. Just if it is
     * marked as finished.
     *
     * @param tuple_id  index of tuples[] vector
     * @return true
     * @return false
     */
    bool inline isFinished(size_t tuple_id) const { return cursers[tuple_id][0] < 0; }

    /**
     * @brief Checks if a single integral might be relevant for calculation.
     *
     * This might be performed by using estimates, heuristics or already solved
     * integrals that have a relation to the current integral.
     *
     * This method is called BEFORE an integral is calculated, i.e. in getNextTasks().
     *
     * This method needs to be implemented by derived classes depending on the
     * type of integral.
     *
     * @param a         the integral to test
     * @return true     integral is relevant and should be calculated
     * @return false    integral is not relevant and can be omitted
     */
    virtual bool estimateRelevance(Integral const& a) const {
        return true;
    }

    /**
     * @brief Checks if an already calculated integral is relevant.
     *
     * This method is called AFTER the integral has been calculated, i.e. in update().
     * The return of this method is used to update the status of the shell and helps
     * to determine if a shell can be marked as finished or not.
     *
     * This method needs to be implemented by derived classes depending on the
     * type of integral.
     *
     * @param a        the calculated integral
     * @return true
     * @return false
     */
    virtual bool isRelevant(Integral const& a) =0;


    /**
     * @brief Criterion to mark a tuple as finished.
     *
     * Depending on the values of the integrals in the last complete shell it might
     * be appropriate to mark the tuple as finished. E.g. when all values are below
     * a certain threshold.
     *
     * This method is called after Coulomb integrals are calculated, e.g. in update().
     *
     * This method needs to be implemented by derived classes depending on the
     * type of integral.
     *
     * @param tuple_id
     * @param lastCompleteShell
     * @return true     tuple should be marked as finished, no more calculations for this tuple
     * @return false
     */
    // virtual bool shouldTupleMarkedAsFinished(size_t tuple_id, size_t lastCompleteShell) {
    //     return false;
    // }

    // virtual bool shouldTupleMarkedAsFinished(size_t tuple_id)  // TODO: rename in something like updateTuple
    // {
    //     int lastCompleteShell = cursers[tuple_id][0] -1;
    //     if (lastCompleteShell<0) {  // no shell completed yet
    //         return false;
    //     }
    //     return shouldTupleMarkedAsFinished(tuple_id, lastCompleteShell);
    // }

    /**
     * @brief Initializes curser, status and worker2Tuple after the tuple is created.
     *
     * This needs to be called by after the tuples vector is created (and populated).
     *
     */
    void resetInternalState() {

        // tuples are already initialized at this point (usually by constructors of derived classes.)

        // initialize cursers
        cursers = vector<vector<int>>(tuples.size());
        for (size_t i=0; i<tuples.size(); i++) cursers[i] = vector<int>{0,0};

        // initialize status
        status = vector<vector<vector<int>>>(tuples.size());
        for (size_t i=0; i<tuples.size(); i++) {
            status[i] = vector<vector<int>>(shells.size());
            for (size_t j=0; j<shells.size(); j++) status[i][j] = vector<int>{0,0};
        }

        // initialize worker2Tuple
        this->worker2Tuple = vector<int>(num_worker);
        size_t j = 0;
        for (int i=0; i<num_worker; i++) {
            worker2Tuple[i] = j;
            j++;
            if (j>=tuples.size()) j=0;
        }
        nextTuple = j;
    }

public:

    CoulombScheduler(int batch_size, int num_worker, vector<vector<vector<int>>> const& shells_l)
    : Scheduler(batch_size, num_worker), shells{shells_l.size()}, nextTuple{0}
    {
        // copy all shells
        int i = 0;
        for (auto const& shell: shells_l) {
            shells[i] = vector<vector<int>>{ std::begin(shell), std::end(shell) };
            i++;
        }

        // tuples will be initialized by derived classes depending on the type of integral
        // resetInternalState() needs to be called by derived classes to initialize all data structures
        // at the end
    }

    vector<vector<Integral>> getNextTasks() override
    {
        vector<vector<Integral>> tasks(num_worker);
        lastTasks = vector<vector<vector<int>>>(num_worker);

        // initialize first
        for (int worker=0; worker< num_worker; worker++) {
            tasks[worker] = vector<Integral>(batch_size);
            lastTasks[worker] = vector<vector<int>>(batch_size);
            for( int i=0; i<batch_size; i++) lastTasks[worker][i] = vector<int>{0,0,0};
        }

        // assign integrals to every worker
        for (int worker=0; worker< num_worker; worker++) {
            for (int j=0; j<batch_size; j++) {

                Integral a;
                bool isConsidered = false;
                int shell, element_in_shell;
                do {
                    a = getNextIntegral(worker2Tuple[worker], shell, element_in_shell);  // get new integral for current tuple

                    // integral could be empty if e.g. tuple is finished
                    if (a.isEmpty()){
                        worker2Tuple[worker] = getNextTuple();  // go to next tuple and get new integral
                        a = getNextIntegral(worker2Tuple[worker], shell, element_in_shell);

                        // if integral for new tuple remains empty that means that there are no integrals
                        // anymore, i.e. this mpi-process is finished.
                        if (a.isEmpty()) break;
                    }

                    isConsidered = estimateRelevance(a);

                    // if integral is not considered it should be marked as not relevant in the status
                    if (!isConsidered) { updateStatus(worker2Tuple[worker], shell, false); }

                }while(!isConsidered);  // integrals that are not relevant (by some heuristics) should be skipped.

                // assign new integral
                tasks[worker][j] = a;
                lastTasks[worker][j][0] = worker2Tuple[worker];  // save indexes for later reference
                lastTasks[worker][j][1] = shell;
                lastTasks[worker][j][2] = element_in_shell;
            }
        }
        return tasks;
    }

    void update( vector<vector<Integral>> const& solved_integrals, bool verbose=true ) override
    {
        // add value to solved_ints map
        double maxValue = 0.0;
        for (size_t i=0; i<solved_integrals.size(); i++){
            for (size_t j=0; j<solved_integrals[i].size(); j++){
                Integral a = solved_integrals[i][j];
                if (a.isEmpty()) continue;  // skip No-Integrals

                solved_ints.insert({a.indexes, a});
                // the hermitian conjugate integral is only added in the very last step
                // otherwise we would have to deal with it in the Yukawa integral scheduler (FixedListScheduler)

                // update status
                int tuple_id = lastTasks[i][j][0];
                int shell_id = lastTasks[i][j][1];
                // int element_in_shell = lastTasks[i][j][2];
                updateStatus(tuple_id, shell_id, isRelevant(a));

                // some statistics and user output
                maxValue = max(maxValue, abs(a.value));  // max value of mini-batch (just for user output)
                //if (verbose) a.print();
                //if (verbose) conj_int.print();
            }
        }
        // cout << "# Max value of this batch: " << maxValue << endl;
        // double time_spend  =chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - startTime).count();
        // cout << "# Total number of skiped integrals (Indicator) = " << skiped_integrals << endl;
        // cout << "# Avg time per integral            = " << time_spend / solved_ints.size()  << " ms\n";

        // cout << "# Cursers:" << endl;
        // for (size_t i=0; i<cursers.size(); i++){
        //     if (cursers[i][0] < 0){
        //         cout << "#   " << i << " :  done." << endl;
        //     }else {
        //         cout << "#   " << i << " :  " << cursers[i][0] << ", " << cursers[i][1] << endl;
        //     }
        // }

        if (verbose) {
            cout << "\r"
                << "Solved integrals: " << solved_ints.size() << ";  "
                << "max value of batch: " << maxValue << "eV;  "
                << "next tuple: " << nextTuple << " of " << tuples.size()  // TODO: maybe print number of relevant integrals
                << flush;
        }
    }

    bool isFinished() const override {
        for (size_t i=0; i<tuples.size(); i++) {
            if (!isFinished(i)) return false;
        }
        return true;
    }

    // TODO: getAllCalculatedIntegrals() const
    //       getAllRelevantIntegrals() const { copy_if(... isRelevant )}

    bool saveState(string state_file) const override {
        ofstream file(state_file);
        if (!file.is_open()) {
            cout << "Cannot open file!" << endl;
            return false;
        }

        file << "CoulombScheduler\n";
        file << this->tuples.size() << endl;
        file << "Tuples:" << endl;
        for (auto& e: this->tuples) {
            file << e[0].toString() << "\t" << e[1].toString() << endl;
        }
        file << "Cursers:" << endl;
        for (auto& e: this->cursers) {
            file << e[0] << "\t" << e[1] << endl;
        }
        file << "Solved integrals:" << endl;
        file << this->solved_ints.size() << endl;
        // for (auto& e: this->solved_ints) {
        //     file << e.toString() << endl;
        // }

        file.close();
        return true;
    }
    bool loadState(string state_file, string coulomb_file) override {

        // open Coulomb file if exits
        if (!filesystem::exists(coulomb_file)) {
            cerr << "Cannot open file: " << coulomb_file << " because it does not exist!"<< endl;
            return false;
        }
        cout << "Read previous calculated integrals from " << coulomb_file << endl;
        auto f = CoulombFileReader(coulomb_file);
        if (! f.wasSuccessful()) {
            cerr << "Cannot read integrals from file.\n";
            return false;
        }
        const vector<Integral>& integrals = f.getElements();

        if (!filesystem::exists(state_file)) {
            cerr << "Cannot open file: " << state_file << " because it does not exist!"<< endl;
            return false;
        }

        ifstream file(state_file);
        if (!file.is_open()) {
            cerr << "Cannot open file: " << state_file << endl;
            return false;
        }

        string line;
        getline(file, line);
        if (line != "CoulombScheduler") {
            cerr << "Cannot load state of the scheduler of type " << line << " (should be CoulombScheduler)\n";
            file.close();
            return false;
        }
        getline(file, line);
        int Num = stoi(line); // get number of tuples
        cout << "Number of tuples: " << Num << endl;

        // read tubles
        getline(file, line);
        if (line != "Tuples:") {
            cerr << "Error while parsing tuples (cannot find keyword)\n";
            file.close();
            return false;
        }
        vector<vector<Density_descr>> new_tuples(Num);
        for (int i=0; i<Num; i++) {
            getline(file, line);
            vector<string> v = splitBySpacesAndTabs(line);//split(line, string("\t"));
            if (v.size() != 10) {
                cerr << "Error while parsing tuples\n";
                file.close();
                return false;
            }
            new_tuples[i] = vector<Density_descr>{
                Density_descr(stoi(v[0]), stoi(v[1]), vector<int>{stoi(v[2]),stoi(v[3]),stoi(v[4])}),
                Density_descr(stoi(v[5]), stoi(v[6]), vector<int>{stoi(v[7]),stoi(v[8]),stoi(v[9])})
            };
        }

        // read cursers
        getline(file, line);
        if (line != "Cursers:") {
            cerr << "Error while parsing cursers (cannot find keyword)\n";
            file.close();
            return false;
        }
        vector<vector<int>> new_cursers(Num);
        for (int i=0; i<Num; i++) {
            getline(file, line);
            vector<string> v = splitBySpacesAndTabs(line); //split(line, string("\t"));
            if (v.size() != 2) {
                cerr << "Error while parsing cursers\n";
                file.close();
                return false;
            }
            new_cursers[i] = vector<int>{stoi(v[0]), stoi(v[1])};
        }

        // read solved ints
        getline(file, line);
        if (line != "Solved integrals:") {
            cerr << "Error while parsing solved integrals (cannot find keyword)\n";
            file.close();
            return false;
        }
        getline(file, line);
        int Num_ints = stoi(line);
        cout << "Number of solved ints: " << Num_ints << endl;
        file.close();

        if (Num_ints != int(integrals.size())) {
            cerr << "Cannot load state of the scheduler: Expected "<< Num_ints << " calculated integrals but found " << integrals.size() << ".\n";
            return false;
        }

        // actually set state
        this->tuples= new_tuples;
        this->setKnownIntegrals(integrals);
        this->resetInternalState();
        this->cursers= new_cursers;

        return true;
    }

};


class DensityDensityScheduler : public CoulombScheduler  //!< density-density integrals. maxDistance is automatically obtained
{
private:
    vector<vector<double>> pos_c, pos_v;    //!< positions of the charge centers
    vector<vector<double>> unitcell;        //!< unitcell of each Wannier function
    const CoulombPotential pot;             //!< potential to evaluate monopole interaction with
    const double DIFF_MONO_THRESHOLD;       //!< threshold for the maximal absolute relative difference between actual value and monopole contribution
    double maxDistance = 0;                 //!< maximal distance that needs to be calculated
    map<int,int> cWannierId_to_id;          //!< maps the conduction Wannier function id to the position in pos_c (Wannier ids can be arbitrary)
    map<int,int> vWannierId_to_id;          //!< maps the valence Wannier function id to the position in pos_c (Wannier ids can be arbitrary)
    const double MAX_ALLOWED_DISTANCE;

    vector<double> monopole_estimate(int id_c, int id_v, vector<int> const& RD) const {
        vector<double> vec_shift = matVecMul3x3(transpose3x3(unitcell), RD);
        vector<double> values(2);

        auto itr = cWannierId_to_id.find(id_c);
        if (itr == cWannierId_to_id.end()) throw runtime_error("In monopole_estimate of scheduler: cannot find position of given conduction Wannier function id.");
        int c = itr->second;

        itr = vWannierId_to_id.find(id_v);
        if (itr == vWannierId_to_id.end()) throw runtime_error("In monopole_estimate of scheduler: cannot find position of given valuence Wannier function id.");
        int v = itr->second;

        // monopole contribution in eV
        values[0] = pot.realCart(
            pos_c[c][0]-pos_v[v][0]-vec_shift[0],
            pos_c[c][1]-pos_v[v][1]-vec_shift[1],
            pos_c[c][2]-pos_v[v][2]-vec_shift[2]);

        // distance in Ang.
        values[1] = L2Norm(vector<double>{
            pos_c[c][0]-pos_v[v][0]-vec_shift[0],
            pos_c[c][1]-pos_v[v][1]-vec_shift[1],
            pos_c[c][2]-pos_v[v][2]-vec_shift[2]
        });

        return values;
    }

protected:

    bool isRelevant(Integral const& a) override
    {
        if(a.isEmpty()) return false;

        vector<double> mono = monopole_estimate(a.getConDensity().id1, a.getValDensity().id1, a.getRD());
        double rel_diff = abs((a.value - mono[0]) / a.value);

        if (rel_diff >= DIFF_MONO_THRESHOLD) {
            // update maximal electron-hole distance
            maxDistance = max(maxDistance, mono[1]);
            return true;
        }

        return false;
    }


public:
    DensityDensityScheduler(int batch_size, int num_worker, vector<int> const& id_c, vector<int> const& id_v,
        vector<vector<vector<int>>> const& shells_l, vector<vector<double>> const& pos_c_, vector<vector<double>> const& pos_v_,
        vector<vector<double>> const& unitcell_, const double DIFF_MONO_THRESHOLD_, const double MAX_ALLOWED_DISTANCE_,
        string const& restart_file="", string const& coulomb_file="")
        : CoulombScheduler(batch_size, num_worker, shells_l), pos_c(pos_c_), pos_v(pos_v_),
         unitcell(unitcell_), DIFF_MONO_THRESHOLD(DIFF_MONO_THRESHOLD_), MAX_ALLOWED_DISTANCE(MAX_ALLOWED_DISTANCE_)
    {
        if ((id_c.size() != pos_c.size()) || (id_v.size() != pos_v.size())) {
            throw runtime_error("ids and positons arrays need to have the same length!");
        }
        for (size_t i=0; i<id_c.size(); i++) {
            cWannierId_to_id.insert({id_c[i], i});
        }
        for (size_t i=0; i<id_v.size(); i++) {
            vWannierId_to_id.insert({id_v[i], i});
        }

        // load from previous calculation if possible
        if (tryRestartScheduler(restart_file, coulomb_file)) return;

        // create vector of all possible combinations (c,v)
        tuples = vector<vector<Density_descr>>(id_c.size()*id_v.size());
        int i = 0;
        for (size_t c=0; c<id_c.size(); c++) {
            for (size_t v=0; v<id_v.size(); v++) {
                tuples[i] = vector<Density_descr>{
                    Density_descr(id_c[c],id_c[c],vector<int>{0,0,0}),
                    Density_descr(id_v[v],id_v[v],vector<int>{0,0,0})
                };
                i++;
            }
        }

        // initialize curser and worker2Tuple
        this->resetInternalState();
    }

    double getMaxDist() const override
    {
        if (maxDistance > MAX_ALLOWED_DISTANCE) {
            cerr << "\n\n-------------------------------------------------------------\n\n";
            cerr << "▗▖ ▗▖ ▗▄▖ ▗▄▄▖ ▗▖  ▗▖▗▄▄▄▖▗▖  ▗▖ ▗▄▄▖\n";
            cerr << "▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌▐▛▚▖▐▌  █  ▐▛▚▖▐▌▐▌   \n";
            cerr << "▐▌ ▐▌▐▛▀▜▌▐▛▀▚▖▐▌ ▝▜▌  █  ▐▌ ▝▜▌▐▌▝▜▌\n";
            cerr << "▐▙█▟▌▐▌ ▐▌▐▌ ▐▌▐▌  ▐▌▗▄█▄▖▐▌  ▐▌▝▚▄▞▘\n";
            cerr << "\n";
            cerr << "Maximum electron-hole distance reached before meeting the\n";
            cerr << "MONOPOLE_RELATIVE_ERROR_THRESHOLD criterion.\n";
            cerr << "Density-density integrals may be incomplete.\n";
            cerr << "To resolve this, increase the supercell size (SUPERCELL_DIM_* in the config file)\n";
            cerr << "or adjust the DIFF_MONO_THRESHOLD value.\n\n";
            cerr << "MAX_ALLOWED_DISTANCE = "<< MAX_ALLOWED_DISTANCE << " Angström" << endl;
            cerr << "current distance     = " << maxDistance << " Angström" << endl;
            cerr << "\n\n-------------------------------------------------------------\n";
        }

        return min(maxDistance, MAX_ALLOWED_DISTANCE);
    }

    void printStartupSummary() const override
    {
        cout << "- For all (c, v), calculate integrals with increasing electron-hole distance (RD)\n";
        cout << "  until the relative difference to the monopole-monopole approximation is smaller\n";
        cout << "  than DIFF_MONO_THRESHOLD = "<< DIFF_MONO_THRESHOLD*100. <<" %.\n";
        cout << "- RD is varied in shells, and the criterion is checked only for complete shells.\n";
        cout << "- Other threshold parameters (DISTANCE_THRESHOLD, ENERGY_THRESHOLD,\n";
        cout << "  ABSCHARGE_THRESHOLD) are ignored.\n";
    }
};


/**
 * @brief Concatenation of the DensityDensityScheduler and the DensityFixedRScheduler
 *
 * This scheduler calculates density-density Coulomb integrals only!
 *
 * The maximal electron-hole distance obtained from DensityDensityScheduler (first phase)
 * is used as an input for the DensityFixedRScheduler (second phase). This makes sure that
 * the COULOMB file contains all integrals up to a certain distance and interactions larger
 * than this electron-hole distance can be approximated by monopole-monopole interaction.
 *
 * With the DIFF_MONO_THRESHOLD parameter the user can control how much the relative error
 * compared to monopole-monopole approximation is allowed to be. This determines the maximal
 * electron-hole distance where integrals have to calculated in full detail. The second phase
 * is needed in order to make sure that no integral with a smaller distance is forgotten even
 * though it may satisfy relative error condition.
 *
 *
 */
class CombinedDensityDensityFixedRScheduler : public Scheduler
{
private:
    vector<int> id_c, id_v;                 //!< conduction and valence Wannier function ids
    vector<vector<double>> pos_c, pos_v;    //!< positions of the charge centers
    vector<vector<vector<int>>> shells;     //!< shells for RD
    vector<vector<double>> unitcell;        //!< unitcell of each Wannier function

    unique_ptr<DensityDensityScheduler> firstPhase{};  //!< scheduler for the 1st phase
    unique_ptr<DensityFixedRScheduler> secondPhase{};  //!< scheduler for the 2nd phase
    Scheduler* currentPhase{};              //!< scheduler of the current phase, should be either firstPhase or secondPhase

public:
    CombinedDensityDensityFixedRScheduler(int batch_size, int num_worker, vector<int> id_c_, vector<int> id_v_,
        vector<vector<vector<int>>> const& shells_l, vector<vector<double>>const& pos_c_, vector<vector<double>>const& pos_v_,
        vector<vector<double>>const& unitcell_, const double DIFF_MONO_THRESHOLD_, const double MAX_ALLOWED_DISTANCE_,
        string const& restart_file="", string const& coulomb_file="")
        : Scheduler(batch_size, num_worker), id_c(id_c_), id_v(id_v_), pos_c(pos_c_), pos_v(pos_v_), shells{shells_l.size()},
          unitcell(unitcell_)
    {
        int i = 0;
        for (auto const& shell: shells_l) {
            shells[i] = vector<vector<int>>{ std::begin(shell), std::end(shell) };
            i++;
        }

        // create scheduler for first phase
        firstPhase = make_unique<DensityDensityScheduler>(batch_size, num_worker,id_c_, id_v_, shells_l, pos_c_,pos_v_,
            unitcell_, DIFF_MONO_THRESHOLD_, MAX_ALLOWED_DISTANCE_);
        secondPhase = nullptr;
        currentPhase = firstPhase.get();

        this->tryRestartScheduler(restart_file, coulomb_file);
    }

    vector<vector<Integral>> getNextTasks() override { return currentPhase->getNextTasks(); }

    void update( vector<vector<Integral>> const& solved_integrals, bool verbose=true ) override
    {
        currentPhase->update(solved_integrals, verbose);

        // check if current phase is finished and second phase has not started
        if (! secondPhase) {
            if (firstPhase->isFinished()) {
                double maxDist = firstPhase->getMaxDist();
                cout << "\nFirst phase finished.\n";
                cout << "--> electron-hole distances larger " << maxDist << "A will be approximated with monopole-monopole approx." << endl;
                cout << "Start second phase of the scheduler...\n";
                secondPhase = make_unique<DensityFixedRScheduler>(batch_size, num_worker, maxDist, id_c, id_v,
                                    shells, pos_c, pos_v, unitcell, firstPhase->getKnownIntegrals_asmap());
                currentPhase = secondPhase.get();

                // free memory
                firstPhase.reset();
            }
        }
    }

    bool isFinished() const override {
        // cannot be finished if second phase has not even started
        if (secondPhase == nullptr) return false;

        return secondPhase->isFinished();
    }

    const list<Integral> getKnownIntegrals() const override { return currentPhase->getKnownIntegrals(); }
    const map<vector<int>, Integral> getKnownIntegrals_asmap() const override {return currentPhase->getKnownIntegrals_asmap(); }
    void setKnownIntegrals(vector<Integral> const& knownInts) override { currentPhase->setKnownIntegrals(knownInts); }
    void setKnownIntegrals(map<vector<int>, Integral> const& knownInts) override { currentPhase->setKnownIntegrals(knownInts); }
    double getMaxDist() const override {return currentPhase->getMaxDist();}

    void printStartupSummary() const override
    {
        cout << "You are calculating density-density (2-center) integrals with only classical\n";
        cout << "charge densities (v1=v2=v, c1=c2=c, Rv=Rc=0). These integrals can be specified\n";
        cout << "by (c, v, RD). The calculation is performed in two phases:\n";
        cout << "\nPhase 1:\n";
        if (! firstPhase) {
            cout << "(skipped)\n";
        }else{
            firstPhase->printStartupSummary();
        }
        cout << "\nPhase 2:\n";
        cout << "- The scheduler determines the maximum el-hole distance needed for the\n";  // TODO do something like DensityFixedRScheduler::printStartupSummary()
        cout << "  monopole-monopole approximation to be appropriate.\n";
        cout << "- All integrals with smaller el-hole distances that were not calculated in Phase 1\n";
        cout << "  will be computed.\n";
    }

    bool saveState(string restart_file) const override { return currentPhase->saveState(restart_file); }

    bool loadState(string restart_file, string coulomb_file) override {

        if (!filesystem::exists(restart_file)) {
            cerr << "Cannot open file: " << restart_file << " because it does not exist!"<< endl;
            return false;
        }

        ifstream file(restart_file);
        if (!file.is_open()) {
            cerr << "Cannot open file: " << restart_file << endl;
            file.close();
            return false;
        }

        string line;
        getline(file, line);
        file.close();

        cout << "Found restart file for " << line << "\n";

        if (line == "FixedListScheduler") {
            cout << "Try to load second phase of the scheduler.\n";

            secondPhase = make_unique<DensityFixedRScheduler>(batch_size, num_worker, firstPhase->getMaxDist(), id_c, id_v,
                                    shells, pos_c, pos_v, unitcell, firstPhase->getKnownIntegrals_asmap());
            currentPhase = secondPhase.get();

            // free memory
            firstPhase.reset();
            return secondPhase->loadState(restart_file, coulomb_file);
        }else if (line == "CoulombScheduler") {
            cout << "Try to load first phase of the scheduler.\n";
            return firstPhase->loadState(restart_file, coulomb_file);
        }else {
            cerr << "Cannot load state form restart file.\n";
            return false;
        }
    }

};


/**
 * @brief Scheduler to calculate overlap-density and overlap-overlap integrals
 *
 * This scheduler uses estimates and relations between Coulomb integrals.
 *
 */
class OverlapDensityScheduler : public CoulombScheduler
{

private:
    map<Density_descr,Density_indicator> const& vIndicators;
    map<Density_descr,Density_indicator> const& cIndicators;

    vector<vector<double>> unitcell;        //!< unitcell of each Wannier function

    const double DISTANCE_THRESHOLD;        //!< maximal electron-hole distance in Anstrom (RD)
    const double ENERGY_THRESHOLD;          //!< minimal energy (in eV) that is relevant

    const string mode;

    double getDistance(Integral const& a) const {
        // get Indicator-parameters for densities
        auto val_itr = vIndicators.find(a.getValDensity());
        if (val_itr == vIndicators.end()) {
            cout << "Density: " << a.getValDensity().toString() << endl;
            throw runtime_error("[FATAL ERROR] Cannot find Indicator-parameters for valence density.");
        }

        auto con_itr = cIndicators.find(a.getConDensity());
        if (con_itr == cIndicators.end()) {
            cout << "Density: " << a.getConDensity().toString() << endl;
            throw runtime_error("[FATAL ERROR] Cannot find Indicator-parameters for conduction density.");
        }

        // get values and distances
        auto RD = a.getRD();
        vector<double> pos_c = con_itr->second.center();
        vector<double> pos_v = val_itr->second.center();
        vector<double> vec_shift = matVecMul3x3(transpose3x3(unitcell), RD);

        // distance in Ang.
        double distance = L2Norm(vector<double>{
            pos_c[0]-pos_v[0]-vec_shift[0],
            pos_c[1]-pos_v[1]-vec_shift[1],
            pos_c[2]-pos_v[2]-vec_shift[2]
        });

        return distance;
    }

protected:


    bool isRelevant(Integral const& a) override
    {
        if (a.isEmpty()) return false;
        return a.value > ENERGY_THRESHOLD;
    }

    bool estimateRelevance(Integral const& a) const override
    {
        if (a.isEmpty()) return false;
        if (getDistance(a) > DISTANCE_THRESHOLD) return false;

        // Indicator estimate in eV
        // double mbie0 = abs_mono_c * abs_mono_v * pot->real(abs(distance - extend_c - extend_v));

        // if (mbie0 < ENERGY_THRESHOLD) {
        //     skiped_integrals++;
        //     // cout << "skipt the follwing integral (due to Indicator value):\n";
        //     // a.print();
        //     // cout << "distance = " << distance << endl;
        //     // cout << "distance - extend_c - extend_v = " << distance - extend_c - extend_v << endl;
        //     // cout << "mbie = " << mbie0 << endl;
        //     return false;
        // }

        // TODO use Cauchy-Schwarz inequality
        return true;
    }

public:
    OverlapDensityScheduler(
        int batch_size, int num_worker,
        map<Density_descr,Density_indicator> const& cIndicators_,
        map<Density_descr,Density_indicator> const& vIndicators_,
        vector<vector<vector<int>>> const& shells_l, vector<vector<double>> const& unitcell_,
        const double ABSCHARGE_THRESHOLD, const double DISTANCE_THRESHOLD_, const double ENERGY_THRESHOLD_,
        const string& mode_="OverlapDensity", string  const& restart_file="", string const& coulomb_file=""
    )
    : CoulombScheduler(batch_size, num_worker, shells_l), vIndicators(vIndicators_), cIndicators(cIndicators_),
      unitcell(unitcell_), DISTANCE_THRESHOLD(DISTANCE_THRESHOLD_), ENERGY_THRESHOLD(ENERGY_THRESHOLD_), mode(mode_) {

        if (vIndicators.size() == 0) throw runtime_error("vIndicators is empty! This should not happen!");
        if (cIndicators.size() == 0) throw runtime_error("cIndicators is empty! This should not happen!");

        tuples = vector<vector<Density_descr>>();  // clear tuples

        // load from previous calculation if possible
        if (tryRestartScheduler(restart_file, coulomb_file)) return;

        // create list of all possible combinations of densities (ignoring RD)
        if (mode == "OverlapDensity") {
            for(auto const& cmap: cIndicators) {
                if (cmap.second.absCharge < ABSCHARGE_THRESHOLD) continue;

                for(auto const& vmap: vIndicators) {

                    if (vmap.second.absCharge < ABSCHARGE_THRESHOLD) continue;

                    // we only want Density-Overlap integrals in this mode
                    if (cmap.first.isOverlapDensity() && vmap.first.isOverlapDensity()) continue;  // skip overlap-overlap
                    if (cmap.first.isClassicalDensity() && vmap.first.isClassicalDensity() ) continue;  // skip density-density

                    // hermiticity (only calculate matrix elements that are on the upper or lower triangular part of the Hamiltonian)
                    vector<int> Rc = cmap.first.R;
                    vector<int> Rv = vmap.first.R;
                    vector<int> id1 = vector<int>{ cmap.first.id1, vmap.first.id1, 0,0,0 };
                    vector<int> id2 = vector<int>{ cmap.first.id2, vmap.first.id2, Rc[0]-Rv[0], Rc[1]-Rv[1], Rc[2]-Rv[2] };

                    // we want to use the Hermitian property of the Hamiltonian, which means that we only
                    // need to calculate half of the elements. For this we can use the lexicographic ordering
                    // as implemented for vectors in C++
                    if (id1 < id2) continue;

                    // tuple satisfies all criteria and is added to the list
                    tuples.push_back(vector<Density_descr>{cmap.first, vmap.first});
                }
            }
        }else if (mode=="OverlapOverlap") {
            for(auto const& cmap: cIndicators) {
                if (cmap.second.absCharge < ABSCHARGE_THRESHOLD) continue;

                for(auto const& vmap: vIndicators) {
                    if (vmap.second.absCharge < ABSCHARGE_THRESHOLD) continue;

                    // we only want Overlap-Overlap integrals in this mode
                    if (cmap.first.isOverlapDensity() && vmap.first.isOverlapDensity()) {

                        // hermiticity (only calculate matrix elements that are on the upper or lower triangular part of the Hamiltonian)
                        vector<int> Rc = cmap.first.R;
                        vector<int> Rv = vmap.first.R;
                        vector<int> id1 = vector<int>{ cmap.first.id1, vmap.first.id1, 0,0,0 };
                        vector<int> id2 = vector<int>{ cmap.first.id2, vmap.first.id2, Rc[0]-Rv[0], Rc[1]-Rv[1], Rc[2]-Rv[2] };

                        // we want to use the Hermitian property of the Hamiltonian, which means that we only
                        // need to calculate half of the elements. For this we can use the lexicographic ordering
                        // as implemented for vectors in C++
                        if (id1 < id2) continue;

                        // tuple satisfies all criteria and is added to the list
                        tuples.push_back(vector<Density_descr>{cmap.first, vmap.first});
                    }
                }
            }
        }else {
            throw runtime_error("Mode of OverlapDensityScheduler should be either 'overlap-density' or 'overlap-overlap'.");
        }

        cout << "tuples.size() = " << tuples.size() << endl;

        // initialize worker2Tuple and cursers
        this->resetInternalState();
    }

    void printStartupSummary() const override {
        if (mode=="OverlapDensity") {
            cout << "Calculating 3-center integrals (overlap density and classical density)\n";
            cout << "while ignoring 2-center (density-density) integrals.\n";
        }else{
            cout << "Calculating 4-center integrals (overlap density and overlap density)\n";
            cout << "while ignoring any integrals that contain classical densites (2-center\n";
            cout << "and 3-center).\n";
        }

        cout << "\n";
        cout << " - Only densities with absolute charge > ABSCHARGE_THRESHOLD are used.\n";
        cout << " - Only integrals with electron-hole distance < DISTANCE_THRESHOLD are considered.\n";
        cout << " - Only integrals potentially > ENERGY_THRESHOLD are selected using estimates\n";
        cout << "   and upper bounds.\n";
        cout << " - DIFF_MONO_THRESHOLD is ignored.\n";

        // cout << "\nShells:\n";
        // for (size_t i=0; i<shells.size(); i++) {
        //     for (size_t j=0; j<shells[i].size(); j++) {
        //         cout << i << " : r_"<< j << " = (" << shells[i][j][0] << ", "<< shells[i][j][1]<< ", "<< shells[i][j][2]<< ")\n ";
        //     }
        // }
    }
};

// TODO: use integration test with wannier optics --> make sure that we use the correct notation/minus signs and so on!!!
class LocalFieldEffectsScheduler : public Scheduler
{
    vector<Density_descr> tuples;         //!< all non-zero overlap densities (c1,v1,S1)
                                    //!< Please note that these densities are defined with a negative shift
                                    //!< rho(x) = w_c1,0(x) w_v1,-S1(x)
                                    //!< for the integral object we have S1=Rc (the minus is intrinsically in the solver routine)
                                    //!< with this we make sure that the output files have columns of S1 instead of -S1

    vector<int> worker2Tuple;       //!< maps each worker to a tuple that he is working on.
                                    //!< This makes sure that workers can reuse the same density
                                    //!< over and over again.

    size_t nextTuple;               //!< Candidate for the next tuple that might be picked if a worker finished
                                    //!< his current tuple. However, this tuple might already be finished! Please use
                                    //!< getNextTuple() to get a proper next tuple.

    vector<int> cursers;            //!< next tuple2 for given tuple1 (size = tuples.size)

    chrono::steady_clock::time_point startTime;  //!< to measure performance

    /**
     * @brief Returns next integral for a given tuple1 (Density) and increments the curser
     *
     * If no integral can be found (tuple is finished) then the curser is set to (-1,0)
     * and a NO-integral is returned.
     *
     * @param tuple     tuple of Wannier indexes (conduction WF, valence WF)
     * @return Integral
     */
    Integral setNextElement(int tuple) {

        while (cursers[tuple]>=0) {

            // we want to use the Hermitian property of the Hamiltonian, which means that we only
            // need to calculate half of the elements. For this we can use the lexicographic ordering
            // as implemented for vectors in C++
            if (tuples[cursers[tuple]] <= tuples[tuple]) {

                // create integral data structure
                Integral i(
                    tuples[tuple].id1,   // c1
                    tuples[cursers[tuple]].id1,  // c2
                    tuples[tuple].id2,   // v1
                    tuples[cursers[tuple]].id2,  // v2
                    vector<int>{0,0,0},  // RD=0
                    tuples[tuple].R,     // Rc=S1
                    tuples[cursers[tuple]].R     // Rv=S2
                );

                // update curser
                cursers[tuple]++;
                if (cursers[tuple] >= int(tuples.size())) {
                    cursers[tuple] = -1;  // mark tuple as finished
                }

                return i;
            }

            // update curser
            cursers[tuple]++;
            if (cursers[tuple] >= int(tuples.size())) {
                cursers[tuple] = -1;  // mark tuple as finished
            }

        }

        return Integral::CreateEmptyIntegral();  // tuple is marked as finished --> return NO-integral

    }

    /**
     * @brief Returns next tuple that can be assigned to a worker.
     *
     * It uses the private attribute nextTuple but also checks if the tuple is already
     * finished. If all tuples are finished it returns 0.
     * The value of nextTuple is updated (however it might point to a finished tuple).
     *
     * @return int   index of tuples[] vector
     */
    int getNextTuple()
    {
        size_t j=nextTuple;
        for (size_t i=0; i<tuples.size(); i++) {

            if (cursers[j] >= 0) {  // checks if tuple j is not finished

                // update nextTuple
                nextTuple = j+1;
                if (nextTuple>=tuples.size()) nextTuple=0;

                // return not finished tuple
                return j;
            }

            j++;
            if (j>=tuples.size()) j=0;
        }

        // all tuples are finished
        return 0;
    }

    void initializeWorker() {
        // initialize worker2Tuple
        worker2Tuple = vector<int>(num_worker);
        size_t j = 0;
        for (int i=0; i<num_worker; i++) {
            worker2Tuple[i] = j;
            j++;
            if (j>=tuples.size()) j=0;
        }
        nextTuple = j;
    }

public:

    LocalFieldEffectsScheduler(
        int batch_size, int num_worker,
        map<Density_descr,Density_indicator> const& lfe_Indicators,
        const double ABSCHARGE_THRESHOLD, string const& restart_file="", string const& coulomb_file=""
    )
    : Scheduler(batch_size, num_worker) {

        if (lfe_Indicators.size() == 0) throw runtime_error("lfe_Indicators is empty! This should not happen!");

        startTime = chrono::steady_clock::now();

        // get all important tuples from lfe_Indicators
        for(auto const& imap: lfe_Indicators) {
            if (imap.second.absCharge >= ABSCHARGE_THRESHOLD) {

                // convert (c1,v1,R) --> (c1,v1,S1) where S1 = -R
                // we need to do this because the LocalFieldEffectsSolver assumes that Rc=S1, Rv=S2, RD=0
                // with the corresponding overlap densities rho = w_c1,0(x) w_v1,-S1(x)
                // but the Indicator routines use the usual convention for overlap densities rho = w_c1,0(x) w_v1,R(x)
                Density_descr convertedDensity(imap.first);
                convertedDensity.R[0] = - convertedDensity.R[0];
                convertedDensity.R[1] = - convertedDensity.R[1];
                convertedDensity.R[2] = - convertedDensity.R[2];
                tuples.push_back(convertedDensity);
            }
        }

        // load from previous calculation if possible
        if (tryRestartScheduler(restart_file, coulomb_file)) return;

        if (tuples.size() == 0) {
            runtime_error("No overlap densities satisfy the absolute monopole criterion. Please check your parameters!");
        }

        // initialize cursers
        cursers = vector<int>(tuples.size());
        for (size_t i=0; i<tuples.size(); i++) {
            cursers[i] = 0;
        }

        cout << "tuples.size() = " << tuples.size() << endl;

        // initialize worker2tuple array
        this->initializeWorker();
    }

    vector<vector<Integral>> getNextTasks() override {

        vector<vector<Integral>> tasks(num_worker);

        // initialize first
        for (int worker=0; worker< num_worker; worker++) {
            tasks[worker] = vector<Integral>(batch_size);
        }

        // assign integrals
        for (int worker=0; worker< num_worker; worker++) {
            for (int j=0; j<batch_size; j++) {
                Integral a = setNextElement(worker2Tuple[worker]);
                if (a.isEmpty()) {  // no more integrals for this tuple --> start next tuple
                    worker2Tuple[worker] = getNextTuple();
                    a = setNextElement(worker2Tuple[worker]);
                }
                tasks[worker][j] = a;
            }
        }

        return tasks;
    }

    void update( vector<vector<Integral>> const& solved_integrals, bool verbose=true ) override {

        // add value to solved_ints map
        double maxValue = 0.0;
        for (size_t i=0; i<solved_integrals.size(); i++){
            for (size_t j=0; j<solved_integrals[i].size(); j++){
                Integral a = solved_integrals[i][j];
                if (a.isEmpty()) continue;  // skip No-Integrals

                solved_ints.insert({a.indexes, a});

                // add Hermitian conjugate to list
                Density_descr tuple1(a.indexes[0], a.indexes[2], vector<int>{a.indexes[7],a.indexes[8],a.indexes[9]});
                Density_descr tuple2(a.indexes[1], a.indexes[3], vector<int>{a.indexes[10],a.indexes[11],a.indexes[12]});
                if (tuple1 != tuple2){  // check if it is a diagonal element of the H^LFE

                    // create hermitian conjugated integral
                    Integral conj_int(
                        tuple2.id1,   // c1
                        tuple1.id1,  // c2
                        tuple2.id2,   // v1
                        tuple1.id2,  // v2
                        vector<int>{0,0,0},  // RD=0
                        tuple2.R,     // S1
                        tuple1.R     // S2
                    );
                    conj_int.value = a.value;
                    solved_ints.insert({conj_int.indexes, conj_int});
                    // if (verbose) conj_int.print();
                }

                maxValue = max(maxValue, abs(a.value));  // max value of mini-batch (just for user output)
                // if (verbose) a.print();
            }
        }
        // if (verbose) {
            // cout << "# Max value of this batch: " << maxValue << endl;
            // double time_spend  =chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - startTime).count();
            // cout << "# Total number of integrals = " << tuples.size()*tuples.size() << endl;
            // cout << "# Total number of solved integrals = " << solved_ints.size() << " ( "<< solved_ints.size()/double(tuples.size()*tuples.size())*100 <<" % )" << "\n";
            // cout << "# Avg time per integral            = " << time_spend / solved_ints.size()  << " ms\n";


            // cout << "# Cursers:" << endl;
            // for (size_t i=0; i<cursers.size(); i++){
            //     if (cursers[i] == -1){
            //         // cout << "#   " << i << " :  done." << endl;
            //     }else if (cursers[i] > 0) {
            //         cout << "#   " << i << " :  " << cursers[i] << " / " << tuples.size() << endl;
            //     }
            // }
            // cout << "# nextTuple = " << nextTuple << endl;
        // }

        if (verbose) {
            cout << "\r"
                << "Solved integrals: " << solved_ints.size() << ";  "
                << "max value of batch: " << maxValue << "eV;  "
                << "next tuple: " << nextTuple << " of " << tuples.size()
                << flush;
        }
    }

    bool isFinished() const override {
        for (size_t i=0; i<tuples.size(); i++) {
            if (cursers[i] >= 0) return false;
        }

        return true;
    }

    bool saveState(string state_file) const override {
        ofstream file(state_file);
        if (!file.is_open()) {
            cout << "Cannot open file!" << endl;
            return false;
        }

        file << "LocalFieldEffectsScheduler\n";
        file << this->tuples.size() << endl;
        file << "Tuples:" << endl;
        for (auto& e: this->tuples) {
            file << e.toString() << endl;
        }
        file << "Cursers:" << endl;
        for (auto& e: this->cursers) {
            file << e << endl;
        }
        file << "Solved integrals:" << endl;
        file << this->solved_ints.size() << endl;
        // for (auto& e: this->solved_ints) {
        //     file << e.toString() << endl;
        // }

        file.close();
        return true;
    }
    bool loadState(string state_file, string coulomb_file) override {

        // open Coulomb file if exits
        if (!filesystem::exists(coulomb_file)) {
            cerr << "Cannot open file: " << coulomb_file << " because it does not exist!"<< endl;
            return false;
        }
        cout << "Read previous calculated integrals from " << coulomb_file << endl;
        auto f = CoulombFileReader(coulomb_file);
        if (! f.wasSuccessful()) {
            cerr << "Cannot read integrals from file.\n";
            return false;
        }
        const vector<Integral>& integrals = f.getElements();

        if (!filesystem::exists(state_file)) {
            cerr << "Cannot open file: " << state_file << " because it does not exist!"<< endl;
            return false;
        }

        ifstream file(state_file);
        if (!file.is_open()) {
            cerr << "Cannot open file: " << state_file << endl;
            return false;
        }

        string line;
        getline(file, line);
        if (line != "LocalFieldEffectsScheduler") {
            cerr << "Cannot load state of the scheduler of type " << line << " (should be LocalFieldEffectsScheduler)\n";
            file.close();
            return false;
        }
        getline(file, line);
        int Num = stoi(line); // get number of tuples
        cout << "Number of tuples: " << Num << endl;

        // read tuples
        getline(file, line);
        if (line != "Tuples:") {
            cerr << "Error while parsing tuples (cannot find keyword)\n";
            file.close();
            return false;
        }
        vector<Density_descr> new_tuples(Num);
        for (int i=0; i<Num; i++) {
            getline(file, line);
            vector<string> v = splitBySpacesAndTabs(line); //split(line, string("\t"));
            if (v.size() != 5) {
                cerr << "Error while parsing tuples\n";
                file.close();
                return false;
            }
            new_tuples[i] = Density_descr(stoi(v[0]), stoi(v[1]), vector<int>{stoi(v[2]),stoi(v[3]),stoi(v[4])});
        }

        // read cursers
        getline(file, line);
        if (line != "Cursers:") {
            cerr << "Error while parsing cursers (cannot find keyword)\n";
            file.close();
            return false;
        }
        vector<int> new_cursers(Num);
        for (int i=0; i<Num; i++) {
            getline(file, line);
            new_cursers[i] = stoi(line);
        }

        // read solved ints
        getline(file, line);
        if (line != "Solved integrals:") {
            cerr << "Error while parsing solved integrals (cannot find keyword)\n";
            file.close();
            return false;
        }
        getline(file, line);
        int Num_ints = stoi(line);
        cout << "Number of solved ints: " << Num_ints << endl;
        file.close();

        if (Num_ints != int(integrals.size())) {
            cerr << "Cannot load state of the scheduler: Expected "<< Num_ints << " calculated integrals but found " << integrals.size() << ".\n";
            return false;
        }

        // actually set state
        this->tuples= new_tuples;
        this->cursers= new_cursers;
        this->setKnownIntegrals(integrals);
        this->initializeWorker();

        return true;
    }

    void printStartupSummary() const override {
        cout << "The program calculates local field effects integrals, using two overlap\n";
        cout << "densities between electrons and holes (c1, v1, S1) and (c2, v2, S2). This\n";
        cout << "differs from usual Coulomb integrals, which involve overlap densities between\n";
        cout << "either electrons or holes. A short-range Coulomb potential is assumed.\n";
        cout << "\n";
        cout << "Only overlap densities with an absolute charge greater than ABSCHARGE_THRESHOLD\n";
        cout << "are considered. With that all combinations of (c1, v1, S1) and (c2, v2, S2) are\n";
        cout << "calculated. Other thresholds (DIFF_MONO_THRESHOLD, DISTANCE_THRESHOLD,\n";
        cout << "ENERGY_THRESHOLD) are ignored.\n";
    }
};


#endif // SCHEDULER_H