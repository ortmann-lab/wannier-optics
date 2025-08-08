#include "welcome.h"
#include "mpi_tags.h"
#include "confGenerator.h"
#include "coulombIntegral.h"
#include "wannierfunction.h"
#include "density.h"
#include "potential.h"
#include "solver.h"
#include "filehandler.h"
#include "scheduler.h"
#include "parallel_computation.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <chrono>
#include <filesystem>
#include <mpi.h>
#include <cmath>

#include "external/simpleini/SimpleIni.h"

using namespace std;

/**
 * @brief Prints a little helper message on the screen.
 *
 * @param name Name of the program
 */
void printHelp(string name = "wo-coulomb.x", int rank=0) {
    if (rank!=0) return;
    cout << "\nThis programs calculates Coulomb integrals for Wannier functions on a real space grid.\n";
    cout << "To generate all necessary configuration files please run:\n";
    cout << "\n  " << name << " -g config_file.ini\n";
    cout << "\nYou can than edit the generated files and adapt to your desires.\n";
    cout << "To run actual calculations use:\n";
    cout << "\n  " << name << " config_file.ini\n\n\n";
    cout << "Have fun :D\n";
}


void calcScreening(Scheduler const* myScheduler, map< int,WannierFunction > const& vWannMap, map< int,WannierFunction > const& cWannMap,
                map<int,double> const& vMeanDensity, map<int,double> const& cMeanDensity,
                const double SCREENING_RELATIVE_PERMITTIVITY, const double SCREENING_ALPHA,
                string const& outfile, string const& restartfile, const int NUM_OMP_THREADS)
{
    MPI_Barrier(MPI_COMM_WORLD);
    int rank, num_worker;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);


    if (rank==0) cout << "Calculate Yukawa integrals for model screening...\n";
    string yukawa_outfile = outfile + string("_YUK_") + to_string(SCREENING_RELATIVE_PERMITTIVITY) + string("_") + to_string(SCREENING_ALPHA);
    string yukawa_restartfile = restartfile + string("_YUK_") + to_string(SCREENING_RELATIVE_PERMITTIVITY) + string("_") + to_string(SCREENING_ALPHA);

    // setup solver
    YukawaSolver impl(vWannMap, cWannMap, vMeanDensity, cMeanDensity, SCREENING_RELATIVE_PERMITTIVITY, SCREENING_ALPHA);

    // setup scheduler
    unique_ptr<FixedListScheduler> yukawa_scheduler{};
    if (rank==0) {  // master
        int BATCH_SIZE = myScheduler->getBatchSize();
        yukawa_scheduler = make_unique<FixedListScheduler>(BATCH_SIZE, num_worker);

        // load from previous calculation (plan file is not used)
        if (! yukawa_scheduler->tryRestartScheduler(yukawa_restartfile, yukawa_outfile)) { // TODO: make sure that the loaded data is with the same screening parameters

            cout << "\nCalculating list of integrals from previous scheduler... " << endl;
            list<Integral> plan_l = myScheduler->getKnownIntegrals();
            vector<Integral> plan{ begin(plan_l), end(plan_l) };

            yukawa_scheduler->setIntegrals(plan);
        }
    }

    run_parallel_calculations(impl, yukawa_scheduler.get(), NUM_OMP_THREADS, yukawa_outfile, yukawa_restartfile);

    if (rank==0) {
        cout << "Create superposition of Coulomb and Yukawa integrals...\n";

        // make superposition of Coulomb and Yukawa potential
        const list<Integral> coulomb_ints{ myScheduler->getKnownIntegrals() };
        const list<Integral> yukawa_ints{ yukawa_scheduler->getKnownIntegrals() };

        // Check if both lists are of the same size
        if (coulomb_ints.size() != yukawa_ints.size()) {
            throw runtime_error("Different number of Coulomb and Yukawa integrals found!");
        }

        // Create a third list to hold the sums
        list<Integral> superposition;

        // Use iterators to iterate through both lists
        auto it1 = coulomb_ints.begin();
        auto it2 = yukawa_ints.begin();

        while (it1 != coulomb_ints.end() && it2 != yukawa_ints.end()) {
            // cout << "COULOMB: " << (*it1).toString() << endl;
            // cout << "YUKAWA : " << (*it2).toString() << endl;

            // check if elements have same index
            auto ids1 = (*it1).indexes;
            auto ids2 = (*it2).indexes;

            for (size_t i=0; i< ids1.size(); i++) {
                if (ids1[i] != ids2[i]) throw runtime_error("Ids of Coulomb and Yukawa integrals are not the same!!!");
            }

            Integral newInt(*it1);
            newInt.value = (newInt.value/SCREENING_RELATIVE_PERMITTIVITY + (1.-1./SCREENING_RELATIVE_PERMITTIVITY) * (*it2).value)*SCREENING_RELATIVE_PERMITTIVITY;

            // cout << "SUPERPOS: " << newInt.toString() << endl << endl;

            superposition.push_back(newInt);
            ++it1;
            ++it2;
        }

        OutputService::writeResults(outfile + string("_model_screening"), superposition,
            string("model screening (coulomb + yukawa) with parameters: epsilon = ") +
            to_string(SCREENING_RELATIVE_PERMITTIVITY) + string("  alpha = ") + to_string(SCREENING_ALPHA),
            myScheduler->getMaxDist(), 1.0/SCREENING_RELATIVE_PERMITTIVITY);
    }
}



void calc_transition(shared_ptr<RealMeshgrid> const& mesh, map< int,WannierFunction > const& vWannMap, map< int,WannierFunction > const& cWannMap,
                    const vector<vector<double>>& pos_c, const vector<vector<double>>& pos_v, const map<int,int>& wId_to_arrayId_c,
                    const map<int,int>& wId_to_arrayId_v, const int NUM_OMP_THREADS, bool ENABLE_TRANSITION_CORRECTIONS=true)
{
    const string outfile = "TRANSITION";
    int rank, num_worker;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);

    if (filesystem::exists(outfile)) {
        if (rank==0) cout << "Found TRANSITION file. Skip evaluating optical dipoles.\n";
        return;
    }

    int N = mesh->getNumDataPoints();
    double dV = mesh->getdV();
    const double* XX = mesh->getXX();
    const double* YY = mesh->getYY();
    const double* ZZ = mesh->getZZ();

    vector<double> supercell_d = cWannMap.begin()->second.getLatticeInUnitcellBasis();
    vector<int> supercell = vector<int>(3);
    for (int i=0; i<3; i++) {
        supercell[i] = round(supercell_d[i]);
    }
    vector<double> origin_d = cWannMap.begin()->second.getOriginInUnitcellBasis();
    vector<int> origin = vector<int>(3);
    for (int i=0; i<3; i++) {
        origin[i] = int(round(origin_d[i]));
    }

    // prepare everything
    if (rank==0) {
        cout << "Supercell\t: ";
        for (int i=0; i<3; i++) {
            cout <<  supercell[i] << "  ";
        }
        cout << endl << "Origin\t\t: ";
        for (int i=0; i<3; i++) {
            cout <<  origin[i] << "  ";
        }
        cout << endl;
    }

    // create vector of dipoles (depending on rank of process)
    list<OpticalDipole> dipoles = list<OpticalDipole>();
    int counter = 0;
    for (auto v = vWannMap.begin(); v != vWannMap.end(); ++v) {
        for (auto c = cWannMap.begin(); c != cWannMap.end(); ++c) {
            for (int i=origin[0]; i<supercell[0]+origin[0]; i++){
                for (int j=origin[1]; j<supercell[1]+origin[1]; j++){
                    for (int k=origin[2]; k<supercell[2]+origin[2]; k++){
                        if (counter==rank)
                            dipoles.push_back(OpticalDipole(c->second.getId(), v->second.getId(), vector<int>{i,j,k}));
                        counter++;
                        if (counter>=num_worker) counter=0;
                    }
                }
            }
        }
    }

    /**
     *  Actual calculations
     **/
    if (rank==0) {
        cout << "\n\nStart Calculations...\n";
        cout << setprecision(2);
        //cout << "# c	v	shift vector		dipole (x,y,z)\t\t\t\t\t\t <-- DATA\n";
    }
    chrono::steady_clock::time_point beginCalc = chrono::steady_clock::now();

    counter = 0;
    int progress = 0;
    vector<vector<double>> unitcell_t = transpose3x3(cWannMap.begin()->second.getUnitcell());
    for (auto& dipole: dipoles){

        auto itr = vWannMap.find(dipole.getIdv());
        if (itr == vWannMap.end()) throw runtime_error("valence Wannier function not found!");
        WannierFunction const& valence = itr->second;

        itr = cWannMap.find(dipole.getIdc());
        if (itr == cWannMap.end()) throw runtime_error("conduction Wannier function not found!");
        WannierFunction const& cond = itr->second;

        vector<int> shift = dipole.getShift();
        vector<int> negativeShift = {-shift[0],-shift[1],-shift[2]};
        vector<double> vec_shift = matVecMul3x3(unitcell_t, negativeShift);

        int v = wId_to_arrayId_v.at(dipole.getIdv());
        int c = wId_to_arrayId_c.at(dipole.getIdc());

        const vector<double>& mono_val = pos_v[v];
        const vector<double>& mono_cond= pos_c[c];

        vector<double> value={0.0,0.0,0.0};

        unique_ptr<double[], free_deleter>density{ joinedDensity(cond, valence, negativeShift) };  // can also be negative

        // if (isZero(density,N)) {
        //     delete density;
        //     dipoles[n]->setDipole(dipole);   // zero
        //     printResult(dipoles[n]);
        //     continue;
        // }

        // calculate usual dipole moment
        double dx=0.0;
        double dy=0.0;
        double dz=0.0;
        omp_set_num_threads(NUM_OMP_THREADS);
        #pragma omp parallel for shared(density, XX, YY, ZZ) firstprivate(N,dV) reduction(+:dx) reduction(+:dy) reduction(+:dz)
        for (int i=0; i<N; i++) {
            dx += density[i]*XX[i]*dV;
            dy += density[i]*YY[i]*dV;
            dz += density[i]*ZZ[i]*dV;
        }
        value[0] = dx;
        value[1] = dy;
        value[2] = dz;

        // use correction of dipole operator for slightly not orthogonal WF
        if (ENABLE_TRANSITION_CORRECTIONS)
        {
            // overlap of valence and conduction WF (should be very small)
            double overlap = 0.0;
            for (int i=0.0; i<N; i++){
                overlap += density[i] * dV;
            }

            // make corrections due to non-orthogonality
            for (size_t j=0; j<3;j++) {
                value[j] -= 0.5*overlap*( vec_shift[j] + mono_val[j] + mono_cond[j]);
            }
        }

        dipole.setDipole(value);
        counter++;
        progress++;

        if (rank==0) {
            if (counter>=10) {
                counter=0;
                cout << "Progress: " << progress*100.0/dipoles.size() << " %\t\t\r" << flush;
            }
        }
    }

    // collect results from all mpi-processes
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) {
        cout << "Done calculation.\t\t\t\n";
        cout << "Collect data from workers.\n";
    }
    gather_results(dipoles);

    MPI_Barrier(MPI_COMM_WORLD);

    chrono::steady_clock::time_point endCalc = chrono::steady_clock::now();
    if (rank==0) {
        cout << endl<< dipoles.size() << " transition dipole moments have been calculated." << endl;
        cout << "Time spend for calculation = " << chrono::duration_cast<chrono::milliseconds>(endCalc - beginCalc).count() << "ms";
        cout << " = " << chrono::duration_cast<chrono::seconds>(endCalc - beginCalc).count() << "s";
        cout << " = " << chrono::duration_cast<chrono::minutes>(endCalc - beginCalc).count() << "min";
        cout << " = " << chrono::duration_cast<chrono::hours>(endCalc - beginCalc).count() << "h\n";
        cout << "\nWrite data in " << outfile << endl;

        // sort results
        dipoles.sort();
        OutputService::writeResults(outfile, dipoles);
    }

}


void calc_custom_file(int BATCH_SIZE, map< int,string > const& vMapping, map< int,string > const& cMapping,
                        map< int,WannierFunction > const& vWannMap, map< int,WannierFunction > const& cWannMap,
                        map<int,double> const& vMeanDensity, map<int,double> const& cMeanDensity,
                        const int NUM_OMP_THREADS, const bool ENABLE_SCREENING_MODEL, const double SCREENING_RELATIVE_PERMITTIVITY, const double SCREENING_ALPHA)
{
    int rank, num_worker;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);

    // File names:
    const string planfile = "CUSTOM";
    const string outfile = "RESULTS_custom";
    const string restartfile = "RESTART_custom";

    if (! filesystem::exists(planfile)) {
        throw runtime_error("Cannot find CUSTOM file.");
    }

    // setup solver
    CoulombSolver impl(vWannMap, cWannMap);

    // setup scheduler
    unique_ptr<FixedListScheduler> file_scheduler{};
    if (rank==0) {
        file_scheduler = make_unique<FixedListScheduler>(BATCH_SIZE, num_worker);

        // load from previous calculation (plan file is not used)
        if (! file_scheduler->tryRestartScheduler(restartfile, outfile)) {

            // start from scratch (plan file is used)
            cout << "Start calculation from scratch using custom integrals from " << planfile << endl;
            vector<Integral> plan;
            readPlan(plan, planfile);

            cout << "\nCheck compatibility of custom integrals and mapping files..." << endl;
            if (!checkPlan(plan, vMapping, cMapping)) {
                throw runtime_error("Given integrals are not compatible with Wannier functions.");
            }

            file_scheduler->setIntegrals(plan);
        }

        // show details of the scheduler
        cout << "\n---------------------------------------\n";
        cout << "\nStart-up Summary of the scheduler:\n\n";
        file_scheduler->printStartupSummary();
        cout << "\n---------------------------------------\n\n";
    }

    run_parallel_calculations(impl, file_scheduler.get(), NUM_OMP_THREADS, outfile, restartfile);
    if (ENABLE_SCREENING_MODEL)
        calcScreening(file_scheduler.get(), vWannMap, cWannMap, vMeanDensity, cMeanDensity, SCREENING_RELATIVE_PERMITTIVITY, SCREENING_ALPHA, outfile, restartfile, NUM_OMP_THREADS);
}


void calc_local_field_effects(int BATCH_SIZE, map< int,WannierFunction > const& vWannMap, map< int,WannierFunction > const& cWannMap, vector<vector<vector<int>>> const& shells, const int NUM_OMP_THREADS, const double ABSOLUTE_CHARGE_THRESHOLD, filesystem::path const& DATA_DIR)
{
    int rank, num_worker;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);

    // File names:
    const string outfile = DATA_DIR / "LOCALFIELDEFFECTS";
    const string restartfile = DATA_DIR / "RESTART_LFE";

    // setup solver
    LocalFieldEffectsSolver impl(vWannMap, cWannMap);


    // get Indicator parameters for valence WF
    map<Density_descr,Density_indicator> lfe_indicators{ readIndicator_estimates(DATA_DIR / "lfe_indicators.dat", ABSOLUTE_CHARGE_THRESHOLD)};
    if (lfe_indicators.size() == 0) {
        if (rank==0) cout << "Prepare Indicator for LFE (this may take a while)\n";
        chrono::steady_clock::time_point time_calcIndicator_start = chrono::steady_clock::now();
        lfe_indicators = calcLFE_estimates_parallel(cWannMap, vWannMap, shells, ABSOLUTE_CHARGE_THRESHOLD);

        if (rank==0) {
            chrono::steady_clock::time_point time_calcIndicator_stop = chrono::steady_clock::now();
            cout << "Done. Time spent (min): " << chrono::duration_cast<chrono::minutes>(time_calcIndicator_stop - time_calcIndicator_start).count() << "\n";
            cout << "Save lfe_indicators.dat\n";
            saveIndicator_estimates(DATA_DIR / "lfe_indicators.dat", lfe_indicators, ABSOLUTE_CHARGE_THRESHOLD);
        }
    }

    // setup scheduler
    unique_ptr<LocalFieldEffectsScheduler> lfe_scheduler{};

    // check if all Wannier functions that are used in the Indicator files are actually loaded. (otherwise we get a SegFault later)
    if (rank==0) {
        cout << "Check compatibility for lfe_indicators data\n";
        checkCompatibilityIndicator_lfe(lfe_indicators, vWannMap, cWannMap);

        // create scheduler
        lfe_scheduler = make_unique<LocalFieldEffectsScheduler>(
            BATCH_SIZE, num_worker, lfe_indicators, ABSOLUTE_CHARGE_THRESHOLD, restartfile, outfile);
    }

    run_parallel_calculations(impl, lfe_scheduler.get(), NUM_OMP_THREADS, outfile, restartfile);
}


int main(int argc, char *argv[]) {
    int rank, num_worker;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);

    printWelcomeCoulombIntegral(rank);
    cout << fixed;
    cout << setprecision(12);

    if (argc < 2) {
        printHelp(argv[0], rank);
        MPI_Finalize();
        exit(0);
    }
    string parameter = argv[1];
    if (parameter == "-g") {
        if (rank==0) {
            cout << "\nGenerate configuration files for you...\n";
            if (argc > 2) ConfGenerator::writeConfIni(argv[2]);
            else ConfGenerator::writeConfIni();
            ConfGenerator::writeEmptyPlan();
            ConfGenerator::writeMappings();
            ConfGenerator::createPOSFILE();
        }
        MPI_Finalize();
        exit(0);
    }else if(parameter == "-h") {
        printHelp(argv[0], rank);
        MPI_Finalize();
        exit(0);
    }

    /**
     *  read configuration file
     **/

    CSimpleIniA ini;
    if (rank==0) cout << "Reading input file \t: " << argv[1] << endl;
	SI_Error rc = ini.LoadFile(argv[1]);
	if (rc < 0) {
        if (rank==0) cerr << "Error while reading configuration file. Error code: " << rc << endl;
        MPI_Finalize();
        exit(1);
     };


    // general settings
    const string vMappingFile = "vmapping.txt";
    const string cMappingFile = "cmapping.txt";
    const filesystem::path DATA_DIR = "DATA";

    const bool ENABLE_TRANSITION_MATRIX    = bool(stoi(ini.GetValue("general", "ENABLE_TRANSITION_MATRIX", "1")));
    const bool ENABLE_TWO_CENTER_INTEGRALS   = bool(stoi(ini.GetValue("general", "ENABLE_TWO_CENTER_INTEGRALS", "1")));
    const bool ENABLE_THREE_CENTER_INTEGRALS = bool(stoi(ini.GetValue("general", "ENABLE_THREE_CENTER_INTEGRALS", "1")));
    const bool ENABLE_FOUR_CENTER_INTEGRALS  = bool(stoi(ini.GetValue("general", "ENABLE_FOUR_CENTER_INTEGRALS", "1")));
    const bool ENABLE_LOCAL_FIELD_EFFECTS = bool(stoi(ini.GetValue("general", "ENABLE_LOCAL_FIELD_EFFECTS", "1")));
    const bool USE_CUSTOM_COULOMB_FILE  = bool(stoi(ini.GetValue("general", "USE_CUSTOM_COULOMB_FILE", "0")));

    const vector<bool> CRYSTAL_PERIODIC{
        bool(stoi(ini.GetValue("general", "CRYSTAL_PERIODIC_X", "1"))),
        bool(stoi(ini.GetValue("general", "CRYSTAL_PERIODIC_Y", "1"))),
        bool(stoi(ini.GetValue("general", "CRYSTAL_PERIODIC_Z", "1")))
    };

    // screening settings
    const bool ENABLE_SCREENING_MODEL = bool(stoi(ini.GetValue("general", "ENABLE_SCREENING_MODEL", "1")));
    const double SCREENING_RELATIVE_PERMITTIVITY = stod(ini.GetValue("general", "SCREENING_RELATIVE_PERMITTIVITY", "1.0"));
    const double SCREENING_ALPHA = stod(ini.GetValue("general", "SCREENING_ALPHA", "1.0"));
    const double NUM_VALENCE_ELECTRONS = stod(ini.GetValue("general", "NUM_VALENCE_ELECTRONS", "1.0"));

    // advanced settings
    const bool ENABLE_TRANSITION_CORRECTIONS = bool(stoi(ini.GetValue("general", "ENABLE_TRANSITION_CORRECTIONS", "1")));
    const double MONOPOLE_RELATIVE_ERROR_THRESHOLD = stod(ini.GetValue("general", "MONOPOLE_RELATIVE_ERROR_THRESHOLD", "0.05"));
    const double ELECTRON_HOLE_DISTANCE_THRESHOLD = stod(ini.GetValue("general", "ELECTRON_HOLE_DISTANCE_THRESHOLD", "10.0"));
    const double ENERGY_THRESHOLD = stod(ini.GetValue("general", "ENERGY_THRESHOLD", "0.001"));
    const double ABSOLUTE_CHARGE_THRESHOLD = stod(ini.GetValue("general", "ABSOLUTE_CHARGE_THRESHOLD", "0.1"));


    //const uint BATCH_SIZE = stoi(ini.GetValue("general", "BATCH_SIZE", "10"));
    const int NUM_OMP_THREADS = stoi(ini.GetValue("general", "NUM_OMP_THREADS", "1"));
    const vector<int> Nsupercell = configureSupercell(ini);

    // hidden options:
    const bool NORMALIZE = bool(stoi(ini.GetValue("general", "NORMALIZE", "1")));
    const double criterion_extend = stod(ini.GetValue("general", "criterion_extend", "1e-2"));  // TODO: better name and better documentation


    omp_set_num_threads(NUM_OMP_THREADS);

    if (rank==0) {  // display settings
        cout << "You are going to calculate:" << endl;
        if (ENABLE_TRANSITION_MATRIX)      cout << " - transition matrix elements (optical dipole elements)"<< endl;
        if (USE_CUSTOM_COULOMB_FILE)       cout << " - custom list of integrals from file." << endl;
        if (ENABLE_TWO_CENTER_INTEGRALS)   cout << " - 2-center integrals (density-density)"<< endl;
        if (ENABLE_THREE_CENTER_INTEGRALS) cout << " - 3-center integrals (overlap-density) excluding 2-center integrals" << endl;
        if (ENABLE_FOUR_CENTER_INTEGRALS)  cout << " - 4-center integrals (overlap-overlap) excluding 2-center and 3-center integrals" << endl;
        if (ENABLE_LOCAL_FIELD_EFFECTS)    cout << " - local field effects" << endl;
        cout << endl;

        if (ENABLE_TWO_CENTER_INTEGRALS || ENABLE_THREE_CENTER_INTEGRALS || ENABLE_FOUR_CENTER_INTEGRALS || USE_CUSTOM_COULOMB_FILE) { // TODO: maybe write down screening formula and explain calculations
            cout << "Use screening model for Coulomb integrals:\n";
            cout << "ENABLE_SCREENING_MODEL \t\t: " << ENABLE_SCREENING_MODEL << endl;
            cout << "SCREENING_RELATIVE_PERMITTIVITY : " << SCREENING_RELATIVE_PERMITTIVITY << endl;
            cout << "SCREENING_ALPHA \t\t: " << SCREENING_ALPHA << endl;
            cout << "NUM_VALENCE_ELECTRONS    \t: " << NUM_VALENCE_ELECTRONS << endl;
        }else{
            cout << "No need to apply any screening function.\n";
            cout << "Use static screening with\n";
            cout << "SCREENING_RELATIVE_PERMITTIVITY \t: " << SCREENING_RELATIVE_PERMITTIVITY << endl;
        }
        cout << "CRYSTAL_PERIODIC \t\t: ( " << CRYSTAL_PERIODIC[0] << ", " << CRYSTAL_PERIODIC[1] << ", " << CRYSTAL_PERIODIC[2] << " )" << endl;

        cout << endl;
        cout << "Technical / advanced settings: " << endl;
        cout << "MPI-Processes \t\t\t : " << num_worker << endl;
        cout << "NUM_OMP_THREADS \t\t : " << NUM_OMP_THREADS << endl;
        //cout << "BATCH_SIZE \t\t\t : " << BATCH_SIZE << endl;
        cout << "MONOPOLE_RELATIVE_ERROR_THRESHOLD: " << MONOPOLE_RELATIVE_ERROR_THRESHOLD << endl;
        cout << "ABSOLUTE_CHARGE_THRESHOLD \t : " << ABSOLUTE_CHARGE_THRESHOLD << endl;
        cout << "ELECTRON_HOLE_DISTANCE_THRESHOLD : " << ELECTRON_HOLE_DISTANCE_THRESHOLD << endl;
        cout << "ENERGY_THRESHOLD \t\t : " << ENERGY_THRESHOLD << endl;
        cout << "criterion_extend \t\t : " << criterion_extend << endl;
        cout << "Supercell \t\t\t : ( " << Nsupercell[0] << ", " << Nsupercell[1] << ", "<< Nsupercell[2] << " )\n";
        cout << "Normalize WF \t\t\t : " << NORMALIZE << endl;

    }

    if ((ABSOLUTE_CHARGE_THRESHOLD>=1.0) || (ABSOLUTE_CHARGE_THRESHOLD<0)) {
        if (rank==0) {
            cout << "\n\n##############################################################################################\n\n";
            cout << "[ERROR] Hear your master's voice:\n";
            cout << "[ERROR]    ABSOLUTE_CHARGE_THRESHOLD needs to be between 0 and 1\n";
            cout << "\n##############################################################################################\n\n";
        }
        MPI_Finalize();
        exit(1);
    }

    if (ELECTRON_HOLE_DISTANCE_THRESHOLD<=0) {
        if (rank==0) {
            cout << "\n\n##############################################################################################\n\n";
            cout << "[ERROR] Hear your master's voice:\n";
            cout << "[ERROR]    ELECTRON_HOLE_DISTANCE_THRESHOLD needs to be larger than 0\n";
            cout << "\n##############################################################################################\n\n";
        }
        MPI_Finalize();
        exit(1);
    }

    if (ENERGY_THRESHOLD<=0) {
        if (rank==0) {
            cout << "\n\n##############################################################################################\n\n";
            cout << "[ERROR] Hear your master's voice:\n";
            cout << "[ERROR]    ENERGY_THRESHOLD needs to be larger than 0\n";
            cout << "\n##############################################################################################\n\n";
        }
        MPI_Finalize();
        exit(1);
    }

    bool useSupercell = (Nsupercell[0] > 1) || (Nsupercell[1] > 1) || (Nsupercell[2] > 1);
    if ((useSupercell && ENABLE_TRANSITION_MATRIX) || (useSupercell && ENABLE_LOCAL_FIELD_EFFECTS)) {
        if (rank==0) {
            cout << "\n\n##############################################################################################\n\n";
            cout << "[ERROR] Hear your master's voice:\n";
            cout << "[ERROR]    Using a supercell >1 is not compatible with calculation of transitions and local field effects.\n";
            cout << "[ERROR]    Please do calculations with ENABLE_TRANSITION_MATRIX=1 and ENABLE_LOCAL_FIELD_EFFECTS=1 in a separate calculation\n";
            cout << "[ERROR]    without supercell.\n";
            cout << "[ERROR]    For Coulomb integrals: Only use supercells if necessary to avoid aliasing!\n";
            cout << "\n##############################################################################################\n\n";
        }
        MPI_Finalize();
        exit(1);
    }

    if ((! CRYSTAL_PERIODIC[0])&&(! CRYSTAL_PERIODIC[1])&&(! CRYSTAL_PERIODIC[2])) {
        if (rank==0) {
            cout << "\n\n##############################################################################################\n\n";
            cout << "[ERROR] Hear your master's voice:\n";
            cout << "[ERROR]    You need at least one periodic direction of the crystal.\n";
            cout << "[ERROR]    Please set CRYSTAL_PERIODIC_X=1, CRYSTAL_PERIODIC_Y=1 or CRYSTAL_PERIODIC_Z=1\n";
            cout << "\n##############################################################################################\n\n";
        }
        MPI_Finalize();
        exit(1);
    }

    if ((SCREENING_RELATIVE_PERMITTIVITY<=1.0)&&(ENABLE_SCREENING_MODEL)) {
        if (rank==0) {
            cout << "\n\n##############################################################################################\n\n";
            cout << "[ERROR] Hear your master's voice:\n";
            cout << "[ERROR]    If you want to use the screening model SCREENING_RELATIVE_PERMITTIVITY needs to be larger than 1\n";
            cout << "\n##############################################################################################\n\n";
        }
        MPI_Finalize();
        exit(1);
    }

    /**
     * Create DATA_DIR
     **/
    if (! filesystem::exists(DATA_DIR)) filesystem::create_directory(DATA_DIR);


    /**
     *  Read Mapping
     **/
    if (rank==0) cout << "\nReading mapping files... " << endl;
    // map from Wannier function id to path of xsf-file for valence WF
    map< int,string > vMapping;
    map< int,string > cMapping;
    readMapping(vMapping, vMappingFile, "-->", rank);
    readMapping(cMapping, cMappingFile, "-->", rank);

    // print mappings
    if (rank==0) {
        cout << "Valence Wannier functions:\n";
        for (auto itr = vMapping.begin(); itr != vMapping.end(); ++itr) {
            cout << "    " << itr->first << " --> " << itr->second << "\n";
        }
        cout << "Conduction Wannier functions:\n";
        for (auto itr = cMapping.begin(); itr != cMapping.end(); ++itr) {
            cout << "    " << itr->first << " --> " << itr->second << "\n";
        }
        cout << endl;
    }


    /**
     *  Read and check Wannier functions
     **/
    map< int,WannierFunction > vWannMap{ openAllWannierFunctions(vMapping, Nsupercell, rank) };   // map from Wannier function id to WannierFunction object for valence WF
    map< int,WannierFunction > cWannMap{ openAllWannierFunctions(cMapping, Nsupercell, rank) };
    if (rank==0) cout << "Done reading all files.\n\n";

    // check if WF are compatible and make them use the same meshgrid object to save memory
    if ( ! tryToShareMeshgrid(vWannMap, cWannMap)) {
        throw runtime_error("Wannier functions are not compatible to each other!");
    }
    shared_ptr<RealMeshgrid> mesh = cWannMap.begin()->second.getSharedMeshgridPtr();
    if (rank==0) cout << "The meshgrid is used " << mesh.use_count() << " times." << endl;

    /**
     *  Preparation steps for meshgrids, Wannier functions, centers and shells
     **/
    mesh->createMeshgridArrays();

    if (NORMALIZE) {
        if (rank==0) cout << "Normalize all Wannier functions ... ";
        for_each(vWannMap.begin(), vWannMap.end(), [](auto& p){ p.second.normalize();});
        for_each(cWannMap.begin(), cWannMap.end(), [](auto& p){ p.second.normalize();});
        if (rank==0) cout << "done.\n";
    }

    // prepare Wannier centers
    if (rank == 0) cout << "Calculate Wannier centers ...";
    vector<int> id_v, id_c;                           // array position to wannierId
    map<int,int> wId_to_arrayId_v, wId_to_arrayId_c;  // wannierId to array position
    vector<vector<double>> pos_c, pos_v;
    id_v = vector<int>(vWannMap.size());
    int i=0;
    for (auto const& itr : vWannMap) {
        id_v[i] = itr.first;
        wId_to_arrayId_v.insert({itr.first, i});
        i++;
    }

    id_c = vector<int>(cWannMap.size());
    i=0;
    for (auto const& itr : cWannMap) {
        id_c[i] = itr.first;
        wId_to_arrayId_c.insert({itr.first, i});
        i++;
    }
    pos_c = calcWannierCenter(cWannMap, id_c);
    pos_v = calcWannierCenter(vWannMap, id_v);
    if (rank == 0) cout << " done\n";

    // prepare shells and cell information
    vector<double> origin_wf = cWannMap.begin()->second.getOriginInUnitcellBasis();
    vector<double> supercell = cWannMap.begin()->second.getLatticeInUnitcellBasis();
    vector<vector<double>> unitcell = cWannMap.begin()->second.getUnitcell();
    vector<vector<vector<int>>> shells{ createShells(origin_wf, supercell, unitcell, vector<double>{0,0,0}, vector<double>{0,0,0}, CRYSTAL_PERIODIC) };
    double MAX_ALLOWED_DISTANCE = getMaxRadiusInsideCell(mesh->getLattice(), CRYSTAL_PERIODIC, pos_c, pos_v, supercell, origin_wf, unitcell);
    int BATCH_SIZE = min(shells[0].size(), size_t(100));

    if (rank==0) {
        cout << "BATCH_SIZE           = " << BATCH_SIZE << endl;
        cout << "MAX_ALLOWED_DISTANCE = " << MAX_ALLOWED_DISTANCE << " Angström\n\n";

        cout << "Notice: MAX_ALLOWED_DISTANCE is determined by the supercell and Wannier centers.\n";
        cout << "    Integrals beyond this limit are skipped to avoid aliasing errors.\n";
        cout << "    To support larger distances, increase the supercell through wannier90\n";
        cout << "    or using SUPERCELL_DIM_* in the config file.\n\n";

        cout << "Write POSFILE.\n";
        writePOSFILE("POSFILE", unitcell, pos_c, pos_v, CRYSTAL_PERIODIC);

        // save shells for debugging
        cout << "Write "<< DATA_DIR / "shells.txt" << " for debugging.\n";
        ofstream file(DATA_DIR / "shells.txt");
        if (!file.is_open()) {
            cout << "Cannot save shells!" << endl;
            return 1;
        }

        i=0;
        for (const auto& shell: shells) {
            for (const auto& R: shell) {
                file << i << "\t" << R[0] << "  " << R[1] << "  " << R[2] << "\n";
            }
            i++;
        }
        file.close();
    }


    /**
     *  Preparation model screening, i.e. calculate expectation values of the ground state density for every WF
     **/
    map<int,double> vMeanDensity{};
    map<int,double> cMeanDensity{};
    if (filesystem::exists("CHGCAR")) {
        if (rank==0) cout << "\nFound CHGCAR. Calculate ground state expecation values for screening model (experimental!).\n";
        CHGCAR chg{ read_CHGCAR("CHGCAR", rank) };
        if (rank==0) cout << "Number of valence electrons in CHGCAR: " << chg.getNumElectrons() << endl;

        chg.makeCompatible(vWannMap.begin()->second);


        if (rank==0) cout << "\nMean densities for valence WF used for model screening:\n";
        for (auto itr = vWannMap.begin(); itr != vWannMap.end(); ++itr) {
            double value = calcExpectationDensity(chg, itr->second);
            if (rank==0) cout << itr->first << " :\t" << value << endl;
            vMeanDensity.insert({itr->first, value});
        }

        if (rank==0) cout << "Mean densities for conduction WF used for model screening:\n";
        for (auto itr = cWannMap.begin(); itr != cWannMap.end(); ++itr) {
            double value = calcExpectationDensity(chg, itr->second);
            if (rank==0) cout << itr->first << " :\t" << value << endl;
            cMeanDensity.insert({itr->first, value});
        }

        chg.reset();  // not needed anymore (free memory)
        if (rank==0) cout << endl;
    }else{
        if (rank==0) cout << "\nNo CHGCAR was found. Use average electron density in screening model.\n";
        double value = NUM_VALENCE_ELECTRONS / det3x3(unitcell);

        if (rank==0) cout << "Mean densities for valence WF used for model screening:\n";
        for (auto itr = vWannMap.begin(); itr != vWannMap.end(); ++itr) {
            if (rank==0) cout << itr->first << " :\t" << value << endl;
            vMeanDensity.insert({itr->first, value});
        }

        if (rank==0) cout << "Mean densities for conduction WF used for model screening:\n";
        for (auto itr = cWannMap.begin(); itr != cWannMap.end(); ++itr) {
            if (rank==0) cout << itr->first << " :\t" << value << endl;
            cMeanDensity.insert({itr->first, value});
        }
    }
    // calculate yukawa parameter
    // if (rank==0) cout << "\nCalculate screening parameters:\n";
    // const double YUK_PARAMETER = YukawaPotential::calc_yukawa_screening_factor(NUM_VALENCE_ELECTRONS,det3x3(unitcell),SCREENING_RELATIVE_PERMITTIVITY, SCREENING_ALPHA);



    /**
     * Summary user output:
     */

    if (rank == 0) {
        cout << "\nSummary of first conduction Wannier function:\n";
        cWannMap.begin()->second.printSummary();

        cout << "\n\nValence Wannier center:\n";
        for (size_t i=0; i<pos_v.size(); i++) {
            cout << id_v[i] << "\t" << pos_v[i][0] << " "<< pos_v[i][1] << " "<< pos_v[i][2] << endl;
        }

        cout << "\n\nConduction Wannier center:\n";
        for (size_t i=0; i<pos_c.size(); i++) {
            cout << id_c[i] << "\t" << pos_c[i][0] << " "<< pos_c[i][1] << " "<< pos_c[i][2] << endl;
        }
        cout << endl;
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) cout << "\nDone preparation. Start actual calculations.\n\n";

    /**
     *  Start calculations
     **/

    if (ENABLE_TRANSITION_MATRIX) { // calculate optical dipoles
        if (rank==0) {
            cout << "\n----------------------------------------------------------------------\n";
            cout << "\nCalculate electron-hole transition elements (optical dipoles):\n";
            cout << "\n----------------------------------------------------------------------\n";
        }
        calc_transition(mesh, vWannMap, cWannMap, pos_c, pos_v, wId_to_arrayId_c, wId_to_arrayId_v, NUM_OMP_THREADS, ENABLE_TRANSITION_CORRECTIONS);
    }



    MPI_Barrier(MPI_COMM_WORLD);
    if (USE_CUSTOM_COULOMB_FILE) { // calculate custom integrals
        if (rank==0) {
            cout << "\n----------------------------------------------------------------------\n";
            cout << "\nCalculate individual Coulomb integrals from file:\n";
            cout << "\n----------------------------------------------------------------------\n";
        }
        calc_custom_file(BATCH_SIZE, vMapping, cMapping, vWannMap, cWannMap, vMeanDensity, cMeanDensity, NUM_OMP_THREADS, ENABLE_SCREENING_MODEL, SCREENING_RELATIVE_PERMITTIVITY, SCREENING_ALPHA);
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if (ENABLE_TWO_CENTER_INTEGRALS)  // calculate density-density integrals
    {
        if (rank==0) {
            cout << "\n----------------------------------------------------------------------\n";
            cout << "\nCalculate 2-center Coulomb integrals (density-density):\n";
            cout << "\n----------------------------------------------------------------------\n";
        }

        const string outfile = DATA_DIR / "DENSITY_DENSITY";
        const string restartfile = DATA_DIR / "RESTART_2";

        // setup solver
        CoulombSolver impl(vWannMap, cWannMap);

        // setup scheduler
        unique_ptr<CombinedDensityDensityFixedRScheduler> density_scheduler{};
        if (rank==0) {  // master
            density_scheduler = make_unique<CombinedDensityDensityFixedRScheduler>(
                BATCH_SIZE, num_worker, id_c, id_v, shells, pos_c, pos_v,
                unitcell, MONOPOLE_RELATIVE_ERROR_THRESHOLD, MAX_ALLOWED_DISTANCE, restartfile, outfile);
        }

        // run all calculations
        run_parallel_calculations(impl, density_scheduler.get(), NUM_OMP_THREADS, outfile, restartfile);
        if (ENABLE_SCREENING_MODEL)
            calcScreening(density_scheduler.get(), vWannMap, cWannMap, vMeanDensity, cMeanDensity, SCREENING_RELATIVE_PERMITTIVITY, SCREENING_ALPHA, outfile, restartfile, NUM_OMP_THREADS);

    }


    MPI_Barrier(MPI_COMM_WORLD);
    map<Density_descr,Density_indicator> vIndicators{};
    map<Density_descr,Density_indicator> cIndicators{};
    if (ENABLE_THREE_CENTER_INTEGRALS || ENABLE_FOUR_CENTER_INTEGRALS)  // prepare indicator data for 3- and 4-center integrals
    {
        if (rank==0) {
            cout << "\n----------------------------------------------------------------------\n";
            cout << "\nPrepare indicators for valence and conduction densities:\n";
            cout << "\n----------------------------------------------------------------------\n";
        }

        // get Indicator parameter for valence WF
        vIndicators = readIndicator_estimates(DATA_DIR /"vIndicators.dat", criterion_extend);
        if (vIndicators.size() == 0) {
            if (rank==0) cout << "Prepare Indicator for valence WF (this may take a while)\n";
            chrono::steady_clock::time_point time_calcIndicator_start = chrono::steady_clock::now();
            vIndicators = calcIndicator_estimates_parallel(vWannMap, shells, criterion_extend, ABSOLUTE_CHARGE_THRESHOLD);
            if (rank==0) {
                chrono::steady_clock::time_point time_calcIndicator_stop = chrono::steady_clock::now();
                cout << "Done. Time spent (min): " << chrono::duration_cast<chrono::minutes>(time_calcIndicator_stop - time_calcIndicator_start).count() << "\n";
                cout << "Save vIndicators.dat\n";
                saveIndicator_estimates(DATA_DIR / "vIndicators.dat", vIndicators, criterion_extend);
            }
        }

        // get Indicator parameter for conduction WF
        cIndicators = readIndicator_estimates(DATA_DIR / "cIndicators.dat", criterion_extend);  // TODO only read if ABSOLUTE_CHARGE_THRESHOLD is the same
        if (cIndicators.size() == 0) {
            if (rank==0) cout << "Prepare Indicator for conduction WF (this may take a while)\n";
            chrono::steady_clock::time_point time_calcIndicator_start = chrono::steady_clock::now();
            cIndicators = calcIndicator_estimates_parallel(cWannMap, shells,criterion_extend, ABSOLUTE_CHARGE_THRESHOLD);
            if (rank==0) {
                chrono::steady_clock::time_point time_calcIndicator_stop = chrono::steady_clock::now();
                cout << "Done. Time spent (min): " << chrono::duration_cast<chrono::minutes>(time_calcIndicator_stop - time_calcIndicator_start).count() << "\n";
                cout << "Save cIndicators.dat\n";
                saveIndicator_estimates(DATA_DIR /"cIndicators.dat", cIndicators, criterion_extend); // TODO save also ABSOLUTE_CHARGE_THRESHOLD
            }
        }

        // check if all Wannier functions that are used in the Indicator files are actually loaded. (otherwise we get a SegFault later)
        if (rank==0) {
            cout << "Check compatibility for vIndicators data\n";
            checkCompatibilityIndicator(vIndicators, vWannMap);
            cout << "Check compatibility for cIndicators data\n";
            checkCompatibilityIndicator(cIndicators, cWannMap);
        }
    }

    if (ENABLE_THREE_CENTER_INTEGRALS) // calculate overlap-overlap integrals
    {
        if (rank==0) {
            cout << "\n----------------------------------------------------------------------\n";
            cout << "\nCalculate 3-center Coulomb integrals (overlap-density):\n";
            cout << "\n----------------------------------------------------------------------\n";
        }

        const string outfile = DATA_DIR / "OVERLAP_DENSITY";
        const string restartfile = DATA_DIR / "RESTART_3";

        // setup solver
        CoulombSolver impl(vWannMap, cWannMap);

        // create scheduler
        unique_ptr<OverlapDensityScheduler> overlap_scheduler{};
        if (rank==0) {
            overlap_scheduler = make_unique<OverlapDensityScheduler>(
                BATCH_SIZE, num_worker, cIndicators, vIndicators, shells, unitcell,
                ABSOLUTE_CHARGE_THRESHOLD, ELECTRON_HOLE_DISTANCE_THRESHOLD, ENERGY_THRESHOLD * SCREENING_RELATIVE_PERMITTIVITY, "OverlapDensity",
                restartfile, outfile);
        }

        // run all calculations
        run_parallel_calculations(impl, overlap_scheduler.get(), NUM_OMP_THREADS, outfile, restartfile);
        if (ENABLE_SCREENING_MODEL)
            calcScreening(overlap_scheduler.get(), vWannMap, cWannMap, vMeanDensity, cMeanDensity, SCREENING_RELATIVE_PERMITTIVITY, SCREENING_ALPHA, outfile, restartfile, NUM_OMP_THREADS);

    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (ENABLE_FOUR_CENTER_INTEGRALS) // calculate overlap-overlap integrals
    {
        if (rank==0) {
            cout << "\n----------------------------------------------------------------------\n";
            cout << "\nCalculate 4-center Coulomb integrals (overlap-overlap):\n";
            cout << "\n----------------------------------------------------------------------\n";
        }

        const string outfile = DATA_DIR / "OVERLAP_OVERLAP";
        const string restartfile = DATA_DIR / "RESTART_4";

        // setup solver
        CoulombSolver impl(vWannMap, cWannMap);

        // setup scheduler
        unique_ptr<OverlapDensityScheduler> overlap_overlap_scheduler{};
        if (rank==0) {
            overlap_overlap_scheduler = make_unique<OverlapDensityScheduler>(
                BATCH_SIZE, num_worker, cIndicators, vIndicators, shells, unitcell,
                ABSOLUTE_CHARGE_THRESHOLD, ELECTRON_HOLE_DISTANCE_THRESHOLD, ENERGY_THRESHOLD * SCREENING_RELATIVE_PERMITTIVITY, "OverlapOverlap",
                restartfile, outfile);
        }

        // run all calculations
        run_parallel_calculations(impl, overlap_overlap_scheduler.get(), NUM_OMP_THREADS, outfile, restartfile);
        if (ENABLE_SCREENING_MODEL)
            calcScreening(overlap_overlap_scheduler.get(), vWannMap, cWannMap, vMeanDensity, cMeanDensity, SCREENING_RELATIVE_PERMITTIVITY, SCREENING_ALPHA, outfile, restartfile, NUM_OMP_THREADS);

    }

    // filter and merge Coulomb integrals into one single file
    cout << setprecision(12);
    if (rank==0) {
        cout << "\n----------------------------------------------------------------------\n";
        cout << "\nReopen, filter and merge all Coulomb integrals into a single file:\n";
        cout << "\n----------------------------------------------------------------------\n";

        vector<string> coulomb_files;
        string comment;
        if (ENABLE_SCREENING_MODEL) {
            coulomb_files = vector<string>();
            if (ENABLE_TWO_CENTER_INTEGRALS) coulomb_files.push_back(DATA_DIR / "DENSITY_DENSITY_model_screening");
            if (ENABLE_THREE_CENTER_INTEGRALS) coulomb_files.push_back(DATA_DIR / "OVERLAP_DENSITY_model_screening");
            if (ENABLE_FOUR_CENTER_INTEGRALS) coulomb_files.push_back(DATA_DIR / "OVERLAP_OVERLAP_model_screening");
            comment = string("Model screening with epsilon = ") + to_string(SCREENING_RELATIVE_PERMITTIVITY)
                + string("  alpha = ") + to_string(SCREENING_ALPHA) + string("  num_electron = ") + to_string(NUM_VALENCE_ELECTRONS) +
                string("; Energy threshold for 3- and 4-center integrals = ") + to_string(ENERGY_THRESHOLD) + string("eV");
        }else{
            coulomb_files = vector<string>();
            if (ENABLE_TWO_CENTER_INTEGRALS) coulomb_files.push_back(DATA_DIR / "DENSITY_DENSITY");
            if (ENABLE_THREE_CENTER_INTEGRALS) coulomb_files.push_back(DATA_DIR / "OVERLAP_DENSITY");
            if (ENABLE_FOUR_CENTER_INTEGRALS) coulomb_files.push_back(DATA_DIR / "OVERLAP_OVERLAP");
            comment = string("Constant screening with epsilon = ") + to_string(SCREENING_RELATIVE_PERMITTIVITY) +
                string("; Energy threshold for 3- and 4-center integrals = ") + to_string(ENERGY_THRESHOLD)+ string("eV");
        }

        list<Integral> all;
        double max_dist = -1;
        for (size_t i=0; i<coulomb_files.size(); i++){
            string f = coulomb_files[i];
            cout << "open " << f << endl;
            auto file = CoulombFileReader(f);
            const vector<Integral>& tmp = file.getElements();
            auto unitcell_T = transpose3x3(unitcell);

            if (i==00 && ENABLE_TWO_CENTER_INTEGRALS) {
                max_dist = file.getDist();
                if (max_dist > 0) {
                    cout << "Filter values electron-hole distance <= " << max_dist << " A" << endl;
                    for(const auto& e: tmp) {

                        vector<double> vec_shift = matVecMul3x3(unitcell_T, e.getRD());
                        int c = wId_to_arrayId_c.at(e.indexes[0]);
                        int v = wId_to_arrayId_v.at(e.indexes[2]);

                        double dist = L2Norm(vector<double>{
                                pos_c[c][0]-pos_v[v][0]-vec_shift[0],
                                pos_c[c][1]-pos_v[v][1]-vec_shift[1],
                                pos_c[c][2]-pos_v[v][2]-vec_shift[2]
                            });

                        if (dist<=max_dist)
                            all.push_back(e);
                    }
                }else{
                    cout << "Don't filter electron-hole distance because no distance is given." << endl;
                    for(const auto& e: tmp) {
                        all.push_back(e);
                    }
                }
            }else{
                cout << "Only values >= " << ENERGY_THRESHOLD << " eV" << endl;
                for(const auto& e: tmp) {
                    if (abs(e.value / SCREENING_RELATIVE_PERMITTIVITY)>= ENERGY_THRESHOLD) {  // filter by energy
                        all.push_back(e);
                    }
                }
            }

            // add hermitian conjugate integrals to the list
            list<Integral> all_conj;
            for (const auto& a: all) {
                vector<int> RD = a.getRD();
                vector<int> Rc = a.getRc();
                vector<int> Rv = a.getRv();

                // add hermitian conjugated integral
                vector<int> RD_conj = vector<int>{
                    RD[0]-Rc[0]+Rv[0], RD[1]-Rc[1]+Rv[1], RD[2]-Rc[2]+Rv[2]
                };
                vector<int> Rc_conj = vector<int>{-Rc[0], -Rc[1], -Rc[2]};
                vector<int> Rv_conj = vector<int>{-Rv[0], -Rv[1], -Rv[2]};
                Integral conj_int(
                    a.indexes[1],  // c1
                    a.indexes[0],  // c2
                    a.indexes[3],  // v1
                    a.indexes[2],  // v2
                    RD_conj,
                    Rc_conj,
                    Rv_conj
                );
                conj_int.value = a.value;

                if (a != conj_int) {
                    all_conj.push_back(conj_int);
                }
            }

            // merge all integrals and their hermitian conjugates
            all.splice(all.end(), all_conj);

            cout << "Total number of Coulomb integrals (+ hermitian conjugate): " << all.size() << endl;
        }

        cout << "Save COULOMB file." << endl;
        OutputService::writeResults("COULOMB", all,comment, max_dist,1.0/SCREENING_RELATIVE_PERMITTIVITY);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (ENABLE_LOCAL_FIELD_EFFECTS)
    {
        if (rank==0) {
            cout << "\n----------------------------------------------------------------------\n";
            cout << "\nCalculate local field effects integrals:\n";
            cout << "\n----------------------------------------------------------------------\n";
        }

        calc_local_field_effects(BATCH_SIZE, vWannMap, cWannMap, shells, NUM_OMP_THREADS, ABSOLUTE_CHARGE_THRESHOLD, DATA_DIR);

        if (rank==0) {
            cout << "\n----------------------------------------------------------------------\n";
            cout << "\nReopen and filter local field effects integrals:\n";
            cout << "\n----------------------------------------------------------------------\n";

            // filter by energy
            string f = DATA_DIR / "LOCALFIELDEFFECTS";
            cout << "open " << f << endl;
            auto file = CoulombFileReader(f);
            const vector<Integral>& tmp = file.getElements();

            list<Integral> lfe;
            cout << "Only values >= " << ENERGY_THRESHOLD << endl;
            for(const auto& e: tmp) {
                if (abs(e.value)>= ENERGY_THRESHOLD) {  // filter by energy
                    lfe.push_back(e);
                }
            }

            // the hermitian conjugate is already added in the scheduler

            cout << "Save LOCALFIELDEFFECTS file." << endl;
            OutputService::writeResults("LOCALFIELDEFFECTS", lfe,"local field effects, only values >= " + to_string(ENERGY_THRESHOLD) + string("eV") ,0.0,1.0);
        }

    }else{
        // write empty LFE
        list<Integral> lfe;
        if(rank==0) cout << "Save empty LOCALFIELDEFFECTS file." << endl;
        OutputService::writeResults("LOCALFIELDEFFECTS", lfe,"local field effects",0.0,1.0);
    }


    // clean up
    MPI_Finalize();
    return 0;
}