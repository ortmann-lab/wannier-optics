#ifndef CONFGENERATOR_H
#define CONFGENERATOR_H

/**
 * @file confGenerator.h
 * @author Konrad Merkel
 * @brief Generates all possible config files to be more user friendly.
 *
 */

#include "filehandler.h"
#include "algebra.h"
#include "wannierfunction.h"

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <regex>

#include "external/simpleini/SimpleIni.h"


using namespace std;


/**
 * @brief Creates a suitable supercell vector from the input
 */
vector<int> configureSupercell(CSimpleIniA const & ini) {
    string Nx = string(ini.GetValue("general", "SUPERCELL_DIM_X", "1"));
    string Ny = string(ini.GetValue("general", "SUPERCELL_DIM_Y", "1"));
    string Nz = string(ini.GetValue("general", "SUPERCELL_DIM_Z", "1"));

    vector<int> N{stoi(Nx), stoi(Ny), stoi(Nz)};
    return N;
}



/**
 * @brief Template input file for CoulombIntegral.x
 */

// to generate string use:   cat config.ini | awk '{print $0 "\\n\\"}'
// TODO: update!
const string conf_ini = "\
[general]\n\
\n\
#-------------------------------------------------------------------------------\n\
#                        Calculation Options (0=disabled, 1=enabled)\n\
#-------------------------------------------------------------------------------\n\
\n\
# Calculate electron-hole transition matrix elements (optical dipole moments)\n\
# (required for optical calculations)\n\
ENABLE_TRANSITION_MATRIX = 1\n\
\n\
# Calculate Coulomb integrals (required for optical calculations)\n\
ENABLE_TWO_CENTER_INTEGRALS = 1\n\
ENABLE_THREE_CENTER_INTEGRALS = 1\n\
ENABLE_FOUR_CENTER_INTEGRALS = 1\n\
\n\
# Calculate local field effect integrals (recommended for optical calculations)\n\
ENABLE_LOCAL_FIELD_EFFECTS = 1\n\
\n\
# Load custom Coulomb integrals from an external file (optional)\n\
USE_CUSTOM_COULOMB_FILE = 0\n\
\n\
\n\
#-------------------------------------------------------------------------------\n\
#                       Screening Model Parameters\n\
#-------------------------------------------------------------------------------\n\
\n\
# Enable a screening model for Coulomb integrals (0=disabled, 1=enabled)\n\
# For details, see Cappellini et al., Phys. Rev. B 47, 9892.\n\
ENABLE_SCREENING_MODEL = 1\n\
\n\
# Relative permittivity of the material (example: Silicon = 11.68)\n\
SCREENING_RELATIVE_PERMITTIVITY = 11.68\n\
\n\
# Universal screening parameter (example: Silicon = 1.563)\n\
SCREENING_ALPHA = 1.563\n\
\n\
# Number of valence electrons in the material used for calculating the average\n\
# electron density. It is independent of the number of valence Wannier function\n\
# that are actual used. (example: Silicon = 8)\n\
NUM_VALENCE_ELECTRONS = 8\n\
\n\
# Periodicity of the crystal for each lattice direction. Only relevant for 2D\n\
# and 1D materials. (1=periodic, 0=not periodic)\n\
CRYSTAL_PERIODIC_X = 1\n\
CRYSTAL_PERIODIC_Y = 1\n\
CRYSTAL_PERIODIC_Z = 1\n\
\n\
\n\
#-------------------------------------------------------------------------------\n\
#                    Technical and Advanced Parameters\n\
#-------------------------------------------------------------------------------\n\
\n\
# Apply correction scheme to transition matrix elements\n\
# (0=disabled, 1=enabled)\n\
ENABLE_TRANSITION_CORRECTIONS = 1\n\
\n\
# Maximum allowable relative difference between monopole-monopole approximation\n\
# and value of 2-center integrals.\n\
# Acts as a stopping criterion for 2-center integrals. Only integrals exceeding\n\
# this threshold will be calculated in full detail.\n\
# Recommended range: 0.01 - 0.05\n\
# Unit: None\n\
MONOPOLE_RELATIVE_ERROR_THRESHOLD = 0.05\n\
\n\
# Maximum electron-hole distance considered in 3- and 4-center integrals.\n\
# Acts as stopping criterion for 3- and 4-center integrals.\n\
# Integrals with larger electron-hole distances are ignored.\n\
# Recommended range: 6.0 - 10.0\n\
# Unit: Angstrom\n\
ELECTRON_HOLE_DISTANCE_THRESHOLD = 10.0\n\
\n\
# Minimum energy threshold for calculating 3- and 4-center integrals and\n\
# local field effect integrals.\n\
# Acts as a stopping and filter criterion for 3- and 4-center integrals.\n\
# For local field effects, this acts as a filter only, not as a stoping\n\
# criterion.\n\
# Integral values below this energy are excluded.\n\
# Recommended range: 0.001 - 0.01\n\
# Unit: eV\n\
ENERGY_THRESHOLD = 0.001\n\
\n\
# Minimum absolute charge of an overlap density to be considered for 3- and\n\
# 4-center integrals and local field effects.\n\
# Recommended range: 0.01 - 0.1\n\
# Unit: None\n\
ABSOLUTE_CHARGE_THRESHOLD = 0.1\n\
\n\
\n\
#-------------------------------------------------------------------------------\n\
#                        Performance and Parallelization\n\
#-------------------------------------------------------------------------------\n\
\n\
# Number of OpenMP threads allocated per MPI process (for nested parallelization).\n\
# Increase if additional threads are available for improved parallel performance.\n\
NUM_OMP_THREADS = 1\n\
\n\
# Supercell dimensions for Fourier-space calculations (applies only to\n\
# Coulomb integrals). Helps prevent aliasing in Fourier space but should\n\
# only be adjusted if necessary, as larger values significantly increase\n\
# computation time. A supercell >1 is incompatible with calculations other\n\
# than Coulomb integrals.\n\
SUPERCELL_DIM_X = 1\n\
SUPERCELL_DIM_Y = 1\n\
SUPERCELL_DIM_Z = 1\n\
";

const string mapping_header = "\
# Mapping of Wannier functions (valence or conduction)\n\
#\n\
# This file contains the mapping Wannier functions given by their\n\
# quantum number (id) to the path of the associated XSF file (XCrysDen Structure\n\
# File) that is obtained from wannier90.\n\
# The id also specifies the c1,c2,v1 and v2 variables in the plan and output file.\n\
#\n\
# Format: id --> path/to/xsf-file\n\
#\n\
";


/**
 * @brief Interactive generator of all different config files that are needed to run calculations
 *
 */
class ConfGenerator
{
public:

    /**
     * @brief Generates a input.ini file with default inputs for CoulombIntegral.x
     */
    static bool writeConfIni(string filename="input.ini") {
        if (filesystem::exists(filename)) {
            cout << "[-] Skip configuration ini-file (already exists).\n";
            return false;
        }

        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Cannot open file! filename = " << filename << endl;
            return false;
        }

        cout << "[+] Write configuration file " << filename << endl;
        file << conf_ini << endl;
        file.close();
        return true;
    }

    /**
     * @brief Writes empty plan file.
     *
     * This should be a help for the user to see the format of the plan file.
     */
    static bool writeEmptyPlan() {
        string filename = "CUSTOM";  // is fixed to be compatible with input.ini

        if (filesystem::exists(filename)) {
            cout << "[-] Skip plan file (already exists).\n";
            return false; // do not overwrite any file ever!
        }

        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Cannot open file! filename = " << filename << endl;
            return false;
        }
        cout << "[+] Write plan file " << filename << endl;
        file << "4\n";
        file << "#\n";
        file << "# You can use this file to specify a collection of Coulomb integrals that you want to\n";
        file << "# calculate. To use this file you need to set \n";
        file << "#      CUSTOM_FILE = 1\n";
        file << "# in the configuration file. \n";
        file << "# \n";
        file << "# c1, c2, v1, v2,   RD,         Rc,         Rv\n";
        file << "  1   1   1   1    0 0 0      0 0 0      0 0 0  # this is an example\n";
        file << "  1   1   1   1    1 0 0      0 0 0      0 0 0\n";
        file << "  1   1   1   1    2 0 0      0 0 0      0 0 0\n";
        file << "  1   1   1   1    3 0 0      0 0 0      0 0 0\n";
        file << "#\n";
        file << "# TODO: create your own list of Coulomb integrals\n";
        file << "# TODO: update first line of this file to the number of Coulomb integrals\n";

        file.close();
        return true;
    }

    static bool writeMappingFile(string filename, map<int,string> filenameMapping) {
        // write files
        cout << "[+] Write " << filename << endl;
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Cannot open file!" << endl;
            return false;
        }

        // write file header
        file << mapping_header;

        // write actual mapping
        for (auto itr = filenameMapping.begin(); itr != filenameMapping.end(); ++itr) {
            cout << "\t" << itr->first << " --> " << itr->second << endl;
            file << itr->first << " --> " << itr->second << endl;
        }
        file.close();
        return true;
    }

    static map<int,string> searchForXsfFiles(string path, string seedname) {
        map<int,string> filenameMapping;

        string pattern =  string("^") + seedname + string("_([0-9]+)\\.xsf");
        regex expr(pattern);
        smatch m;
        string filename;

        for (const auto & entry : filesystem::directory_iterator(path)) {

            filename = entry.path().filename();
            if (regex_search (filename,m,expr)) {
                if (m.size()==2) {
                    //cout << "\t" << stoi(m[1]) << " --> " << string(entry.path()) << endl;
                    filenameMapping.insert({stoi(m[1]), string(entry.path())});
                }
            }
        }
        return filenameMapping;
    }

    /**
     * @brief Generates vmapping.txt and cmapping.txt from given paths.
     *
     * The method looks into the paths and tries to find the Wannier funcitons as xsf-files
     * using regular expressions. Needed data like paths or wannier90 seedname will be ask
     * interactively.
     */
    static bool writeMappings() {

        string path, seedname="wannier90";

        if (filesystem::exists("vmapping.txt")) {
            cout << "[-] Skip vmapping.txt (already exists).\n";
        }else {
            cout << "Path to valence Wannier functions? ";
            cin >> path;

            cout << "wannier90 seedname for valence WF? [default=wannier90] ";
            cin.ignore();
            if (std::cin.peek() == '\n') {          //check if next character is newline
                seedname = "wannier90";             //and assign the default
            } else if (!(std::cin >> seedname)) {   //be sure to handle invalid input
                std::cout << "Invalid input. Use default value!\n";
                seedname = "wannier90";
            }

            map<int,string> vMapping = searchForXsfFiles(path, seedname);
            writeMappingFile("vmapping.txt",vMapping);


        }

        if (filesystem::exists("cmapping.txt")) {
            cout << "[-] Skip cmapping.txt (already exists).\n";
        }else {
            cout << "Path to conduction Wannier functions? ";
            cin >> path;

            cout << "wannier90 seedname for conduction WF? [default=wannier90] ";
            cin.ignore();
            if (std::cin.peek() == '\n') {          //check if next character is newline
                seedname = "wannier90";             //and assign the default
            } else if (!(std::cin >> seedname)) {   //be sure to handle invalid input
                std::cout << "Invalid input. Use default value!\n";
                seedname = "wannier90";
            }

            map<int,string> cMapping = searchForXsfFiles(path, seedname);
            writeMappingFile("cmapping.txt",cMapping);
        }

        return true;
    }


    /**
     * @brief Shifts a group of Wannier functions to their individual home cell and saves them as xsf-files.
     *
     * Corresponding monopole, wannier and filename maps are updated such that they can be used later for generating
     * all kinds of input files or for calculations. The transformation is like a rotate algorithm. That means
     * no information lost.
     *
     * All input references are modified during the process.
     *
     * @param MonoMap
     * @param wannierMap
     * @param filenameMapping
     */
    static void rotateToHomeCell(map< int,Monopole >& MonoMap, map< int,WannierFunction >& wannierMap, map< int,string >& filenameMapping)
    {
        vector<double> pos_d{0.,0.,0.};
        vector<int> R_shift{0,0,0};
        vector< vector<double> > invUnitcell = invMat3x3(transpose3x3(wannierMap.begin()->second.getUnitcell()));

        for (auto itr = MonoMap.begin(); itr != MonoMap.end(); ++itr) {

            // calcuate shift
            pos_d = matVecMul3x3(invUnitcell, vector<double>{itr->second.x, itr->second.y, itr->second.z});
            R_shift[0] = -floor(pos_d[0]+1e-6);
            R_shift[1] = -floor(pos_d[1]+1e-6);
            R_shift[2] = -floor(pos_d[2]+1e-6);


            if ( (R_shift[0]==0) && (R_shift[1]==0) && (R_shift[2]==0) ) {
                cout << "\t" << itr->first<< "\t: already in home cell\n";
                continue;
            }
            cout << "\t" << itr->first<< "\t: apply shift to home cell (";
            cout << R_shift[0] << ", " << R_shift[1] << ", " << R_shift[2] << ")" << endl;

            auto wann_itr = wannierMap.find(itr->first);
            if (wann_itr == wannierMap.end()) {  // not found
                cerr << "Cannot find Wannier function for index "<< itr->first <<" ... skip "<< endl;
                continue;
            }

            // update maps for WF and Monopoles
            WannierFunction& wann{wann_itr->second};
            rotateWannierFunction(wann, R_shift);


            unique_ptr<double[], free_deleter>density{ joinedDensity(wann, wann, vector<int>{0,0,0}) };
            MonoMap[itr->first] = getMonopole(density.get(), wann.getMeshgrid());

            string filename = filenameMapping.find(itr->first)->second;

            string ending = filename.substr(filename.size()-4, filename.size()-1);
            if (ending==".xsf") {
                filename = filename.substr(0, filename.size()-4) + "_shifted.xsf";
            }else{
                filename += "_shifted.xsf";
            }

            cout << "write shifted WF to file: " << filename << endl;
            if (XSF_controller::save_file(filename, wann)) {  // save new file
                filenameMapping[itr->first] = filename;  // update file name
            }
        }
    }

    static bool createPOSFILE() {  // TODO make a better interface and use writePOSFILE(...)

        cout << "[+] Read mappings from file\n";

        string vMappingFile = "vmapping.txt";
        string cMappingFile = "cmapping.txt";
        map< int,string > vMapping;
        map< int,string > cMapping;
        readMapping(vMapping, vMappingFile, "-->");
        readMapping(cMapping, cMappingFile, "-->");

        // print mappings
        cout << "Valence Wannier functions:\n";
        for (auto itr = vMapping.begin(); itr != vMapping.end(); ++itr) {
            cout << itr->first << " --> " << itr->second << "\n";
        }
        cout << "Conduction Wannier functions:\n";
        for (auto itr = cMapping.begin(); itr != cMapping.end(); ++itr) {
            cout << itr->first << " --> " << itr->second << "\n";
        }
        cout << endl;


        // read Wannier functions and calculate monopols
        cout << "[+] Read Wannier functions (xsf-files):\n";
        map< int,WannierFunction > vWannMap{ openAllWannierFunctions(vMapping) };
        map< int,WannierFunction > cWannMap{ openAllWannierFunctions(cMapping) };
        cout << "[+] Done reading all files.\n";


        /**
         * Calculate monopoles
         **/
        map< int,Monopole > vMonoMap;
        map< int,Monopole > cMonoMap;

        // read Wannier functions and calculate monopols
        cout << "[+] Calculate Wannier centers ... ";
        for (auto itr = vWannMap.begin(); itr != vWannMap.end(); ++itr) {
            unique_ptr<double[], free_deleter>density{ joinedDensity(itr->second, itr->second, vector<int>{0,0,0}) };
            vMonoMap.insert({ itr->first ,getMonopole(density.get(), itr->second.getMeshgrid()) });
        }

        for (auto itr = cWannMap.begin(); itr != cWannMap.end(); ++itr) {
            unique_ptr<double[], free_deleter>density{ joinedDensity(itr->second, itr->second, vector<int>{0,0,0})};
            cMonoMap.insert({ itr->first ,getMonopole(density.get(), itr->second.getMeshgrid()) });
        }
        cout << "done.\n";

        cout << "Wannier centers for valence WF (without shift):\n";
        for (auto itr = vMonoMap.begin(); itr != vMonoMap.end(); ++itr) {
            cout  << "\t"<< itr->first<< "\t: r=(";
            cout << itr->second.x << ", " << itr->second.y << ", " << itr->second.z << ")" << endl;
        }
        cout << "Wannier centers for conduction WF (without shift):\n";
        for (auto itr = cMonoMap.begin(); itr != cMonoMap.end(); ++itr) {
            cout  << "\t"<< itr->first<< "\t: r=(";
            cout << itr->second.x << ", " << itr->second.y << ", " << itr->second.z << ")" << endl;
        }


        /**
         * shift to home cells
         **/
        cout << "\n[+] Shift valence Wannier functions to home cell...\n";
        rotateToHomeCell(vMonoMap, vWannMap, vMapping);

        // update mapping file
        writeMappingFile(vMappingFile, vMapping);

        cout << "\n[+] Shift conduction Wannier functions to home cell...\n";
        rotateToHomeCell(cMonoMap, cWannMap, cMapping);

        // update mapping file
        writeMappingFile(cMappingFile, cMapping);


        // show monopoles again
        cout << "\n(shifted) Wannier centers for valence WF:\n";
        for (auto itr = vMonoMap.begin(); itr != vMonoMap.end(); ++itr) {
            cout  << "\t"<< itr->first<< "\t: r=(";
            cout << itr->second.x << ", " << itr->second.y << ", " << itr->second.z << ")" << endl;
        }
        cout << "(shifted) Wannier centers for conduction WF:\n";
        for (auto itr = cMonoMap.begin(); itr != cMonoMap.end(); ++itr) {
            cout  << "\t"<< itr->first<< "\t: r=(";
            cout << itr->second.x << ", " << itr->second.y << ", " << itr->second.z << ")" << endl;
        }
        cout << endl;


        /**
         * create POSFILE
         **/
        cout << "[+] Create POSFILE\n";
        vector<vector<double>> unitcell = vWannMap.begin()->second.getUnitcell();

        ofstream file("POSFILE");
        if (!file.is_open()) {
            cerr << "Cannot write POSFILE!" << endl;
            return false;
        }
        file << fixed;
        file << setprecision(12);

        file << " 20 20 20 y y y" << endl;
        file << "model" << endl;
        for (int i=0; i<3; i++) {
            file << " ";
            for (int j=0; j<3; j++) {
                file << unitcell[i][j] << '\t';
            }
            file << endl;
        }
        file << vMonoMap.size() << endl;
        file << "C" << endl;
        for (auto itr = vMonoMap.begin(); itr != vMonoMap.end(); ++itr) {
            file << " " << itr->second.x << "\t" << itr->second.y << "\t" << itr->second.z << endl;
        }
        file << cMonoMap.size() << endl;
        file << "C" << endl;
        for (auto itr = cMonoMap.begin(); itr != cMonoMap.end(); ++itr) {
            file << " " << itr->second.x << "\t" << itr->second.y << "\t" << itr->second.z << endl;
        }
        file.close();

        cout << "[+] Done generating POSFILE.\n\n";
        cout << "[NOTICE]   If you want to use WannierOptics later, you need to make sure\n";
        cout << "[NOTICE]   that the generated POSFILE is in agreement with your tight-binding\n";
        cout << "[NOTICE]   model for valence and conduction band structure!\n";

        return true;

    }
};

bool writePOSFILE(const string filename, const vector<vector<double>>& unitcell, const vector<vector<double>>& pos_c, const vector<vector<double>>& pos_v,
                    vector<bool> const& CRYSTAL_PERIODIC = vector<bool>{true, true, true})
{

    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot write " << filename << endl;
        return false;
    }
    file << fixed;
    file << setprecision(12);

    // dimensions of the supercell
    for (size_t i=0; i<3; i++) {
        if (CRYSTAL_PERIODIC[i]) file << " 50";  // TODO better heuristics
        else file << " 1";
    }

    // periodic boundary conditions
    for (size_t i=0; i<3; i++) {
        if (CRYSTAL_PERIODIC[i]) file << " y";
        else file << " n";
    }
    file << endl;

    // primitive unit cell
    file << "model" << endl;
    for (int i=0; i<3; i++) {
        file << " ";
        for (int j=0; j<3; j++) {
            file << unitcell[i][j] << '\t';
        }
        file << endl;
    }

    // valence Wannier center
    file << pos_v.size() << endl;
    file << "C" << endl;
    for (const auto& pos : pos_v) {
        file << " " << pos[0] << "\t" << pos[1] << "\t" << pos[2] << endl;
    }

    // conduction Wannier center
    file << pos_c.size() << endl;
    file << "C" << endl;
    for (const auto& pos : pos_c) {
        file << " " << pos[0] << "\t" << pos[1] << "\t" << pos[2] << endl;
    }
    file.close();
    return true;
}

#endif // CONFGENERATOR_H