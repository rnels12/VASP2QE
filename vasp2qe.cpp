/*
 *  Created by Ryky Nelson July 2019
 */

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include <sstream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <cstring>
#include <cstdlib>
#include "elements.h"

using namespace std;
using namespace Eigen;
using namespace boost::filesystem;

void Usage(string exe_name){    
    cerr << "please specify your POSCAR or CONTCAR;\n"
	 << "e.g.: " << exe_name << " POSCAR\n";
}

void author(){    
    cerr << "Program is created by Ryky Nelson.\n"
	 << "The use of the program is subject to the author's permission!\n";
}

string GetEnv( const string & var ) {
     const char * val = ::getenv( var.c_str() );
     if ( val == 0 ) {
         return "";
     }
     else {
         return val;
     }
}

int main(int argc, char* argv[]){

    directory_iterator end_itr; // default construction yields past-the-end
    string vaspin = "";
    if (argc < 2 ) {
    	try {
    	    for ( directory_iterator itr( "." ) ; itr != end_itr; ++itr ) {
    		string fname = itr->path().filename().generic_string();
		// cerr << "fname = " << fname << endl;
		if (fname.compare("POSCAR") == 0) {
		    vaspin="POSCAR";
		    break;
		}
		else if (fname.compare("CONTCAR") == 0) {
		    vaspin="CONTCAR";
		    break;
		}
    	    }
	    if (vaspin.empty()) throw runtime_error("ERROR: can't find POSCAR or CONTCAR!\n");
    	}
    	catch (const exception& e){
	    cerr <<  e.what();
    	    Usage(argv[0]);
    	    return 1;
    	}

    }
    else { 
	string flag(argv[1]);
	if (flag == "-c" || flag == "-C") {
	    author();
	    return 0;
	}

	vaspin = argv[1]; 
    }

    std::ifstream vaspFile( vaspin.c_str() );
    if( !vaspFile || (vaspFile.peek() == std::ifstream::traits_type::eof()) ) {
    	cerr << "Unable to open " << vaspin << " or it is empty!" << endl;
    	return 0;
    }

    string line; stringstream ss;
    string trash, sval;
    double dval;

    // string nspin("1");
    string NSpin = GetEnv("NSPIN");
    if ( NSpin.empty() ) NSpin = '1';
    ss.str(string()); ss.clear(); ss.str(NSpin);
    int nspin;
    ss >> nspin;

    getline(vaspFile,line); // ignore the header

    getline(vaspFile,line);
    ss.str(string()); ss.clear(); ss.str(line);
    double acell;
    ss >> acell;

    // get the primitive vector from POSCAR 
    Matrix3d primVec = Matrix3d::Zero();
    for (int idummy = 0; idummy < 3; ++idummy) { // the elements' order is kept!!!
	getline(vaspFile,line);
	ss.str(string()); ss.clear(); ss.str(line);
	ss >> primVec(idummy,0) >> primVec(idummy,1) >> primVec(idummy,2);
    }
    double maxEl = primVec.array().abs().maxCoeff();
    acell *= maxEl;
    primVec /= maxEl;

    // get the atoms and position-vectors from POSCAR
    size_t natom(0), ntyp(0);
    string atyp;
    vector<string> atomType;   
    vector<double> magmomType;   
    vector<size_t> numPerType;

    getline(vaspFile,line);

    if (nspin == 1) {
	ss.str(string()); ss.clear(); ss.str(line);
	while(ss>>atyp) atomType.push_back(atyp);
	ntyp = atomType.size();
    }
    else if (nspin == 2) {
    	string magmom = GetEnv("MAGMOM");
	ss.str(string()); ss.clear(); ss.str(magmom);

	size_t tdummy;
	while(ss>>tdummy) {
	    numPerType.push_back(tdummy);
	    double magmomV;
	    ss >> magmomV;
	    magmomType.push_back(magmomV);
	    natom += tdummy;
	}
	ntyp = numPerType.size();	
    }
	

    getline(vaspFile,line);
    if (nspin == 1) {
	ss.str(string()); ss.clear(); ss.str(line);
	for(vector<string>::iterator it = atomType.begin(); it != atomType.end(); ++it) {
	    size_t tdummy; ss >> tdummy;
	    natom += tdummy;
	    numPerType.push_back(tdummy);
	}
    }
    
    getline(vaspFile,line); // kind of coordinates: D or C, for now assuming D
    vector<Vector3d> apos;
    for(size_t idummy = 0; idummy < ntyp; ++idummy) {
	size_t npTyp = numPerType[idummy];
	for(size_t jdummy = 0; jdummy < npTyp; ++jdummy) {
	    getline(vaspFile,line);
	    ss.str(string()); ss.clear(); ss.str(line);
	    Vector3d pos;
	    ss >> pos(0) >> pos(1) >> pos(2) >> atyp;
	    apos.push_back(pos);
	    
	    if (nspin == 2 && jdummy == 0) atomType.push_back(atyp);

	}
    }

    path currDir(boost::filesystem::current_path());
    string cwd = currDir.filename().string();
    string potDir = GetEnv("QEPP");
    if (potDir == "") potDir="./";

    //checking PP files
    vector<string> PPFiles;
    for (vector<string>::iterator it = atomType.begin(); it != atomType.end(); ++it) {
    	string Ename = *it + '.';
    	size_t strSize = Ename.size();
    	    bool found = false;
    	    for ( directory_iterator itr( potDir ) ; itr != end_itr; ++itr ) {
    	    	string fname = itr->path().filename().generic_string();
    	    	// cerr << "fname = " << fname << endl;
    	    	int fsize = fname.size();
    	    	if ( fsize >= strSize && fname.substr(0,strSize).find(Ename) != string::npos ) {
    	    	    string PPFile = itr->path().filename().generic_string();
    	    	    PPFiles.push_back(PPFile);
    	    	    found = true;
    	    	    break;
    	    	}
    	    }
    	    if ( found ) continue;
	    PPFiles.push_back(" ");
    }

    cout << "&control\n" 
	 << "prefix = '" << cwd << "',\n"
	 << "pseudo_dir = './',\n" 
	 << "outdir = './',\n" 
	 << "wf_collect = .true.\n"
	 << "verbosity = 'high'\n"
	 << "/\n";

    string occ = GetEnv("OCCUP");
    if ( occ.empty() ) occ = "smearing";
    string smearing = GetEnv("SMEARING");
    if ( smearing.empty() ) smearing = "gauss";
    if ( occ.compare("tetrahedra") == 0 ) smearing = "";


    // print system part
    cout << "&system\n"
    	 << "ibrav = 0 \nA = " << acell << endl
    	 << "nat = " << natom << "\nntyp = " << ntyp << endl
    	 << "ecutwfc = XXECUTWFCXX" << endl 
	 << "nbnd = XXNBNDXX" << "\nnspin = " << nspin << endl
    	 << "nosym = .TRUE." << endl
    	 << "occupations = '" << occ << "' " << endl;
    if (occ.compare("tetrahedra") != 0) cout << "smearing = '" << smearing << "' " << endl;

    if (nspin == 2) {
    	for (int iTyp = 0; iTyp < ntyp; ++iTyp) {
    	    cout << "starting_magnetization(" << (iTyp+1) << ") = "
		 << magmomType[iTyp] << endl;
    	}
    }
    cout << "/\n";

    // electron part
    cout << "&electrons\n" << "conv_thr =  1.0d-9\n" << "/\n";

    // CELL_PARAMETERS 
    cout << fixed << setprecision(8);
    cout << "CELL_PARAMETERS {alat}\n" << primVec << endl;

    // ATOMIC_SPECIES
    cout << fixed << setprecision(3);
    cout << "ATOMIC_SPECIES\n";
    size_t PPidx(0);
    if (nspin == 1) {
	for (vector<string>::iterator it = atomType.begin(); it != atomType.end(); ++it) {
	    double massOfE = elementsMass.at(*it);
	    cout << *it << ' ' << massOfE << ' ' << PPFiles[PPidx] << endl;
	    ++PPidx;
	}
    }
    else if (nspin == 2) {
	size_t atmIdx(0);
	string prev_atm;
	for (vector<string>::iterator it = atomType.begin(); it != atomType.end(); ++it) {
	    if ( *it != prev_atm ) atmIdx = 0;
	    ++atmIdx;
	    double massOfE = elementsMass.at(*it);
	    cout << *it << atmIdx << ' ' << massOfE << ' ' << PPFiles[PPidx] << endl;
	    ++PPidx;
	    prev_atm = *it;
	}
	
    }

    // print ATOMIC_POSITIONS
    cout << fixed << setprecision(8);
    cout << "ATOMIC_POSITIONS {crystal}\n";
    if (nspin == 1) {
	size_t iterdummy(0);
	for(size_t idummy = 0; idummy < ntyp; ++idummy) {
	    string aname = atomType[idummy];
	    size_t npTyp = numPerType[idummy];
	    for(size_t jdummy = 0; jdummy < npTyp; ++jdummy) {
		cout << ' ' << aname << ' ' << apos[iterdummy].transpose() << endl;
		++iterdummy;
	    }
	}
    }
    else if (nspin == 2) {
	size_t atmIdx(0);
	string prev_atm;
	size_t iterdummy(0);
	for(size_t idummy = 0; idummy < ntyp; ++idummy) {
	    string aname = atomType[idummy];
	    if ( aname != prev_atm ) atmIdx = 0;
	    ++atmIdx;
	    size_t npTyp = numPerType[idummy];
	    for(size_t jdummy = 0; jdummy < npTyp; ++jdummy) {
		cout << ' ' << aname << atmIdx << ' ' << apos[iterdummy].transpose() << endl;
		++iterdummy;
	    }
	    prev_atm = aname;
	}

    }

    // print K_POINTS {automatic}
    cout << "K_POINTS {automatic}\n"
	 << "  XXKPOINTXX 0 0 0\n";

    return 0;

}
