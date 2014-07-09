#include "TGUtils.hpp"
#include "TG.hpp"

using namespace std;

namespace TG {
	int TGSystem::readResiduesDihed(istream &ifs, unsigned int &lineCount) {
		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			TGResidue res; stringstream ss(line); ss >> res.dihedral.name;
			if (std::any_of(POSITIVE_RESIDUES.begin(), POSITIVE_RESIDUES.end(), [&](const string &s) {return s == res.dihedral.name;})) {
				res.dihedral.name.append("+");
			} else if (std::any_of(NEGATIVE_RESIDUES.begin(), NEGATIVE_RESIDUES.end(), [&](const string &s) {return s == res.dihedral.name;})) {
				res.dihedral.name.append("-");
			} else if(res.dihedral.name == "HISE" or res.dihedral.name == "HIS+" == 0) {
			}

			res.name = res.dihedral.name;
			res.num = localCount+1;

			if (molecules.empty()) molecules.push_back( TGMol() );
			molecules.back().residues.push_back( res );
		} return 0;
	}

	int TGSystem::readSSBondDihed(istream &ifs, unsigned int &lineCount) {
		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;
			ssbonds.push_back( TGSSBond::readInTPLDihed(line, lineCount) );
		} return 0;
	}

	int TGSystem::readAnglesDihed(istream &ifs, unsigned int &lineCount) {
		if (molecules.empty()) molecules.push_back( TGMol() );

		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;
			if (localCount < molecules.back().residues.size()) molecules.back().residues.push_back( TGResidue() );
			molecules.back().residues.back().readInTPLDihed(line, lineCount);
		} return 0;
	}

	int TGSystem::loadFromDihed(istream &ifs, bool output_flag) {
		// need to add code to clear data!!!
		ifs.clear(); ifs.seekg(0, ios::beg);
		if (not ifs) { cerr << " ERROR> TGSystem::loadFromDihed\nstream is invalid; exiting!\n"; std::exit(1); }

		vector<string> fileLines;
		std::copy( istream_iterator<FileLine>(ifs), istream_iterator<FileLine>(), std::back_inserter(fileLines) );
		modec = 1;

		// scan file, and flag ss_flg and circ_flg if "PRE>SSBOND" or "PRE>CIRC" exist, if BOTH exist, then exit DIHED_SS_CIRC_ERROR
		bool ssFlag=false, circFlag=false;
		for (const auto &l : fileLines) {
			if (l.compare(0, 10, "PRE>SSBOND") == 0) ssFlag=true;
			if (l.compare(0, 8, "PRE>CIRC") == 0) circFlag=true;
			if (ssFlag and circFlag) {
				cerr << "\n ERROR> tgReadDihedral\n"
				<< "       File Format Error!\n"
				<< "       Two keywords \"PRE>SSBOND\" and \"PRE>CIRC\" were found,\n"
				<< "       Cannot specify both keywords!\n"; std::exit(1);
			}
		}

		// scan file, and upon "PRE>SSBOND", and count number of lines until "PRE>" - NOT NEEDED

		// read into TGSystem
		ifs.clear(); ifs.seekg(0, ios::beg);
		unsigned int lineCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); ) {
			if (line == "PRE>SEQUENCE") readResiduesDihed(ifs, lineCount);
			if (line == "PRE>SSBOND") readSSBondDihed(ifs, lineCount);
			if (line == "PRE>CIRCLE") { if (not molecules.empty()) molecules.back().icrclf = true; }
			if (line == "PRE>DIHEDRAL-ANGLES") { modec = 2; readAnglesDihed(ifs, lineCount); }
		}

		// SSBOND checks
		for (auto &molecule : molecules) {
			if (ssbonds.size() == 0) {
				std::fill(molecule.icnres.begin(), molecule.icnres.end(), TGMol::ISSOFS); continue;
			}

			for (unsigned int i=0; i < molecule.residues.size(); i++) {
				auto &residue = molecule.residues[i];
				for (const auto &ssbond : ssbonds) {

					if (i+1 == ssbond.seqNum1_org or i+1 == ssbond.seqNum2_org) {
						if (residue.name == "CYS") { residue.name += "S"; residue.dihedral.name = residue.name; }
						else {
							cerr << "\n ERROR> tgReadDihedral\n"
								<< "       Residue Number Error!\n"
								<< "       Residue Number in SSBOND is wrong.\n"
								<< "       Residue num  = " << i+1 << "\n"
								<< "       Residue name = " << residue.dihedral.name << "\n";
						}

						if (i+1 == ssbond.seqNum1_org) {
							molecule.icnres[i] = (i+1 > ssbond.seqNum2_org) ? -ssbond.seqNum2_org : ssbond.seqNum2_org;
						} else if (i+1 == ssbond.seqNum2_org) {
							molecule.icnres[i] = (i+1 > ssbond.seqNum1_org) ? -ssbond.seqNum1_org : ssbond.seqNum1_org;
						}
					}
				}
			}
		}

		// Adjust terminal ends
		for (auto &molecule : molecules) {
			if (not molecule.icrclf) {
				// N+ terminal
				auto &frontres = molecule.residues.front();
				if (frontres.dihedral.name == "CYSS" or (frontres.dihedral.name.compare(3, 2, "N+") != 0 and frontres.dihedral.name.compare(4, 2, "N+") != 0)) {
					frontres.dihedral.name += "N+";
					frontres.name = frontres.dihedral.name;
				}

				// C- terminal
				auto &backres = molecule.residues.back();
				if (backres.dihedral.name == "CYSS" or (backres.dihedral.name.compare(3, 2, "C-") != 0 and frontres.dihedral.name.compare(4, 2, "C-") != 0)) {
					backres.dihedral.name += "C-";
					backres.name = backres.dihedral.name;
				}
			} else {
				molecule.circles.push_back( TGCircle(1, molecule.residues.size()) );
			}
		}

		if (output_flag) printAsDihed();
		return 0;
	}

	int TGSystem::printAsDihed() {
		cout << "\n INFORMATION> tgReadDihedral\n             Amino acid Sequence of the protein\n";
		cout.width(7);
		auto &residues = molecules.back().residues;

		for (unsigned int i=0; i < residues.size(); i++) {
			cout << std::left << residues[i].name;
			if (i % RES_DISPLAY == 0) cout << "\n";
		} cout << "\n"; cout.clear();

		if (modec == 2) {
			cout << "\n INFORMATION> tgReadDihedral\n             The following Dihedral Angle are read\n";
			for (unsigned int i=0; i < residues.size(); i++) cout << i << " " << residues[i].dihedral.angleStr() << endl;
		} return 0;
	}
}
