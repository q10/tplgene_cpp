#include "TGUtils.hpp"
#include "TG.hpp"

using namespace std;

namespace TG {
	int TGSystem::loadFromDNA(istream &ifs, bool output_flag) {
		ifs.clear(); ifs.seekg(0, ios::beg);
		if (not ifs) { cerr << " ERROR> TGSystem::loadFromDNA\nstream is invalid; exiting!\n"; std::exit(1); }

		// need to implement clear all data prior!!!

		NAType nucleicAcidType = NAType::NONE;

		// check if U or T (RNA or DNA)
		unsigned int lineCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); ) {
			if (line.compare(0, 4, "ATOM") == 0) {
				TGAtom atom = TGAtom::readInPDB(line, lineCount, true);
				if (atom.residue_name == "U") {
					nucleicAcidType = NAType::RNA; break;
				} else if (atom.residue_name == "T") {
					nucleicAcidType = NAType::DNA; break;
				}
			}
		}

		if (nucleicAcidType == NAType::NONE) { cerr << "ERROR: cannot deduce whether molecule is DNA or RNA; exiting!\n"; std::exit(1); }
		ifs.clear(); ifs.seekg(0, ios::beg);

		lineCount = 0;
		int resNumOffset = -1, lastResNum = -1;
		bool prevTerminate = true, firstAtomInResidue = false;
		for (string line; TPLGetLine( ifs, line, lineCount ); ) {

            /******************************
            0123456789012345678901234567890
            ATOM      1  O5*  +C A   1
            ATOM      2  C5*  +C A   1
            ******************************/
			if (line.compare(0, 4, "ATOM") == 0) {
				// create new molecule if necessary
				if (prevTerminate) { prevTerminate = false; molecules.push_back( TGMol() ); }

				// read in atom
				if (line[18] == '+') line[18] = ' ';
				TGAtom atom = TGAtom::readInPDB(line, lineCount, true);

				// set the resnum offset if not set yet
				if (resNumOffset < 0) resNumOffset = atom.residue_num_mod - 1;


				// if we are the first atom in the residue, set flag
				auto &residues = molecules.back().residues;
				if (residues.empty() or lastResNum != atom.residue_num_org-resNumOffset) firstAtomInResidue = true;

				// fix atom names
				atom.fixResNamesForDNA(nucleicAcidType, firstAtomInResidue);

				// if this is the first atom in residue, push a new residue into molecule first before pushing atom in
				if (firstAtomInResidue) {
					residues.push_back( TGResidue() );
					firstAtomInResidue = false;
				}

				residues.back().atoms.push_back( atom );

				auto &residue = residues.back();
				residue.iCode = residue.atoms.back().iCode;

				// update resinfo and atominfo

			} else if (line.compare(0, 3, "TER") == 0 or line.compare(0, 3, "END") == 0) {
				prevTerminate = true;
				cout << "\n INFORMATION> tgReadDNASequence\n"
					 <<	"             Molecule Number           :" << molecules.size()+1 << "\n"
					 << "             Total number of residues  :" << molecules.back().residues.size() << "\n";
				if (molecules.back().icrclf) cerr << "             \nThis peptide is a circular peptide.\n";
			}
		}

		for (auto &molecule : molecules) {
			molecule.labelCTerminalForDNA();
			molecule.ensureUniqueResidueNumbers();
		}

		if (output_flag) printAsDNA();
		return 0;
	}

	int TGSystem::printAsDNA() {
		cout << "\n INFORMATION> tgReadDNASequence\n             Nucleotide Sequence\n";
		unsigned int molCount = 0;
		for (auto &molecule : molecules) {
			cout << "\n MOLECULE NUMBER : " << ++molCount << "\n\n";
			cout.width(7);
			for (unsigned int j=0; j < molecule.residues.size(); j++) {
				cout << std::left << molecule.residues[j].name;
				if (j % RES_DISPLAY == 0) cout << "\n";
			} cout << "\n"; cout.clear();
		} return 0;
	}
}
