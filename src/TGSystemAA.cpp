#include <sstream>
#include <numeric>
#include "TGUtils.hpp"
#include "TG.hpp"

using namespace std;

namespace TG {
	int TGSystem::loadFromAA(istream &ifs, bool output_flag) {
		ifs.clear(); ifs.seekg(0, ios::beg);
		if (not ifs) { cerr << " ERROR> TGSystem::loadFromAA\nstream is invalid; exiting!\n"; std::exit(1); }

		// need to implement clear all data prior!!!
		chainGroups.clear(); ssbonds.clear(); molecules.clear();

		const int ATOM=1, HETATM=2;
		TGAtom prevAtom;
		int inter_ss=0, ss_circ_flg=0, atom_type=0;
		bool prevTerminate = true;
		vector<char> orderedChainIDs; // the order of chainIDs

		unsigned int lineCount = 0;
		for (string line; PDB::PDBGetLine( ifs, line, lineCount ); ) {
			// create new molecule if necessary
			if (prevTerminate) { prevTerminate = false; molecules.push_back( TGMol() ); }

			auto &molecule = molecules.back();
			auto &residues = molecule.residues;
			auto &circles = molecule.circles;

			if (line.compare(0, 6, "SSBOND") == 0) {
				ssbonds.push_back( TGSSBond::readInPDB(line, lineCount) );

			} else if (line.compare(0, 6, "CIRCLE") == 0) {
				molecule.icrclf = true;
				ss_circ_flg = 1;

			} else if (line.compare(0, 6, "HETATM") == 0) {
				atom_type = HETATM;

			} else if (line.compare(0, 4, "ATOM") == 0) {
				atom_type = ATOM;

				TGAtom atom = TGAtom::readInPDB(line, lineCount);

				// If we need a new residue, push a new one back residue and auto-number it
				if ( residues.empty() or prevAtom.residue_name != atom.residue_name or
					prevAtom.iCode != atom.iCode or prevAtom.residue_num_org != atom.residue_num_org ) {
						residues.push_back( TGResidue() );
						residues.back().num = residues.size();
				}

				// If the current residue is not the first, then check the atom type of the first atom of the second-to-last residue
				if (residues.size() > 1 and residues[ residues.size()-2 ].atoms.front().type == "ATOM  " and atom.type == "HETATM") {
					cerr << "THROWING HETATM ERROR\n";
					std::exit(1);
				}

				// set the correct residue number
				atom.residue_num_mod = residues.size();

				// append charge appropriately
				if (std::any_of(POSITIVE_RESIDUES.begin(), POSITIVE_RESIDUES.end(), [&](const string &s) {return s == atom.residue_name;})) {
					atom.residue_name.append("+");
				} else if (std::any_of(NEGATIVE_RESIDUES.begin(), NEGATIVE_RESIDUES.end(), [&](const string &s) {return s == atom.residue_name;})) {
					atom.residue_name.append("-");
				}

				// if the current residue is the first in the chain and is not circular, append N+ terminus
				if (residues.size() == 1 and not molecule.icrclf) atom.residue_name.append("N+");

				// if circle flag is on, then set mol->circ_start
				if (ss_circ_flg == 1) {
					circles.push_back( TGCircle(residues.size()) );
					ss_circ_flg = 2;
				}

				// in this case we also add N+
				if (inter_ss == 1 && ss_circ_flg != 2) atom.residue_name.append("N+");

				// add atom's chainID to list of chains
				if (atom.chainID != ' ' and (orderedChainIDs.size() == 0 or orderedChainIDs.back() != atom.chainID))
					orderedChainIDs.push_back(atom.chainID);

				// assign name back to parent, and update current temp variables
				auto &residue = residues.back();
				residue.name = atom.residue_name;

				// push atom back and set prevAtom to be current atom
				residue.atoms.push_back( atom );
				prevAtom = atom;

			} else if (line.compare(0, 3, "TER") == 0 or line.compare(0, 3, "END") == 0) {
				auto &residue = residues.back();

				if (ss_circ_flg == 2) {
					circles.back().end = residues.size();
					ss_circ_flg = 0;

				} else if (ss_circ_flg == 0 and atom_type == ATOM) {
					residue.appendToResidueName("C-");
				}

				if (atom_type == ATOM) {
					cout << "\n INFORMATION> tgReadAminoSequence\n"
						<< "             Molecule Number           : " << molecules.size() << endl
						<< "             Total number of residues  : " << molecule.residues.size() << endl;
					if (molecule.icrclf) cout << "             This peptide is a circular peptide.\n";
				}
				prevTerminate = true;

			} else if (line.compare(0, 4, "INTE") == 0) { // Intermolecular SS-Bond Calculation
				auto &residue = residues.back();

				if (ss_circ_flg == 2) {
					circles.back().end = residues.size();
					ss_circ_flg = 0;
				}

				array<string, 2> positive_residues {{ "LYS", "ARG" }};
				array<string, 2> negative_residues {{ "ASP", "GLU" }};

				// iterate through all atoms in current residue (not future)
				for (auto &atom : residue.atoms) {
					if (std::any_of(positive_residues.begin(), positive_residues.end(), [&](const string &s) {return s == atom.residue_name;})) {
						atom.residue_name.append("+");
					} else if (std::any_of(negative_residues.begin(), negative_residues.end(), [&](const string &s) {return s == atom.residue_name;})) {
						atom.residue_name.append("-");
					}

					if (not molecule.icrclf or circles.size() < 1 or circles.back().end != residue.num) {
						atom.residue_name.append("C-");
						residue.name = atom.residue_name;
					}
				}
				inter_ss = 1;
			}
		}

		for (auto &molecule : molecules) molecule.ensureUniqueResidueNumbers();

		// copy ssbond info over from args
		// TGSSBond &ssbond = ssbond;
		// ssbond.chainID1_mod = cID1;
		// ssbond.chainID2_mod = cID2;
		// ssbond.seqNum1_mod = sNum1_mod;
		// ssbond.seqNum2_mod = sNum2_mod;

		// sort SSBONDs
		for (auto &ssbond : ssbonds) ssbond.sortInternal(orderedChainIDs);

		// renumber residue numbers for ssbond info
		for (auto &molecule : molecules) {
			unsigned int resNum = 1;
			for (auto &residue : molecule.residues) {
				for (auto &atom : residue.atoms) {
					for (auto ssbond : ssbonds) {
						if (ssbond.chainID1_mod == residue.chainID or ssbond.chainID2_mod == residue.chainID) {
							if (ssbond.seqNum1_org == atom.residue_num_org) ssbond.seqNum1_org = resNum;
							if (ssbond.seqNum2_org == atom.residue_num_org) ssbond.seqNum2_org = resNum;
						}
					}
				} resNum++;
			}
		}

		processSSBondsAA();
		//if (output_flag)
		return 0;
	}

	int TGSystem::processSSBondsAA() {
		for (auto &molecule : molecules) {
			if (ssbonds.size() == 0) {
				std::fill(molecule.icnres.begin(), molecule.icnres.end(), TGMol::ISSOFS); continue;
			}

			for (unsigned int i=0; i < molecule.residues.size(); i++) {
				auto &residue = molecule.residues[i];
				auto &atom = residue.atoms.front();
				auto ssbond_loc = std::find_if(ssbonds.begin(), ssbonds.end(),
					[&](TGSSBond &ssbond) { return atom.residue_num_mod == ssbond.seqNum1_mod or atom.residue_num_mod == ssbond.seqNum2_mod; }
				);

				if (ssbond_loc == ssbonds.end()) {
					molecule.icnres[i] = TGMol::ISSOFS;

				} else {
					auto &ssbond = *ssbond_loc;
					if (residue.chainID == ssbond.chainID1_mod or residue.chainID == ssbond.chainID1_org or residue.chainID == ssbond.chainID2_org) {
						if (atom.residue_name == "CYS") {
							residue.appendToResidueName("S");

						} else if (atom.residue_name == "CYSC-") {
							residue.setResidueName("CYSSC-");

						} else if (atom.residue_name == "CYSN+") {
							residue.setResidueName("CYSSN+");

						} else if (atom.residue_name == "CYSS") {

						} else {
							cerr << "\n ERROR> tgReadAminoSequence\n"
								 << "       Residue Number Error!\n"
								 << "       Residue Number in SSBOND is wrong.\n"
								 << "       Mol num     = " << mol_num << "\n"
								 << "       Residue num = " << atom.residue_num << "\n"
								 << "       Residue name= " << residue.name << "\n";
							std::exit(1);
						}

						// change seqNum1_org, seqNum2_org to seqNum1_mod, seqNum2_mod
						if (i+1 == ssbond.seqNum1_mod) {
							molecule.icnres[i] = (i+1 > ssbond.seqNum2_mod) ? -ssbond.seqNum2_mod : ssbond.seqNum2_mod;
						} else if (i+1 == ssbond.seqNum2_mod){
							molecule.icnres[i] = (i+1 > ssbond.seqNum1_mod) ? -ssbond.seqNum1_mod : ssbond.seqNum1_mod;
						} else {
							molecule.icnres[i] = TGMol::ISSOFS;
						}
					}
				}
			}
		} return 0;
	}

	int TGSystem::printAsAA() {
		cout << "\n INFORMATION> tgReadAminoSequence\n" << "             Amino acid Sequence of the protein\n";
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









/*
		// review/output what we have
		cout << "\n INFORMATION> tgReadAminoSequence\n" << "             Amino acid Sequence of the protein\n";
		unsigned int mol_count = 0;
		for (auto &molecule : molecules) {
			cout << "\n Molecule Number : " << mol_count + 1 << "\n ";

			unsigned int res_count = 0;
			for (auto &residue : molecule.residues) {
				printf("%-7s", residue.atoms.front().residue_name.c_str());
				if (++res_count % RES_DISPLAY == 0) cout << "\n ";
			} cout << endl;

			if (molecules.size() != sdata.mol_num) {
				if (sdata.input_short != 2 && sdata.method != 2) {
					cerr << "\n\n WARNING> tgReadAminoSequence\n"
						<< "          Number of Molecules does not match to Number of Molecule Names.\n\n"
						<< " %%%% Please input Molecular Name of Molecule Number " << mol_count + 1 << " . %%%%\n";

					string input;
					while (input.size() > 39) {
						getline(cin, input);
						if (input.size() > 39) cerr << "MOLECULAR NAME IS LONG . \nINPUT LESS THAN 39 CHARS . \n";
					}
					sdata.inpmol[mol_count].mol_name = input;

				} else {
					sdata.inpmol[mol_count].mol_name = string("Pro") + std::to_string(mol_count + 1);
				}
			}
			mol_count++;
		}
*/

