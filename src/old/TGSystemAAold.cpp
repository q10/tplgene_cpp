#include <fstream>
#include <sstream>
#include <array>
#include <numeric>

#include "TGUtils.hpp"
#include "TG.hpp"

using namespace std;

namespace TG {
	int TGSystem::loadFromAA(istream &ifs, bool output_flag) {
		ifstream ifs(pdbfilename, ifstream::in);
		if (not ifs) { cerr << " ERROR> TGSystem::load\n" << "        Cannot open file!"; std::exit(1); }

		// Count the number of lines with SSBOND tag
		unsigned int ssbond_num = std::count_if(istream_iterator<fileLine>(ifs), istream_iterator<fileLine>(), [](const string &tmpline) { return tmpline.compare(0, 6, "SSBOND") == 0; } );

		// Go back to beginning of filestream
		ifs.clear(); ifs.seekg(0, ios::beg);

		chainGroups.clear(); ssbonds.clear(); molecules.clear();
		molecules.push_back( TGMol() );

		// Initialize
		int alternate_indicator = 0;

		//SSbond_num = 0; // probably replaced by ssbonds.size()
		//mol_num = 1; // probably replaced by molecules.size()
		//molecules[0].res_num = 1; // probably replaced by molecules[0].residues.size()
		//molecules[0].residues[0].num = 1;
		//molecules[0].residues[0].atom_num = 0; // probably replaced by molecules[0].residues[0].atoms.size()
		//molecules.atom_num = 0; // probably replaced by molecules[0].numAtoms()

		int inter_ss = 0;
		int circ_flg = 0;
		int ss_circ_flg = 0;
		int circ_num = 0; // this is the index of circ_end and circ_start in molecule class - can be replaced by push_back
		int flag = 0;
		int atom_type = 0;
		int work_lig_flg = 0;
		int work_flg;
		vector<char> orderedChainIDs; // the order of chainIDs
		char previousICode = 0, previousChainID;
		string previousResidueName;
		int lastResidueNum;

		unsigned int lineCount = 0;
		for (string line; PDBGetLine( ifs, line, lineCount ); ) {

			if (line.compare(0, 6, "SSBOND") == 0) {
				ssbonds.push_back( TGSSBond::readInPDB(line, lineCount) );

			} else if (line.compare(0, 6, "CIRCLE") == 0) {
				TGMol &mol = molecules.back();
				mol.icrclf = 1;
				ss_circ_flg = 1;

			} else if (line.compare(0, 6, "HETATM") == 0) {
				atom_type = 2;

			} else if (line.compare(0, 4, "ATOM") == 0) {
				atom_type = 1;
				work_lig_flg = 0;

				// get current molecule
				TGMol &molecule = molecules.back();

				// Get all data from line
				char currentICode = isalpha(line[26]) ? line[26] : 0;
				char currentChainID = isalnum(line[21]) ? line[21] : ' ';
				line[21] = line[26] = ' ';

				int tmp_fileatomcount, tmp_resnumorg, currentResidueNum;
				double tmp_x, tmp_y, tmp_z;
				string tmp_type, tmp_atomname, currentResidueName;
				stringstream linestream(line);
				if (not (linestream >> tmp_type >> tmp_fileatomcount >> tmp_atomname >> currentResidueName >> currentResidueNum >> tmp_x >> tmp_y >> tmp_z)) {
					cerr << "bad line\n"; std::exit(1);
				}

				if ( molecule.residues.empty() or
					(previousResidueName != currentResidueName) or
					(previousICode != currentICode) or
					(lastResidueNum != currentResidueNum) or
					(previousICode == 0 and currentICode == 0 and lastResidueNum != currentResidueNum) ) {

						molecule.residues.push_back( TGResidue() );
						molecule.residues.back().num = molecule.residues.size();
				}

				// If the current residue is not the first, then check the atom type of the first atom of the second-to-last residue
				if (molecule.residues.size() > 1 and
					molecule.residues[ molecule.residues.size()-2 ].atoms[0].type == "ATOM" and
					tmp_type == "HETATM") {
						cerr << "THROWING HETATM ERROR\n";
						std::exit(1);
				}

				// set the correct references
				TGResidue &residue = molecules.back().residues.back();
				residue.atoms.push_back( TGAtom() );
				TGAtom &atom = residue.atoms.back();

				// copy data over to current atom and residue
				atom.type = tmp_type;
				atom.atom_name = tmp_atomname;
				atom.residue_name = currentResidueName;
				atom.residue_num_org = currentResidueNum;
				atom.coord.x = tmp_x; atom.coord.y = tmp_y; atom.coord.z = tmp_z;
				atom.iCode = residue.iCode = currentICode;
				atom.chainID = atom.chainID_mod = residue.chainID = currentChainID;
				atom.residue_num_mod = molecule.residues.size();

				if (currentChainID != ' ' and (orderedChainIDs.size() == 0 or orderedChainIDs.back() != currentChainID))
					orderedChainIDs.push_back(currentChainID);

				// if one of the following residues, append with appropriate charge
				if (std::any_of(POSITIVE_RESIDUES.begin(), POSITIVE_RESIDUES.end(), [&](const string &s) {return s == atom.residue_name;})) {
					atom.residue_name.append("+");
				} else if (std::any_of(NEGATIVE_RESIDUES.begin(), NEGATIVE_RESIDUES.end(), [&](const string &s) {return s == atom.residue_name;})) {
					atom.residue_name.append("-");
				}

				// if the current residue is the first in the chain and is not circular, append N+ terminus
				if ((molecule.residues.size() == 1) and (molecule.icrclf != 1)) {
					atom.residue_name.append("N+");
				}

				// if circle flag is on, then set mol->circ_start
				if (ss_circ_flg == 1) {
					molecule.circ_start[circ_num] = atom.residue_num_mod;
					ss_circ_flg = 2;
				}

				// in this case we also add N+
				if (inter_ss == 1 && ss_circ_flg != 2) {
					atom.residue_name.append("N+");
				}

				// assign name back to parent, and update current temp variables
				residue.name = atom.residue_name;
				previousICode = currentICode;
				previousChainID = currentChainID;
				previousResidueName = currentResidueName;
				lastResidueNum = lastResidueNum;


			} else if (line.compare(0, 3, "TER") == 0 or line.compare(0, 3, "END") == 0) {
				if (work_lig_flg == 0) {

					TGMol &molecule = molecules.back();
					TGResidue &residue = molecule.residues.back();

					if (ss_circ_flg == 2) {
						molecule.circ_end[circ_num] = molecule.residues.size();
						circ_num++;
						ss_circ_flg = 0;

					} else if (ss_circ_flg == 0 and atom_type == 1) { // atom_type 1 == ATOM, 2 == HETATM
						residue.name.append("C-");
						for (auto &atom : residue.atoms) { atom.residue_name.append("C-"); }
					}

					if (alternate_indicator == 0 and atom_type == ATOM) {
						cout << "\n INFORMATION> tgReadAminoSequence\n"
							<< "             Molecule Number           : " << molecules.size() << endl
							<< "             Total number of residues  : " << molecule.residues.size() << endl;
						if (molecule.icrclf == 1) {
							cout << "             This peptide is a circular peptide.\n";
						}
					}
				}


			// Intermolecular SS-Bond Calculation
			} else if (line.compare(0, 4, "INTE") == 0) {
				auto &molecule = molecules.back();
				auto &residue = molecule.residues.back();
				if (ss_circ_flg == 2) {
					molecule.circ_end[circ_num] = residue.num;
					circ_num++;
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

					if (molecule.icrclf == 0 or molecule.circ_end.at(circ_num - 1) != residue.num) {
						atom.residue_name.append("C-");
						residue.name = atom.residue_name;
					}
				}
				inter_ss = 1;
			}
		}

		// copy ssbond info over from args
		TGSSBond &ssbond = ssbond;
		ssbond.chainID1_mod = cID1;
		ssbond.chainID2_mod = cID2;
		ssbond.seqNum1_mod = sNum1_mod;
		ssbond.seqNum2_mod = sNum2_mod;
		// NEED TO CHECK all vectors in ssbond are same length!

		/******************************
		 * Sorting of SSBOND
		 ******************************/
		for (auto &ssbond : ssbonds) ssbond.sortInternal(orderedChainIDs);



		/******************************
		 * Renumber of Residue Number
		 ******************************/
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
				}
				resNum++;
			}
		}


		// sanity checks with chainID
		unsigned int errors = 0;
		for (auto &molecule : molecules) {
			for (auto &residue : molecule.residues) {
				TGAtom *previousAtom = nullptr;
				for (auto &atom : residue.atoms) {
					if (not previousAtom) previousAtom = &atom;

					if (atom.chainID == previousAtom->chainID and atom.iCode == previousAtom->iCode and atom.residue_num_org == previousAtom->residue_num_org) {
						cerr <<"\n ERROR> tgReadAminoSequence\n" << endl
								<< "     The same residue sequence number is defined.\n"
								<< "     Please check pdb file.\n"
								<< "     chainID: " << previousAtom->chainID << endl
								<< "     residue sequence number: " << previousAtom->residue_num_org << previousAtom->iCode << endl
								<< "     residue names: " << previousAtom->residue_name << ", " << atom.residue_name << endl << endl;
						errors++;
					}
					previousAtom = &atom;
				}
			}
		}
		if (errors > 0) { std::exit(1); }


		// process ssbonds
		/******************************
		 * SS-Bond Process
		 ******************************/
/*		unsigned int resnum = 0;

		for (auto &molecule : molecules) {
			// When These Molecules Have No SS-Bond
			if (SSbond_num == 0) {
				molecule.icnres = vector<int>(molecule.residues.size(), ISSOFS);

			} else {
				unsigned int i = 0;
				for (auto &residue : molecule.residues) {
					auto &atom = residue.atoms.front();

					if (resnum != atom.residue_num_mod) {
						vector<unsigned int> all_js(SSbond_num);
						std::iota(all_js.begin(), all_js.end(), 0);
						auto id_compare = [&](unsigned int j) {
							return ( atom.residue_num_mod == ssbond.seqNum1_mod[j] or atom.residue_num_mod == ssbond.seqNum2_mod[j] ) and
									( residue.chainID == ssbond.chainID1_mod[j] or residue.chainID == ssbond.chainID1_org[j] or residue.chainID == ssbond.chainID2_org[j] );
						};
						auto it = std::find_if(all_js.begin(), all_js.end(), id_compare);

						if (it == all_js.end()) {
							molecule.icnres[i] = ISSOFS;

						} else {
							unsigned int j = *it;

							if (atom.residue_name == "CYS") {
								residue.name.append("S");
								for (auto &t_atom : residue.atoms) t_atom.residue_name.append("S");

							} else if (atom.residue_name == "CYSC-") {
								residue.name = "CYSSC-";
								for (auto &t_atom : residue.atoms) t_atom.residue_name = "CYSSC-";

							} else if (atom.residue_name == "CYSN+") {
								residue.name = "CYSSN+";
								for (auto &t_atom : residue.atoms) t_atom.residue_name = "CYSSN+";

							} else if (atom.residue_name == "CYSS") {
								// do nothing

							} else {
								//tgError(SSBOND_WRONG_RES_NUM, atom->residue_name, NULL, NULL, amber->mol_num, atom->residue_num_mod);
							}

							// change seqNum1_org, seqNum2_org to seqNum1_mod, seqNum2_mod
							if (i + 1 == ssbond.seqNum1_mod[j]) {
								if(i + 1 > ssbond.seqNum2_mod[j]) {
									molecule.icnres[i] = -ssbond.seqNum2_mod[j];
								} else {
									molecule.icnres[i] = ssbond.seqNum2_mod[j];
								}
							}
							else if(i + 1 == ssbond.seqNum2_mod[j]){
								if(i + 1 > ssbond.seqNum1_mod[j]){
									molecule.icnres[i] = -ssbond.seqNum1_mod[j];
								} else {
									molecule.icnres[i] = ssbond.seqNum1_mod[j];
								}
							} else {
								molecule.icnres[i] = ISSOFS;
							}
						}
					}
					resnum = atom.residue_num_mod;
					i++;
				}
			}
		}
*/


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
	}
}

/*


	int ReadSSBond(TGSystem *amberSystem, const string &line) {
		stringstream linestream(line);
		TGSSBond &ssbond = ssbond;

		// get chainID1
		ssbond.chainID1_org.push_back( string(&line[15], 1) );
		ssbond.chainID1_mod.push_back( string(&line[15], 1) );

		// get seqNum1
		linestream.clear(); linestream.seekg(17, ios::beg);
		int temp; linestream >> temp;
		ssbond.seqNum1_org.push_back( temp );

		// get iCode1
		ssbond.iCode1.push_back( string(&line[21], 1) );

		// get chainID2
		ssbond.chainID2_org.push_back( string(&line[29], 1) );
		ssbond.chainID2_mod.push_back( string(&line[29], 1) );

		// get seqNum2
		linestream.clear(); linestream.seekg(31, ios::beg);
		int temp2; linestream >> temp2;
		ssbond.seqNum2_org.push_back( temp2 );

		// get iCode2
		ssbond.iCode2.push_back( string(&line[35], 1) );

		SSbond_num++;
		return 0;
	}


*/