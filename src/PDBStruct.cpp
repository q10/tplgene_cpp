#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include "TGUtils.hpp"
#include "PDB.hpp"

using namespace std;
using namespace TG;

namespace PDB {
	ostream& operator<<(ostream &os, const PDBStruct &p) {
		return os;
	}

	PDBStruct PDBStruct::readIn(const string &filenamepath, ChainSpecies chainType) {
		ifstream ifs(filenamepath, ifstream::in);
		PDBStruct pdb;

		bool circleFlag = false;
		unsigned int lineCount = 0;
		for (string line; PDBGetLine( ifs, line, lineCount ); ) {

			if (line.compare(0, 4, "ATOM") == 0 or line.compare(0, 6, "HETATM") == 0) {
				pdb.addAtom(line, chainType, lineCount);

				if (circleFlag) {
					pdb.residues.back().circ_flg = true;
					circleFlag = false;
				}

			} else if (line.compare(0, 4, "LINK") == 0) {
				Link link = Link::readIn(line, lineCount);
				if (chainType == ChainSpecies::PEPTIDE and (link.name1 != " N  " or link.name1 != " C  ")) {
					pdb.links.push_back( Link() ); continue;

				} else if (chainType == ChainSpecies::NUCLEOTIDE and (link.name1 != " P  " or link.name1 != " O3'")) {
					pdb.links.push_back( Link() ); continue;
				}
				pdb.links.push_back( link );

			} else if (line.compare(0, 6, "SSBOND") == 0) {
				pdb.ssbonds.push_back( SSBond::readIn(line, lineCount) );

			} else if (line.compare(0, 4, "CIRC") == 0) {
				circleFlag = true;
			}

		} ifs.close(); return pdb;
	}

	void PDBStruct::addAtom(const string &line, ChainSpecies chainType, unsigned int lineCount) {
		Atom t_atom = Atom::readIn(line, lineCount);

		// append new residue if necessary
		if (residues.size() == 0) {
			residues.push_back( Residue() );
		} else {
			auto &prevAtom = residues.back().atoms.back();
			if (t_atom.resName != prevAtom.resName or t_atom.resSeq != prevAtom.resSeq or t_atom.iCode != prevAtom.iCode) {
				residues.push_back( Residue() );
			}
		}
		auto &residue = residues.back();
		residue.atoms.push_back( t_atom );
		auto &atom = residue.atoms.back();

		residue.resName = atom.resName;
		residue.chainID = atom.chainID;
		residue.resSeq = atom.resSeq;
		residue.iCode = atom.iCode;
		residue.chainID_mod = atom.chainID;

		availableChainIDs.remove_if( [&](const char &cID) {return cID == residue.chainID;} );

		if (chainType == ChainSpecies::PEPTIDE) {
			if (atom.atom_name == "N") {
				residue.atom_name_Nterm = atom.atom_name;
				residue.coord_Nterm = atom.coord;
				residue.atom_num_Nterm = atom.atom_serial_num;

			} else if (atom.atom_name == "C") {
				residue.atom_name_Cterm = atom.atom_name;
				residue.coord_Cterm = atom.coord;
				residue.atom_num_Cterm = atom.atom_serial_num;
			}
		} else if (chainType == ChainSpecies::NUCLEOTIDE) {
			if (atom.atom_name == "P") {
				residue.atom_name_Nterm = atom.atom_name;
				residue.coord_Nterm = atom.coord;
				residue.atom_num_Nterm = atom.atom_serial_num;

			} else if (atom.atom_name == "O3") {
				residue.atom_name_Cterm = atom.atom_name;
				residue.coord_Cterm = atom.coord;
				residue.atom_num_Cterm = atom.atom_serial_num;
			}
		}
	}

	void PDBStruct::setToDeleteAltLocs() {
		for (auto &residue : residues) {
			for (unsigned int i=0; i < residue.atoms.size()-1; i++) {
				auto &atom = residue.atoms[i];

				if (not isspace(atom.altLoc)) {
					for (unsigned int j=i+1; j < residue.atoms.size(); j++) {
						auto &other_atom = residue.atoms[j];
							if (atom.atom_name != other_atom.atom_name or atom.resName != other_atom.resName or
								atom.resSeq != other_atom.resSeq or atom.iCode != other_atom.iCode) break;
							other_atom.del_flg = true;
					}
				}
			}
		}
	}

	void PDBStruct::setModifiedChainIDSsbond() {
		for (auto &residue : residues) {
			if (residue.resName != "CYS") continue;

			for (auto &ssbond : ssbonds) {
				if (residue.chainID == ssbond.chainID1 and residue.resSeq == ssbond.seqNum1 and residue.iCode == ssbond.iCode1)
					ssbond.chainID1_mod = residue.chainID_mod;
				if (residue.chainID == ssbond.chainID2 and residue.resSeq == ssbond.seqNum2 and residue.iCode == ssbond.iCode2)
					ssbond.chainID2_mod = residue.chainID_mod;
			}
		}
	}

	void PDBStruct::setModifiedChainID(unsigned int startingFrom) {
		if (availableChainIDs.empty()) { cerr << "ERROR> setModifiedChainID : NO MORE AVAILABLE CHAIN IDs!\n"; std::exit(1); }

		char newChainID = availableChainIDs.front();
		availableChainIDs.pop_front();

		for (unsigned int i=startingFrom; i < residues.size(); i++) {
			auto &residue = residues[i];
			residue.chainID_mod = newChainID;
			for (auto &atom : residue.atoms) atom.chainID_mod = residue.chainID_mod;
		}
	}

	void PDBStruct::setSameChainID(unsigned int currentIndex, ChainSpecies chainType) {
		if (currentIndex >= residues.size()-1) { cerr << "setSameChainID: currentIndex too large\n"; return; }
		auto &constant_residue = residues[ currentIndex ];

		for (unsigned int i=currentIndex; i < residues.size()-1; i++) {
			auto &residue = residues[i];
			auto &nextResidue = residues[i+1];

			if (residue.chainID != nextResidue.chainID) break;

			double dist = vectorDistance(residue.coord_Cterm, nextResidue.coord_Nterm);
			if ((chainType == ChainSpecies::PEPTIDE && dist > 2.0) or (chainType == ChainSpecies::NUCLEOTIDE and dist > 3.0))
				break;

			nextResidue.chainID_mod = constant_residue.chainID_mod;
			for (auto &atom : nextResidue.atoms) atom.chainID_mod = constant_residue.chainID_mod;
		}
	}

	void PDBStruct::tgModify4LackedRes(ChainSpecies chainType) {
		for (unsigned int i=0; i < residues.size()-1; i++) {
			auto &residue = residues[i];
			auto &nextResidue = residues[i+1];

			// If residue has no chainID, assign a new one (based on list of available unique chain IDs)
			if (residue.chainID == ' ' and residue.chainID_mod == ' ') {
				if (availableChainIDs.empty()) { cerr << "NO MORE AVAILABLE CHAIN IDs!\n"; std::exit(1); }
				residue.chainID_mod = availableChainIDs.front();
				for (auto &atom : residue.atoms) atom.chainID_mod = residue.chainID_mod;
				availableChainIDs.pop_front();
				break;
			}

			// check on the metals
			ResAtomSpecies ret_cod = CheckAtomNameOfMetals(residue.atoms[0].atom_name, residue.atoms[0].residue_name, &AAResiduesList, &MEResiduesList);

			// 1: AA or NA residues
			if (ret_cod == 0) {

				if (residue.chainID != nextResidue.chainID) continue;

				double dist = vectorDistance(residue.coord_Cterm, nextResidue.coord_Nterm);
				if ((chainType == ChainSpecies::PEPTIDE and dist > 2.0) or (chainType == ChainSpecies::NUCLEOTIDE and dist > 3.0)) {
					puts("\n WARNING> tgModify4LackedRes");
					if (chainType == ChainSpecies::PEPTIDE) {
						printf("          It is too away between N and C atoms in same chain.");
						printf("(%.1lf angstrom)\n", dist);

					} else if (chainType == ChainSpecies::NUCLEOTIDE) {
						printf("          It is too away between P and O3 atoms in same chain.");
						printf("(%.1lf angstrom)\n", dist);
					}
					puts("          It is regarded as different chain.");
					printf("ATOM  %5d  %-3s %-4s%c%4d    %8.3lf%8.3lf%8.3lf\n", residue.atom_num_Cterm, residue.atom_name_Cterm.c_str(), residue.resName.c_str(),
						residue.chainID, residue.resSeq, residue.coord_Cterm[0], residue.coord_Cterm[1], residue.coord_Cterm[2]);
					printf("ATOM  %5d  %-3s %-4s%c%4d    %8.3lf%8.3lf%8.3lf\n", nextResidue.atom_num_Nterm, nextResidue.atom_name_Nterm.c_str(), nextResidue.resName.c_str(),
						nextResidue.chainID, nextResidue.resSeq, nextResidue.coord_Nterm[0], nextResidue.coord_Nterm[1], nextResidue.coord_Nterm[2]);

					setModifiedChainID(i);
				} else {
					setSameChainID(i, chainType);
				}

			} else {
				setSameChainID(i, chainType);
			}
		}
		setModifiedChainIDSsbond();
	}

	void PDBStruct::checkSsbondAndResInfo() {
		for (const auto &ssbond : ssbonds) {
			bool hit_flg1 = false, hit_flg2 = false;

			for (const auto &residue : residues) {
				if(ssbond.resName1.compare(0, 4, residue.resName) == 0 and ssbond.chainID1 == residue.chainID and
					ssbond.seqNum1 == residue.resSeq and ssbond.iCode1 == residue.iCode) {
					hit_flg1 = true;
				}
				if(ssbond.resName2.compare(0, 4, residue.resName) == 0 and ssbond.chainID2 == residue.chainID and
					ssbond.seqNum2 == residue.resSeq and ssbond.iCode2 == residue.iCode) {
					hit_flg2 = true;
				}
				if (hit_flg1 and hit_flg2) break;
			}

			if (not (hit_flg1 and hit_flg2)) {
				printf(" ERROR> checkSsbondAndResInfo\n");
				puts  ("        SSBOND line is wrong.");
				puts  ("        Please check SSBOND line and related ATOM lines.");
				printf("SSBOND %3d %3s %1c %4d%1c   %3s %1c %4d%1c\n",
					ssbond.serNum, ssbond.resName1.c_str(), ssbond.chainID1, ssbond.seqNum1, ssbond.iCode1,
					ssbond.resName2.c_str(), ssbond.chainID2, ssbond.seqNum2, ssbond.iCode2);
				exit(1);
			}
		}
	}

	void PDBStruct::autoSSDetection() {
    	double ssbond_thr2 = SSBOND_THRESHOLD * SSBOND_THRESHOLD;
    	double s_metal_thr2 = S_METAL_DIST_THRESHOLD * S_METAL_DIST_THRESHOLD;


		vector<char> cys_ChainID, cys_ChainID_mod, cys_iCode;
		vector<int> cys_SeqNum;
		vector< vector<double> > cys_SgCoord, metal_coord;
		vector<string> metal_atom_name;

		// Collect relevant info from SG(SYS) and HETATM atoms
		for (auto &residue : residues) {
			if (residue.resName.compare(0, 3, "CYS")) {
				for (auto &atom : residue.atoms) {
					if (atom.atom_name == "SG") {
						cys_SeqNum.push_back(atom.resSeq);
						cys_ChainID.push_back(atom.chainID);
						cys_ChainID_mod.push_back(atom.chainID_mod);
						cys_iCode.push_back(atom.iCode);
						cys_SgCoord.push_back(atom.coord);
					}
				}
			} else {
				auto &atom = residue.atoms.front();
				if(atom.record_name.compare(0, 6, "HETATM") == 0) {
					if (CheckAtomNameOfMetals(atom.atom_name, atom.resName, "ZZZ") == ResAtomSpecies::ION) {
						metal_atom_name.push_back(atom.atom_name.substr(0, 4));
						metal_coord.push_back(atom.coord);
					}
				}
			}
		}

		// Flag SG(CYS)-HETATM pairs whose distance is below the threshold;
		vector<bool> cys_valid(cys_SgCoord.size(), true);
		bool print_metal_header = true;
		unsigned int sg_i=0;
		for (auto &sg : cys_SgCoord) {
			unsigned int mt_j=0;
			for (auto &mt : metal_coord) {
				if( vectorDistance(sg, mt) < s_metal_thr2) {
					char buf[81]; sprintf(buf, "%s %8.3lf %8.3lf %8.3lf", metal_atom_name[mt_j].c_str(), mt[0], mt[1], mt[2]);
					//tgError(SSBOND_AUTODETECT_ERROR, buf, NULL, NULL, DMY_INT, DMY_INT);
					cys_valid[sg_i] = false; // invalid

					if (print_metal_header) {
						cout << "\n INFORMATION> tgCheckSS\n             The distance of CYS and metal atom is checked\n";
						print_metal_header = false;
					}
					printf("\n The following atom was found near the SG(CYS) atom.\n");
					printf(" %s\n", buf);
					printf(" SG(CYS) Information\n");
					printf("      chain ID       : %c\n", cys_ChainID[sg_i]);
					printf("      residue number : %d\n", cys_SeqNum[sg_i]);
				} mt_j++;
			} sg_i++;
		}

		// For SG(CYS)-SG(CYS) pairs whose distance is below the threshold and are unflagged, add corresponding ssbond info if not already in ssbonds
		for (unsigned int i=0; i < cys_SgCoord.size()-1; i++) {
			for (unsigned int j=i+1; j < cys_SgCoord.size(); j++) {
				double distSquared = squaredVectorDistance(cys_SgCoord[i], cys_SgCoord[j]);
				if(cys_valid[i] and distSquared < ssbond_thr2) {
					auto ssbonds_iter = std::find_if(ssbonds.begin(), ssbonds.end(), [&](const SSBond& ss) {
						return (cys_SeqNum[i] == ss.seqNum1 and cys_ChainID[i] == ss.chainID1 and cys_iCode[i] == ss.iCode1) or
								(cys_SeqNum[j] == ss.seqNum2 and cys_ChainID[j] == ss.chainID2 and cys_iCode[j] == ss.iCode2) or
								(cys_SeqNum[i] == ss.seqNum2 and cys_ChainID[i] == ss.chainID2 and cys_iCode[i] == ss.iCode2) or
								(cys_SeqNum[j] == ss.seqNum1 and cys_ChainID[j] == ss.chainID1 and cys_iCode[j] == ss.iCode1);
					});
					if (ssbonds_iter == ssbonds.end()) {
						SSBond ss;
						ss.iCode1 = cys_iCode[i]; ss.iCode2 = cys_iCode[j];
						ss.chainID1 = cys_ChainID[i]; ss.chainID2 = cys_ChainID[j];
						ss.chainID1_mod = cys_ChainID_mod[i]; ss.chainID2_mod = cys_ChainID_mod[j];
						ss.seqNum1 = cys_SeqNum[i]; ss.seqNum2 = cys_SeqNum[j];
						ss.resName1 = ss.resName2 = "CYS";
						ss.serNum = ssbonds.size() + 2;
						ss.length = sqrt(distSquared);
						ssbonds.push_back(ss);
					}
				}
			}
		}
	}

	bool PDBStruct::writeToFile(const string &filenamepath) {
		ofstream ofs(filenamepath, ofstream::out);
		if (not ofs) { cerr << "Cannot open file " << filenamepath << "!"; return false; }

		/* link */

		unsigned int ss_count=0;
		for (auto &ssbond : ssbonds) {
			char c_line[81];
			snprintf(c_line, 80, "SSBOND %3d CYS %1c %4d%1c   CYS %1c %4d%1c\n", ++ss_count,
				ssbond.chainID1_mod, ssbond.seqNum1, ssbond.iCode1,
				ssbond.chainID2_mod, ssbond.seqNum2, ssbond.iCode2);
			string line(c_line); line.resize(80, ' ');
			ofs << line << endl;
		}

		char cID_bef=0;
		string res_bef;
		bool atFirstAtomOFFirstRes = true;
		unsigned int atom_count = 0;
		for (auto &residue : residues) {
			for (auto &atom : residue.atoms) {

				if(++atom_count > 99999) atom_count = 0;

				// determine output format
				string line_format = (atom.atom_name.length() > 3) ? "%6s%5d %4s %-4s%1c%4d%1c   %8.3lf%8.3lf%8.3lf\n" : "%6s%5d  %-3s %-4s%1c%4d%1c   %8.3lf%8.3lf%8.3lf\n";

				// add TER and CIRCs
				if (cID_bef != 0 and cID_bef != atom.chainID_mod) {
					ofs << "TER\n";
					cID_bef = atom.chainID_mod;
					if (residue.circ_flg) ofs << "CIRC\n";

				} else if(atom.record_name == "HETATM" and res_bef != atom.resName) {
					ofs << "TER\n";
					if (residue.circ_flg) ofs << "CIRC\n";

				} else if (atFirstAtomOFFirstRes) {
					atFirstAtomOFFirstRes = false;
					if (residue.circ_flg) ofs << "CIRC\n";
				}
				res_bef = atom.resName;

				char c_line[81];
				snprintf(c_line, 80, line_format.c_str(), atom.record_name.c_str(), atom_count, atom.atom_name.c_str(), atom.resName.c_str(),
						atom.chainID_mod, atom.resSeq, atom.iCode, atom.coord[0], atom.coord[1], atom.coord[2]);
				string line(c_line); line.resize(80, ' ');
				ofs << line << endl;
			}
		}
		ofs << "TER\n"; ofs.close();
	}

	int PDBStruct::cleanPDBStruct(SDATA &sdata, vector<ResChainIDConv> &residueConversionTable) {
		PDBFile::checkNonStandardResidues(sdata.coordFile);

		// Read in PDB
		PDBStruct pdb = PDBStruct::readIn("/Users/bm/BPWORK/sample.pdb", sdata.chainType);

		// Flag ALTLOC atoms for deletion
		pdb.setToDeleteAltLocs();

		// Set/correct chainIDs based on distances between the C- and N-term atoms of each successive pair of residues
		pdb.tgModify4LackedRes(sdata.chainType);

	    if (sdata.chainType == ChainSpecies::PEPTIDE) {
	    	// Checks SSBond info corresponds to residue info
	    	pdb.checkSsbondAndResInfo();

	    	// If flag is activated, automatically determine undocumented SSBonds by checking SG(CYS) and HETATM atoms' pairwise distances
	        if (sdata.autoSSDetect) pdb.autoSSDetection();
	    }

	    residueConversionTable.clear();
	    for (auto &res : residues) {
			residueConversionTable.push_back( ResChainIDConv(res) );
		}

		sdata.coordFile = "changed_file0";
		pdb.writeToFile(sdata.coordFile);
	}
}
