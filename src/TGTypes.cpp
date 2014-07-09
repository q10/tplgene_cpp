#include <sstream>
#include "TGUtils.hpp"
#include "TG.hpp"
#include "cpplinq.hpp"

using namespace std;
namespace LINQ = cpplinq;

namespace TG {
	const array<string, 8> TGSystem::POSITIVE_RESIDUES = array<string, 8> {{ "LYS", "ARG", "MML", "DML", "TML", "SDA", "ADA", "MMA" }};
	const array<string, 2> TGSystem::NEGATIVE_RESIDUES = array<string, 2> {{ "ASP", "GLU" }};
	const int TGMol::ISSOFS = 200000;
	const unsigned int TGSystem::RES_DISPLAY = 10;

	bool TPLGetLine(istream &ifs, string &line, unsigned int &lineCount, unsigned int *bufferpos) {
		if (bufferpos) *bufferpos = ifs.tellg();
		while (std::getline( ifs, line )) { // getline implicitly removes the \n
			lineCount++;
			removeWhitespaces(line);
			if (line[0] == ';') {
				if (bufferpos) *bufferpos = ifs.tellg(); continue; // loop past comment lines
			}
			return true; // else it is a valid line so we should return
		} return false;
	}

	ostream& operator<<(ostream& os, const TGPair& tg) {
		return os << JSON( vector<string> {JSONPair("pair", tg.pair), JSONPair("tri", tg.tri), JSONPair("tetra", tg.tetra)} );
	}
	ostream& operator<<(ostream &os, const TGZmat &tg) { return os; }
	ostream& operator<<(ostream &os, const TGCoord &tg) { return os; }
	ostream& operator<<(ostream &os, const TGAtom &tg) { return os; }
	ostream& operator<<(ostream &os, const TGBond &tg) { return os; }
	ostream& operator<<(ostream &os, const TGAngle &tg) { return os; }
	ostream& operator<<(ostream &os, const TGTParts &tg) { return os; }
	ostream& operator<<(ostream &os, const TGTorsion &tg) { return os; }
	ostream& operator<<(ostream &os, const TGDihedRes &tg) { return os; }
	ostream& operator<<(ostream &os, const TGResidue &tg) { return os; }
	ostream& operator<<(ostream &os, const TGFunction &tg) { return os; }
	ostream& operator<<(ostream &os, const TGNonBond &tg) { return os; }
	ostream& operator<<(ostream &os, const TGMol &tg) { return os; }
	ostream& operator<<(ostream &os, const TGSSBond &tg) { return os; }

/*
	TGAtom TGAtom::readInTPL(const string &line, unsigned int lineNum) {

	}
*/
	TGAtom TGAtom::readInPDB(const string &line, unsigned int lineNum, bool clearWhitespacesForDNA) {
		if (line.length() < 80) {
			cerr << "WARNING: LINE " << lineNum << " has less than 80 characters; ATOM may be parsed incorrectly!\n";
		}
		TGAtom atom;
		atom.type = line.substr(0, 6);

		string temp = line.substr(6, 5);
		atom.num = atoi(temp.c_str());
		atom.atom_name = line.substr(12, 4); removeWhitespaces(atom.atom_name);
		atom.altLoc = line[16];
		atom.residue_name = line.substr(17, 3);
		atom.chainID = line[21];

		temp = line.substr(22, 4);
		atom.residue_num_org = atoi(temp.c_str());
		atom.iCode = line[26];

		temp = line.substr(30, 8);
		atom.coord.x = atof(temp.c_str());
		temp = line.substr(38, 8);
		atom.coord.y = atof(temp.c_str());
		temp = line.substr(46, 8);
		atom.coord.z = atof(temp.c_str());

		atom.chainID_mod = atom.chainID;

		if (clearWhitespacesForDNA) {
			removeWhitespaces(atom.residue_name);
		}

		cerr << "LINE " << lineNum << ": " << atom << endl;
		return atom;

	}

	void TGAtom::fixResNamesForDNA(NAType nucleicAcidType, bool firstAtomInResidue) {
		if (nucleicAcidType == NAType::RNA) {
			if (residue_name == "C") residue_name = "CYT";
			else if (residue_name == "G") residue_name = "GUA";
			else if (residue_name == "A") residue_name = "ADE";
			else if (residue_name == "U") residue_name = "URA";

		} else if (nucleicAcidType == NAType::RNA) {
			if (firstAtomInResidue) {
				if (residue_name == "C") residue_name = "DCY";
				else if (residue_name == "G") residue_name = "DGU";
				else if (residue_name == "A") residue_name = "DAD";
				else if (residue_name == "T") residue_name = "DTH";

			} else {
				if (residue_name == "C") residue_name = "DC";
				else if (residue_name == "G") residue_name = "DG";
				else if (residue_name == "A") residue_name = "DA";
				else if (residue_name == "T") residue_name = "DT";

			}
		}

		if (atom_name.back() == '*') atom_name.back() = '\'';
	}

/*
	TGBond TGBond::readInTPL(const string &line, unsigned int lineNum) {

	}

	TGAngle TGAngle::readInTPL(const string &line, unsigned int lineNum) {

	}

	TGTorsion TGTorsion::readInTPL(const string &line, unsigned int lineNum) {

	}
*/
	int TGResidue::readInTPLDihed(const string &line, unsigned int lineNum) {
		stringstream linestream(line);
		if (not (linestream >> dihedral.DihedAngle[0] >> dihedral.DihedAngle[1] >> dihedral.DihedAngle[2]
				>> dihedral.DihedAngle[3] >> dihedral.DihedAngle[4] >> dihedral.DihedAngle[5] >> dihedral.DihedAngle[6]
				>> dihedral.DihedAngle[7] >> dihedral.DihedAngle[8] >> dihedral.DihedAngle[9])) {
			cerr << "TGResidue::readInDihed failed to read line " << lineNum << "; exiting!!!\n"; std::exit(1);
		} return 0;
	}

	int TGResidue::appendToResidueName(const string &s) {
		name.append(s);
		for (auto &atom : atoms) atom.residue_name.append(s);
		return 0;
	}

	int TGResidue::setResidueName(const string &s) {
		name = s;
		for (auto &atom : atoms) atom.residue_name = s;
		return 0;
	}

	string TGDihedRes::angleStr() {
		char line[81];
		snprintf(line, 80, "%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f", DihedAngle[0], DihedAngle[1], DihedAngle[2],
				DihedAngle[3], DihedAngle[4], DihedAngle[5], DihedAngle[6], DihedAngle[7], DihedAngle[8], DihedAngle[9]);
		return string(line);
	}

/*
	TGFunction TGFunction::readInTPL(const string &line, unsigned int lineNum) {

	}

	TGNonBond TGNonBond::readInTPL(const string &line, unsigned int lineNum) {

	}
*/

	TGSSBond TGSSBond::readInPDB(const string &line, unsigned int lineNum) {
		if (line.length() < 80) {
			cerr << "WARNING: LINE " << lineNum << " has less than 80 characters; SSBOND may be parsed incorrectly!\n";
		}
		stringstream linestream(line);
		TGSSBond ssbond;

		ssbond.chainID1_org = ssbond.chainID1_mod = line[15];
		ssbond.chainID2_org = ssbond.chainID2_mod = line[29];
		ssbond.iCode1 = line[21];
		ssbond.iCode2 = line[35];

		linestream.clear(); linestream.seekg(17, ios::beg);
		linestream >> ssbond.seqNum1_org;

		linestream.clear(); linestream.seekg(31, ios::beg);
		linestream >> ssbond.seqNum2_org;
		return ssbond;
	}

	TGSSBond TGSSBond::readInTPLDihed(const string &line, unsigned int lineNum) {
		stringstream linestream(line);
		TGSSBond ssbond;
		if (not (linestream >> ssbond.seqNum1_org >> ssbond.seqNum2_org)) { cerr << "TGSSBond::readInDihed failed to read line " << lineNum << "; exiting!!!\n"; std::exit(1); }
		ssbond.seqNum1_mod = ssbond.seqNum1_org; ssbond.seqNum2_mod = ssbond.seqNum2_org;
		return ssbond;
	}

	void TGSSBond::sortInternal(const vector<char> &orderedChainIDs) {
		if (seqNum1_org > seqNum2_org) {
			auto tmp = seqNum1_org;
			seqNum1_org = seqNum2_org; seqNum2_org = tmp;

		} else {
			auto index_chainID1 = std::find(orderedChainIDs.begin(), orderedChainIDs.end(), chainID1_mod);
			auto index_chainID2 = std::find(orderedChainIDs.begin(), orderedChainIDs.end(), chainID2_mod);
			if (index_chainID1 == orderedChainIDs.end()) {
				cerr << "\n ERROR> tgReadAminoSequence\n"
						"       ChainID Error!\n"
						"       ChainID \" " << chainID1_mod << " \" specified in the SSBOND line does not exist.\n"
						"       Please check Inputdata.\n";
			} else if (index_chainID2 == orderedChainIDs.end()) {
				cerr << "\n ERROR> tgReadAminoSequence\n"
						"       ChainID Error!\n"
						"       ChainID \" " << chainID2_mod << " \" specified in the SSBOND line does not exist.\n"
						"       Please check Inputdata.\n";
			}
			if (index_chainID1 > index_chainID2) {
				auto tmp = chainID1_mod;
				chainID1_mod = chainID2_mod; chainID2_mod = tmp;
			}
		}
	}

	int TGPair::numInteractions() { return pair_num + tri_num + tetra_num; }

	unsigned int TGMol::numAtoms() const {
		return LINQ::from(residues) >> LINQ::select( [](const TGResidue &res) {return res.atoms.size();} ) >> LINQ::sum();
		/* unsigned int count = 0;
			for (const auto &res : residues) count += res.atoms.size();
			return count; */
		}

	void TGMol::labelCTerminalForDNA() {
		auto &residue = residues.back();
		for (auto &atom : residue.atoms) atom.residue_name += "3*";
		residue.name = residue.atoms.back().residue_name;
	}

	void TGMol::ensureUniqueResidueNumbers() {
		// check if adjacent residues have same resnum, using first atom of each residue
		for (unsigned int i=0; i < residues.size()-1; i++) {
			auto &currentResAtom = residues[i].atoms.front(); auto &nextResAtom = residues[i+1].atoms.front();
			if (currentResAtom.chainID == nextResAtom.chainID and currentResAtom.iCode == nextResAtom.iCode
				and currentResAtom.residue_num_org == nextResAtom.residue_num_org) {
				cerr << "\n ERROR> TGMol::ensureUniqueResidueNumbers\n"
					 << "     The same residue sequence number is defined.\n"
					 << "     Please check pdb file.\n"
					 << "     chainID: " << currentResAtom.chainID << "\n"
					 << "     residue sequence number: " << currentResAtom.residue_num_org << currentResAtom.iCode << "\n"
					 << "     residue names: " <<  currentResAtom.residue_name << nextResAtom.residue_name << "\n\n";
				std::exit(1);
			}
		}
	}

	ResChainIDConversion::ResChainIDConversion(const PDB::Residue &res) {
		chainID_org = res.chainID;
		chainID_mod = res.chainID_mod;
		resSeq_org = res.resSeq;
		resSeq_mod = res.resSeq;
		iCode_org = res.iCode;
		iCode_mod = res.iCode_mod;
		res_name = res.resName;
	}
}



/*

			if(amber->mol_num != sdata.mol_num) {
				if(sdata.input_short != 2 && sdata.method != 2){
					printf("\n");
					printf("\n WARNING> tgReadDNASequence\n");
					printf("          Number of Molecules does not match ");
					printf("to Number of Molecule Names .\n");
					printf("\n %%%% Please input Molecular Name of ");
					printf("Molecule Number %d . %%%%\n", j + 1);
					strcpy(textbuf, "");
					do {
						if(strlen(textbuf) > 40) {
							printf("Molecular name is long . \n");
							printf("Input less than 39 chars . \n");
						}
						scanf("%s", textbuf);
					} while(strlen(textbuf) > 39);
					strcpy(sdata.inpmol[j].mol_name, textbuf);
				}
				else {
					sprintf(sdata.inpmol[j].mol_name, "Nuc%d",j + 1);
				}
			}

			if(j == amber->mol_num - 1)
				break;
			rescount = 0;
			amber->mol = amber->mol->next;
			amber->mol->residue = amber->mol->residuebegin;
			amber->mol->residue->atom = amber->mol->residue->atombegin;
		}
		printf("\n");
	}
*/

