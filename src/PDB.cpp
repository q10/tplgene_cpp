#include <cstdio>
#include <fstream>
#include "TGUtils.hpp"
#include "PDB.hpp"

using namespace std;
using namespace PDB;

namespace PDB {
	ostream& operator<<(ostream &os, const Coord &p) {
		return os;
	}

	ostream& operator<<(ostream &os, const Atom &p) {
		return os << JSON( vector<string> {JSONPair("record_name", p.record_name), JSONPair("atom_serial_num", p.atom_serial_num), JSONPair("atom_name", p.atom_name),
					JSONPair("altLoc", p.altLoc), JSONPair("resName", p.resName), JSONPair("chainID", p.chainID), JSONPair("resSeq", p.resSeq),
					JSONPair("iCode", p.iCode), JSONPair("coord", p.coord)} );
	}

	ostream& operator<<(ostream &os, const Link &p) {
		return os << JSON( vector<string> {JSONPair("record_name", p.record_name), JSONPair("name1", p.name1), JSONPair("altLoc1", p.altLoc1),
					JSONPair("resName1", p.resName1), JSONPair("chainID1", p.chainID1), JSONPair("resSeq1", p.resSeq1), JSONPair("iCode1", p.iCode1),
					JSONPair("name2", p.name1), JSONPair("altLoc2", p.altLoc1), JSONPair("resName2", p.resName1), JSONPair("chainID2", p.chainID1),
					JSONPair("resSeq2", p.resSeq1), JSONPair("iCode2", p.iCode1), JSONPair("sym1", p.sym1), JSONPair("sym2", p.sym2),
					JSONPair("length", p.length)} );
	}

	ostream& operator<<(ostream &os, const SSBond &p) {
		return os << JSON( vector<string> {JSONPair("record_name", p.record_name), JSONPair("serNum", p.serNum), JSONPair("resName1", p.resName1),
					JSONPair("chainID1", p.chainID1), JSONPair("seqNum1", p.seqNum1), JSONPair("iCode1", p.iCode1), JSONPair("resName2", p.resName2),
					JSONPair("chainID2", p.chainID2), JSONPair("seqNum2", p.seqNum2), JSONPair("iCode2", p.iCode2),
					JSONPair("sym1", p.sym1), JSONPair("sym2", p.sym2), JSONPair("length", p.length)} );
	}

	ostream& operator<<(ostream &os, const Residue &p) {
		return os;
	}

	ostream& operator<<(ostream &os, const PDBFile &p) {
		return os;
	}

	bool PDBGetLine(istream &ifs, string &line, unsigned int &lineCount) {
		while (std::getline( ifs, line )) { // getline implicitly removes the \n
			lineCount++;
			return true; // else it is a valid line so we should return
		} return false;
	}

	ResAtomSpecies isMatchResName(const string &residueName, vector< vector<string> > &AAResiduesTable, vector< vector<string> > &MEResiduesTable) {
		if (residueName.empty()) return ResAtomSpecies::NOMATCH;
		for (const auto &lst : AAResiduesTable) {
			if (std::any_of(lst.begin(), lst.end(), [&](const string &r) {
				return residueName == r or ((residueName[3] == '+' or residueName[3] == '-') and residueName.compare(0, 3, r) == 0);
			})) return ResAtomSpecies::AA;
		}
		for (const auto &lst : MEResiduesTable) {
			if (std::any_of(lst.begin(), lst.end(), [&](const string &r) {return r == residueName;})) return ResAtomSpecies::METALS_WATER;
		}
		return ResAtomSpecies::LIGAND;
	}

	/*
	http://www.wwpdb.org/documentation/format23/sect6.html
	COLUMNS      DATA TYPE        FIELD        DEFINITION
	-----------------------------------------------------------------
	 1 - 6       Record name      "LINK   "
	13 - 16      Atom             name1        Atom name.
	17           Character        altLoc1      Alternate location indicator.
	18 - 20      Residue name     resName1     Residue name.
	22           Character        chainID1     Chain identifier.
	23 - 26      Integer          resSeq1      Residue sequence number.
	27           AChar            iCode1       Insertion code.
	43 - 46      Atom             name2        Atom name.
	47           Character        altLoc2      Alternate location indicator.
	48 - 50      Residue name     resName2     Residue name.
	52           Character        chainID2     Chain identifier.
	53 - 56      Integer          resSeq2      Residue sequence number.
	57           AChar            iCode2       Insertion code.
	60 - 65      SymOP            sym1         Symmetry operator for 1st atom.
	67 - 72      SymOP            sym2         Symmetry operator for 2nd atom.

	74 - 78                       length       link distance (PDB format extension)
	*/
	Link Link::readIn(const string &line, unsigned int lineNum) {
		if (line.length() < 80) {
			cerr << "WARNING: LINE " << lineNum << " has less than 80 characters; LINK may be parsed incorrectly!\n";
		}
		Link link;
		link.record_name = line.substr(0, 6);
		link.name1 = line.substr(12, 4);
		link.altLoc1 = line[16];
		link.resName1 = line.substr(17, 3);
		link.chainID1 = line[21];

		string temp = line.substr(22, 4);
		link.resSeq1 = atoi(temp.c_str());

		link.iCode1 = line[26];
		link.name2 = line.substr(42, 4);
		link.altLoc2 = line[46];
		link.resName2 = line.substr(47, 3);
		link.chainID2 = line[51];

		temp = line.substr(52, 4);
		link.resSeq2 = atoi(temp.c_str());

		link.iCode2 = line[56];
		link.sym1 = line.substr(59, 6);
		link.sym2 = line.substr(66, 6);

		temp = line.substr(73, 5);
		link.length = atof(temp.c_str());
		cerr << "LINE " << lineNum << ": " << link << endl;
		return link;
	}

	/*
	http://www.wwpdb.org/documentation/format23/sect6.html
	COLUMNS        DATA TYPE       FIELD         DEFINITION
	-------------------------------------------------------------------
	 1 -  6        Record name     "SSBOND"
	 8 - 10        Integer         serNum       Serial number.
	12 - 14        LString(3)      "CYS"        Residue name.
	16             Character       chainID1     Chain identifier.
	18 - 21        Integer         seqNum1      Residue sequence number.
	22             AChar           icode1       Insertion code.
	26 - 28        LString(3)      "CYS"        Residue name.
	30             Character       chainID2     Chain identifier.
	32 - 35        Integer         seqNum2      Residue sequence number.
	36             AChar           icode2       Insertion code.
	60 - 65        SymOP           sym1         Symmetry oper for 1st resid
	67 - 72        SymOP           sym2         Symmetry oper for 2nd resid

	74 - 78                        length       link distance (PDB format extension)
	*/
	SSBond SSBond::readIn(const string &line, unsigned int lineNum) {
		if (line.length() < 80) {
			cerr << "WARNING: LINE " << lineNum << " has less than 80 characters; SSBOND may be parsed incorrectly!\n";
		}
		SSBond ssbond;
		ssbond.record_name = line.substr(0, 6);

		string temp = line.substr(7, 3);
		ssbond.serNum = atoi(temp.c_str());

		ssbond.resName1 = line.substr(11, 3);
		ssbond.chainID1 = line[15];

		temp = line.substr(17, 4);
		ssbond.seqNum1 = atoi(temp.c_str());

		ssbond.iCode1 = line[21];
		ssbond.resName2 = line.substr(25, 3);
		ssbond.chainID2 = line[29];

		temp = line.substr(31, 4);
		ssbond.seqNum2 = atoi(temp.c_str());

		ssbond.iCode2 = line[35];
		ssbond.sym1 = line.substr(59, 6);
		ssbond.sym2 = line.substr(66, 6);

		temp = line.substr(73, 5);
		ssbond.length = atof(temp.c_str());

		cerr << "LINE " << lineNum << ": " << ssbond << endl;
		return ssbond;
	}

	/*
	http://www.wwpdb.org/documentation/format23/sect9.html
	COLUMNS      DATA TYPE        FIELD      DEFINITION
	------------------------------------------------------
	 1 -  6      Record name      "ATOM    "
	 7 - 11      Integer          serial     Atom serial number.
	13 - 16      Atom             name       Atom name.
	17           Character        altLoc     Alternate location indicator.
	18 - 20      Residue name     resName    Residue name.
	22           Character        chainID    Chain identifier.
	23 - 26      Integer          resSeq     Residue sequence number.
	27           AChar            iCode      Code for insertion of residues.
	31 - 38      Real(8.3)        x          Orthogonal coordinates for X in
											 Angstroms
	39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in
											 Angstroms
	47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in
											 Angstroms
	55 - 60      Real(6.2)        occupancy  Occupancy.
	61 - 66      Real(6.2)        tempFactor Temperature factor.
	77 - 78      LString(2)       element    Element symbol, right-justified.
	79 - 80      LString(2)       charge     Charge on the atom.
	*/
	Atom Atom::readIn(const string &line, unsigned int lineNum) {
		if (line.length() < 80) {
			cerr << "WARNING: LINE " << lineNum << " has less than 80 characters; ATOM may be parsed incorrectly!\n";
		}
		Atom atom;
		atom.record_name = line.substr(0, 6);

		string temp = line.substr(6, 5);
		atom.atom_serial_num = atoi(temp.c_str());
		atom.atom_name = line.substr(12, 4); removeWhitespaces(atom.atom_name);
		atom.altLoc = line[16];
		atom.resName = line.substr(17, 3);
		atom.chainID = line[21];

		temp = line.substr(22, 4);
		atom.resSeq = atoi(temp.c_str());
		atom.iCode = line[26];

		temp = line.substr(30, 8);
		atom.coord[0] = atof(temp.c_str());
		temp = line.substr(38, 8);
		atom.coord[1] = atof(temp.c_str());
		temp = line.substr(46, 8);
		atom.coord[2] = atof(temp.c_str());

		atom.chainID_mod = atom.chainID;

		cerr << "LINE " << lineNum << ": " << atom << endl;
		return atom;
	}
}
