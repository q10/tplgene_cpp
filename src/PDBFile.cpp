#include <fstream>
#include <array>
#include "cpplinq.hpp"
#include "TGUtils.hpp"
#include "PDB.hpp"

using namespace std;
namespace LINQ = cpplinq;

namespace PDB {
	PDBFile PDBFile::readIn(const string &filenamepath) {
		ifstream ifs(filenamepath, ifstream::in);
		PDBFile pdb;

		unsigned int lineCount = 0;
		for (string line; PDBGetLine( ifs, line, lineCount ); ) {
			if (line.compare(0, 6, "SSBOND") == 0) {
				pdb.ssbonds.push_back( SSBond::readIn(line, lineCount) );

			} else if (line.compare(0, 4, "LINK") == 0) {
				pdb.links.push_back( Link::readIn(line, lineCount) );

			} else if (line.compare(0, 4, "ATOM") == 0 or line.compare(0, 6, "HETATM") == 0) {
				pdb.atoms.push_back( Atom::readIn(line, lineCount) );

			} else {
				cerr << "SKIPPING LINE " << lineCount << ": \"" << line << "\"\n";
			}

		}
		cout << "PDB File \'" << filenamepath <<"\':"
				<< "\n    Number of atoms:   " << pdb.atoms.size()
				<< "\n    Number of ssbonds: " << pdb.ssbonds.size()
				<< "\n    Number of links:   " << pdb.links.size() << endl;
		return pdb;
	}

	int PDBFile::checkNonStandardResidues(const string &filenamepath) {
		ifstream ifs(filenamepath, ifstream::in);
		if (not ifs) { cerr << " ERROR> tgCheckNonStdRes\n" << "        Can not open file"; std::exit(1); }

		const int ATOM = 1, HETATM = 2, TER = 3;
		int previousAtomType = 0;
		unsigned int lineCount = 0;
		for (string line; PDBGetLine( ifs, line, lineCount ); ) {
			if (line.compare(0, 4, "ATOM") == 0) {
				if (previousAtomType == HETATM){
					cerr << "\n ERROR> tgCheckNonStdRes\n"
						<< "        The non-standard residue was found in protein, nucleic acid.\n"
						<< "        If it is ligand molecule, you should insert TER line in PDB file.\n";
					std::exit(1);
				} previousAtomType = ATOM;

			} else if (line.compare(0, 3, "TER") == 0) {
				previousAtomType = TER;

			} else if (line.compare(0, 6, "HETATM") == 0) {
				if (previousAtomType == ATOM) {
					cerr << "\n ERROR> tgCheckNonStdRes"
						<< "        The non-standard residue was found in protein, nucleic acid.\n"
						<< "        If it is ligand molecule, you should insert TER line in PDB file.\n";
					std::exit(1);
				} previousAtomType = HETATM;
			}
		} return 0;
	}

	int PDBFile::extractRemarks(istream &ifs, vector<string> &remarks) {
		ifs.clear(); ifs.seekg(0, ios::beg);
		if (not ifs) { cerr << " ERROR> extractRemarks\nstream is invalid; exiting!\n"; std::exit(1); }

		vector<string> lines;
		std::copy( istream_iterator<FileLine>(ifs), istream_iterator<FileLine>(), std::back_inserter(lines) );
		remarks = LINQ::from(lines)
				>> LINQ::where( [](const string &l){ return l.compare(0, 6, "REMARK") == 0; } )
				>> LINQ::to_vector();
		return 0;
	}

	int PDBFile::cleanPDBFile(istream &ifs, ostream &ofs, vector< vector<string> > &AAResiduesTable, vector< vector<string> > &MEResiduesTable) {
		ifs.clear(); ifs.seekg(0, ios::beg);
		if (not ifs) { cerr << " ERROR> cleanPDBFile\nstream is invalid; exiting!\n"; std::exit(1); }


		array<string, 30> ignore6Letters {{ "ANISOU", "HETNAM", "MODRES", "HEADER", "COMPND", "SOURCE", "KEYWDS", "EXTDTA", "AUTHOR", "SEQRES",
											"FORMUL", "CRYST1", "EXPDTA", "REVDAT", "CONECT", "MASTER", "CISPEP", "NUMMDL", "SPRSDE", "SEQADV" }};
		array<string, 6> ignore5Letters {{ "HELIX", "SHEET", "SCALE", "DBREF", "TITLE", "MODEL", }};
		array<string, 4> ignore4Letters {{ "ORIG", "JRNL", "SITE", "HET" }};

		const int NONE=0, ATOM=1, HETATM=2;

		bool prevTerminate = true;
		unsigned int lineCount = 0;
		unsigned int prevAtomType = NONE;
		vector<Link> links;
		Atom prevAtom;
		for (string line; PDB::PDBGetLine( ifs, line, lineCount ); ) {

			if (line.compare(0, 4, "ATOM") == 0 or line.compare(0, 6, "HETATM") == 0) {
				fixAtomRecordTag(line, AAResiduesTable, MEResiduesTable);
				Atom atom = Atom::readIn(line, lineCount);

				if (isCircular(atom, links)) ofs << "CIRCULAR\n";

				// if current atom is the first non-OXT atom, it means we have a new chain
				if ( prevAtom.atom_name == "OXT" and
					(atom.resName != prevAtom.resName or atom.chainID != prevAtom.chainID or atom.resSeq != prevAtom.resSeq) and
					not prevAtom.resName.empty() ) {
					prevTerminate = true; ofs << "TER\n";
				}

				// so previous atom is not OXT, BUT chainID's are different, it means we have a new chain
				else if (prevAtomType != NONE and atom.chainID != prevAtom.chainID and not prevTerminate) {
					prevTerminate = true; ofs << "TER\n";
				}

				// if previous atom is HETATM, and resname is different, then add TER (signify end of non-AA molecule)
				else if (prevAtomType == HETATM and atom.resName != prevAtom.resName and not prevTerminate) {
					prevAtomType = NONE; prevTerminate = true; ofs << "TER\n";
				}

				prevAtomType = (line.compare(0, 4, "ATOM") == 0) ? ATOM : (line.compare(0, 6, "HETATM", 6) == 0) ? HETATM : NONE;

				if (isalnum(atom.altLoc)) ofs << line << endl;
				else {
					line[16] = ' ';
					if (atom.atom_name == prevAtom.atom_name and atom.resName == prevAtom.resName and atom.resSeq == prevAtom.resSeq) {
					} else if (atom.resName != prevAtom.resName and atom.resSeq != prevAtom.resSeq) {
					} else {
						ofs << line << endl;
					}
				}

				prevTerminate = false;
				prevAtom = atom;

			} else if (line.compare(0, 4, "LINK")) {
				prevAtomType = NONE;
				if (line.compare(12, 4, " N  ") != 0 or line.compare(42, 4, " C  ") != 0) continue;
				links.push_back( Link::readIn(line, lineCount) );
				prevTerminate = false; ofs << line << endl;

			} else if ( std::find_if(ignore6Letters.begin(), ignore6Letters.end(), [&](const string& s) {return line.compare(0, 6, s) == 0;}) == ignore6Letters.end() or
						std::find_if(ignore5Letters.begin(), ignore5Letters.end(), [&](const string& s) {return line.compare(0, 5, s) == 0;}) == ignore5Letters.end() or
						std::find_if(ignore4Letters.begin(), ignore4Letters.end(), [&](const string& s) {return line.compare(0, 4, s) == 0;}) == ignore4Letters.end() ) {
				continue;

			} else if ((line.compare(0, 3, "TER") or line.compare(0, 3, "END")) and not prevTerminate) {
				prevTerminate = true; ofs << "TER\n";

			} else {
				prevAtomType = NONE; ofs << line << endl;
			}
		}

		return 0;
	}

	bool PDBFile::isCircular(const Atom &atom, vector<Link> &links) {
		return std::any_of(links.begin(), links.end(), [&](const PDB::Link &l) {
			return atom.atom_name == l.name1 and atom.resName == l.resName1 and atom.chainID == l.chainID1 and
					atom.resSeq == l.resSeq1 and atom.iCode == l.iCode1;
		});
	}

	int PDBFile::fixAtomRecordTag(string &pdb_line, vector< vector<string> > &AAResiduesTable, vector< vector<string> > &MEResiduesTable) {
		if (pdb_line.empty()) return -1;
		string residue_name = pdb_line.substr(17, 4);
		auto res_atom_species = PDB::isMatchResName(residue_name, AAResiduesTable, MEResiduesTable);
		if (res_atom_species == ResAtomSpecies::AA) {
			if (pdb_line.compare(0, 6, "HETATM") == 0) {
				pdb_line[0] = 'A'; pdb_line[1] = 'T'; pdb_line[2] = 'O'; pdb_line[3] = 'M'; pdb_line[4] = ' '; pdb_line[5] = ' ';
			}
		} else {
			if (pdb_line.compare(0, 4, "ATOM") == 0) {
				pdb_line[0] = 'H'; pdb_line[1] = 'E'; pdb_line[2] = 'T'; pdb_line[3] = 'A'; pdb_line[4] = 'T'; pdb_line[5] = 'M';
			}
		} return 0;
	}
}

