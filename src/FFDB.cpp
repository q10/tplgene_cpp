#include <cstdlib>
#include <fstream>
#include <array>
#include "TGUtils.hpp"
#include "PDB.hpp"
#include "TG.hpp"

using namespace std;

namespace TG {
	int FFDB::readMolecule(ifstream &ifs, unsigned int &lineCount) {
		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			stringstream thestream(line); TGResidue res;
			if (thestream >> res.name) {
				cout << "LINE " << lineCount << ": " << res.name << endl;
				if (res.name.empty()) continue;
				residues.push_back(res);
			} else {
				return -1;
			}
		} return 0;
	}

	int FFDB::readAtoms(ifstream &ifs, unsigned int &lineCount) {
		string residue; TPLGetLine(ifs, residue, lineCount);
		cout << "\n\n\n\n" << residue << "\n";

		auto residue_iter = std::find_if(residues.begin(), residues.end(), [&](const TGResidue& res) -> bool { return res.name == residue; });
		if (residue_iter == residues.end()) { cerr << "ERROR: CANNOT FIND RESIDUE " << residue << endl; }

		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			TGAtom atom;
			stringstream thestream(line); string dmy;
			if (thestream >> atom.atom_name >> atom.atom_name_type >> atom.atom_type >> atom.residue_name >> atom.residue_num >> atom.mass
							>> atom.radius >> atom.charge >> atom.pair.pair_num >> atom.pair.tri_num >> atom.pair.tetra_num >> dmy >> atom.num) {
				cerr << "OK" << endl;
			} else {
				return -1;
			}
			cout  << atom.atom_name << " " << atom.atom_name_type << " " << atom.atom_type << " " << atom.residue_name << " " << atom.residue_num << " "
					<< atom.mass << " " << atom.radius << " " << atom.charge << " " << atom.pair.pair_num << " " << atom.pair.tri_num << " "
					<< atom.pair.tetra_num << " " << dmy << " " << atom.num << endl;

			if (atom.pair.numInteractions()) {
				TPLGetLine(ifs, line, lineCount); localCount++;
				stringstream thestream2(line);
				string dmy2; thestream2 >> dmy2;

				vector<int> sec_line( (istream_iterator<int>(thestream2)), istream_iterator<int>() );
				for (const auto &i : sec_line) cout << i << " "; cout << endl;

				if (sec_line.size() != atom.pair.numInteractions()) cerr << "ERROR!!!" << endl;

				int begin = 0, end = atom.pair.pair_num;
				atom.pair.pair = vector<int>( &sec_line[begin], &sec_line[end] );

				begin = end; end += atom.pair.tri_num;
				atom.pair.tri = vector<int>( &sec_line[begin], &sec_line[end] );

				begin = end; end += atom.pair.tetra_num;
				atom.pair.tetra = vector<int>( &sec_line[begin], &sec_line[end] );
			}

			// third line
			TPLGetLine(ifs, line, lineCount); localCount++;
			stringstream thestream3(line); string dmy3;
			if (thestream3 >> dmy3 >> atom.zmat.ith >> atom.zmat.jth >> atom.zmat.kth >> atom.zmat.lth >> atom.zmat.bond
							>> atom.zmat.angle >> atom.zmat.torsion) {
				cerr << "OK\n";
			} else {
				return -1;
			}
			// cout << dmy3 << " " << atom.zmat.ith << " " << atom.zmat.jth << " " << atom.zmat.kth << " " << atom.zmat.lth << " "
			// 		<< atom.zmat.bond << " " << atom.zmat.angle << " " << atom.zmat.torsion << endl;

			residue_iter->atoms.push_back(atom);
		} return 0;
	}

	int FFDB::readBonds(ifstream &ifs, unsigned int &lineCount) {
		string residue; TPLGetLine(ifs, residue, lineCount);
		// cout << "\n\n\n\n" << residue << "\n";

		auto residue_iter = std::find_if(residues.begin(), residues.end(), [&](const TGResidue& res) { return res.name == residue; });
		if (residue_iter == residues.end()) { cerr << "ERROR: CANNOT FIND RESIDUE " << residue << endl; }

		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			TGBond bond;
			stringstream thestream(line); string dmy;
			if (thestream >> bond.ith >> bond.jth >> bond.force >> bond.distance >> dmy >> bond.num) {
				cout << "OK" << endl;
			} else {
				cout << "FAIL" << endl; return -1;
			}
			// cout << bond.ith << " " << bond.jth << " " << bond.force << " " << bond.distance << " " << dmy << " " << bond.num << endl;

			if (bond.ith == 0 or bond.jth == 0) continue;
			residue_iter->bonds.push_back(bond);
		}
		return 0;
	}

	int FFDB::readAngles(ifstream &ifs, unsigned int &lineCount) {
		double pi = 3.141592653589793, pai = 180.0;
		double rad = (pi / pai);

		string residue; TPLGetLine(ifs, residue, lineCount);
		cout << "\n\n\n\n" << residue << "\n";

		auto residue_iter = std::find_if(residues.begin(), residues.end(), [&](const TGResidue& res) { return res.name == residue; });
		if (residue_iter == residues.end()) { cerr << "ERROR: CANNOT FIND RESIDUE " << residue << endl; }

		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			TGAngle angle;
			stringstream thestream(line); string dmy;

			if (forceFieldType == FFType::AMBER) {
				thestream >> angle.ith >> angle.jth >> angle.kth >> angle.force >> angle.angle >> dmy >> angle.num;
			} else if (forceFieldType == FFType::CHARMM19 or forceFieldType == FFType::CHARMM22) {
				thestream >> angle.ith >> angle.jth >> angle.kth >> angle.force >> angle.fyfub >> angle.fyqub >> angle.angle >> dmy >> angle.num;
			} else { cerr << "FFType is unknown; don't know how to parse! Exiting..."; std::exit(5); }

			cerr << (thestream ? "OK" : "FAIL") << endl;
			if (not thestream) return -1;

			if (angle.ith == 0 or angle.jth == 0 or angle.kth == 0) continue;
			// cout << angle.ith << " " << angle.jth << " " << angle.kth << " " << angle.force << " " << angle.angle << " " << dmy << " " << angle.num << endl;
			angle.angle *= rad;
			residue_iter->angles.push_back(angle);
		}
		return 0;
	}

	int FFDB::readTorsions(ifstream &ifs, unsigned int &lineCount) {
		string residue; TPLGetLine(ifs, residue, lineCount);
		cout << "\n\n\n\n" << residue << "\n";

		auto residue_iter = std::find_if(residues.begin(), residues.end(), [&](const TGResidue& res) { return res.name == residue; });
		if (residue_iter == residues.end()) { cerr << "ERROR: CANNOT FIND RESIDUE " << residue << endl; }

		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			TGTorsion torsion;
			stringstream thestream(line); string dmy;
			if (thestream >> torsion.ith >> torsion.jth >> torsion.kth >> torsion.lth >> torsion.parts.force >> torsion.parts.tors_num
							>> torsion.parts.symmetry >> torsion.parts.phase >> torsion.parts.type) {
				cout << "OK" << endl;
			} else {
				cout << "FAIL" << endl; return -1;
			}
			if (torsion.ith == 0 or torsion.jth == 0 or torsion.kth == 0 or torsion.lth == 0) continue;
			torsion.num = localCount+1;

			// cout << torsion.ith << " " << torsion.jth << " " << torsion.kth << " " << torsion.lth << " " << torsion.parts.force << " "
			// 		<< torsion.parts.tors_num << " " << torsion.parts.symmetry << " " << torsion.parts.phase << " " << torsion.parts.type << endl;
			residue_iter->torsions.push_back(torsion);
		}
		return 0;
	}

	int FFDB::readImproperTorsions(ifstream &ifs, unsigned int &lineCount) {
		string residue; TPLGetLine(ifs, residue, lineCount);
		cout << "\n\n\n\n" << residue << "\n";

		auto residue_iter = std::find_if(residues.begin(), residues.end(), [&](const TGResidue& res) { return res.name == residue; });
		if (residue_iter == residues.end()) { cerr << "ERROR: CANNOT FIND RESIDUE " << residue << endl; }

		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			TGTorsion torsion;
			stringstream thestream(line); string dmy;
			if (thestream >> torsion.ith >> torsion.jth >> torsion.kth >> torsion.lth >> torsion.parts.force >> torsion.parts.tors_num
							>> torsion.parts.symmetry >> torsion.parts.phase >> torsion.parts.type >> dmy >> torsion.num) {
				cout << "OK" << endl;
			} else {
				cout << "FAIL" << endl; return -1;
			}
			if (torsion.ith == 0 or torsion.jth == 0 or torsion.kth == 0 or torsion.lth == 0) continue;
			// cout << torsion.ith << " " << torsion.jth << " " << torsion.kth << " " << torsion.lth << " " << torsion.parts.force << " "
			// 		<< torsion.parts.tors_num << " " << torsion.parts.symmetry << " " << torsion.parts.phase << " " << torsion.parts.type << endl;

			residue_iter->impropers.push_back(torsion);
		}
		return 0;
	}

	int FFDB::readFunctions(ifstream &ifs, unsigned int &lineCount) {
		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			TGFunction func;
			stringstream thestream(line);
			if (thestream >> func.num >> func.parm_num >> func.comment) {
				cout << "OK" << endl;
			} else {
				cout << "FAIL" << endl; return -1;
			}
			cout << func.num << " " << func.parm_num << " " << func.comment << endl;
			functions.push_back(func);
		} return 0;
	}

	int FFDB::readNonBonds(ifstream &ifs, unsigned int &lineCount) {
		unsigned int localCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); localCount++) {
			if (line[0] == 0 or line.compare(0, 4, "PRE>") == 0) break;

			TGNonBond nonbond;
			stringstream thestream(line);

			if (forceFieldType == FFType::AMBER) {
				thestream >> nonbond.atomtypei >> nonbond.atomtypej >> nonbond.functype >> nonbond.radius >> nonbond.depth >> nonbond.electro >> nonbond.force;
			} else if (forceFieldType == FFType::CHARMM19 or forceFieldType == FFType::CHARMM22) {
				thestream >> nonbond.atomtype >> nonbond.radius >> nonbond.depth >> nonbond.electro >> nonbond.force;
			} else { cerr << "FFType is unknown; don't know how to parse! Exiting..."; std::exit(5); }

			cerr << (thestream ? "OK" : "FAIL") << endl;
			if (not thestream) return -1;

			// cout << nonbond.atomtypei << " " << nonbond.atomtypej << " " << nonbond.functype << " " << nonbond.radius << " "
			//		<< nonbond.depth << " " << nonbond.electro << " " << nonbond.force << endl;
			nonbonds.push_back(nonbond);
		} return 0;
	}

	void FFDB::parseInFFFile(const string &filenamepath) {
		ifstream ifs(filenamepath, ifstream::in);
		if (not ifs) { cerr << "\n ERROR> tgReadAmber\n       File Open Error!\n       Can not open Force Field DB File \"" << filenamepath << "\"\n"; std::exit(5); }
		residues.clear(); functions.clear(); nonbonds.clear();

		unsigned int lineCount = 0; string line;

		TPLGetLine( ifs, line, lineCount );
		if (line.compare(0, 9, ";CHARMm19") == 0) forceFieldType = FFType::CHARMM19;
		else if (line.compare(0, 9, ";CHARMm22") == 0) forceFieldType = FFType::CHARMM22;
		else if (line.compare(0, 6, ";Amber") == 0) forceFieldType = FFType::AMBER;
		else forceFieldType = FFType::UNKNOWN;

		if (forceFieldType == FFType::AMBER) cout << "\n INFORMATION> tgReadInputTopology\n             Amber Type Topology Database File is read\n";
		else if (forceFieldType == FFType::CHARMM19 or forceFieldType == FFType::CHARMM22) cout << "\n INFORMATION> tgReadInputTopology\n             Charmm Type Topology Database File is read\n";
		else {cerr << "First line of Force Field DB does not contain descriptor of the force field; don't know how to parse! Exiting..."; std::exit(5); }

		for (; TPLGetLine( ifs, line, lineCount ); ) {
			if (line == "PRE>MOLECULES") readMolecule(ifs, lineCount);
			if (line == "PRE>ATOMS") readAtoms(ifs, lineCount);
			if (line == "PRE>BONDS") readBonds(ifs, lineCount);
			if (line == "PRE>ANGLES") readAngles(ifs, lineCount);
			if (line == "PRE>TORSIONS") readTorsions(ifs, lineCount);
			if (line == "PRE>IMPROPER-TORSIONS") readImproperTorsions(ifs, lineCount);
			if (line == "PRE>FUNCTIONS") readFunctions(ifs, lineCount);
			if (line == "PRE>NONBONDS") readNonBonds(ifs, lineCount);
		}
	}

	unsigned int FFDB::loadResiduesTable(const string &filenamepath, vector< vector<string> > &residuesTable) {
		ifstream ifs(filenamepath, ifstream::in);
		if (not ifs) { cerr << " ERROR> getResLists\n        Can not open topology file\"" << filenamepath << "\"\n"; std::exit(5); }

		residuesTable.clear();
		unsigned int lineCount = 0, flg=0;
		for (string line; TPLGetLine( ifs, line, lineCount ); ) {
			if (line[0] == 0) {
				continue;
			} else if (line == "PRE>MOLECULES") {
				flg=1; continue;
			} else if (line == "PRE>ATOMS") {
				break;
			}
			if(flg == 0) continue;

			int pos;
			if ((pos = line.find("N+")) != string::npos) line = line.substr(0, pos);
			else if ((pos = line.find("C-")) != string::npos) line = line.substr(0, pos);
			std::replace(line.begin(), line.end(), '-', ' ');
			for (unsigned int i=0; i < line.length(); i++) {
				if (line[i] == '+' and i >= 3 and strncmp(&line[i-3], "HIS+", 4) != 0) line[i] = ' ';
			}

			stringstream s(line);
			vector<string> res_patterns( (istream_iterator<string>(s)), istream_iterator<string>() );
			res_patterns.erase(std::remove(res_patterns.begin(), res_patterns.end(), "1"), res_patterns.end());
			residuesTable.push_back( res_patterns );
		} return residuesTable.size();
	}

	ResAtomSpecies FFDB::checkAtomNameOfMetals(const string &atom_name, const string &residue_name) {

		array<string, 4> water_residue_names {{ "HOH", "WAT", "H2O", "TIP" }};
		array<string, 5> h_atom_names {{ "H", "H1", "H2", "1H", "2H" }};

	    ResAtomSpecies ret = isMatchResName(residue_name);

	    if (ret == ResAtomSpecies::AA) { // AA/NA
	        //return ResAtomSpecies::NOMATCH;
	        return ret;

	    } else if (ret == ResAtomSpecies::LIGAND) { // Ligand
	        return ret;

	    } else if (ret == ResAtomSpecies::METALS_WATER) { // water or metals
	    	if (std::any_of(water_residue_names.begin(), water_residue_names.end(), [&](const string &a) {return a == residue_name;})) {
	            if (atom_name == "O") return ResAtomSpecies::WATER;
	            if (std::any_of(h_atom_names.begin(), h_atom_names.end(), [&](const string &a) {return a == atom_name;})) {
	                water_hflg = 1; // formerly sdata.water_hflg
			        return ResAtomSpecies::WATER;
	            }
	            if (atom_name == "M") {
	                water_mflg = 1; // formerly sdata.water_mflg
					return ResAtomSpecies::WATER;
	            }
	        } else return ResAtomSpecies::ION;
	    } return ResAtomSpecies::NOMATCH;
	}

	ResAtomSpecies FFDB::isMatchResName(const string &residueName) {
		return PDB::isMatchResName(residueName, AAResiduesTable, MEResiduesTable);
	}

	bool FFDB::checkResidueName(const string &pdb_line) {
		PDB::Atom a = PDB::Atom::readIn(pdb_line);
		if (isMatchResName(a.resName) == ResAtomSpecies::AA) return true;
		else if (a.resName == "I" and a.atom_name != "I") return true; // inosine residue
		return false;
	}

	FFDB::FFDB(const string &dbfilename) { load(dbfilename); }

	bool FFDB::load(const string &dbfilename) {
		// create list of AA, NA residues
		loadResiduesTable(dbfilename, AAResiduesTable);

		// create list of ions
		string DB_PATH = std::getenv("TPL_DB_PATH");
		if (DB_PATH.empty()) DB_PATH = "./";
		else if (DB_PATH.back() != '/') DB_PATH += "/";
		DB_PATH += "metals.tpl";
		loadResiduesTable(DB_PATH, MEResiduesTable);

		// load the DB file into FFDB
		parseInFFFile(dbfilename);
		return true;
	}










	int checkResidueNumber(FFDB &DB, const string &coordfile) {
		ifstream ifs(coordfile, ifstream::in);
		if (not ifs) { cerr << "\nCannot open PDB file \"" << coordfile << "\"!\n"; std::exit(5); }

		unsigned int lineCount = 0;
		for (string line; TPLGetLine( ifs, line, lineCount ); ) {
			if (line.compare(0, 4, "PRE>") == 0) { cerr << "\"" << coordfile << "\" is not a PDB file!\n"; std::exit(5); }
			if (line.compare(0, 4, "ATOM") == 0 or line.compare(0, 4, "HETA") == 0) {
				if (not DB.checkResidueName(line)) {
					//sdata.hetatms = 1;
					PDB::Atom a = PDB::Atom::readIn(line);

					if (not TPLGetLine( ifs, line, lineCount )) break;
					if (not DB.checkResidueName(line)) continue;
					else {
						PDB::Atom b = PDB::Atom::readIn(line);
						if (a.resSeq == b.resSeq - 1) {
							cerr << "\n ERROR> tgCheckResidueNumber\n"
								 <<	"        The non-standard residue was found in protein, nucleic acid.\n"
								 <<	"        If it is ligand molecule, you should insert TER line in PDB file.\n";
						}
					}
				}
			}
		} return 0;
	}

}




/*
int tgCheckResidueNumber(TGADBPtr DB)
{
    int heta_res_num, atom_res_num;

    int ck_code;

    char heta_res_name[BUFSIZE];

    char buf[BUFSIZE];
    char textbuf[BUFSIZE];
    char temp_num[NUM_MAX];
    char temp_name[NUM_MAX];

    FILE *fp;

    strcpy(buf, "\0");
    strcat(buf, sdata.coordfile);

    if(strcmp(buf, "\n") == 0 || strncmp(buf, " ", 1) == 0) {
        tgError(PDB_FILE_NAME_NOT_EXIST, buf, NULL, NULL, DMY_INT, DMY_INT);
    }
    if((fp = fopen(buf, "r+")) == NULL) {
        tgError(PDB_FILE_OPEN_ERROR, buf, NULL, NULL, DMY_INT, DMY_INT);
    }
    if(fgets(textbuf, BUFSIZE, fp) == 0) {
        tgError(PDB_FILE_READ_ERROR, buf, NULL, NULL, DMY_INT, DMY_INT);
    }
    while(WHILE_LOOP) {
        if(strncmp(textbuf, "PRE>", 4) == 0) {
            tgError(NOT_PDB_FILE, buf, NULL, NULL, DMY_INT, DMY_INT);
        }

        if(strncmp(textbuf, "ATOM", PRE_WORD_NUM) == 0 ||
           strncmp(textbuf, "HETA", PRE_WORD_NUM) == 0){
            ck_code = 0;

            if(tgCheckResidueName(textbuf, DB) == 0){
              // set hetatm flg on
                sdata.hetatms = 1;

              // Get residue name and number
                strncpy(temp_name,&textbuf[17],3);
                temp_name[3] = '\0';
                sscanf(temp_name,"%s",heta_res_name);
                strncpy(temp_num,&textbuf[22],4);
                temp_num[4] = '\0';
                sscanf(temp_num,"%d",&heta_res_num);

              // Check next line
                if(fgets(textbuf, BUFSIZE, fp) == 0) break;

                ck_code = tgCheckResidueName(textbuf,DB);

                if(ck_code == 0){
                    continue;
                }
                else if(ck_code == 1){
                    strncpy(temp_num,&textbuf[22],4);
                    temp_num[4] = '\0';
                    sscanf(temp_num,"%d",&atom_res_num);

                    if(heta_res_num == atom_res_num - 1){
                        tgError(WRONG_RES_NUM, heta_res_name, NULL, NULL,
                                heta_res_num, DMY_INT);
                    }
                }
            }
        }

        if(fgets(textbuf, BUFSIZE, fp) == 0) break;
    }

    fclose(fp);

    return 0;

}

*/

