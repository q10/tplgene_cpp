#ifndef PDB_HPP
#define PDB_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>

namespace PDB {
	const double SSBOND_THRESHOLD 		= 2.5;
	const double S_METAL_DIST_THRESHOLD = 2.5;

	enum ChainSpecies { PEPTIDE, NUCLEOTIDE, COMPOUND };
	enum ResAtomSpecies { NOMATCH, AA, ION, WATER, METALS_WATER, LIGAND };

	typedef struct Coord {
		double       x;
		double       y;
		double       z;
		friend std::ostream& operator<<(std::ostream &os, const Coord &tg);
	} Coord, *CoordPtr;

	typedef struct Atom {
		std::string         record_name;        // "ATOM  " or "HETATM"
		int                 atom_serial_num;    // Atom serial number
		std::string         atom_name;          // Atom name
		char         		altLoc = ' ';		// Alternate location indicator
		std::string         resName;            // Residue name
		char				chainID = ' ';		// Chain identifier
		int                 resSeq;             // Residue sequence number
		char				iCode = ' ';		// Code for insertion of residues
		std::vector<double> coord = std::vector<double>(3);             // Orthogonal coordinates for X,Y,Z in Angstroms

		char         		iCode_mod = ' ';	// Code for insertion of residues
		char         		chainID_mod = ' ';	// modified std::stringinID
		int                 resSeq_mod;			// modified resSeq
		char         		iCode_mod2 = ' ';	// Code for insertion of residues
		char         		chainID_mod2 = ' ';	// modified std::stringinID
		int                 resSeq_mod2;        // modified resSeq
		bool  				del_flg = false;	// delete flag for altLoc

		friend std::ostream& operator<<(std::ostream &os, const Atom &p);
		static Atom readIn(const std::string &line, unsigned int lineNum=0);
	} Atom, *AtomPtr;

	typedef struct Link {
		std::string 	record_name;    // "LINK  "
		std::string 	name1;          // Atom name 1
		char			altLoc1 = ' ';	// Alternate location indicator
		std::string 	resName1;       // Residue name 1
		char			chainID1 = ' ';	// Chain identifier
		int  			resSeq1;		// Residue sequence number
		char			iCode1 = ' ';	// Insersion code
		std::string 	name2;          // Atom name 2
		char			altLoc2 = ' ';	// Alternate location indicator
		std::string 	resName2;		// Residue name 2
		char			chainID2 = ' ';	// Chain identifier
		int  			resSeq2;		// Residue sequence number
		char 			iCode2 = ' ';	// Insersion code
		std::string 	sym1;           // Symmetry operator atom 1
		std::string 	sym2;           // Symmetry operator atom 2
		double 			length;         // Link distance

		char 			iCode1_mod = ' ';		// Insersion code
		char 			iCode2_mod = ' ';		// Insersion code
		char 			chainID1_mod = ' ';		// modified chainID1
		char 			chainID2_mod = ' ';		// modified chainID2
		int  			resSeq1_mod;			// modified resSeq1
		int  			resSeq2_mod;			// modified resSeq2
		char 			iCode1_mod2 = ' ';		// Insersion code
		char 			iCode2_mod2 = ' ';		// Insersion code
		char 			chainID1_mod2 = ' ';	// modified chainID1
		char 			chainID2_mod2 = ' ';	// modified chainID2
		int  			resSeq1_mod2;			// modified resSeq1
		int  			resSeq2_mod2;			// modified resSeq2

		friend std::ostream& operator<<(std::ostream &os, const Link &p);
		static Link readIn(const std::string &line, unsigned int lineNum=0);
	} Link, *LinkPtr;

	typedef struct SSBond {
		std::string 	record_name;    // "SSBOND"
		int  			serNum;         // Serial number
		std::string 	resName1;       // "CYS"
		char 			chainID1 = ' ';	// Chain identifier
		int  			seqNum1;        // Residue sequence number
		char 			iCode1 = ' ';	// Insertion code
		std::string 	resName2;       // "CYS"
		char 			chainID2 = ' ';	// Chain identifier
		int  			seqNum2;        // Residue sequence number
		char 			iCode2 = ' ';	// Insertion code
		std::string 	sym1;           // Symmetry operator atom 1
		std::string 	sym2;           // Symmetry operator atom 2
		double 			length;         // Link distance

		char 			iCode1_mod = ' ';	// Insersion code
		char 			iCode2_mod = ' ';	// Insersion code
		char 			chainID1_mod = ' ';	// modified chainID1
		char 			chainID2_mod = ' ';	// modified chainID2
		int  			seqNum1_mod;		// modified resSeq1
		int  			seqNum2_mod;		// modified resSeq2

		char 			iCode1_mod2 = ' ';		// Insersion code
		char		 	iCode2_mod2 = ' ';		// Insersion code
		char 			chainID1_mod2 = ' ';	// modified chainID1
		char		 	chainID2_mod2 = ' ';	// modified chainID2
		int  			resSeq1_mod2;			// modified resSeq1
		int  			resSeq2_mod2;			// modified resSeq2

		friend std::ostream& operator<<(std::ostream &os, const SSBond &p);
		static SSBond readIn(const std::string &line, unsigned int lineNum=0);
	} SSBond, *SSBondPtr;

	typedef struct Residue {
		int  					atom_num;
		char 					chainID = ' ';
		int  					resSeq;
		char 					iCode = ' ';
		std::string 			resName;

		int 					atom_num_bef;

		std::string 			atom_name_Nterm;
		std::string 			atom_name_Cterm;
		std::vector<double> 	coord_Nterm = std::vector<double>(3);
		std::vector<double> 	coord_Cterm = std::vector<double>(3);
		int  					atom_num_Nterm;
		int  					atom_num_Cterm;

		char 					iCode_mod = ' ';
		char 					chainID_mod = ' ';
		int  					resSeq_mod;
		char 					iCode_mod2 = ' ';
		char 					chainID_mod2 = ' ';
		int  					resSeq_mod2;

		bool 					circ_flg = false;
		std::vector<Atom>		atoms;

		// void update();
		friend std::ostream& operator<<(std::ostream &os, const Residue &p);
	} Residue, *ResiduePtr;

	typedef struct PDBStruct {
		int atom_num; // replace with atomnum()
		int ret_num; // no idea

		std::vector<Residue> 	residues;
		std::vector<Link> 		links;
		std::vector<SSBond> 	ssbonds;

		static PDBStruct readIn(const std::string &filenamepath, ChainSpecies chainType);
		bool writeToFile(const std::string &filenamepath);
		friend std::ostream& operator<<(std::ostream &os, const PDBStruct &p);

		void addAtom(const std::string &line, ChainSpecies chainType, unsigned int lineCount=0);
		void setToDeleteAltLocs();
		void tgModify4LackedRes(ChainSpecies chainType);

		PDBStruct() { std::copy( chainIDList.begin(), chainIDList.end(), std::back_inserter(availableChainIDs) ); }
		// int generateCleanedPDBFile(TG::SDATA &sdata, std::vector<TG::ResChainIDConv> &residueConversionTable);

	private:
		std::string chainIDList = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
		std::list<char> availableChainIDs;
		void setModifiedChainIDSsbond();
		void setModifiedChainID(unsigned int startingFrom);
		void setSameChainID(unsigned int currentIndex, ChainSpecies chainType);
		void checkSsbondAndResInfo();
		void autoSSDetection();


	} PDBStruct, *PDBStructPtr;

	typedef struct PDBFile {
		std::vector<Atom> 		atoms;
		std::vector<Link> 		links;
		std::vector<SSBond> 	ssbonds;

		static PDBFile readIn(const std::string &filenamepath);
		static int checkNonStandardResidues(const std::string &filenamepath);
		static int extractRemarks(std::istream &ifs, std::vector<std::string> &remarks);
		static int cleanPDBFile(std::istream &ifs, std::ostream &ofs,
								std::vector< std::vector<std::string> > &AAResiduesTable,
								std::vector< std::vector<std::string> > &MEResiduesTable);
		static bool isCircular(const Atom &atom, std::vector<Link> &links);
		static int fixAtomRecordTag(std::string &pdb_line,
									std::vector< std::vector<std::string> > &AAResiduesTable,
									std::vector< std::vector<std::string> > &MEResiduesTable);

		friend std::ostream& operator<<(std::ostream &os, const PDBFile &p);
	} PDBFile, *PDBFilePtr;


	std::ostream& operator<<(std::ostream &os, const Coord &tg);
	std::ostream& operator<<(std::ostream &os, const Atom &p);
	std::ostream& operator<<(std::ostream &os, const Link &p);
	std::ostream& operator<<(std::ostream &os, const SSBond &p);
	std::ostream& operator<<(std::ostream &os, const Residue &p);
	std::ostream& operator<<(std::ostream &os, const PDBStruct &p);
	std::ostream& operator<<(std::ostream &os, const PDBFile &p);

	bool PDBGetLine(std::istream &ifs, std::string &line, unsigned int &lineCount);

	ResAtomSpecies isMatchResName(const std::string &residueName,
								std::vector< std::vector<std::string> > &AAResiduesTable,
								std::vector< std::vector<std::string> > &MEResiduesTable);

};


#endif
