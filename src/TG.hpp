#ifndef TG_HPP
#define TG_HPP

#include <string>
#include <vector>
#include <array>
#include <iostream>
#include "PDB.hpp"

namespace PDB {
	struct Residue;
};

namespace TG {
	enum FileFormat { PDB, TPL, DIHED };
	enum FFType { UNKNOWN, AMBER, CHARMM19, CHARMM22 };
	enum PDBFormat { STANDARD, PRESTO };
	enum WaterModel { TIP3P, TIP4P };
	enum NAType { NONE, DNA, RNA };

	bool TPLGetLine(std::istream &ifs, std::string &line, unsigned int &lineCount, unsigned int *bufferpos=nullptr);


	/*************************************************************/
	/* Struct of Pair Information                                */
	/*************************************************************/
	typedef struct TGPair {
		int    pair_num = 0;        /* Number of 1-2 Interaction */
		int    tri_num = 0;         /* Number of 1-3 Interaction */
		int    tetra_num = 0;       /* Number of 1-4 Interaction */
		std::vector<int>    pair;           /* Arrangement of 1-2 Interaction */
		std::vector<int>    tri;            /* Arrangement of 1-3 Interaction */
		std::vector<int>    tetra;          /* Arrangement of 1-4 Interaction */
		int numInteractions();
		friend std::ostream& operator<<(std::ostream &os, const TGPair &tg);
	} TGPair, *TGPairPtr;

	/*************************************************************/
	/* Struct of Zmat Informations                               */
	/*************************************************************/
	typedef struct TGZmat {
		int    ith = 0;             /* Relative number of 1-2 Interaction */
		int    jth = 0;             /* Relative number of 1-3 Interaction */
		int    kth = 0;             /* Relative number of 1-4 Interaction */
		int    lth = 0;
		int    mth = 0;             /* ON/OFF COORDINATE */
		int    iflgtr = 0;     /* Torsion Flag */
		double    bond = 0;         /* Bond Length */
		double    angle = 0;        /* Bond Angle*/
		double    torsion = 0;      /* Dihedral Angle */
		friend std::ostream& operator<<(std::ostream &os, const TGZmat &tg);
	} TGZmat, *TGZmatPtr;

	/*************************************************************/
	/* Struct of Coord.                                          */
	/*************************************************************/
	typedef struct TGCoord {
		double       x = 0;
		double       y = 0;
		double       z = 0;
		friend std::ostream& operator<<(std::ostream &os, const TGCoord &tg);
	} TGCoord, *TGCoordPtr;

	/*************************************************************/
	/* Struct of Atom                                            */
	/*************************************************************/
	typedef struct TGAtom {
		int 			num = 0;              // Serial Number of All Atoms
		int    			atom_type = 0;        // Atomic Type Number
		int    			residue_num = 0;      // Residue Number (for amberDB)
	    int    			residue_num_org = 0;      // Residue Number (for TGSys)
		int    			residue_num_mod = 0; // Modified Residue Number
		int    			chain_num = 0;        // Chain Number
		int    			mol_num = 0;          // Number of Molecules
		std::string 	type;
		std::string 	atom_name;     // Atomic Name
		std::string 	atom_name_type;// Atomic Type Name
		std::string 	residue_name;  // Residue Name
		std::string 	elem_symbol;
		double 			mass = 0;          // weight
		double 			radius = 0;        // vdW radius
		double 			charge = 0;        // Partial Charge
		char 			chainID = 0;    // Chain ID
		char 			chainID_mod = 0;   // Modified Chain ID
		char 			chainID_mod2 = 0;   // Modified Chain ID
		char 			iCode = 0;        // insertion code
		char 			altLoc = 0;       // alternate location indicator
		TGPair 			pair;
		TGZmat 			zmat;
		TGCoord 		coord;
		friend std::ostream& operator<<(std::ostream &os, const TGAtom &tg);
		static TGAtom readInTPL(const std::string &line, unsigned int lineNum=0);
		static TGAtom readInPDB(const std::string &line, unsigned int lineNum=0, bool clearWhitespacesForDNA=false);

		void fixResNamesForDNA(NAType nucleicAcidType, bool firstAtomInResidue);
	} TGAtom, *TGAtomPtr;

	/*************************************************************/
	/* Struct of Bond                                            */
	/*************************************************************/
	typedef struct TGBond {
		int    num = 0;             /* Serial Number of Bond */
		int    ith = 0;             /* The Number of Atom Which Constitutes Bond */
		int    jth = 0;             /* The Number of Atom Which Constitutes Bond */
		double    force = 0;        /* Force Constant */
		double    distance = 0;     /* Distance */
		friend std::ostream& operator<<(std::ostream &os, const TGBond &tg);
		static TGBond readInTPL(const std::string &line, unsigned int lineNum=0);
	} TGBond, *TGBondPtr;

	/*************************************************************/
	/* Struct of Angle                                           */
	/*************************************************************/
	typedef struct TGAngle {
		int    num = 0;             /* Serial Number of Angle */
		int    ith = 0;             /* The Number of Atom Which Constitutes Angle*/
		int    jth = 0;             /* The Number of Atom Which Constitutes Angle */
		int    kth = 0;             /* The Number of Atom Which Constitutes Angle */
		double    force = 0;        /* Force Constant */
		double    angle = 0;        /* Angle */
		double    fyfub = 0;
		double    fyqub = 0;
		friend std::ostream& operator<<(std::ostream &os, const TGAngle &tg);
		static TGAngle readInTPL(const std::string &line, unsigned int lineNum=0);
	} TGAngle, *TGAnglePtr;

	/*************************************************************/
	/* Struct of Torsion, Improper_Torsion Parts                 */
	/*************************************************************/
	typedef struct TGTParts {
		int 		tors_num = 0;
		int 		type = 0;
		double 		force = 0;
		double 		symmetry = 0;
		double 		phase = 0;
		double 		fyvwrd = 0;
		double 		fyvwme = 0;
		friend std::ostream& operator<<(std::ostream &os, const TGTParts &tg);
	} TGTParts, *TGTPartsPtr;

	/*************************************************************/
	/* Struct of Torsions                                        */
	/*************************************************************/
	typedef struct TGTorsion {
		int 		num = 0;
		int 		ith = 0;
		int 		jth = 0;
		int 		kth = 0;
		int 		lth = 0;
		double 		vdw = 0;
		double 		ess = 0;
		TGTParts 	parts;
		friend std::ostream& operator<<(std::ostream &os, const TGTorsion &tg);
		static TGTorsion readInTPL(const std::string &line, unsigned int lineNum=0);
	} TGTorsion, *TGTorsionPtr;

	/*************************************************************/
	/* Struct of Dihedrals                                       */
	/*************************************************************/
	typedef struct TGDihedRes {
		int 					num = 0;
		std::string 			name;
		std::vector<double> 	DihedAngle = std::vector<double>(11); // size 11
		std::string angleStr();
		friend std::ostream& operator<<(std::ostream &os, const TGDihedRes &tg);
	} TGDihedRes, *TGDihedResPtr;

	/*************************************************************/
	/* Struct of Residues                                        */
	/*************************************************************/
	typedef struct TGResidue {
		std::string   name;
		std::string   name_line;
		char   chainID;              /* Chain Number */
		int    num = 0;
		int    atom_num = 0; // replace with atoms.size()
		int    bond_num = 0; // replace with bonds.size()
		int    angle_num = 0; // replace with angles.size()
		int    torsion_num = 0; // replace with torsions.size()
		int    improper_num = 0; // replace with impropers.size()
		int    total_atom_num = 0;  /* use in tgMkTopologyOfMetals */
		int    join_flg = 0;       /* Flag for SS Bond */
		bool   ss_circ_flg = false;    /* Flag for "Inter SS + Circ" */
		int    resSeq_org = 0;     /* residue number of input pdb */
		char   iCode = 0;        /* insertion code */ //2


		std::vector<TGAtom> atoms;
		std::vector<TGBond> bonds;
		std::vector<TGAngle> angles;
		std::vector<TGTorsion> torsions;
		std::vector<TGTorsion> impropers;
		TGDihedRes dihedral;

		friend std::ostream& operator<<(std::ostream &os, const TGResidue &tg);
		int readInTPLDihed(const std::string &line, unsigned int lineNum=0);
		int appendToResidueName(const std::string &s);
		int setResidueName(const std::string &s);
	} TGResidue, *TGResiduePtr;

	/*************************************************************/
	/* Struct of Function                                        */
	/*************************************************************/
	typedef struct TGFunction {
		int    num = 0;
		int    parm_num = 0;
		std::string    comment;
		friend std::ostream& operator<<(std::ostream &os, const TGFunction &tg);
		static TGFunction readInTPL(const std::string &line, unsigned int lineNum=0);
	} TGFunction, *TGFunctionPtr;

	/*************************************************************/
	/* Struct of NonBond                                         */
	/*************************************************************/
	typedef struct TGNonBond {
		int    atomtype = 0;          /* Atom Type Number */
		int    atomtypei = 0;         /* Atom Type Number */
		int    atomtypej = 0;         /* Atom Type Number */
		int    functype = 0;          /* Function Type */
		double    radius = 0;         /* Radius */
		double    depth = 0;          /* Ei */
		double    force = 0;          /* Force Const */
		double    electro = 0;        /* Electrostatic Interaction Parameter */
		friend std::ostream& operator<<(std::ostream &os, const TGNonBond &tg);
		static TGNonBond readInTPL(const std::string &line, unsigned int lineNum=0);
	} TGNonBond, *TGNonBondPtr;

	/*************************************************************/
	/* Struct of Circle                                          */
	/*************************************************************/
	typedef struct TGCircle {
		int start = 0;
		int end = 0;
		friend std::ostream& operator<<(std::ostream &os, const TGCircle &tg);
		TGCircle(int t_start) : start(t_start) {}
		TGCircle(int t_start, int t_end) : start(t_start), end(t_end) {}
	} TGCircle, *TGCirclePtr;

	// TGSystem
	/*************************************************************/
	/* Struct of Molecules                                       */
	/*************************************************************/
	typedef struct TGMol {
		static const int 	ISSOFS;
		std::string    		mol_name;							// Molecular Name
		int 				atom_num = 0;							// Number of Atoms
		int 				res_num = 0;							// Number of Residues
		bool				icrclf = false;						// defaults to zero
		std::vector<int>    icnres = std::vector<int>(1000);	// 1000
		int 				inter_ss;
		std::vector<int>    circ_start = std::vector<int>(10);	// 10
		std::vector<int>    circ_end = std::vector<int>(10);	// 10

		TGAtom    **atom;
		std::vector<TGResidue> residues;
		std::vector<TGCircle> circles;

		double    **fxcord1;
		double    **fxcord2;
		double    **fxcord3;
		int    **ixincd1;
		int    **ixincd2;
		int    **ixincd3;
		int    **ixincd4;
		int    **ixincd5;
		double    **fxincd1;
		double    **fxincd2;
		double    **fxincd3;
		int    *r_num;

		std::string    chainID; //2
		std::string    chainID_mod; //2
		bool    ss_circ_flg = false;

		unsigned int numAtoms() const;

		TGMol() {}
		friend std::ostream& operator<<(std::ostream &os, const TGMol &tg);

		void labelCTerminalForDNA();
		void ensureUniqueResidueNumbers();
	} TGMol, *TGMolPtr;

	/*************************************************************/
	/* Struct of SS Bond                                         */
	/*************************************************************/
	// I think it is TGSSBondInfo, since it describes all ssbonds between and within a chaingroup
	typedef struct TGSSBond {
		int 	seqNum1_org = 0;
		int 	seqNum2_org = 0;
		int 	seqNum1_mod = 0;
		int 	seqNum2_mod = 0;
		char 	chainID1_mod = 0;
		char 	chainID2_mod = 0;
		char 	chainID1_org = 0;
		char 	chainID2_org = 0;
		char 	iCode1 = 0;
		char	iCode2 = 0;

		friend std::ostream& operator<<(std::ostream &os, const TGSSBond &tg);
		static TGSSBond readInPDB(const std::string &line, unsigned int lineNum=0);
		static TGSSBond readInTPLDihed(const std::string &line, unsigned int lineNum=0);
		void sortInternal(const std::vector<char> &orderedChainIDs);
	} TGSSBond, *TGSSBondPtr;

	/*************************************************************/
	/* Struct of Molecular System                                */
	/*************************************************************/
	typedef struct TGSystem {
		const static unsigned int RES_DISPLAY;
		const static std::array<std::string, 8> POSITIVE_RESIDUES;
		const static std::array<std::string, 2> NEGATIVE_RESIDUES;

		int    modec = 0;
		int    modeo = 0;
		int    mol_num = 0;        /* Molecular Number */ // deprecated?
		int    num_metals = 0;     /* Molecular Number of Metals */
		int    SSbond_num = 0;     /* Number of SSbond */ // deprecated?

		std::vector<char>			chainGroups;
		std::vector<TGSSBond>		ssbonds;
		std::vector<TGMol>			molecules;

		unsigned int numAtoms() {
			unsigned int count = 0;
			for (const auto &mol : molecules) count += mol.numAtoms();
			return count;
		}
		// TGSystem() : SSbond_num(0) {}

		TGSystem() {}
		int loadFromDihed(std::istream &ifs, bool output_flag=true);
		int loadFromDNA(std::istream &ifs, bool output_flag=true);
		int loadFromAA(std::istream &ifs, bool output_flag=true);
		int printAsDihed();
		int printAsDNA();
		int printAsAA();
		friend std::ostream& operator<<(std::ostream &os, const TGSystem &tg);

	private:
		int readResiduesDihed(std::istream &ifs, unsigned int &lineCount);
		int readSSBondDihed(std::istream &ifs, unsigned int &lineCount);
		int readAnglesDihed(std::istream &ifs, unsigned int &lineCount);
		int processSSBondsAA();

	} TGSystem, *TGSystemPtr;

	/*************************************************************/
	/* Struct of DB Information                                  */
	/*************************************************************/
	typedef struct FFDB {
		int    atom_type_num = 0;          /* Number of Atom Type */
		int    mol_num = 0;                /* Number of Molecules */
		int    total_mol_num = 0;          /* Total Number of Molecules */
		int    func_num = 0;               /* Number of Functions */
		int    nonbond_num = 0;            /* Number of Nonbond */
		int    SSbond_num = 0;             /* Number of S-S Bridges */

		int 	water_hflg = 0;            /* Check Hydrogen of water */
		int 	water_mflg = 0;            /* Check Mass point of water */

		FFType			forceFieldType;

		std::vector<TGResidue> residues;
		std::vector<TGFunction> functions;
		std::vector<TGNonBond> nonbonds;

		std::vector< std::vector<std::string> > AAResiduesTable;	// residue name list in DB file
		std::vector< std::vector<std::string> > MEResiduesTable;	// residue name list in metals DB file

		FFDB() {}
		FFDB(const std::string &dbfilename);
		bool load(const std::string &dbfilename);
		PDB::ResAtomSpecies checkAtomNameOfMetals(const std::string &atom_name, const std::string &residue_name);
		PDB::ResAtomSpecies isMatchResName(const std::string &residueName);
		bool checkResidueName(const std::string &pdb_line);
		friend std::ostream& operator<<(std::ostream &os, const FFDB &tg);

	private:
		void parseInFFFile(const std::string &filenamepath);
		unsigned int loadResiduesTable(const std::string &filenamepath, std::vector< std::vector<std::string> > &residuesTable);
		int readMolecule(std::ifstream &ifs, unsigned int &lineCount);
		int readAtoms(std::ifstream &ifs, unsigned int &lineCount);
		int readBonds(std::ifstream &ifs, unsigned int &lineCount);
		int readAngles(std::ifstream &ifs, unsigned int &lineCount);
		int readTorsions(std::ifstream &ifs, unsigned int &lineCount);
		int readImproperTorsions(std::ifstream &ifs, unsigned int &lineCount);
		int readFunctions(std::ifstream &ifs, unsigned int &lineCount);
		int readNonBonds(std::ifstream &ifs, unsigned int &lineCount);
	} FFDB, *FFDBPtr;

	int checkResidueNumber(FFDB &DB, const std::string &coordfile);

	/******************************
	 * Struct "SDATA"
	 ******************************/
	typedef struct SDATA {
		std::string   	inputPDBFile; 		// Input pdb file
		std::string   	inputCoordFile; 	// Input Coord. File Name
		std::string   	FFDBFile; 			// Input FF DB File Name
		std::string   	outputCoordFile; 	// Output Coord. File Name
		std::string   	outputFFFile; 		// Output TPL File Name
		std::string   	envFile; 			// Env. File Name

		PDB::ChainSpecies	chainType = PDB::ChainSpecies::PEPTIDE;
		FileFormat			inputCoordFormat;
		PDBFormat			outputCoordFormat = PDBFormat::PRESTO;
		WaterModel 			waterModel = WaterModel::TIP3P;

		bool 			autoSSDetect = false;		// SS bond detection mode
		bool 			autoTerminalCap = false;	// N, C terminal cap mode

		int 			hetatms = 0;
		int 			metals = 0;
		int 			titleline = 0;               /* Number of Title Lines */
		int 			mol_num = 0;                 /* Number of Molecules*/

		std::vector<std::string>   ligandNames;		// ligand names -  LGA, LGB, ..., LGZ
		std::vector<std::string>   ligandFFFiles;	// ligand tpl filenames

		std::vector<std::string> titles;
		std::vector<std::string> moleculeNames;

		SDATA() {}
		SDATA(int argc, char *argv[]);
		friend std::ostream& operator<<(std::ostream &os, const SDATA &sdata);
	} SDATA, *SDATAPtr;


	/******************************
	 * Struct of chainid conversion table for Residues
	 ******************************/
	typedef struct ResChainIDConversion {
		char 			chainID_org = 0;
		char 			chainID_mod = 0;
		char 			iCode_org = 0;
		char 			iCode_mod = 0;
		std::string		res_name;
		int 			resSeq_org = 0;
		int 			resSeq_mod = 0;
		ResChainIDConversion() {}
		ResChainIDConversion(const PDB::Residue &res);
	} ResChainIDConversion, *ResChainIDConversionPtr;

	/* add start by tahara 2013/01/04 */
	/******************************
	 * Struct "RemarkMessage"
	 ******************************/
	typedef struct RemarkMessage{
	//    char    remark_message[5000][BUFSIZE_2];
	//    long    remark_position[5000];
	//    char    remark_position_name[THOUSAND][COMMENT_MAX];
	//    long    remark_ss_position;
	} REMMES, *REMMESPtr;

	std::ostream& operator<<(std::ostream &os, const TGPair &tg);
	std::ostream& operator<<(std::ostream &os, const TGZmat &tg);
	std::ostream& operator<<(std::ostream &os, const TGCoord &tg);
	std::ostream& operator<<(std::ostream &os, const TGAtom &tg);
	std::ostream& operator<<(std::ostream &os, const TGBond &tg);
	std::ostream& operator<<(std::ostream &os, const TGAngle &tg);
	std::ostream& operator<<(std::ostream &os, const SDATA &sdata);
};

#endif
