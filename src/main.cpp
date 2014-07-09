#include "TGTypes.hpp"
#include "TGUtils.hpp"

using namespace std;
using namespace PDB;
using namespace TG;


/*
int readInputSequence() {
	if (sdata.chainType == ChainSpecies::PEPTIDE) {

		if (sdata.inputCoordFormat == FileFormat::DIHED) {
			// read in dihed
		} else {
			PDBFile::cleanPDBFile(ifs, ofs, sdataAAResiduesTable, sdata.MEResiduesTable);
			PDBFile::checkNonStandardResidues(ofs);

		}
	} else if (sdata.chainType == ChainSpecies::NUCLEOTIDE) {

	} else if (sdata.chainType == ChainSpecies::COMPOUND) {
		PDBFile::cleanPDBFile(ifs, ofs, sdataAAResiduesTable, sdata.MEResiduesTable);

	}
}
*/


int main(int argc, char *argv[]) {
	SDATA sdata(argc, argv);

	FFDB DB(sdata.FFDBFile);

	if (sdata.inputCoordFormat == FileFormat::PDB and sdata.chainType == ChainSpecies::PEPTIDE)
		TG::checkResidueNumber(DB, sdata.inputCoordFile);



	/******************************
	* Read Amino Acid Sequence
	******************************/
	ret_code = tgReadInputSequence(amber, pdb, resCnvTbl);

	/* add 2011.11.27, start */
	/******************************
	* Check non-standard AA, NA
	******************************/
	ret_code = tgCheckResNames(amber, DB);
	/* add end */

	/******************************
	* Read  PDB File of Metals
	******************************/
	if(sdata.hetatms == 1){
		/* Modify 20080522 */
		ret_code = tgReadPDBFileOfMetals(pdb_metals,
		res_list, &res_list_count, res_num);
		/* Modify end */

		if(ret_code == 1){
			sdata.metals = 1;
		}
	}



	cout << sdata;
	return 0;
}


/*
tgReadInputSequence(TGASystemPtr amber, TGASystemPtr pdb, resConvTblPtr resCnvTbl)

if sdata.chnspecies is peptide
	if pdb file

	else if dihedral file
        amber->modec = 1;
		read dihedral info into amberdb and into pdbstruct

else if species is nucleotide
	tgDeleteAltLocIndicator
	tgReadPdb
	tgCheckNonStdRes
	if pdb file
	    amber->modec = 3;
	    ired = tgCheckInputPDB();
		read in to amber and by using tgReadDNASequence
	else
		error since we cant do DNA and DIHED

else if species is compound
	tgDeleteAltLocIndicator
	tgModify4LackedRes1
	amber->modec = 3;
	tgCheckInputPDB
*/









/*


int main() {
	TGAmberDB *amberDB = ParseIntoForceFieldDB("/Users/bm/Downloads/tplgene/tplgene/DB/C96_aa.tpl");

	for (const auto &res : amberDB->residues) {
		cout << res.name << endl;
		for (const auto &torsion : res.torsions) cout << torsion.parts.force << " ";
			cout << endl;
	}
	cout << "\n\n" << amberDB->nonbonds.size() <<  "\n\n";
	for (const auto &nb : amberDB->nonbonds) {
		cout << nb.depth << endl;
	}

	// lambda that returns a lambda to generate Fibonacci numbers
	cout << "\n\n\nFIBONACCI\n";
	auto fib = [](int n1, int n2) {
		return [=]() mutable {
			auto n = n1;
			n1 = n2;
			n2 += n;
			return n;
		};
	};

	TGSystem* n = ParsePDBFile("/Users/bm/Downloads/tplgene/sample/tplgene_sample/ACE-ALA-ALA-NME.pdb");

	// Output portion of Fibonacci sequence
	auto fibGen = fib(0, 1);
	for (int i = 0; i < 10; ++i)
		std::cout << fibGen() << " ";
	std::cout << "\n";
	std::generate_n(std::ostream_iterator<int>(std::cout, "\n"), 5, fibGen);


	TGPair p;
	p.pair.push_back(3);
	p.tri.push_back(4);
	p.tri.push_back(8);
	cout << p << endl;

}
*/
