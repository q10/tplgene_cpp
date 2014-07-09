#include <cstdlib>
#include <iostream>
#include <array>
#include "OptionParser.h"
#include "TGTypes.hpp"

using namespace std;
using namespace optparse;

namespace po = boost::program_options;

namespace TG {

	SDATA::SDATA(int argc, char *argv[]) {
		OptionParser parser = OptionParser()
			.usage("Usage: %prog [OPTION]...")
			.version("\n%prog 1.0\nCopyright (C) 2014 Benson Ma (マ・ベンソン)\n")
			//.description("Please read instructions")
			.epilog(" ");

		parser.set_defaults("stdpdb", "0");
		parser.set_defaults("ss", "0");
		parser.set_defaults("term", "0");
		parser.set_defaults("water", "TIP3P");

		parser.add_option("--stdpdb") .action("store_true") .dest("stdpdb") .help("Ouput to standard PDB instead of PRESTO PDB");
		parser.add_option("--ss") .action("store_true") .dest("ss") .help("Enable automatic SSBond detection");
		parser.add_option("--term") .action("store_true") .dest("term") .help("Enable automatic N, C terminal capping");
		parser.add_option("--name") .action("append") .help("Molecule name (can invoke multiple times for list of names)");
		parser.add_option("--inputcoord") .help("Input coords filename (*.PDB/*.DIH)");
		parser.add_option("--db") .help("Input topology DB filename (*.TPL)");
		parser.add_option("--outputcoord") .help("Output coords filename");
		parser.add_option("--outputff") .help("Output topology filename");
		array<string, 3> chaintypes {{ "PEP", "NUC", "COM" }};
		parser.add_option("--chain") .choices(chaintypes.begin(), chaintypes.end()).help("Chain species");
		array<string, 2> inputtypes {{ "PDB", "DIHED" }};
		parser.add_option("--inputformat") .choices(inputtypes.begin(), inputtypes.end()).help("Input coord format (defaults to PDB)");
		array<string, 2> watertypes {{ "TIP3P", "TIP4P" }};
		parser.add_option("--water") .choices(watertypes.begin(), watertypes.end()) .help("Water model (defaults to TIP3P)");

		Values& options = parser.parse_args(argc, argv);

		cerr << "stdpdb: " << (options.get("stdpdb") ? "true" : "false") << endl
			<< "ss: " << (options.get("ss") ? "true" : "false") << endl
			<< "term: " << (options.get("term") ? "true" : "false") << endl
			<< "inputcoord: " << options["inputcoord"] << endl
			<< "db: " << options["db"] << endl
			<< "outputcoord: " << options["outputcoord"] << endl
			<< "outputff: " << options["outputff"] << endl
			<< "chain: " << options["chain"] << endl
			<< "inputformat: " << options["inputformat"] << endl
			<< "water: " << options["water"] << endl;
			cerr << "molname: [ "; for (auto it = options.all("name").begin(); it != options.all("name").end(); ++it) cerr << *it << ", "; cerr << "]\n";

		autoSSDetect = options.get("ss");
		autoTerminalCap = options.get("term");
		waterModel = (options["water"] == "TIP3P") ? WaterModel::TIP3P : WaterModel::TIP4P;

		inputFileFormat = (options["inputformat"] == "DIHED") ? FileFormat::DIHED : FileFormat::PDB;
		// inputPDBFile = options["inputcoord"];
		chainType = (options["chain"] == "COM") ? ChainSpecies::COMPOUND : (options["chain"] == "NUC") ? ChainSpecies::NUCLEOTIDE : ChainSpecies::PEPTIDE;
		FFDBFile = options["db"];

		outputCoordFormat = (options.get("stdpdb")) ? PDBFormat::STANDARD : PDBFormat::PRESTO;
		outputCoordFile = options["outputcoord"];
		outputFFFile = options["outputff"];

	}

}

using namespace TG;

int main(int argc, char *argv[]) {




	//SDATA sdata(argc, argv);

	//cout << "sdfdsf " << sdata.chainType << endl;
	return 0;
}
