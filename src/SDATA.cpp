// c++ -std=c++11 sample.cpp -lboost_program_options -O4 && strip a.out
#include <fstream>
#include <boost/program_options.hpp>
#include "TG.hpp"
#include "TGUtils.hpp"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

namespace TG {
	template<class T>
	ostream& operator<<(ostream& os, const vector<T>& v) {
		std::copy(v.begin(), v.end(), ostream_iterator<T>(os, ", "));  return os;
	}

	void toUpper(string &s) { std::transform(s.begin(), s.end(), s.begin(), ::toupper); }

	time_t getTimeNow() {
		time_t rawtime;
		time(&rawtime);
		return rawtime;
	}

	string timeStamp(time_t &rawtime, const string &format) {
		struct tm * timeinfo;
		char buffer[80];
		timeinfo = localtime (&rawtime);
		std::strftime(buffer, 80, format.c_str(), timeinfo);
		return string(buffer);
	}

	string extensionChopped(const string &s) { return s.substr(0, s.find_last_of(".")); }

	void fileTypeErrorExit(const string &type1, const string &type2) {
		cerr << "\n ERROR> tgReadEnv\n"
			 << "       Selection Error!\n"
			 << "       File Type for " << type1 << " was selected.\n"
			 << "       Please select " << type2 << " Type.\n";
		 std::exit(2);
	}

	void showMismatchFileTypeError(const string &filetype, const string &inputfile) {
        cerr << "\n ERROR> tgCheckInputPDB\n"
			 << "       filetype and inputfile are no matching.\n"
			 << "           filetype  = " << filetype << "\n"
			 << "           inputfile = " << inputfile << "\n\n";
	}

	SDATA::SDATA(int argc, char *argv[]) {
		try {
			time_t now = getTimeNow();
			bool use_stdpdb = false;
			string chain_type, input_format, water_model, currenTime = timeStamp(now, "%Y%m%d_%H%M%S");
			vector<string> chaintypes { "PEP", "NUC", "COM" };
			stringstream chaintypesdesc; chaintypesdesc << "Chain species [ " << chaintypes << "]";
			vector<string> inputtypes { "PDB", "DIHED" };
			stringstream inputtypesdesc; inputtypesdesc << "Input coord format [ " << inputtypes << "]";
			vector<string> watertypes { "TIP3P", "TIP4P" };
			stringstream watertypesdesc; watertypesdesc << "Water model [ " << watertypes << "]";

			po::options_description desc("Allowed options");
			desc.add_options()
				("help,h", "Display this help message")
				("version,v", "Display the version number")
				("interact", "Interactive mode")
				("stdpdb", po::bool_switch(&use_stdpdb), "Ouput to standard PDB instead of PRESTO PDB")
				("ss", po::bool_switch(&autoSSDetect), "Enable automatic SSBond detection")
				("term", po::bool_switch(&autoTerminalCap), "Enable automatic N, C terminal capping")
				("names", po::value< vector<string> >(&moleculeNames)->multitoken(), "Names of the molecules")
				("titles", po::value< vector<string> >(&titles)->multitoken()->default_value(vector<string>(), string("tplgene")+timeStamp(now, "%Y%m%d")), "Titles information")
				("inputcoord", po::value<string>(&inputCoordFile)->required(), "Input coords filename (*.PDB/*.DIH)")
				("db", po::value<string>(&FFDBFile)->required(), "Input topology DB filename (*.TPL)")
				("outputcoord", po::value<string>(&outputCoordFile), "Output coords filename")
				("outputff", po::value<string>(&outputFFFile), "Output topology filename")
				("chain", po::value<string>(&chain_type)->default_value("PEP"), chaintypesdesc.str().c_str())
				("inputformat", po::value<string>(&input_format)->default_value("PDB"), inputtypesdesc.str().c_str())
				("water", po::value<string>(&water_model)->default_value("TIP3P"), watertypesdesc.str().c_str())
			;

			po::variables_map vm;
			po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
			if (vm.count("help")) {
				cerr << "Usage: options_description [options]\n" << desc; std::exit(0);
			} else if (vm.count("version")) {
				cerr << "\n" << argv[0] << " 0.01\nCopyright (C) 2014 Benson Ma (マ・ベンソン)\n\n"; std::exit(0);
			} else if (vm.count("interact")) {
				cerr << "\nInteractive mode still under development; exiting!\n"; std::exit(0);
			}
			po::notify(vm); // Notify throws error messages about missing parameters, etc; so call this AFTER checking for help

			toUpper(chain_type); toUpper(input_format); toUpper(water_model);
			if (std::find_if(chaintypes.begin(), chaintypes.end(), [&](string &s) {return s == chain_type;}) == chaintypes.end()) {
				cerr << "\nChain species must be one of the following - [ " << chaintypes << "] !\n"; std::exit(-1);

			} else if (std::find_if(inputtypes.begin(), inputtypes.end(), [&](string &s) {return s == input_format;}) == inputtypes.end()) {
				cerr << "\nInput coord format must be one of the following - [ " << inputtypes << "] !\n"; std::exit(-1);

			} else if (std::find_if(watertypes.begin(), watertypes.end(), [&](string &s) {return s == water_model;}) == watertypes.end()) {
				cerr << "\nWater model must be one of the following - [ " << watertypes << "] !\n"; std::exit(-1);
			}

			if (outputCoordFile.empty()) outputCoordFile = extensionChopped(inputCoordFile) + "_tpl_" + currenTime + ".pdb";
			if (outputFFFile.empty()) outputFFFile = extensionChopped(inputCoordFile) + currenTime + ".tpl";

			// Convert choices to enums
			outputCoordFormat = use_stdpdb ? PDBFormat::STANDARD : PDBFormat::PRESTO;
			waterModel = (water_model == "TIP4P") ? WaterModel::TIP4P : WaterModel::TIP3P;
			inputCoordFormat = (input_format == "DIHED") ? FileFormat::DIHED : FileFormat::PDB;
			chainType = (chain_type == "NUC") ? PDB::ChainSpecies::NUCLEOTIDE : (chain_type == "COM") ? PDB::ChainSpecies::COMPOUND : PDB::ChainSpecies::PEPTIDE;

			if (chainType == PDB::ChainSpecies::COMPOUND) inputCoordFormat = FileFormat::PDB;

			// defaults C99_aa.tpl C99_na.tpl
			// read envfile instead of reading in terminal .data.inp

			ifstream fffile(FFDBFile); FFType forceFieldType;
			if (not fffile) {
				cerr << "\nCannot open force field file at \'" << FFDBFile << "\'! Exiting...\n"; std::exit(1);
			} else {
				string line; std::getline( fffile, line );
				if (line.compare(0, 9, ";Amber_aa") == 0) {
					forceFieldType = FFType::AMBER;
					if (chainType == PDB::ChainSpecies::NUCLEOTIDE) fileTypeErrorExit("Amino Acid", "Nucleotide");

				} else if (line.compare(0, 9, ";Amber_na") == 0) {
					forceFieldType = FFType::AMBER;
					if (chainType != PDB::ChainSpecies::NUCLEOTIDE) fileTypeErrorExit("Nucleotide", "Amino Acid");

				} else if (line.compare(0, 12, ";CHARMm19_aa") == 0 or line.compare(0, 12, ";CHARMm22_aa") == 0) {
					forceFieldType = (line.compare(0, 12, ";CHARMm19_aa") == 0) ? FFType::CHARMM19 : FFType::CHARMM22;
					if (chainType == PDB::ChainSpecies::NUCLEOTIDE) fileTypeErrorExit("Amino Acid", "Nucleotide");

				} else if (line.compare(0, 12, ";CHARMm19_na") == 0 or line.compare(0, 12, ";CHARMm22_na") == 0) {
					forceFieldType = (line.compare(0, 12, ";CHARMm19_aa") == 0) ? FFType::CHARMM19 : FFType::CHARMM22;
					if (chainType != PDB::ChainSpecies::NUCLEOTIDE) fileTypeErrorExit("Nucleotide", "Amino Acid");

				} else {
					cerr << "\nFirst line of Force Field DB does not contain descriptor of the force field! Exiting...\n"; std::exit(1);
				}
			}

			if (chainType == PDB::ChainSpecies::NUCLEOTIDE and (forceFieldType == FFType::CHARMM19 or forceFieldType == FFType::CHARMM22)) {
				cerr << "\n ERROR> tgReadEnv\n"
					 << "       Selection Error!\n"
					 << "       Nucleotide Calculation for CHARMM type under construction now.\n"; std::exit(3);
			}

			// Checks that the input files are of proper format
			ifstream pdbfile(inputCoordFile); bool fileLooksOK = false;
			if (not pdbfile) {
				cerr << "\nCannot open coord file at \'" << inputCoordFile << "\'! Exiting...\n"; std::exit(1);
			} else {
				string line;
				while (std::getline( pdbfile, line )) {
					if (line.compare(0, 4, "ATOM") == 0) {
						if (inputCoordFormat != FileFormat::PDB) showMismatchFileTypeError(input_format, inputCoordFile);
						else fileLooksOK = true;
						break;

					} else if (line.compare(0, 14, "PRE>SEQUENCE") == 0) {
						if (inputCoordFormat != FileFormat::DIHED) showMismatchFileTypeError(input_format, inputCoordFile);
						else fileLooksOK = true;
						break;
			        }
				}
			}
			if (not fileLooksOK) { cerr << "\nInput data looks incorrectly formatted; are you sure you are pointing to the correct file?\n"; std::exit(4); }

		} catch(std::exception& e) {
			cerr << e.what() << "\n"; std::exit(1);
		}
	}

	ostream& operator<<(ostream &os, const SDATA &sdata) {
		return os << "\nSDATA : {\n\t" << JSONPair("inputPDBFile", sdata.inputPDBFile)
					<< "\n\t" << JSONPair("inputCoordFile", sdata.inputCoordFile)
					<< "\n\t" << JSONPair("FFDBFile", sdata.FFDBFile)
					<< "\n\t" << JSONPair("outputCoordFile", sdata.outputCoordFile)
					<< "\n\t" << JSONPair("outputFFFile", sdata.outputFFFile)
					<< "\n\t" << JSONPair("chainType", sdata.chainType)
					// << "\n\t" << JSONPair("forceFieldType", sdata.forceFieldType)
					<< "\n\t" << JSONPair("inputCoordFormat", sdata.inputCoordFormat)
					<< "\n\t" << JSONPair("autoSSDetect", sdata.autoSSDetect)
					<< "\n\t" << JSONPair("autoTerminalCap", sdata.autoTerminalCap)
					<< "\n\t" << JSONPair("autoTerminalCap", sdata.autoTerminalCap)
					<< "\n}\n\n";
	}
}
