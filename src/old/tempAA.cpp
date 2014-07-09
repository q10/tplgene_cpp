for (auto &molecule : molecules) {
    if (ssbonds.size() == 0) {
        std::fill(molecule.icnres.begin(), molecule.icnres.end(), TGMol::ISSOFS);
        continue;
    }

    // find the first SSBond where atom.residue_num_mod == ssbond.seqNum1_mod or atom.residue_num_mod == ssbond.seqNum2_mod
    // if not found,         molecule.icnres[i] = ISSOFS; continue;
    // else if inner loop
    // do bulk of stuff

    for (unsigned int i=0; i < molecule.residues.size(); i++) {
        auto &residue = molecule.residues[i];
        auto &atom = residue.atoms.front();
        auto ssbond_loc = std::find_if(ssbonds.begin(), ssbonds.end(),
            [&](TGSSBond &ssbond) { return atom.residue_num_mod == ssbond.seqNum1_mod or atom.residue_num_mod == ssbond.seqNum2_mod; }
        );

        if (ssbond_loc == ssbonds.end()) {
            molecule.icnres[i] = ISSOFS;

        } else {
            if (residue.chainID == ssbond.chainID1_mod or residue.chainID == ssbond.chainID1_org or residue.chainID == ssbond.chainID2_org) {
                if (atom.residue_name == "CYS") {
                    residue.appendToResidueName("S");

                } else if (atom.residue_name == "CYSC-") {
                    residue.setResidueName("CYSSC-");

                } else if (atom.residue_name == "CYSN+") {
                    residue.setResidueName("CYSSN+");

                } else if (atom.residue_name == "CYSS") {

                } else {
                    cerr << "\n ERROR> tgReadAminoSequence\n"
                         << "       Residue Number Error!\n"
                         << "       Residue Number in SSBOND is wrong.\n"
                         << "       Mol num     = " << mol_num << "\n"
                         << "       Residue num = " << atom.residue_num << "\n"
                         << "       Residue name= " << residue.name << "\n";
                    std::exit(1);
                }

                // change seqNum1_org, seqNum2_org to seqNum1_mod, seqNum2_mod
                if (i+1 == ssbond.seqNum1_mod) {
                    molecule.icnres[i] = (i+1 > ssbond.seqNum2_mod) ? -ssbond.seqNum2_mod : ssbond.seqNum2_mod;
                } else if (i+1 == ssbond.seqNum2_mod){
                    molecule.icnres[i] = (i+1 > ssbond.seqNum1_mod) ? -ssbond.seqNum1_mod : ssbond.seqNum1_mod;
                } else {
                    molecule.icnres[i] = ISSOFS;
                }
            }
        }
    }
}
