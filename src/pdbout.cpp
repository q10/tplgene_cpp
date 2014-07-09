#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "TP_define.h"
#include "TP_struct.h"
#include "TP_AmberDB.h"
#include "TP_Env.h"
#include "TP_error.h"

int tgError(int, char*, char*, char*, int, int);
void *Malloc(int);
int tgCheckAtomNameOfMetals(char *, char *,char *);


extern SDATA sdata;
/* add start by tahara 2013/01/04 */
extern REMMES remmes;
/* add end by tahara 2013/01/04 */

#define PDB_FORMAT_SSBOND_LENGTH 2.04
#define PDB_FORMAT_LINK_LENGTH 1.34
#define PDB_FORMAT_SYMMETRY "1555"

#define PRESTO_CHAINIDS "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

/*--------------------------------------------------------*/
/* int tgOutputRemarkLine()                               */
/*--------------------------------------------------------*/
int tgOutputRemarkLine(FILE* fp, TGASystemPtr amber)
{
    int chainID_idx = 0;
    int i, j, k;

    char format[BUFSIZE + 1];
    char header[BUFSIZE + 1];
    char chainID_old[2];
    char chainID_mod_old[2];

    char presto_chainIDs[] = PRESTO_CHAINIDS;

    int print_flg = 1;
    int header_print_flg = 0;

    TGAMolPtr mol;
    TGAResiduePtr residue;
    TGAAtomPtr atom;

    chainID_idx = 0;
    chainID_old[0] = '\0';

    /* Header */
    strcpy(header, "REMARK     CHAIN ID WAS MODIFIED BY TPLGENE.\n");
    strcpy(format, "REMARK      INPUT:'%1s', OUTPUT:'%c'\n");

    /* reset chain id.  set A, B, C, ... */
    mol = amber->begin;
    for(i = 0; i < amber->mol_num; i++, chainID_idx++) {
        residue = mol->residuebegin;
        for(j = 0; j < mol->res_num; j++) {
            atom = residue->atombegin;

            if(print_flg == 1){
                if(atom->chainID_mod[0] != presto_chainIDs[chainID_idx % 26]){
                    if(header_print_flg == 0){
                        fprintf(fp, "%s", header);
                        header_print_flg = 1;
                    }
                }

                print_flg = 0;
            }
            /* add 2013.11.29 */
            atom->chainID_mod2[0] = presto_chainIDs[chainID_idx % 26];
            atom->chainID_mod2[1] = '\0';
            /* add end */

            if(j != 0 &&
                strcmp(chainID_old, atom->chainID_mod) != 0){
                chainID_idx++;
                print_flg = 1;
            }
            strcpy(chainID_old, atom->chainID_mod);

            /* add 2013.11.29 */
            for(k=0; k < residue->atom_num; k++){
                if(k != 0) atom = atom->next;
                strcpy(atom->chainID_mod2, residue->atombegin->chainID_mod2);
            }
            /* add end */

            residue = residue->next;
        }

        mol = mol->next;
        print_flg = 1;
    }

    /* add 2013.11.28 */
    /* reset chain id, which connect ssbond (intermolecular ssbond) */
    print_flg = 1;
    chainID_old[0] = '\0';
    mol = amber->begin;
    for(i = 0; i < amber->mol_num; i++) {
        if(i != 0) mol = mol->next;


        residue = mol->residuebegin;
        for(j=0; j < mol->res_num; j++){
            if(j != 0) residue = residue->next;

            atom = residue->atombegin;

            if(strcmp(chainID_old, atom->chainID) != 0 ||
               strcmp(chainID_mod_old, atom->chainID_mod) != 0){
                print_flg = 1;
            }
            else print_flg = 0;

            if(print_flg == 1){
                if(atom->chainID[0] != atom->chainID_mod2[0]){
                    if(header_print_flg == 0){
                        fprintf(fp, "%s", header);
                        header_print_flg = 1;
                    }
                    fprintf(fp, "REMARK      INPUT:'%1s', OUTPUT:'%1s'\n",
                        atom->chainID, atom->chainID_mod2);
                    print_flg = 0;
                }
            }


            strcpy(chainID_old, atom->chainID);
            strcpy(chainID_mod_old, atom->chainID_mod);
        }
    }
    /* add end */

    return chainID_idx;
}

/* add start by tahara 2013/01/08 */
/*--------------------------------------------------------*/
/* int tgWriteRemarkLine()                                */
/*--------------------------------------------------------*/
int tgWriteRemarkLine(FILE* fp,
                      char* r_name,
                      int remark_counter,
                      int remark_name_counter)
{

    long remark_temp = 0;
    int remark_length = 0;

    if(r_name[0] == '\0'){
        return remark_counter;
    }
    else if((!strncmp(r_name,remmes.remark_position_name[remark_name_counter],3)) ||
            (!strncmp(r_name,"OTHER",5))) {
        fprintf(fp, "%s", remmes.remark_message[remark_counter]);
        while(1) {
            remark_temp = remmes.remark_position[remark_counter+1]
                        - remmes.remark_position[remark_counter];

            remark_length = strlen(remmes.remark_message[remark_counter+1]);
            remark_counter++;
            if((int)remark_temp == remark_length) {
                fprintf(fp, "%s", remmes.remark_message[remark_counter]);
             } else {
                break;
             }
        }
    }
    return remark_counter;
}
/* add end by tahara 2013/01/08 */

/* add, 2013.11.28 */
/*--------------------------------------------------------*/
/* void resetChainIDs(TGASystemPtr amber, resConvTblPtr resCnvTbl) */
/*--------------------------------------------------------*/
void resetChainIDs(TGASystemPtr amber, resConvTblPtr resCnvTbl)
{
    int i, j, k, l;
    TGAAtomPtr atom;
    TGAResiduePtr res;
    TGAMolPtr mol;
    TGASSBondPtr ss;
    resConvPtr resCnv;

    /* reset chain id information for ATOM struvt */
    mol = amber->begin;
    for(i=0; i < amber->mol_num; i++){
        if(i != 0) mol = mol->next;

        res = mol->residuebegin;
        for(j=0; j < mol->res_num; j++){
            if(j != 0) res = res->next;

            atom = res->atombegin;
            for(k=0; k < res->atom_num; k++){
                if(k != 0) atom = atom->next;

                resCnv = resCnvTbl->resCnvbegin;
                for(l=0; l < resCnvTbl->res_num; l++){
                    if(l != 0) resCnv = resCnv->next;

                    if(resCnv->chainID_mod[0] == atom->chainID[0] &&
                       resCnv->resSeq_mod == atom->residue_num_org &&
                       strncmp(resCnv->res_name, atom->residue_name, 3) == 0 &&
                       resCnv->iCode_mod[0] == atom->iCode[0]){

                        strcpy(atom->chainID, resCnv->chainID_org);
                        atom->residue_num_org = resCnv->resSeq_org;
                        strcpy(atom->iCode, resCnv->iCode_org);

                        break;
                    }
                }
            }
        }
    }

    /* reset chain id information fro SSBOND struct */
    ss = amber->ssbond;
    for(i=0; i < amber->SSbond_num; i++){
        resCnv = resCnvTbl->resCnvbegin;
        for(l=0; l < resCnvTbl->res_num; l++){
            if(l != 0) resCnv = resCnv->next;
            if(strncmp(resCnv->res_name, "CYS", 3) != 0) continue;

            if(resCnv->chainID_mod[0] == ss->chainID1_org[i][0] &&
               resCnv->resSeq_mod == ss->seqNum1_org[i] &&
               (resCnv->iCode_mod[0] == ss->iCode1[i][0] ||
               (resCnv->iCode_mod[0] == '\0' && ss->iCode1[i][0] == ' '))){

                strcpy(ss->chainID1_org[i], resCnv->chainID_org);
                ss->seqNum1_org[i] = resCnv->resSeq_org;
                strcpy(ss->iCode1[i], resCnv->iCode_org);
            }
            else if(resCnv->chainID_mod[0] == ss->chainID2_org[i][0] &&
               resCnv->resSeq_mod == ss->seqNum2_org[i] &&
               (resCnv->iCode_mod[0] == ss->iCode2[i][0] ||
               (resCnv->iCode_mod[0] == '\0' && ss->iCode2[i][0] == ' '))){

                strcpy(ss->chainID2_org[i], resCnv->chainID_org);
                ss->seqNum2_org[i] = resCnv->resSeq_org;
                strcpy(ss->iCode2[i], resCnv->iCode_org);
            }
        }
    }

}
/* add end */

/*--------------------------------------------------------*/
/* int tgOutputCoordinate(TGASystemPtr, TGASystemPtr)     */
/*--------------------------------------------------------*/
int tgOutputCoordinate(TGASystemPtr amber, TGASystemPtr pdb_metals,
        resConvTblPtr resCnvTbl)
{
    FILE *fp;
    char buf[BUFSIZE + 1];
    char buf1[BUFSIZE + 1];
    char res_name[5], res_name2[5];
    FILE *fplist = NULL;
    char log_file[] = "Original_aminoacid_No.dat";
    int resnum_org = 0;
    char chainID_org[2], chainID_mod[2];
    char chainID_old[2];
    int i, j, k, l, ll;
    int atomnum, resnum;

    double pre_c_coord[3];
    int pre_c_flg;
    int m;
    /* add 2011.6.24 */
    TGAAtomPtr atom1, atom2;
    /* add 2011.10.21 */
    int temp_res_num_org = 0;
    int temp_res_num_modify = 0;
    char temp_chain_id_modify[2];
    char temp_chain_id[2];
    int chainID_idx = 0;

    char chk_symbol_dmy[ATOM_NAME_MAX + 1];

    int pre_resnum = 0;  /* output for TER */
    char pre_resname[5]; /* output for TER */
    char pre_chainID[2]; /* output for TER */
    char pre_iCode[2];   /* output for TER */


    char presto_chainIDs[] = PRESTO_CHAINIDS;
    /* add start by tahara 2013/01/04 */
    int old_remark_counter = 0;
    int new_remark_counter = 0;
    int remark_name_counter = 0;
    /* add end by tahara 2013/01/04 */

    puts("\n INFORMATION> tgOutputCoordinate");
    /* modify 2011.5.26, start */
    /* output pdb which is presto format */
    if(sdata.out_pdb_style == 1){
        puts("             Output pdb file which is presto style.");
    }
    /* output pdb which is standard pdb format */
    else{
        puts("             Output pdb file which is standard style.");
    }
    /* modify end */
    /* add 2011.6.18, start */
    if(sdata.out_pdb_style == 1){
        resnum_org = 0;
        strcpy(chainID_org, " \0");
        strcpy(chainID_mod, " \0");
    }
    /* add end */

    pre_resname[0] = '\0';
    pre_chainID[0] = '\0';
    pre_iCode[0] = '\0';

    amber->modeo = 1;

    /* reset chainids */
    resetChainIDs(amber, resCnvTbl);

    /******************************
     * Output PDB File
     ******************************/
    strcpy(buf, sdata.outcoordfile);
    if((fp = fopen(buf, "w")) == NULL) {
        tgError(OUTPDB_FILE_OPEN_ERROR, buf, NULL, NULL, DMY_INT,
            DMY_INT);
    }
    /* add 2011.6.18, start */
    if(sdata.out_pdb_style == 1){
        if((fplist = fopen(log_file, "wb")) == NULL) {
            tgError(BINARY_FILE_OPEN_ERROR, buf1, NULL, NULL, DMY_INT,
            DMY_INT);
        }
        fprintf(fplist, "Atom         Residue                 ChainID\n");
        fprintf(fplist, "  Num  Name  Name  serialNum orgNum  org  mod\n");
    }
    amber->mol = amber->begin;
    /* add end */

    /* print SSBOND information, -stdpdb only */
    if(sdata.out_pdb_style == 1){
        /* print REMARK line */
        tgOutputRemarkLine(fp, amber);
    }

















    // handle circles

    amber->mol = amber->begin;
    for(i=0; i < amber->mol_num; i++){
        if(i != 0) amber->mol = amber->mol->next;
        if(amber->mol->icrclf != 1) continue;

        for(k=0; k < TEN; k++){
            atom1 = atom2 = NULL;
            amber->mol->residue = amber->mol->residuebegin;
            for(j=0; j < amber->mol->res_num; j++){
                if(j != 0) amber->mol->residue = amber->mol->residue->next;
                amber->mol->residue->atom = amber->mol->residue->atombegin;

                if (amber->mol->circ_start[k] == amber->mol->residue->atom->residue_num_mod) {
                    atom1 = amber->mol->residue->atom;
                }
                if (amber->mol->circ_end[k] == amber->mol->residue->atom->residue_num_mod){
                    atom2 = amber->mol->residue->atom;
                    temp_res_num_org = atom2->residue_num_org;
                    temp_res_num_modify = atom2->residue_num_mod;
                    strcpy(temp_chain_id, atom2->chainID);
                    strcpy(temp_chain_id_modify, atom2->chainID_mod);

                    while(1) {
                        if(strcmp(atom2->atom_name, "C") == 0) break;
                        atom2 = atom2->next;
                    }
                }
                if(atom1 != NULL && atom2 != NULL) break;
            }
            if(atom1 != NULL && atom2 != NULL){
                strncpy(res_name, atom1->residue_name, 3);
                strncpy(res_name2, atom2->residue_name, 3);
                res_name[3] = '\0';
                res_name2[3] = '\0';

                /* myPresto original format */
                if(sdata.out_pdb_style == 1){

                    fprintf(fp, "%-6s    %4s   %3s %1s%4d%1s             %4s   %3s %1s%4d%1s  %6s %6s  %4.2lf\n",
                        "LINK",
                        atom1->atom_name, res_name,
                        atom1->chainID_mod2,
                        atom1->residue_num_mod, " ",
                        atom2->atom_name, res_name2,
                        atom2->chainID_mod2,
                        temp_res_num_modify, " ",
                        PDB_FORMAT_SYMMETRY, PDB_FORMAT_SYMMETRY,
                        PDB_FORMAT_LINK_LENGTH);
                    /* modified 2012/08/11 */
                }
                /* standard pdb format */
                else if(sdata.out_pdb_style == 2){
                    fprintf(fp, "%-6s    %4s   %3s %1s%4d%1s             %4s   %3s %1s%4d%1s  %6s %6s  %4.2lf\n",
                        "LINK",
                        atom1->atom_name, res_name,
                        atom1->chainID,
                        atom1->residue_num_org, atom1->iCode,
                        atom2->atom_name, res_name2,
                        temp_chain_id,
                        temp_res_num_org, atom2->iCode,
                        PDB_FORMAT_SYMMETRY, PDB_FORMAT_SYMMETRY,
                        PDB_FORMAT_LINK_LENGTH);
                    /* modified 2012/08/11 */
                }
            }
            else{
            }
        }
    }













    // handle ssbond

    for(i=0; i < amber->SSbond_num; i++){
        /* add start by tahara 2013/01/09 */
        while(remmes.remark_ss_position > remmes.remark_position[old_remark_counter]) {
            /* add start by tahara 2013/01/10 */
            if(remmes.remark_position[old_remark_counter] == 0)break;
            /* add end by tahara 2013/01/10 */
            new_remark_counter = tgWriteRemarkLine(fp,
                    "OTHER",
                    old_remark_counter,
                    remark_name_counter);
            if(old_remark_counter != new_remark_counter) {
                old_remark_counter = new_remark_counter;
            }
        }
        /* add end by tahara 2013/01/09 */
        /* myPresto Original format */
        if(sdata.out_pdb_style == 1){

            /* Searching SSBOND pare amust use chainID_org as chainID */
            fprintf(fp,  "SSBOND %3d %3s %1s %4d%1s   %3s %1s %4d%1s                       %6s %6s  %4.2lf\n",
                i+1,
                "CYS", amber->ssbond->chainID1_mod[i],
                amber->ssbond->seqNum1_mod[i],
                " ",
                "CYS", amber->ssbond->chainID2_mod[i],
                amber->ssbond->seqNum2_mod[i],
                " ",
                PDB_FORMAT_SYMMETRY, PDB_FORMAT_SYMMETRY,
                PDB_FORMAT_SSBOND_LENGTH);


            /* modified 2012/08/11 */
        }
        /* standard pdb format */
        else if(sdata.out_pdb_style == 2){
            fprintf(fp,  "SSBOND %3d %3s %1s %4d%1s   %3s %1s %4d%1s                       %6s %6s  %4.2lf\n",
                i+1,
                "CYS", amber->ssbond->chainID1_org[i],
                amber->ssbond->seqNum1_org[i],
                amber->ssbond->iCode1[i],
                "CYS", amber->ssbond->chainID2_org[i],
                amber->ssbond->seqNum2_org[i],
                amber->ssbond->iCode2[i],
                PDB_FORMAT_SYMMETRY, PDB_FORMAT_SYMMETRY,
                PDB_FORMAT_SSBOND_LENGTH);
            /* modified 2012/08/11 */
        }
    }

    amber->mol = amber->begin;

    chainID_idx = 0;























    // do something unknown

    l = 1;
    for(i = 0; i < amber->mol_num; i++, chainID_idx++) {

        pre_c_flg = 0;
        for(m = 0; m < 3; m++) {
            pre_c_coord[m] = 0.0;
        }
        amber->mol->residue = amber->mol->residuebegin;
        for(j = 0; j < amber->mol->res_num; j++) {
            amber->mol->residue->atom = amber->mol->residue->atombegin;

            /* add 2011.6.23, start */
            /* modified 20130327
            if(j != 0 &&
                strcmp(chainID_old, amber->mol->residue->atom->chainID) != 0){
                if(sdata.slctform == 1){
                    atomnum = l % 10000;
                    fprintf(fp, "TER   %5d      %-4s%1s%4d%1s\n", atomnum,
                    pre_resname, pre_chainID,pre_resnum, pre_iCode);
                    l++;
                    chainID_idx++;
                }
            }
            */
            /* add end */
            for(k = 0; k < amber->mol->residue->atom_num; k++) {
                /* add, modify 2011.6.17, start */
                if(sdata.out_pdb_style == 1){
                    if(strlen(amber->mol->residue->atom->atom_name) > 3) {
                        strcpy(buf, "ATOM  %5d %4s %-4s%c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n");
                    } else {
                        strcpy(buf, "ATOM  %5d  %-3s %-4s%c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n");
                    }
                } else{
                    if(strlen(amber->mol->residue->atom->atom_name) > 3) {
                        strcpy(buf, "ATOM  %5d %4s %-4s%1s%4d%1s   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n");
                    } else {
                        strcpy(buf, "ATOM  %5d  %-3s %-4s%1s%4d%1s   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n");
                    }
                }

                if((amber->mol->residue->atom->coord->x > -0.0001)
                    && (amber->mol->residue->atom->coord->x < 0.0001))
                    amber->mol->residue->atom->coord->x = 0.0000;
                if((amber->mol->residue->atom->coord->y > -0.0001)
                    && (amber->mol->residue->atom->coord->y < 0.0001))
                    amber->mol->residue->atom->coord->y = 0.0000;
                if((amber->mol->residue->atom->coord->z > -0.0001)
                    && (amber->mol->residue->atom->coord->z < 0.0001))
                    amber->mol->residue->atom->coord->z = 0.0000;

                /* modify 2011.6.17, start */
                if(sdata.out_pdb_style == 1){
                    if(strcmp(amber->mol->residue->atom->residue_name, "CYT") == 0)
                        strcpy(amber->mol->residue->atom->residue_name, "  C\0");
                    else if(strcmp(amber->mol->residue->atom->residue_name, "DCY") == 0)
                        strcpy(amber->mol->residue->atom->residue_name, " DC\0");
                    else if(strcmp(amber->mol->residue->atom->residue_name, "GUA") == 0)
                        strcpy(amber->mol->residue->atom->residue_name, "  G\0");
                    else if(strcmp(amber->mol->residue->atom->residue_name, "DGU") == 0)
                        strcpy(amber->mol->residue->atom->residue_name, " DG\0");
                    else if(strcmp(amber->mol->residue->atom->residue_name, "ADE") == 0)
                        strcpy(amber->mol->residue->atom->residue_name,"  A\0");
                    else if(strcmp(amber->mol->residue->atom->residue_name,"DAD") == 0)
                        strcpy(amber->mol->residue->atom->residue_name," DA\0");
                    else if(strcmp(amber->mol->residue->atom->residue_name, "URA") == 0)
                        strcpy(amber->mol->residue->atom->residue_name,"  U\0");
                    else if(strcmp(amber->mol->residue->atom->residue_name, "DTH") == 0)
                        strcpy(amber->mol->residue->atom->residue_name, " DT\0");
                }
                /* add 2011.5.26, start */
                /* output standard pdb style */
                if(sdata.out_pdb_style == 2){
                    if(strlen(amber->mol->residue->atom->residue_name) > 3){
                        strncpy(res_name, amber->mol->residue->atom->residue_name, 3);
                        res_name[3] = '\0';
                    } else{
                        strncpy(res_name, amber->mol->residue->atom->residue_name, 4);
                        res_name[4] = '\0';
                    }
                }
                /* output presto pdb style */
                else{
                    strncpy(res_name,
                        amber->mol->residue->atom->residue_name, 4);
                    res_name[4] = '\0';
                }

                atomnum = l % 10000;
                resnum  = amber->mol->residue->atom->residue_num_mod % 10000;


                /* add 2011.6.18, start */
                /* set values of complemented data */
                /* mainly hydrogen atoms */
                if(amber->mol->residue->atom->residue_num_org == 0){
                    if(k != 0){
                        amber->mol->residue->atom->residue_num_org = resnum_org;
                    }
                    /* modified 20130327 */
                    if(amber->mol->residue->atom->chainID[0] == '\0'){
                        strcpy(amber->mol->residue->atom->chainID, chainID_org);
                        strcpy(amber->mol->residue->atom->chainID_mod,
                            chainID_mod);
                    }
                    /* modified end */
                }
                resnum_org = amber->mol->residue->atom->residue_num_org;
                strcpy(chainID_org, amber->mol->residue->atom->chainID);
                strcpy(chainID_mod,
                    amber->mol->residue->atom->chainID_mod);

                for(ll=0; ll < amber->SSbond_num; ll++){
                    if(amber->ssbond->chainID1_mod[ll][0] ==
                        amber->mol->residue->atom->chainID_mod[0] &&
                        amber->ssbond->seqNum1_mod[ll] ==
                        amber->mol->residue->atom->residue_num_mod){
                        strcpy(amber->mol->residue->atom->iCode,
                            amber->ssbond->iCode1[ll]);
                        break;
                    }
                    else if(amber->ssbond->chainID2_mod[ll][0] ==
                        amber->mol->residue->atom->chainID_mod[0] &&
                        amber->ssbond->seqNum2_mod[ll] ==
                        amber->mol->residue->atom->residue_num_mod){
                        strcpy(amber->mol->residue->atom->iCode,
                            amber->ssbond->iCode2[ll]);
                        break;
                    }
                }

                /* add start by tahara 2013/01/08 */
                new_remark_counter = tgWriteRemarkLine(fp,
                        amber->mol->residue->atom->residue_name,
                        old_remark_counter,
                        remark_name_counter);
                if(old_remark_counter != new_remark_counter) {
                    remark_name_counter++;
                    old_remark_counter = new_remark_counter;
                }
                /* add end by tahara 2013/01/08 */

                if(sdata.out_pdb_style == 1){
                  fprintf(fplist, "%5d %5s  %-5s %5d  %5d%1s     %2s   %2s\n",
                      atomnum,
                      amber->mol->residue->atom->atom_name,
                      res_name,
                      amber->mol->residue->atom->residue_num_mod,
                      amber->mol->residue->atom->residue_num_org,
                      amber->mol->residue->iCode,
                      amber->mol->residue->atom->chainID,
                      amber->mol->residue->atom->chainID_mod2);
                }
                /* add end */

                /* modify 2011.5.26, start */
                /* write myPresto format (non standard pdb format) */
                if(sdata.out_pdb_style == 1){
                    fprintf(fp, buf,
                        atomnum,
                        amber->mol->residue->atom->atom_name,
                        res_name,
                        /* added 2012/08/08 */
                        presto_chainIDs[chainID_idx%26],
                        /* added 2012/08/08 */
                        resnum,
                        amber->mol->residue->atom->coord->x,
                        amber->mol->residue->atom->coord->y,
                        amber->mol->residue->atom->coord->z,
                        amber->mol->residue->atom->mass,
                        amber->mol->residue->atom->charge);

                    /* added 2012/08/11 */
                    pre_resnum = resnum;
                    strcpy(pre_resname, res_name);
                    sprintf(pre_chainID, "%c", presto_chainIDs[chainID_idx%26]);
                    strcpy(pre_iCode, " ");
                    /* added 2012/08/11 */

                }
                /* write pdb style format */
                else{
                    fprintf(fp, buf,
                        atomnum,
                        amber->mol->residue->atom->atom_name,
                        res_name,
                        amber->mol->residue->atom->chainID,
                        amber->mol->residue->atom->residue_num_org,
                        amber->mol->residue->iCode,
                        amber->mol->residue->atom->coord->x,
                        amber->mol->residue->atom->coord->y,
                        amber->mol->residue->atom->coord->z,
                        amber->mol->residue->atom->mass,
                        amber->mol->residue->atom->charge);

                    /* added 2012/08/11 */
                    pre_resnum = amber->mol->residue->atom->residue_num_org;
                    strcpy(pre_resname, res_name);
                    strcpy(pre_chainID, amber->mol->residue->atom->chainID);
                    strcpy(pre_iCode, amber->mol->residue->iCode);
                    /* added 2012/08/11 */
                }
                /* modify end */

                if(strcmp(amber->mol->residue->atom->atom_name,
                        "C") == 0) {
                    pre_c_flg = 1;
                    pre_c_coord[0] =
                        amber->mol->residue->atom->coord->x;
                    pre_c_coord[1] =
                        amber->mol->residue->atom->coord->y;
                    pre_c_coord[2] =
                        amber->mol->residue->atom->coord->z;
                }
                /* add 2011.5.26, start */
                if(strcmp(amber->mol->residue->atom->atom_name,
                    "OXT") == 0){
                    pre_c_flg = 0;
                }
                /* add end */
                if(pre_c_flg == 1
                    && strcmp(amber->mol->residue->atom->atom_name,
                        "N") == 0) {
                    if(sqrt((pre_c_coord[0] -
                                amber->mol->residue->atom->coord->x) *
                            (pre_c_coord[0] -
                                amber->mol->residue->atom->coord->x) +
                            (pre_c_coord[1] -
                                amber->mol->residue->atom->coord->y) *
                            (pre_c_coord[1] -
                                amber->mol->residue->atom->coord->y) +
                            (pre_c_coord[2] -
                                amber->mol->residue->atom->coord->z) *
                            (pre_c_coord[2] -
                                amber->mol->residue->atom->coord->z)) >
                        3.0) {
                        puts("");
                        puts(" WARNING> tgOutputCoordinate");
                        printf("          It is too away between the ");
                        /* modify 2011.5.27, start */
                        printf("%dth and %dth residue of %dth molcule!\n",
                            amber->mol->residue->atom->residue_num_mod - 1,
                            amber->mol->residue->atom->residue_num_mod,
                            i + 1);
                        /* modify end */
                    }
                }

                l++;

                amber->mol->residue->atom =
                    amber->mol->residue->atom->next;
            }
            /* add 2011.6.23, start */
            strcpy(chainID_old, chainID_org);
            /* add end */

            amber->mol->residue = amber->mol->residue->next;
        }
        atomnum = l % 10000;
        fprintf(fp, "TER   %5d      %-4s%1s%4d%1s\n", atomnum,
            pre_resname, pre_chainID,pre_resnum, pre_iCode);
        l++;
        amber->mol = amber->mol->next;
    }


























    /*******************************/
    /* Output coordinate of Metals */
    /*******************************/
    if(sdata.metals == 1){

        pdb_metals->mol = pdb_metals->begin;
        for(i=0; i < pdb_metals->mol_num; i++, chainID_idx++){
            pdb_metals->mol->residue = pdb_metals->mol->residuebegin;

            /* add start by tahara 2013/01/08 */
            new_remark_counter = tgWriteRemarkLine(fp,
                    pdb_metals->mol->residue->atom->residue_name,
                    old_remark_counter,
                    remark_name_counter);
            if(old_remark_counter != new_remark_counter) {
                remark_name_counter++;
                old_remark_counter = new_remark_counter;
            }
            /* add end by tahara 2013/01/08 */

            for(j=0; j < pdb_metals->mol->res_num; j++){
                pdb_metals->mol->residue->atom =
                                   pdb_metals->mol->residue->atombegin;

                for(k = 0; k < pdb_metals->mol->residue->atom_num;k++){
                    if(strlen
                      (pdb_metals->mol->residue->atom->atom_name) > 3){
                        strcpy(buf, "HETATM%5d %4s %-4s%1c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n");

                    }
                    else {
                        strcpy(buf, "HETATM%5d  %-3s %-4s%1c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n");

                    }
                    atomnum = l % 10000;
                    resnum  = pdb_metals->mol->residue->num % 10000;

                    if(sdata.out_pdb_style == 1){
                        fprintf(fp, buf,atomnum,
                            pdb_metals->mol->residue->atom->atom_name,
                            pdb_metals->mol->residue->atom->residue_name,
                            presto_chainIDs[chainID_idx%26],
                            resnum,
                            pdb_metals->mol->residue->atom->coord->x,
                            pdb_metals->mol->residue->atom->coord->y,
                            pdb_metals->mol->residue->atom->coord->z,
                            pdb_metals->mol->residue->atom->mass,
                            pdb_metals->mol->residue->atom->charge);

                    }
                    else{
                        fprintf(fp, buf,atomnum,
                            pdb_metals->mol->residue->atom->atom_name,
                            pdb_metals->mol->residue->atom->residue_name,
                            ' ',
                            resnum,
                            pdb_metals->mol->residue->atom->coord->x,
                            pdb_metals->mol->residue->atom->coord->y,
                            pdb_metals->mol->residue->atom->coord->z,
                            pdb_metals->mol->residue->atom->mass,
                            pdb_metals->mol->residue->atom->charge);

                    }
                    /* modified 2012/08/13 */


                    /* added 2012/08/11 */
                    pre_resnum = resnum;
                    strcpy(pre_resname, pdb_metals->mol->residue->atom->residue_name);
                    if(sdata.out_pdb_style == 1){
                        sprintf(pre_chainID, "%c", presto_chainIDs[chainID_idx%26]);
                    }
                    else{
                        strcpy(pre_chainID, " ");
                    }
                    strcpy(pre_iCode, " ");
                    /* added 2012/08/11 */

                    l++;
                    if(pdb_metals->mol->residue->atom->next == NULL)
                        break;

                    if(sdata.out_pdb_style == 1){
                       if(pdb_metals->mol->residue->next != NULL){
                            /* if next atom exists. */
                           if(tgCheckAtomNameOfMetals(pdb_metals->mol->residue->atom->atom_name,
                               pdb_metals->mol->residue->atom->residue_name,
                               chk_symbol_dmy) == 1){
                               if(strcmp(pdb_metals->mol->residue->atom->atom_name, "Na") != 0 &&
                                  strcmp(pdb_metals->mol->residue->atom->atom_name, "NA") != 0 &&
                                  strcmp(pdb_metals->mol->residue->atom->atom_name, "Cl") != 0 &&
                                  strcmp(pdb_metals->mol->residue->atom->atom_name, "CL") != 0){
                                   atomnum = l % 10000;
                                   fprintf(fp, "TER   %5d      %-4s%1s%4d%1s\n", atomnum,
                                   pre_resname, pre_chainID,pre_resnum, pre_iCode);

                                   chainID_idx++;
                                   l++;
                               }
                           }
                        }
                    }

                    pdb_metals->mol->residue->atom = pdb_metals->mol->residue->atom->next;

                }
                if(pdb_metals->mol->residue->next == NULL) break;
                pdb_metals->mol->residue = pdb_metals->mol->residue->next;
            }
            if(pdb_metals->mol->res_num > 0){
                atomnum = l % 10000;
                fprintf(fp, "TER   %5d      %-4s%1s%4d%1s\n", atomnum,
                pre_resname, pre_chainID,pre_resnum, pre_iCode);
                l++;
                /* modified 2012/08/11 */
            }
            if(pdb_metals->mol->next == NULL) break;
            pdb_metals->mol = pdb_metals->mol->next;
        }
    }
    fclose(fp);
    /* add 2011.6.18, start */
    if(sdata.out_pdb_style == 1){
        fclose(fplist);
    }
    /* add end */

    return TG_OK;
}
