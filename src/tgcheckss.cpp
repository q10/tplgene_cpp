#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "TP_error.h"
#include "TP_define.h"
#include "TP_struct.h"
#include "TP_AmberDB.h"
#include "TP_Env.h"

int tgError(int, char*, char*, char*, int, int);
void *Malloc(int);
int Free(void*);
int tgNextMalloc(TGASystemPtr, int);
int getChainGroup(TGASystemPtr mod, int ssbond_num);
int tgCheckAtomNameOfMetals(char *, char *,char *);
/* added by Sumita 2012/08/28 */
int tgContainsSameSSAtom(TGASSBondPtr, int, int, int);
/* add 2013.11.14 */
void tgFreeAmberStruct(TGASystemPtr sys);
int InitMemoryAllocate(TGASystemPtr mod);

extern SDATA sdata;


#define SSBOND_THRESHOLD           2.5
#define S_METAL_DIST_THRESHOLD     2.5

/*--------------------------------------------------------*/
/*  double tgCalcDistance                                 */
/*--------------------------------------------------------*/
double tgCalcSquareDistance(double* coord1, double* coord2)
{
    return (coord1[0]-coord2[0])*(coord1[0]-coord2[0]) +
           (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) +
           (coord1[2]-coord2[2])*(coord1[2]-coord2[2]);

}

/*--------------------------------------------------------*/
/*  int tgCheckSS()                                       */
/*--------------------------------------------------------*/
int tgCheckSS(char cID1[MOLNUMMAX][2], char cID2[MOLNUMMAX][2],
              char cID1_org[MOLNUMMAX][2], char cID2_org[MOLNUMMAX][2],
              int sNum1_mod[MOLNUMMAX], int sNum2_mod[MOLNUMMAX],
              int sNum1[MOLNUMMAX], int sNum2[MOLNUMMAX],
              char icod1[MOLNUMMAX][2], char icod2[MOLNUMMAX][2],
              char molCID[MOLNUMMAX][2], char molCID_mod[MOLNUMMAX][2]) {
    int ret_code;
    int ssbond_num;
    char **chainID, **chainID_modify;
    char textbuf[BUFSIZE + 1];
    char buf[BUFSIZE + 1];
    char dmychar1[256];
    char dmychar2[256];
    char dmychar3[256];
    char checkID[2];
    FILE *fp, *fpw;
    TGASystemPtr mod;
    int dmynum1, resnum;
    int atom_num, res_num, mol_num, total_atom_num;
    int *check;
    int *before_res_num;
    int seqNum1, seqNum2;
    int end, flg;
    int i, j, k, l, m;
    int val;
    int chainID_space_num;

    char bef_chainID[2], bef_iCode[2];
    int bef_resnum = 0;
    char temp1, temp2;
    int ter_flg;

    int all_cys_num;
    int cys_SeqNum[MAXSSBOND];
    char cys_ChainID[MAXSSBOND][2];
    char cys_iCode[MAXSSBOND][2];
    double cys_SgCoord[MAXSSBOND][3];
    /* added by Sumita 2012/08/28 */
    int cys_valid[MAXSSBOND]; /* 0:invalid, 1:valid */
    int print_metal_header;
    /* added by Sumita 2012/08/28 */

    double ssbond_thr2;
    double s_metal_thr2;

    char het_atomname[ATOM_NAME_MAX + 1];
    char het_resname[ATOM_NAME_MAX + 1];
    char het_symbol_dmy[ATOM_NAME_MAX + 1];
    double metal_coord[METALS_MAX][3];
    char metal_atom_name[METALS_MAX][ATOM_NAME_MAX + 1];
    int metal_num = 0;

    /******************************
     * Initial setting
     ******************************/
    bef_chainID[0] = '\0'; bef_chainID[1] = '\0';
    bef_iCode[0] = '\0'; bef_iCode[1] = '\0';

    ret_code = 0;
    strcpy(checkID, "\0");

    strcpy(buf, "\0");
    strcat(buf, sdata.coordfile);

    for(i=0; i < METALS_MAX; i++){
        /* modify, 2013.11.14 */
        memset(metal_atom_name[i], '\0', ATOM_NAME_MAX + 1);
        /* mod end */
    }

    /* added by Sumita 2012/08/28 */
    for(i=0; i < MAXSSBOND; i++){
        cys_valid[i] = 1;
    }
    print_metal_header = 0;
    /* added by Sumita 2012/08/28 */

    mod = (TGASystemPtr) Malloc(sizeof(TGASystem));
    mod->mol = (TGAMolPtr) Malloc(sizeof(TGAMol));
    mod->begin = mod->end = mod->mol;

    InitMemoryAllocate(mod);

    /******************************
     * File open
     ******************************/
    if((fp = fopen(buf, "r")) == NULL) {
        tgError(READSEQ_FILE_OPEN_ERROR, buf, NULL, NULL, DMY_INT,
            DMY_INT);
    }

    /******************************
     * Count number of SSBOND
     ******************************/
    if(fgets(textbuf, BUFSIZE, fp) == 0) {
        printf("\n ERROR> tgCheckSS\n");
        tgError(READ_FORMAT_ERROR, buf, NULL, NULL, DMY_INT, DMY_INT);
    }
    ssbond_num = 0;

    // READ in all SSBONDS from file into vector of SSSBonds
/*
    while(WHILE_LOOP) {
        if(strncmp(textbuf, "SSBO", PRE_WORD_NUM) == 0) {
            mod->ssbond->iCode1[ssbond_num][0] = textbuf[21];
            mod->ssbond->iCode1[ssbond_num][1] = '\0';
            textbuf[21] = ' ';

            mod->ssbond->iCode2[ssbond_num][0] = textbuf[35];
            mod->ssbond->iCode2[ssbond_num][1] = '\0';
            textbuf[35] = ' ';

            if(mod->ssbond->iCode2[ssbond_num][0] == '\n' ||
               mod->ssbond->iCode2[ssbond_num][0] == '\0'){
                mod->ssbond->iCode2[ssbond_num][0] = ' ';
                mod->ssbond->iCode2[ssbond_num][1] = '\0';
            }

            mod->ssbond->chainID1_mod[ssbond_num][0] = textbuf[15];
            mod->ssbond->chainID1_mod[ssbond_num][1] = '\0';
            strncpy(&textbuf[15], " ", 1);
            mod->ssbond->chainID2_mod[ssbond_num][0] = textbuf[29];
            mod->ssbond->chainID2_mod[ssbond_num][1] = '\0';
            strncpy(&textbuf[29], " ", 1);

            sscanf(textbuf, "%s%d%s%d%s%d",
                dmychar1, &dmynum1, dmychar2,
                &mod->ssbond->seqNum1_org[ssbond_num],
                dmychar3, &mod->ssbond->seqNum2_org[ssbond_num]);

            if(mod->ssbond->seqNum1_org[ssbond_num] == 0 ||
                mod->ssbond->seqNum2_org[ssbond_num] == 0) {
                tgError(SSBOND_FORMAT_ERROR, textbuf, NULL, NULL, DMY_INT,
                    DMY_INT);
            }
            ssbond_num++;
        }
        if(fgets(textbuf, BUFSIZE, fp) == 0) {
            break;
        }
    }

    if(ssbond_num == 0 && sdata.ss_detect != 2 ) {
        fclose(fp);
        tgFreeAmberStruct(mod);
        return 0;
    }
*/
    rewind(fp);

    /******************************
     * ATOM
     ******************************/
    atom_num = res_num = mol_num = 0;
    ter_flg = 0;

    if(fgets(textbuf, BUFSIZE, fp) == 0) {
        printf("\n ERROR> tgCheckSS\n");
        tgError(READ_FORMAT_ERROR, buf, NULL, NULL, DMY_INT, DMY_INT);
    }
    mod->mol = mod->begin;
    mod->mol->residue = mod->mol->residuebegin;
    mod->mol->residue->atom = mod->mol->residue->atombegin;

    while(WHILE_LOOP) {

        // if atom, build to pdbstruct
        if(strncmp(textbuf, "ATOM", PRE_WORD_NUM) == 0) {

            /* get residue number */
            sscanf(&textbuf[22], "%d", &resnum);

            ter_flg = 0;
            if(atom_num == 0){
                mod->mol->residue->atom =
                    mod->mol->residue->atombegin;
            } else{
                /* next mol */
                if(textbuf[21] != bef_chainID[0]){
                    mod->mol->res_num = res_num;
                    mod->mol->residue->atom_num = atom_num;

                    mod->mol->residue->atomend = mod->mol->residue->atom;
                    mod->mol->residueend = mod->mol->residue;

                    ret_code = tgNextMalloc(mod, OBJECT_MOL);

                    mol_num++;
                    res_num = 0;
                    atom_num = 0;
                }
                /* next residue */
                else if(bef_resnum != resnum ||
                   textbuf[26] != bef_iCode[0]){
                    mod->mol->residue->atom_num = atom_num;

                    mod->mol->residue->atomend = mod->mol->residue->atom;
                    ret_code = tgNextMalloc(mod, OBJECT_RESIDUE);

                    res_num++;
                    atom_num = 0;
                }
                /* next atom */
                else{
                    ret_code = tgNextMalloc(mod, OBJECT_ATOM);
                }
            }

            bef_resnum = resnum;
            bef_chainID[0] = textbuf[21]; bef_chainID[1] = '\0';
            bef_iCode[0] = textbuf[26]; bef_iCode[1] = '\0';

            mod->mol->residue->atom->altLoc[0] = textbuf[16];
            mod->mol->residue->atom->altLoc[1] = '\0';
            mod->mol->residue->atom->chainID[0] = textbuf[21];
            mod->mol->residue->atom->chainID[1] = '\0';
            mod->mol->residue->atom->iCode[0] = textbuf[26];
            mod->mol->residue->atom->iCode[1] = '\0';
            textbuf[16] = ' ';
            textbuf[21] = ' ';
            textbuf[26] = ' ';

            sscanf(textbuf, "%s%d%s%s%d%lf%lf%lf",
                mod->mol->residue->atom->type, &dmynum1,
                mod->mol->residue->atom->atom_name,
                mod->mol->residue->atom->residue_name,
                &mod->mol->residue->atom->residue_num_org,
                &mod->mol->residue->atom->coord->x,
                &mod->mol->residue->atom->coord->y,
                &mod->mol->residue->atom->coord->z);

                mod->mol->residue->atom->residue_num_mod = res_num + 1;

                atom_num++;

        }

        // if CIRC, set the molecule "ss_circ_flag", and unset the ter_flag
        else if(strncmp(textbuf, "CIRC", PRE_WORD_NUM) == 0){
            ter_flg = 0;
            mod->mol->ss_circ_flg = 1;
            mod->mol->residue->ss_circ_flg = 1;
        }

        // if a termination signal, set prevTerm to true; thats it
        else if((strncmp(textbuf, "TER", 3) == 0) ||
            (strncmp(textbuf, "END", 3) == 0) ||
            (strncmp(textbuf, "REMARK ORIG", 11) == 0) ||
            (strcmp(textbuf, "\n") == 0)) {

            if(ter_flg == 0){
                ter_flg = 1;
                strcpy(mod->mol->chainID,
                    mod->mol->residue->atom->chainID);
                /* modify end */

                mod->mol->residue->atom_num = atom_num;
                mod->mol->res_num = res_num + 1;
                mod->mol->residue->atomend = mod->mol->residue->atom;
                mod->mol->residueend = mod->mol->residue;

                ret_code = tgNextMalloc(mod, OBJECT_MOL);

                res_num = 0;
                atom_num = 0;
                mol_num++;
            }
        }

        // if ss_detect is true, and HETA, then get atom_name and resname (12 and 17) and see if it is metal; if so, save to table and coords
        else if(sdata.ss_detect == 2){
            //  read hetero atom information
            if(strncmp(textbuf, "HETA", PRE_WORD_NUM) == 0) {
                strncpy(buf,&textbuf[12],5); buf[5] = '\0';
                sscanf(buf, "%s", het_atomname);
                strncpy(buf,&textbuf[17],4); buf[4] = '\0';
                sscanf(buf, "%s", het_resname);

                if(tgCheckAtomNameOfMetals(het_atomname, het_resname,het_symbol_dmy) == 1){
                    strncpy(metal_atom_name[metal_num], &textbuf[12], 4);
                    metal_atom_name[metal_num][4] = '\0';
                    sscanf(&textbuf[30], "%lf%lf%lf",
                        &metal_coord[metal_num][0],
                        &metal_coord[metal_num][1],
                        &metal_coord[metal_num][2]);

                    if(METALS_MAX > metal_num){
                        metal_num++;
                    }
                }
            }

        }

        if(fgets(textbuf, BUFSIZE, fp) == 0) {
            break;
        }
    }
    if(ter_flg == 0){
        mod->mol_num = mol_num + 1;
        mod->mol->res_num = res_num + 1;
    }
    else{
        mod->mol_num = mol_num;
        mod->mol->res_num = res_num;
    }
    mod->mol->residue->atom_num = atom_num;
    mod->mol->residue->atomend = mod->mol->residue->atom;
    mod->mol->residueend = mod->mol->residue;
    mod->end = mod->mol;

    fclose(fp);
















    // if autodetect
    if(sdata.ss_detect == 2){

        // Detect SS bond
        ssbond_thr2 = SSBOND_THRESHOLD * SSBOND_THRESHOLD;
        s_metal_thr2 = S_METAL_DIST_THRESHOLD * S_METAL_DIST_THRESHOLD;

        // Create Cys SG list
        all_cys_num = 0;

        mod->mol = mod->begin;
        for(i=0; i < mod->mol_num; i++){
            if(i != 0) mod->mol = mod->mol->next;

            mod->mol->residue = mod->mol->residuebegin;

            for(j=0; j < mod->mol->res_num; j++){
                if( j != 0 ){
                    mod->mol->residue = mod->mol->residue->next;
                }

                /* Search atom "SG" */
                mod->mol->residue->atom = mod->mol->residue->atombegin;
                for(k = 0; k < mod->mol->residue->atom_num; k++){
                    if(k != 0){
                        mod->mol->residue->atom = mod->mol->residue->atom->next;
                    }

                    if( (strcmp(mod->mol->residue->atom->residue_name, "CYS") == 0 ||
                         strcmp(mod->mol->residue->atom->residue_name, "CYSN+") == 0 ||
                         strcmp(mod->mol->residue->atom->residue_name, "CYSS") == 0 ||
                         strcmp(mod->mol->residue->atom->residue_name, "CYSSN+") == 0 ||
                         strcmp(mod->mol->residue->atom->residue_name, "CYSC-") == 0 ||
                         strcmp(mod->mol->residue->atom->residue_name, "CYSSC-") == 0) &&
                        strcmp(mod->mol->residue->atom->atom_name, "SG") == 0){
                        cys_SeqNum[all_cys_num] = mod->mol->residue->atom->residue_num_org;

                        /* modify 2013.11.27 *
                        strcpy(cys_ChainID[all_cys_num], mod->mol->residue->atom->chainID);
                        strcpy(cys_iCode[all_cys_num], mod->mol->residue->atom->iCode);
                        */
                        cys_ChainID[all_cys_num][0] = mod->mol->residue->atom->chainID[0];
                        cys_ChainID[all_cys_num][1] = '\0';
                        cys_iCode[all_cys_num][0] = mod->mol->residue->atom->iCode[0];
                        cys_iCode[all_cys_num][1] = '\0';
                        /* modify end */

                        cys_SgCoord[all_cys_num][0] = mod->mol->residue->atom->coord->x;
                        cys_SgCoord[all_cys_num][1] = mod->mol->residue->atom->coord->y;
                        cys_SgCoord[all_cys_num][2] = mod->mol->residue->atom->coord->z;

                        all_cys_num++;
                        break;
                    }
                }
            }
        }

        /* Compare Cys SG-Metal distance */
        for(i = 0; i < all_cys_num; i++){
            for(j = 0; j < metal_num; j++ ){
                if( tgCalcSquareDistance(cys_SgCoord[i], metal_coord[j])
                    < s_metal_thr2){

                    sprintf(buf, "%s %8.3lf %8.3lf %8.3lf",
                        metal_atom_name[j],
                        metal_coord[j][0],
                        metal_coord[j][1],
                        metal_coord[j][2]
                    );
                    cys_valid[i] = 0; /* invalid */

                    if(print_metal_header == 0){

                        printf("\n");
                        printf(" INFORMATION> tgCheckSS\n");
                        printf("             The distance of CYS and metal atom is checked\n");

                        print_metal_header = 1;
                    }

                    printf("\n The following atom was found near the SG(CYS) atom.\n");
                    printf(" %s\n", buf);
                    printf(" SG(CYS) Information\n");
                    printf("      chain ID       : %s\n", cys_ChainID[i]);
                    printf("      residue number : %d\n", cys_SeqNum[i]);
                    /* modified by Sumita 2012/08/28 */
                }
            }
        }

        /* Compare Cys SG-SG distance */
        for(i = 0; i < all_cys_num - 1; i++){
            for(j = i + 1; j < all_cys_num; j++){
                if( cys_valid[i] == 1 &&
                    tgCalcSquareDistance(cys_SgCoord[i], cys_SgCoord[j]) < ssbond_thr2){

                    for(k=0; k < ssbond_num; k++){

                        /* modified by Sumita 2012/08/28 */
                        if( (cys_SeqNum[i] == mod->ssbond->seqNum1_org[k] &&
                             strcmp(cys_ChainID[i], mod->ssbond->chainID1_mod[k]) == 0 &&
                             strcmp(cys_iCode[i], mod->ssbond->iCode1[k]) == 0 ) ||
                            (cys_SeqNum[j] == mod->ssbond->seqNum2_org[k] &&
                             strcmp(cys_ChainID[j], mod->ssbond->chainID2_mod[k]) == 0 &&
                             strcmp(cys_iCode[j], mod->ssbond->iCode2[k]) == 0 ) ||
                            (cys_SeqNum[i] == mod->ssbond->seqNum2_org[k] &&
                             strcmp(cys_ChainID[i], mod->ssbond->chainID2_mod[k]) == 0 &&
                             strcmp(cys_iCode[i], mod->ssbond->iCode2[k]) == 0 ) ||
                            (cys_SeqNum[j] == mod->ssbond->seqNum1_org[k] &&
                             strcmp(cys_ChainID[j], mod->ssbond->chainID1_mod[k]) == 0 &&
                             strcmp(cys_iCode[j], mod->ssbond->iCode1[k]) == 0 ) ){
                        /* modified by Sumita 2012/08/28 */

                            /* Same SSBOND already exists. */
                            break;

                        }

                    }

                    /* Same SSBOND does not exist.  */
                    if(k == ssbond_num){
                        strcpy(mod->ssbond->iCode1[ssbond_num], cys_iCode[i]);
                        strcpy(mod->ssbond->iCode2[ssbond_num], cys_iCode[j]);
                        strcpy(mod->ssbond->chainID1_mod[ssbond_num], cys_ChainID[i]);
                        strcpy(mod->ssbond->chainID2_mod[ssbond_num], cys_ChainID[j]);
                        mod->ssbond->seqNum1_org[ssbond_num] = cys_SeqNum[i];
                        mod->ssbond->seqNum2_org[ssbond_num] = cys_SeqNum[j];
                        ssbond_num++;
                    }

                }
            }
        }
    }

    if(ssbond_num != 0){
        chainID_space_num = 0;
        mod->mol = mod->begin;
        for(i=0; i < mod->mol_num; i++){
            if(i != 0) mod->mol = mod->mol->next;

            mod->mol->residue = mod->mol->residuebegin;
            mod->mol->residue->atom = mod->mol->residue->atombegin;

            if(mod->mol->residue->atom->chainID[0] == ' '){
                chainID_space_num++;
            }
        }
        if(chainID_space_num > 1){
            puts(" ERROR> tgCheckSS");
            puts("        Chain ID is not defined although there are two or");
            puts("        more chains and SSBOND is defined.");
            exit(1);
        }

        /* added by Sumita 2012/08/28 */
        if(mod->ssbond != NULL){
            /* check same ssbond. */
            for(i=0; i < ssbond_num - 1; i++){
                for(j= i + 1; j < ssbond_num; j++){

                    if(tgContainsSameSSAtom(mod->ssbond, i, j, ssbond_num) == 0){
                        puts(" ERROR> tgCheckSS");
                        puts("        SSBOND contains two or more same residues.");
                        exit(1);
                    }
                }
            }
        }
        /* added by Sumita 2012/08/28 */

    }



    /******************************
     * modify residue number
     ******************************/
    chainID = (char **)Malloc(sizeof(char *) * mod->mol_num);
    chainID_modify = (char **)Malloc(sizeof(char *) * mod->mol_num);
    for(i = 0; i < mod->mol_num; i++) {
        chainID[i] = (char *)Malloc(sizeof(char) * 2);
        chainID_modify[i] = (char *)Malloc(sizeof(char) * 2);
    }

    mod->mol = mod->begin;
    for(i = 0; i < mod->mol_num; i++) {
        if(i != 0) mod->mol = mod->mol->next;

        strcpy(chainID[i], mod->mol->chainID);
        strcpy(chainID_modify[i], chainID[i]);
    }

    /******************************
     * modify SSBOND information
     ******************************/
    for(j = 0; j < ssbond_num; j++) {
        strcpy(mod->ssbond->chainID1_org[j], mod->ssbond->chainID1_mod[j]);
        strcpy(mod->ssbond->chainID2_org[j], mod->ssbond->chainID2_mod[j]);
    }

    getChainGroup(mod, ssbond_num);

    // loop for chainGroup
    for(i=0; i < mod->mol_num; i++){
        // loop for chainID
        for(k=0; k < mod->mol_num; k++){
            // loop for ssbond
            for(j=0; j < ssbond_num; j++){
                if(chainID[k][0] == mod->ssbond->chainID1_mod[j][0] ||
                   chainID[k][0] == mod->ssbond->chainID2_mod[j][0]){

                    for(l=0; l < strlen(mod->chainGroup[i]); l++){
                        if(mod->chainGroup[i][l] == chainID[k][0] &&
                           mod->chainGroup[i][l] == chainID_modify[k][0]){
                            for(m=0; m < strlen(mod->chainGroup[i]); m++){
                                temp1 = chainID_modify[k][0];
                                temp2 = mod->chainGroup[i][m];

                                if(temp1 > temp2){
                                    chainID_modify[k][0] =
                                        mod->chainGroup[i][m];
                                    chainID_modify[k][1] = '\0';
                                }
                            }
                        }
                    }
                    if(chainID[k][0] == mod->ssbond->chainID1_mod[j][0]){
                        strcpy(mod->ssbond->chainID1_mod[j],
                            chainID_modify[k]);
                    }
                    else{
                        strcpy(mod->ssbond->chainID2_mod[j],
                            chainID_modify[k]);
                    }
                }
            }
        }
    }

    /******************************
     * correct SSBOND-residue_num
     ******************************/
    before_res_num = (int *)Malloc(sizeof(int) * mod->mol_num);

    /* modify 2011.10.12, start */
    for(i=0; i < mod->mol_num; i++){
        if(i == 0) before_res_num[i] = 0;
        else{

            for(j=0; j < i; j++){
                if(chainID_modify[i][0] == chainID_modify[j][0]) {
                    mod->mol = mod->begin;
                    for(k=0; k < j; k++){
                        mod->mol = mod->mol->next;
                    }
                    before_res_num[i] += mod->mol->res_num;
                }
            }
        }
    }
    /* modify end */


    /******************************/
    val = -999;
    for(k = 0; k < ssbond_num; k++) {
        mod->mol = mod->begin;
        for(m = 0; m < mod->mol_num; m++) {
            if(m != 0) mod->mol = mod->mol->next;
            if(strcmp(mod->ssbond->chainID1_org[k], chainID[m]) == 0){
                val = m;
                break;
            }
        }
        if(val < 0) continue;
        if(strcmp(mod->ssbond->chainID1_org[k], chainID[val]) == 0) {
            mod->mol->residue = mod->mol->residuebegin;
            for(j = 0; j < mod->mol->res_num; j++) {
                if(j != 0) mod->mol->residue =
                    mod->mol->residue->next;

                mod->mol->residue->atom =
                    mod->mol->residue->atombegin;

                if((mod->mol->residue->atom->residue_num_org ==
                        mod->ssbond->seqNum1_org[k])
                    && (mod->mol->residue->atom->iCode[0] ==
                            mod->ssbond->iCode1[k][0])) {

                    mod->ssbond->seqNum1_mod[k] =
                        mod->mol->residue->atom->residue_num_mod;

                    break;

                }
            }
        }
    }
    val = -999;
    for(k = 0; k < ssbond_num; k++) {
        mod->mol = mod->begin;
        for(m = 0; m < mod->mol_num; m++) {
            if(m != 0) mod->mol = mod->mol->next;

            if(strcmp(mod->ssbond->chainID2_org[k], chainID[m]) == 0){
                val = m;
                break;
            }
        }

        if(val < 0) continue;
        if(strcmp(mod->ssbond->chainID2_org[k], chainID[val]) == 0) {
            mod->mol->residue = mod->mol->residuebegin;
            for(j = 0; j < mod->mol->res_num; j++) {
                if(j != 0) mod->mol->residue =
                    mod->mol->residue->next;

                mod->mol->residue->atom =
                    mod->mol->residue->atombegin;
                if(mod->ssbond->iCode2[k][0] == '\n'){
                    mod->ssbond->iCode2[k][0] = ' ';
                    mod->ssbond->iCode2[k][1] = '\0';
                }

                if((mod->mol->residue->atom->residue_num_org ==
                        mod->ssbond->seqNum2_org[k])
                    && (mod->mol->residue->atom->iCode[0] ==
                            mod->ssbond->iCode2[k][0])) {

                    mod->ssbond->seqNum2_mod[k] =
                        mod->mol->residue->atom->residue_num_mod;

                    break;
                }
            }
        }
    }
    /******************************
     * Get Modified chainID
     ******************************/
    mod->mol = mod->begin;
    for(k = 0; k < mod->mol_num; k++) {
        if(k != 0) mod->mol = mod->mol->next;
        for(l = 0; l < mod->mol_num; l++) {
            /* modify 2013.11.27 *
            if(strcmp(chainID[l], mod->mol->chainID) == 0)
                strcpy(mod->mol->chainID_mod, chainID_modify[l]);
            */
            if(chainID[l][0] == mod->mol->chainID[0]) {
                mod->mol->chainID_mod[0] = chainID_modify[l][0];
                mod->mol->chainID_mod[1] = '\0';
            /* modify end */
                break;
            }
        }
    }

    /******************************
     * File Open
     ******************************/
    if((fpw = fopen("changed_file", "w")) == NULL) {
        tgError(FILE_OPEN_ERROR, buf, NULL, NULL, DMY_INT, DMY_INT);
    }

    /******************************
     * Output
     ******************************/
    check = (int *)Malloc(sizeof(int) * mod->mol_num);

    seqNum1 = seqNum2 = 0;
    for(i = 0; i < ssbond_num; i++) {

        if(mod->ssbond->seqNum1_mod[i] != 0
            && mod->ssbond->seqNum2_mod[i] != 0) {

            mod->mol = mod->begin;
            for(j = 0; j < mod->mol_num; j++) {
                if(strcmp(mod->mol->chainID,
                        mod->ssbond->chainID1_org[i]) == 0) {
                    seqNum1 = before_res_num[j];
                }
                if(strcmp(mod->mol->chainID,
                        mod->ssbond->chainID2_org[i]) == 0) {
                    seqNum2 = before_res_num[j];
                }

                if(j == mod->mol_num - 1)
                    break;
                mod->mol = mod->mol->next;
            }
            fprintf(fpw, "SSBOND %3d CYS %1s %4d%1s   CYS %1s %4d%1s\n",
                i + 1,
                mod->ssbond->chainID1_org[i],
                mod->ssbond->seqNum1_org[i], mod->ssbond->iCode1[i],
                mod->ssbond->chainID2_org[i],
                mod->ssbond->seqNum2_org[i], mod->ssbond->iCode2[i]
            );

            /* add 2011.5.26, start */
            sNum1_mod[i] = mod->ssbond->seqNum1_mod[i] + seqNum1;
            sNum2_mod[i] = mod->ssbond->seqNum2_mod[i] + seqNum2;
            /* add end */
        }
    }

    mod->mol = mod->begin;
    end = 0;
    flg = 0;
    total_atom_num = 0;

    strcpy(checkID, mod->mol->chainID_mod);
    while(1){
        for(k = 0; k < mod->mol_num; k++) {
            if(k != 0) mod->mol = mod->mol->next;

            if(check[k] != 0) continue;

            if(strcmp(checkID, mod->mol->chainID_mod) == 0) {

                if(k != 0 && flg == 1 && ssbond_num != 0) {
                    fprintf(fpw, "INTER-SSBOND\n");
                }
                if(k != 0 && flg == 1 && ssbond_num == 0
                    && strcmp(checkID, " ") == 0) {
                    fprintf(fpw, "TER\n");
                }

                mod->mol->residue = mod->mol->residuebegin;
                if(mod->mol->ss_circ_flg == 1 || mod->mol->residue->ss_circ_flg == 1) {
                    fprintf(fpw, "CIRC\n");
                }

                memset(dmychar2, '\0', sizeof(dmychar2));

                for(j = 0; j < mod->mol->res_num; j++) {
                    if(j != 0) mod->mol->residue =
                        mod->mol->residue->next;

                    mod->mol->residue->atom =
                        mod->mol->residue->atombegin;
                    for(i = 0; i < mod->mol->residue->atom_num; i++) {
                        if(i != 0) mod->mol->residue->atom =
                            mod->mol->residue->atom->next;

                        total_atom_num++;

                        if(strlen(mod->mol->residue->atom->atom_name) > 3) {
                          strcpy(buf,
                          "ATOM  %5d %4s%1s%-4s%1s%4d%1s   %8.3lf%8.3lf%8.3lf\n");
                        }
                        else {
                          strcpy(buf,
                          "ATOM  %5d  %-3s%1s%-4s%1s%4d%1s   %8.3lf%8.3lf%8.3lf\n");
                        }

                        strcpy(mod->mol->residue->atom->chainID_mod,
                            mod->mol->chainID_mod);

                        fprintf(fpw, buf,
                            total_atom_num,
                            mod->mol->residue->atom->atom_name,
                            mod->mol->residue->atom->altLoc,
                            mod->mol->residue->atom->residue_name,
                            mod->mol->residue->atom->chainID,
                            mod->mol->residue->atom->residue_num_org,
                            mod->mol->residue->atom->iCode,
                            mod->mol->residue->atom->coord->x,
                            mod->mol->residue->atom->coord->y,
                            mod->mol->residue->atom->coord->z);

                        /* add 20080826 start */
                        strcpy(dmychar2,
                                mod->mol->residue->atom->residue_name);
                        /* add 20080826 end */
                    }
                }
                strcpy(checkID, mod->mol->residue->atom->chainID_mod);
                check[k] = 1;
                flg = 1;
            }

            if(k == mod->mol_num - 1) {
                flg = 0;
                for(l = 0; l < mod->mol_num; l++) {
                    end = 1;
                    if(check[l] == 0) {
                        strcpy(checkID, chainID_modify[l]);
                        end = 0;

                        break;
                    }
                }
                fprintf(fpw, "TER\n");
            }
            if(end == 1) break;
        }

        end = 1;
        for(l=0; l < mod->mol_num; l++){
            if(check[l] == 0){
                end = 0;
                break;
            }
        }

        if(end == 1) break;
        mod->mol = mod->begin;
    }

    fprintf(fpw, "TER\n");
    fclose(fpw);

    strcpy(sdata.coordfile, "changed_file");

    /* add 2011.5.26, start */
    for(i=0; i < ssbond_num; i++){
        strcpy(cID1[i], mod->ssbond->chainID1_mod[i]);
        strcpy(cID2[i], mod->ssbond->chainID2_mod[i]);
        strcpy(cID1_org[i], mod->ssbond->chainID1_org[i]);
        strcpy(cID2_org[i], mod->ssbond->chainID2_org[i]);
        strcpy(icod1[i], mod->ssbond->iCode1[i]);
        strcpy(icod2[i], mod->ssbond->iCode2[i]);
        sNum1[i] = mod->ssbond->seqNum1_org[i];
        sNum2[i] = mod->ssbond->seqNum2_org[i];
    }
    /* add end */

    /* add 2011.11.15, start */
    mod->mol = mod->begin;
    for(i=0; i < mod->mol_num; i++){
	    if(i != 0) mod->mol = mod->mol->next;
        strcpy(molCID[i], mod->mol->chainID);
        strcpy(molCID_mod[i], mod->mol->chainID_mod);
    }
    /* add end */

    for(i = 0; i < mod->mol_num; i++) {
        Free(chainID[i]);
        Free(chainID_modify[i]);
    }
    Free(chainID);
    Free(chainID_modify);
    Free(before_res_num);
    Free(check);

    /* add 2013.11.14 */
    tgFreeAmberStruct(mod);

    return ssbond_num;
}

/*******************************************************************/
/* get ChainGroup() */
/*******************************************************************/
int getChainGroup(TGASystemPtr mod, int ssbond_num) {
    int i, j, k, ii, jj;
    char group[MOLNUMMAX][10];
    short flg, flg2;
    int len, len2;

    /* get ChainID */
    mod->mol = mod->begin;
    for(i=0; i < mod->mol_num; i++){
        if(i != 0) mod->mol = mod->mol->next;
        group[i][0] = mod->mol->chainID[0];
        group[i][1] = '\0';
    }


    for(ii=0; ii < mod->mol_num; ii++){
        flg2 = 0;
        for(i=0; i < mod->mol_num; i++){

            len = strlen(group[i]);
            len2 = len;
            for(jj=0; jj < len; jj++){
                for(j=0; j < ssbond_num; j++){
                    if(mod->ssbond->chainID1_mod[j][0] ==
                       mod->ssbond->chainID2_mod[j][0]) continue;

                    if(group[i][jj] ==
                       mod->ssbond->chainID1_mod[j][0]){

                        flg = 0;
                        for(k=0; k < len2; k++){
                            if(group[i][k] ==
                                mod->ssbond->chainID2_mod[j][0]) flg = 1;
                        }
                        if(flg == 0){
                            group[i][len2] =
                                mod->ssbond->chainID2_mod[j][0];
                            group[i][len2+1] = '\0';
                            len2++;
                            flg2 = 1;
                       }
                    }
                    if(group[i][jj] ==
                       mod->ssbond->chainID2_mod[j][0]){
                        flg = 0;
                        for(k=0; k < len2; k++){
                            if(group[i][k] ==
                                mod->ssbond->chainID1_mod[j][0]) flg = 1;
                        }
                        if(flg == 0){
                            group[i][len2] =
                                mod->ssbond->chainID1_mod[j][0];
                            group[i][len2+1] = '\0';
                            len2++;
                            flg2 = 1;
                        }
                    }
                }
            }

        }
        if(flg2 == 0) break;
    }

    for(i=0; i < mod->mol_num; i++){
        strcpy(mod->chainGroup[i], group[i]);
    }

    return 0;
}

/*******************************************************************/
/* tgContainsSameSSAtom()                                          */
/*******************************************************************/
int tgContainsSameSSAtom(TGASSBondPtr ssbond, int i, int j, int ssbond_num) {
    int retCode = 1;

    if(ssbond == NULL || ssbond_num > MAXSSBOND){
        retCode = -1;
    }

    if( strcmp(ssbond->chainID1_mod[i], ssbond->chainID1_mod[j]) == 0 &&
        ssbond->seqNum1_org[i] == ssbond->seqNum1_org[j] &&
        strcmp(ssbond->iCode1[i], ssbond->iCode1[j]) == 0 ){
        retCode = 0;

    }
    else if(strcmp(ssbond->chainID1_mod[i], ssbond->chainID2_mod[j]) == 0 &&
            ssbond->seqNum1_org[i] == ssbond->seqNum2_org[j] &&
            strcmp(ssbond->iCode1[i], ssbond->iCode2[j]) == 0 ){
        retCode = 0;

    }
    else if(strcmp(ssbond->chainID2_mod[i], ssbond->chainID1_mod[j]) == 0 &&
            ssbond->seqNum2_org[i] == ssbond->seqNum1_org[j] &&
            strcmp(ssbond->iCode2[i], ssbond->iCode1[j]) == 0 ){
        retCode = 0;

    }
    else if(strcmp(ssbond->chainID2_mod[i], ssbond->chainID2_mod[j]) == 0 &&
            ssbond->seqNum2_org[i] == ssbond->seqNum2_org[j] &&
            strcmp(ssbond->iCode2[i], ssbond->iCode2[j]) == 0 ){
        retCode = 0;

    }
    else{
        retCode = 1;

    }

    return retCode;

}

