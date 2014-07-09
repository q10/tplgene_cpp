int tgModify4LackedRes1(char *in_file)
{
    char out_file[256];
    char textbuf[BUFSIZE+1], buf[10], buf1[10], buf2[10], buf3[10];
    FILE *fp, *fpw;
    char work_atom_name[10];
    char atom_name_C[ATOM_NAME_MAX + 1];
    char chainID_C[TWO];
    char res_ins_C[10];
    double cod_C[3], cod_N[3];
    char chainID_N[TWO];
    char res_ins_N[10];
    double dist;
    /* add 2013.10.29, start */
    int interSS_flg = 0;
    /* add end */

    memset(out_file, '\0', sizeof(char) * 256);
    strcpy(out_file, "changed_file2");

    if((fp = fopen(in_file, "r")) == NULL) {
        puts(" ERROR> tgModifyLackedRes");
        puts("        Can not open file");
        exit(1);
    }

    if((fpw = fopen(out_file, "w")) == NULL){
        puts(" ERROR> tgModifyLackedRes");
        puts("        Can not open file");
        exit(1);
    }

    memset(atom_name_C, '\0', sizeof(char) * (ATOM_NAME_MAX + 1));
    while(!feof(fp)){
		if(fgets(textbuf, BUFSIZE, fp) == 0) break;

        if (line.compare(0, 4, "ATOM") == 0 or line.compare(0, 6, "HETATM") == 0) {

            // get atom_name (12, 4), remote whitespace

            // Calculate distance between C and N atoms, or P and C3' atoms
            if((strcmp(work_atom_name, "N") == 0 && sdata.chnspecies == 1) or (strcmp(work_atom_name, "P") == 0 && sdata.chnspecies == 2)){

               if (atom_name_C[0] != '\0') {
                    memset(chainID_N, '\0', sizeof(char) * 2);
                    memset(res_ins_N, '\0', sizeof(char) * 10);
                    strncpy(chainID_N, &textbuf[21], 1);
                    strncpy(res_ins_N, &textbuf[22], 5);

                    if(strcmp(res_ins_N, res_ins_C) != 0){
                        if(strcmp(chainID_N, chainID_C) == 0){
                            strncpy(buf1, &textbuf[30], 8);
                            strncpy(buf2, &textbuf[38], 8);
                            strncpy(buf3, &textbuf[46], 8);

                            sscanf(buf1, "%lf", &cod_N[0]);
                            sscanf(buf2, "%lf", &cod_N[1]);
                            sscanf(buf3, "%lf", &cod_N[2]);

                            dist = tgCalcSquareDistance(cod_N, cod_C);
                            if(sdata.chnspecies == PEPTIDE){
                                if(sqrt(dist) > 2.0) {
                                    // 2013.11.28
                                    //puts(" WARNING> tgModify4LackedRes");
                                    //puts("          It is too away between N and C atoms in same chain..");
                                    //printf("          %sth and %sth residue!\n", res_ins_C, res_ins_N);
                                    //

                                    //puts("          Set TER line.");

                                    //fprintf(fpw, "INTER-SSBOND\n");
                                    if(interSS_flg == 0){
                                        fprintf(fpw, "TER\n");
                                    }

                                    memset(atom_name_C, '\0', sizeof(char) * (ATOM_NAME_MAX + 1));
                                    chainID_N[0] = '\0'; chainID_N[1] = '\0';
                                    chainID_C[0] = '\0'; chainID_C[1] = '\0';
                                    memset(res_ins_N, '\0', sizeof(char) * 10);
                                    memset(res_ins_C, '\0', sizeof(char) * 10);

                                }
                            } else if(sdata.chnspecies == NUCLEOTIDE){
                                if(sqrt(dist) > 3.0){
                                    if(interSS_flg == 0) fprintf(fpw, "TER\n");
                                    memset(atom_name_C, '\0', sizeof(char) * (ATOM_NAME_MAX + 1));
                                    chainID_N[0] = '\0'; chainID_N[1] = '\0';
                                    chainID_C[0] = '\0'; chainID_C[1] = '\0';
                                    memset(res_ins_N, '\0', sizeof(char) * 10);
                                    memset(res_ins_C, '\0', sizeof(char) * 10);
                                }
                            }
                        }
                    }
                }
            }
            /* store atom name, chainID residue number */
            else if((strcmp(work_atom_name, "C") == 0 && sdata.chnspecies == 1)||
               (strcmp(work_atom_name, "C3'") == 0 && sdata.chnspecies == 2)){
                memset(chainID_C, '\0', sizeof(char) * 2);
                memset(res_ins_C, '\0', sizeof(char) * 10);
                strncpy(chainID_C, &textbuf[21], 1);
                strncpy(res_ins_C, &textbuf[22], 5);

                memset(atom_name_C, '\0', sizeof(char) * (ATOM_NAME_MAX + 1));
                strcpy(atom_name_C, work_atom_name);

                memset(buf1, '\0', sizeof(char) * 10);
                memset(buf2, '\0', sizeof(char) * 10);
                memset(buf3, '\0', sizeof(char) * 10);

                strncpy(buf1, &textbuf[30], 8);
                strncpy(buf2, &textbuf[38], 8);
                strncpy(buf3, &textbuf[46], 8);

                sscanf(buf1, "%lf", &cod_C[0]);
                sscanf(buf2, "%lf", &cod_C[1]);
                sscanf(buf3, "%lf", &cod_C[2]);
            }
            fprintf(fpw, "%s", textbuf);

            interSS_flg = 0;
        }

        else if(strncmp(textbuf, "TER", 3) == 0){
            memset(atom_name_C, '\0', sizeof(char) * (ATOM_NAME_MAX + 1));
            chainID_N[0] = '\0'; chainID_N[1] = '\0';
            chainID_C[0] = '\0'; chainID_C[1] = '\0';
            memset(res_ins_N, '\0', sizeof(char) * 10);
            memset(res_ins_C, '\0', sizeof(char) * 10);

            cod_C[0] = cod_C[1] = cod_C[2] = 0.0;

            fprintf(fpw, "%s", textbuf);

            interSS_flg = 0;
        }
        /* add 2013.10.29, start */
        else if(strncmp(textbuf, "INTER-SSBOND", 12) == 0){
            fprintf(fpw, "%s", textbuf);
            interSS_flg = 1;
        }
        /* add end */
        else{
            fprintf(fpw, "%s", textbuf);
            interSS_flg = 0;
        }
    }

    memset(sdata.coordfile, '\0', sizeof(char) * strlen(sdata.coordfile));
    strcpy(sdata.coordfile, out_file);

    fclose(fpw);
    fclose(fp);

    return 0;
}
