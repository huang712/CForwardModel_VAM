//---------------------------------------------------------------------------
//
// Functions to save data into files
// Created by Feixiong Huang on 1/21/18
//
//***************************************************************************

#include "forwardmodel.h"
#include "complex.h"
#include <string.h>

void DDMfm_saveToFile(struct DDMfm ddm_fm, int index, int pathType, char saveDir[1000]) {

    double val;
    char filename[50];
    char temp[1000];
    strcpy(temp,saveDir);

    FILE *outp;
    switch(pathType){
        case 0:
            outp = fopen("DDMfm.dat","wb");
            break;
        case 1:
            sprintf(filename, "DDMfm/DDMfm%d.dat", index);
            outp = fopen(filename, "wb");
            break;
        case 2:
            strcat(temp, "DDMfm.dat");
            outp = fopen(temp,"wb");
            break;
    }

    for (int i = 0; i < ddm_fm.numDelaybins * ddm_fm.numDopplerbins; i++) {
        val = ddm_fm.data[i].power;
        fwrite(&val, 1, sizeof(double), outp);
    }
    fclose(outp);
    printf("save FM DDM into file\n");
}

void Jacobian_saveToFile(struct Jacobian jacob, int index, int pathType, char saveDir[1000]){
    // pathType:
    // 0 for save in current directory;
    // 1 for save in a folder; (not supported)
    // 2 for in VAM folder

    double val;
    char temp[1000];
    strcpy(temp,saveDir);

    FILE *outp;
    switch(pathType){
        case 0:
            outp = fopen("Jacobian.dat", "wb");
            break;
        case 1:
            printf("save Jacobian not available\n");
            break;
        case 2:
            strcat(temp, "Jacobian.dat");
            outp = fopen(temp, "wb");
            break;
    }

    double numDDMbins,numPts_LL;
    numDDMbins = (double)jacob.numDDMbins;
    numPts_LL = (double)jacob.numPts_LL;

    fwrite(&numDDMbins, 1, sizeof(double), outp);
    fwrite(&numPts_LL, 1, sizeof(double), outp);

    for (int i = 0; i < jacob.numDDMbins * jacob.numPts_LL; i++) {
        val = jacob.data[i].value;
        fwrite(&val, 1, sizeof(double), outp);
    }
    fclose(outp);
    printf("save Jacobian into file\n");
}

void indexLL_saveToFile(struct Jacobian jacob, char saveDir[1000]) {
    // indexLL: indices of the grid points in the input wind field

    char temp[1000];
    strcpy(temp,saveDir);
    strcat(temp, "indexLL.dat");
    FILE *outp = fopen(temp, "wb");
    fwrite(jacob.Pts_ind_vec, sizeof(int), jacob.numPts_LL, outp);
    fclose(outp);
    printf("save index of points in LL into file\n");
}
