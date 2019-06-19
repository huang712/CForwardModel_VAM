//
// Created by Feixiong Huang on 1/21/18.
//

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
    double val,lat,lon;
    char temp[1000];
    strcpy(temp,saveDir);
    //complex double valc;

    FILE *outp,*outp1,*outp2;
    switch(pathType){
        case 0:
            outp = fopen("Jacobian.dat", "wb");
            outp1 = fopen("Jacobian_lat.dat", "wb");
            outp2 = fopen("Jacobian_lon.dat", "wb");
            break;
        case 1:
            printf('save Jacobian not available\n');
            break;
        case 2:
            strcat(temp, "Jacobian.dat");
            outp = fopen(temp, "wb");
            //outp1 = fopen("/users/fax/CYGNSS/VAM/MATLAB/Jacobian_lat.dat", "wb");
            //outp2 = fopen("/users/fax/CYGNSS/VAM/MATLAB/Jacobian_lon.dat", "wb");
            break;
    }

    double numDDMbins,numPts_LL;
    numDDMbins = (double)jacob.numDDMbins;
    numPts_LL = (double)jacob.numPts_LL;

    fwrite(&numDDMbins, 1, sizeof(double), outp);
    //fwrite(&numDDMbins, 1, sizeof(double), outp1);
    //fwrite(&numDDMbins, 1, sizeof(double), outp2);
    fwrite(&numPts_LL, 1, sizeof(double), outp);
    //fwrite(&numPts_LL, 1, sizeof(double), outp1);
    //fwrite(&numPts_LL, 1, sizeof(double), outp2);

    for (int i = 0; i < jacob.numDDMbins * jacob.numPts_LL; i++) {
        val = jacob.data[i].value;
        lat = jacob.data[i].lat_deg;
        lon = jacob.data[i].lon_deg;
        fwrite(&val, 1, sizeof(double), outp);
        //fwrite(&lat, 1, sizeof(double), outp1);
        //fwrite(&lon, 1, sizeof(double), outp2);
    }
    fclose(outp);
    //fclose(outp1);fclose(outp2);
    printf("save Jacobian into file\n");
}


void indexLL_saveToFile(struct Jacobian jacob, char saveDir[1000]) {
    //FILE *outp = fopen("PtsVec.dat", "wb");
    char temp[1000];
    strcpy(temp,saveDir);
    strcat(temp, "indexLL.dat");
    FILE *outp = fopen(temp, "wb");
    fwrite(jacob.Pts_ind_vec, sizeof(int), jacob.numPts_LL, outp);
//    double temp;
//    for (int j = 0; j < jacob.numPts_LL; j++) {
//        temp = (double)jacob.Pts_ind_vec[j];
//        //fwrite(&temp, sizeof(double), 1, outp);
//        //fwrite(&jacob.Pts_lat_vec[j], sizeof(double), 1, outp);
//        //fwrite(&jacob.Pts_lon_vec[j], sizeof(double), 1, outp);
//    }
    fclose(outp);
    printf("save index of points in LL into file\n");
}
