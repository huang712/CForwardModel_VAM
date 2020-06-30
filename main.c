// Forward model for DDM assimilation
// Read a config file as input
// Read CYGNSS data and wind field data file to produce a simulated DDM and a Jacobian matrix
// Version 3.01  2019/09/22 by Feixiong Huang

#include <math.h>
#include <time.h>
#include <string.h>
#include "forwardmodel.h"
#include "cygnss.h"

int main(int argc, char *argv[]) {
    double start, end;
    start = clock();

    //read file path, index from a config file (txt)
    char windFilename[1000],L1dataFilename[1000],saveDir[1000], str[100], a[10];
    FILE *fp;
    int len;
    int sampleIndex; //zero based
    int ddmIndex; //ddm channel 0-3

    struct windInfo info;
    struct option opt;
    fp = fopen(argv[1], "r"); // argv[1] is the config file name
    if (fp == NULL){
        printf("Could not open file %s",argv[1]);
        return 1;
    }

    fgets(windFilename, 1000, fp); //read first line (wind file name)  there is a '\n' at the end of the string - need to remove
    fgets(L1dataFilename, 1000, fp); //read second line (CYGNSS file name)
    fgets(saveDir, 1000, fp); //Directory to write DDMfm, Jaocbian and indexLL

    len = strlen(windFilename); windFilename[len-1] = '\0';  //add '\0' at the end of filename which means the end of the string
    len = strlen(L1dataFilename); L1dataFilename[len-1] = '\0';
    len = strlen(saveDir); saveDir[len-1] = '\0';

    char a4[2];
    fgets(str, 100, fp); //read 4th line (DDM index) zeros based 0 - 3
    strncpy(a4, str+14, 1);
    ddmIndex = atoi(a4);

    char a5[10];
    fgets(str, 100, fp); //read 5th line (CYGNSS index)
    strncpy(a5, str+14, 8);  //it will add '\0' at the end of the string
    sampleIndex = atoi(a5);

    char a6[10];
    fgets(str, 100, fp); //read 6th line
    strncpy(a6, str+14, 5);
    info.numPtsLon = atoi(a6);

    fgets(str, 100, fp); //read 7th line
    info.numPtsLat = atoi(strncpy(a, str+14, 5));
    fgets(str, 100, fp); //read 8th line
    info.lon_min_deg = atof(strncpy(a, str+14, 10));
    fgets(str, 100, fp); //read 9th line
    info.lat_min_deg = atof(strncpy(a, str+14, 10));
    fgets(str, 100, fp); //read 10th line
    info.resolution = atof(strncpy(a, str+14, 10));
    fgets(str, 100, fp); //read 11th line
    opt.JacobOnOff = atoi(strncpy(a, str+14, 1));
    fgets(str, 100, fp); //read 12th line
    opt.thermalNoiseOnOff = atoi(strncpy(a, str+20, 1));

    info.lon_max_deg = info.lon_min_deg + info.resolution * (info.numPtsLon-1);
    info.lat_max_deg = info.lat_min_deg + info.resolution * (info.numPtsLat-1);

    //Initilization -----------------------------------------------------------------------------------------------------------------------
    int pathType=2; //0 for save in current directory; 1 for save in folders; 2 for in VAM folder

    struct CYGNSSL1 l1data;
    readL1data(L1dataFilename, sampleIndex, ddmIndex, &l1data);  //read L1 data into the structure l1data

    //quality flag should be checked before running the forward model
    if(l1data.quality_flags != 0){
        printf("Quality flags = %d\n",l1data.quality_flags);
        return 0; //skip data of quality issue
    }

    printf("ddmIndex = %d, sampleIndex = %d, quality_flags = %d\n", ddmIndex, sampleIndex, l1data.quality_flags);
    printf("GPS PRN = %d\n", l1data.prn_code);
    printf("sp delay row = %f, sp doppler col = %f\n", l1data.ddm_sp_delay_row,l1data.ddm_sp_dopp_col);
    printf("sp lat = %f, lon = %f\n",l1data.sp_lat,l1data.sp_lon);
    printf("ant = %d\n",l1data.ddm_ant);

    struct metadata meta;
    struct powerParm pp;
    struct inputWindField iwf;
    struct Geometry geom;
    struct DDMfm ddm_fm;
    struct Jacobian jacob;

    printf("\n");
    printf("Initialize input/output structure...\n");

    init_inputWindField_data(windFilename, &iwf, info);
    init_metadata(l1data, &meta);
    init_powerParm(l1data, &pp);  //cost time ~ 0.3s
    init_Geometry(l1data, &geom);
    init_DDM(l1data, &ddm_fm);
    init_Jacobian(&jacob);

    //Run forward model ---------------------------------------------------------------------------------------------------------------
    forwardModel(meta, pp, iwf, geom, &ddm_fm, &jacob,opt);

    printf("ddm obs = %e\n",l1data.DDM_power[6][5]);
    printf("ddm fm = %e\n",ddm_fm.data[92].power);
    printf("H = %e\n",jacob.data[12247].value);

    //Save simulated DDM and Jacobian -------------------------------------------------------------------------------------------------
    //DDMobs_saveToFile(l1data, sampleIndex,pathType, saveDir);
    DDMfm_saveToFile(ddm_fm, sampleIndex, pathType, saveDir);
    if (opt.JacobOnOff == 1) {
        Jacobian_saveToFile(jacob, sampleIndex, pathType, saveDir);
        indexLL_saveToFile(jacob, saveDir);
    }

    free(pp.data);
    free(iwf.data);
    free(ddm_fm.data);
    free(jacob.data);

    printf("END\n");
    printf("\n");
    end =clock();
    printf("Forward model running time: %f seconds\n", (end-start)/CLOCKS_PER_SEC);
    return 0;
}
