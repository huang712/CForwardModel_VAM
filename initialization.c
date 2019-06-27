//
// Created by Feixiong Huang on 10/19/17.
//
#include <netcdf.h>
#include <math.h>
#include <string.h>
#include "forwardmodel.h"

#define NLAT_core 501   //core: resolution = 0.02deg
#define NLON_core 501
#define NLAT_syno 721   //synopric: resolution = 0.125deg
#define NLON_syno 881

#define ERR(e) {printf("Error: %s\n", nc_strerror(e));}

void init_metadata(struct CYGNSSL1 l1data, struct metadata *meta) {
    //read from CYGNSS data
    meta->meas_ddm_sp_index[0] = l1data.ddm_sp_delay_row;
    meta->meas_ddm_sp_index[1] = l1data.ddm_sp_dopp_col;

    //default numbers
    meta->numDelaybins = 400;
    meta->numDopplerbins = 400;
    meta->delayRez_chips = 0.0510345;
    meta->dopplerRes_Hz = 25;
    meta->resample_startBin[0] = 0;
    meta->resample_startBin[1] = 100;
    meta->resample_resolution_bins[0] = 5;
    meta->resample_resolution_bins[1] = 20;
    meta->resample_numBins[0] = 17;
    meta->resample_numBins[1] = 11;
    meta->specular_delayBinIdx = meta->resample_startBin[0] + (int)round(meta->meas_ddm_sp_index[0] * meta->resample_resolution_bins[0]);
    meta->specular_dopplerBinIdx = meta->resample_startBin[1] + (int)round(meta->resample_resolution_bins[1] * meta->meas_ddm_sp_index[1]);

    meta->temp_K = l1data.ant_temperature_cels+273.15;
    meta->noiseFigure_dB = l1data.noise_figue;
    meta->thermalNoiseOnOff = 0;
    meta->excess_noisefloor_dB = 0;

    meta->grid_resolution_m = 1000;
    meta->numGridPoints[0] = 120; // this must be 120 x 120 km for incidence angle = 60 deg.
    meta->numGridPoints[1] = 120;
    meta->surfaceCurvatureType = 1;//1 = spherical, 2 = flat

    meta->prn_code = l1data.prn_code;
    meta->utc_sec = l1data.utc_sec;

}

void init_powerParm(struct CYGNSSL1 l1data, struct powerParm *pp){
    pp->Rx_antennaGain_dB = 0;
    pp->Tx_antennaGain_dB = 0;
    pp->AntennaType = 0;  //0 for antenna file, 1 for isotropic
    pp->Tx_Power_dB = 10 * log10(l1data.gps_eirp_watt);
    pp->AtmosphericLoss_dB = 0.0;
    pp->Rx_upsamplingFactor = 10;

    //use lidata.sc_num to select antenna patterns
    FILE *file;
    char *Rx_filename = getRxAntenna(l1data.sc_num, l1data.ddm_ant);
    file = fopen(Rx_filename,"rb");
    //file = fopen("../../Data/antennaRx_CYGNSS_Obs1_Nadir01_Port_E2ES_v3.bin","rb");

    if (file == NULL){
        printf("fail to open antenna file\n");
        exit(1);
    }
    else{
        printf("read antenna file\n");
    }
    double a,b;
    fread(&a,sizeof(double),1,file);
    fread(&b,sizeof(double),1,file);
    pp->Rx_numEl = (int) floor(a);
    pp->Rx_numAz = (int) floor(b);
    pp->Rx_numData = pp->Rx_numEl * pp->Rx_numAz;

    //allocate memory
    pp->data = (struct Rx_antennaDataPixel *)calloc(pp->Rx_numData,sizeof(struct Rx_antennaDataPixel));
    if(pp->data == NULL){
        printf("Error: Not enough memory to load antenna data\n");
        exit(1);
    }

    double start1, end1;
    //read from data file into memory
    int N = pp->Rx_numData;
    int M = 3; //num fields
    double *tempBuffer = (double *)calloc(N*M, sizeof(double));
    int numBytesReadFromFile = (int)fread(tempBuffer,sizeof(double),N*M,file);
    if (numBytesReadFromFile != (M*N)) {
        printf("Error: Reading error in antenna data: read %d, expected %d \n", numBytesReadFromFile, N*M);
        exit(3);
    }
    for(int i=0;i<N;i++){
        pp->data[i].gain_dB = tempBuffer[0*N+i];
        pp->data[i].el_deg  = tempBuffer[1*N+i];
        pp->data[i].az_deg  = tempBuffer[2*N+i];
    }
    free(tempBuffer);
    fclose(file);

}
void init_inputWindField_data(char dataFileName[], struct inputWindField *iwf, struct windInfo info){
    //read wind field from a data file
    //lon: 0-360
    //they are hard coded now
    printf("read wind from data\n");
    iwf->numPtsLon=info.numPtsLon;
    iwf->numPtsLat=info.numPtsLat;
    iwf->numPts=iwf->numPtsLon*iwf->numPtsLat;
    iwf->lat_min_deg = info.lat_min_deg;
    iwf->lat_max_deg = info.lat_max_deg;
    iwf->lon_min_deg = info.lon_min_deg;
    iwf->lon_max_deg = info.lon_max_deg;
    iwf->resolution_lat_deg = info.resolution;
    iwf->resolution_lon_deg = info.resolution;
    iwf->data = (struct inputWindFieldPixel *)calloc(iwf->numPts,sizeof(struct inputWindFieldPixel));
    FILE *file;
    file = fopen(dataFileName,"rb");
    if (file==NULL){
        printf("fail to open wind data\n");
        exit(1);
    }

    double *windData = (double *)calloc(iwf->numPts*2, sizeof(double));
    int numBytesReadFromFile = (int)fread(windData, sizeof(double),iwf->numPts*2,file);
    int ind;
    for(int lon = 0; lon < iwf->numPtsLon; lon++){
        for(int lat = 0;lat < iwf->numPtsLat; lat++){
            ind = lat * iwf->numPtsLon + lon;
            iwf->data[ind].windSpeed_U10_ms = windData[lat*iwf->numPtsLon+lon];
            iwf->data[ind].windSpeed_V10_ms = windData[iwf->numPts+lat*iwf->numPtsLon+lon];
            iwf->data[ind].windSpeed_ms=sqrt(iwf->data[ind].windSpeed_U10_ms*iwf->data[ind].windSpeed_U10_ms+
                                             iwf->data[ind].windSpeed_V10_ms*iwf->data[ind].windSpeed_V10_ms);
            iwf->data[ind].rainRate_mmhr = 0;
            iwf->data[ind].freezingHeight_m = 2;
            iwf->data[ind].lat_deg = iwf->lat_min_deg+lat*iwf->resolution_lat_deg;
            iwf->data[ind].lon_deg = iwf->lon_min_deg+lon*iwf->resolution_lon_deg;
        }
    }

    free(windData);
    fclose(file);
}

void init_Geometry(struct CYGNSSL1 l1data, struct Geometry *geom){

    memcpy(geom->rx_position_ecef_m, l1data.rx_position_ecef_m, 3*sizeof(double));
    memcpy(geom->rx_velocity_ecef_ms, l1data.rx_velocity_ecef_ms, 3*sizeof(double));
    memcpy(geom->tx_position_ecef_m, l1data.tx_position_ecef_m, 3*sizeof(double));
    memcpy(geom->tx_velocity_ecef_ms, l1data.tx_velocity_ecef_ms, 3*sizeof(double));
    memcpy(geom->sp_position_ecef_m, l1data.sp_position_ecef_m, 3*sizeof(double));
    memcpy(geom->sc_att_rad, l1data.sc_att_rad, 3*sizeof(double));

}

void init_DDM(struct CYGNSSL1 l1data, struct DDMfm *ddm_fm){
    ddm_fm->numDelaybins = 17;
    ddm_fm->numDopplerbins = 11;
    ddm_fm->delay_min_chip = l1data.ddm_sp_delay_row * (-0.25);
    ddm_fm->delay_max_chip = (17-1-l1data.ddm_sp_delay_row) * 0.25;
    ddm_fm->Doppler_min_Hz = -2500;
    ddm_fm->Doppler_max_Hz = 2500;
    ddm_fm->delay_resolution_chip = 0.25;
    ddm_fm->Doppler_resolution_Hz = 500;

    int numBin = ddm_fm->numDelaybins * ddm_fm->numDopplerbins;
    ddm_fm->data = (struct DDMPixel *)calloc(numBin,sizeof(struct DDMPixel));
    for(int i = 0; i < numBin; i++){
        ddm_fm->data[i].power = 0;
        ddm_fm->data[i].delay = 0;
        ddm_fm->data[i].Doppler = 0;
    }

    ddm_fm->check_sp_align = 1;

}

void init_Jacobian(struct Jacobian *jacob){
    jacob->numDDMbins = 187;
    jacob->numPts_LL = 200; //make it larger

    int numBin = jacob->numDDMbins * jacob->numPts_LL;
    jacob->data = (struct JacobianPixel *)calloc(numBin,sizeof(struct JacobianPixel));
    for(int i = 0; i< numBin; i++){
        jacob->data[i].value = 0;
        jacob->data[i].lat_deg = 0;
        jacob->data[i].lon_deg = 0;
    }
    for(int i = 0; i< jacob->numPts_LL; i++){
        jacob->Pts_lat_vec[i] = 0;
        jacob->Pts_lon_vec[i] = 0;
        jacob->Pts_ind_vec[i] = -1;
    }

}


char* getRxAntenna(int sc_num, int ddm_ant){
    char *filename;
    //char prefix[]="/Data/All_E2ES_antennas/V6";  //must return a pointer, but pointer cannot connect by strcat
    switch(sc_num){

        case 1:
            switch(ddm_ant){
                case 0:
                    printf("No antenna\n");break;
                case 2:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx1_starboard_E2ES_v6.bin"; break;
                case 3:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx1_port_E2ES_v6.bin"; break;
            }
            break;
        case 2:
            switch(ddm_ant){
                case 0:
                    printf("No antenna\n");break;
                case 2:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx2_starboard_E2ES_v6.bin"; break;
                case 3:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx2_port_E2ES_v6.bin"; break;
            }
            break;
        case 3:
            switch(ddm_ant){
                case 0:
                    printf("No antenna\n");break;
                case 2:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx3_starboard_E2ES_v6.bin"; break;
                case 3:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx3_port_E2ES_v6.bin"; break;
            }
            break;

        case 4:
            switch(ddm_ant){
                case 0:
                    printf("No antenna\n");break;
                case 2:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx4_starboard_E2ES_v6.bin"; break;
                case 3:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx4_port_E2ES_v6.bin"; break;

            }
            break;

        case 5:
            switch(ddm_ant){
                case 0:
                    printf("No antenna\n");break;
                case 2:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx5_starboard_E2ES_v6.bin"; break;
                case 3:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx5_port_E2ES_v6.bin"; break;
            }
            break;
        case 6:
            switch(ddm_ant){
                case 0:
                    printf("No antenna\n");break;
                case 2:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx6_starboard_E2ES_v6.bin"; break;
                case 3:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx6_port_E2ES_v6.bin"; break;
            }
            break;
        case 7:
            switch(ddm_ant){
                case 0:
                    printf("No antenna\n");break;
                case 2:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx7_starboard_E2ES_v6.bin"; break;
                case 3:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx7_port_E2ES_v6.bin"; break;
            }
            break;
        case 8:
            switch(ddm_ant){
                case 0:
                    printf("No antenna\n");break;
                case 2:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx8_starboard_E2ES_v6.bin"; break;
                case 3:
                    filename = "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx8_port_E2ES_v6.bin"; break;
            }
            break;

    }

    return filename;

}
