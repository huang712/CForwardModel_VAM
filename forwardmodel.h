//---------------------------------------------------------------------------
//
// Define structures and functions of forward model
// Created by Feixiong Huang on 2017/10/19/
//
//***************************************************************************

#ifndef CFORWARDMODEL_INPUTOUPUT_H
#define CFORWARDMODEL_INPUTOUPUT_H

#include <stdio.h>
#include <stdlib.h>
#include "cygnss.h"

//#define ANTENNA_PATH "/scratch/andoria/a/huang712/FM_data/All_E2ES_antennas/V6/Rx"
//#define PRN_ACF_FILE "/scratch/andoria/a/huang712/FM_data/PRN_ACF.bin"
//#define GMF_PATH "/scratch/andoria/a/huang712/FM_data/GMF_data/"

#define ANTENNA_PATH "/users/fax/CYGNSS/Data/All_E2ES_antennas/V6/Rx"
#define PRN_ACF_FILE "/users/fax/CYGNSS/Data/PRN_ACF.bin"
#define GMF_PATH "/users/fax/CYGNSS/Data/GMF_data/"

#define GMF_OnOff 0  // use modified CYGNSS GMF model unless Katzberg model
#define fastMode_OnOff 1  // 1 for fast mode (in initialization.c)

struct option
{
    int JacobOnOff;
    int thermalNoiseOnOff;
};

struct windInfo  // information of wind field in the data from the config file
{
    int numPtsLat, numPtsLon;
    double lat_min_deg, lat_max_deg, lon_min_deg, lon_max_deg;
    double resolution;
};

struct metadata
{
    // high resolution DDM
    int numDelaybins;
    int numDopplerbins;
    double delayRez_chips;
    double dopplerRes_Hz;

    // specularBin
    int specular_delayBinIdx;
    int specular_dopplerBinIdx;

    // downsampled DDM
    int resample_startBin[2];
    int resample_resolution_bins[2];
    int resample_numBins[2];

    // thermalNoise
    double temp_K;
    double noiseFigure_dB;
    double excess_noisefloor_dB;

    // surfaceGrid
    double grid_resolution_m;
    int numGridPoints[2];
    int surfaceCurvatureType;

    int utc_sec;
    int prn_code;
    double meas_ddm_sp_index[2];
    double fresnel_coeff2;  // square of reflectivity
};


struct powerParm
{
    double Rx_antennaGain_dB;  // excess receiver antenna gain
    double Tx_antennaGain_dB;
    unsigned int AntennaType;
    double Tx_eirp_watt;
    double AtmosphericLoss_dB;
    int Rx_numEl;
    int Rx_numAz;
    int Rx_numData;
    struct Rx_antennaDataPixel *data;
};

struct Rx_antennaDataPixel
{
    double gain_dB, el_deg, az_deg;
};

struct inputWindField
{
    int numPtsLat, numPtsLon, numPts;
    double lat_min_deg, lat_max_deg, lon_min_deg, lon_max_deg;
    double resolution_lat_deg, resolution_lon_deg;
    struct inputWindFieldPixel *data;  // structure in structure
};

struct inputWindFieldPixel
{
    double windSpeed_U10_ms, windSpeed_V10_ms,windSpeed_ms;
    double rainRate_mmhr, freezingHeight_m;
    double lat_deg, lon_deg;
};

struct Geometry
{
    double rx_position_ecef_m[3];
    double rx_velocity_ecef_ms[3];
    double tx_position_ecef_m[3];
    double tx_velocity_ecef_ms[3];
    double sp_position_ecef_m[3];
    double sp_velocity_ecef_m[3];
    double sc_att_rad[3];
    double sp_lat, sp_lon;
};

//**************************** Output Data ****************************

struct DDMfm
{
    int numDelaybins;
    int numDopplerbins;
    double delay_min_chip, delay_max_chip;
    double Doppler_min_Hz, Doppler_max_Hz;
    double delay_resolution_chip, Doppler_resolution_Hz;
    struct DDMPixel *data; //structure in structure
    unsigned int check_sp_align;
};

struct DDMPixel
{
    double power;
    double delay;
    double Doppler;
};

struct Jacobian
{
    int numDDMbins;
    int numPts_LL;
    double Pts_lat_vec[200]; // long enough to larger than numPts_LL
    double Pts_lon_vec[200];
    int Pts_ind_vec[200];
    struct JacobianPixel *data; // structure in structure
};

struct JacobianPixel
{
    double value;
    double lat_deg, lon_deg;
};

//**************************** Functions ****************************

//forward model.c : Forward model main function
void forwardModel(struct metadata meta, struct powerParm pp,
                  struct inputWindField iwf, struct Geometry geom,
                  struct DDMfm *ddm_fm, struct Jacobian *jacob, struct option opt);

// initialization.c : initialization functions
void init_metadata(struct CYGNSSL1 l1data, struct metadata *meta);
void init_powerParm(struct CYGNSSL1 l1data, struct powerParm *pp);
void init_inputWindField_data(char dataName[], struct inputWindField *iwf, struct windInfo);
void init_Geometry(struct CYGNSSL1 l1data, struct Geometry *geom);
void init_DDM(struct CYGNSSL1 l1data, struct DDMfm *ddm_fm);
void init_Jacobian(struct Jacobian *jacob);
char* getRxAntenna(int sc_num, int ddm_ant);

// saveFile.c : file saving functions (for debug)
void DDMfm_saveToFile(struct DDMfm ddm_fm, int index, int pathType, char saveDir[1000]);
void Jacobian_saveToFile(struct Jacobian jacob, int index, int pathType, char saveDir[1000]);
void indexLL_saveToFile(struct Jacobian jacob, char saveDir[1000]);

#endif //CFORWARDMODEL_INPUTOUPUT_H


