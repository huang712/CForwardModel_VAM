//
// Created by Feixiong Huang on 10/19/17.
//

#ifndef CFORWARDMODEL_INPUTOUPUT_H
#define CFORWARDMODEL_INPUTOUPUT_H

#include <stdio.h>
#include <stdlib.h>
#include "cygnss.h"

struct metadata
{
    //high resolution DDM
    int numDelaybins;
    int numDopplerbins;
    double delayRez_chips;
    double dopplerRes_Hz;

    //specularBin
    int specular_delayBinIdx;
    int specular_dopplerBinIdx;

    //downsampled DDM
    int resample_startBin[2];
    int resample_resolution_bins[2];
    int resample_numBins[2];

    //thermalNoise
    double temp_K;
    double noiseFigure_dB;
    unsigned int thermalNoiseOnOff;
    double excess_noisefloor_dB;

    //surfaceGrid
    double grid_resolution_m;
    int numGridPoints[2];
    unsigned int surfaceCurvatureType;

    int utc_sec;
    int prn_code;
    double meas_ddm_sp_index[2];

};


struct powerParm
{
    double Rx_antennaGain_dB; //excess receiver antenna gain
    double Tx_antennaGain_dB;
    unsigned int AntennaType;
    double Tx_Power_dB;
    double AtmosphericLoss_dB;
    double Rx_upsamplingFactor;
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
    struct inputWindFieldPixel *data;  //structure in structure
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
    double sc_att_rad[3];
};

//output data



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
    double Pts_lat_vec[144]; //long enough to larger than numPts_LL
    double Pts_lon_vec[144];
    int Pts_ind_vec[144];
    struct JacobianPixel *data; //structure in structure
};

struct JacobianPixel
{
    double value;
    double lat_deg, lon_deg;
};

//forward model.c : Forward model main function
void forwardModel(struct metadata meta, struct powerParm pp,
                  struct inputWindField iwf, struct Geometry geom,
                  struct DDMfm *ddm_fm, struct Jacobian *jacob, int option);

//initialization.c : initialization functions
void init_metadata(struct CYGNSSL1 l1data, struct metadata *meta);
void init_powerParm(struct CYGNSSL1 l1data, struct powerParm *pp);
void init_inputWindField_core(char windFileName[], struct inputWindField *iwf);
void init_inputWindField_synoptic(char windFileName[], struct inputWindField *iwf);
void init_inputWindField_data(char dataName[], struct inputWindField *iwf);
void init_Geometry(struct CYGNSSL1 l1data, struct Geometry *geom);
void init_DDM(struct CYGNSSL1 l1data, struct DDMfm *ddm_fm);
void init_Jacobian(struct Jacobian *jacob);
char* getRxAntenna(int sc_num, int ddm_ant);

//saveFile.c : file saving functions (for debug)
void DDMfm_saveToFile(struct DDMfm ddm_fm, int index, int pathType);
void Jacobian_saveToFile(struct Jacobian jacob, int index, int pathType);
void PtsVec_saveToFile(struct Jacobian jacob);

#endif //CFORWARDMODEL_INPUTOUPUT_H


