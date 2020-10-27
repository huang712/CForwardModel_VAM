//---------------------------------------------------------------------------
//
// strutures to store CYGNSS L1 data and related functions
// Created by Feixiong Huang on 1/22/18
//
//***************************************************************************

#ifndef CFORWARDMODEL_CYGNSS_H
#define CFORWARDMODEL_CYGNSS_H

struct CYGNSSL1
{
    long int quality_flags;
    double utc_sec;
    int sc_num;  // 1-8 CYGNSS SV number
    int index;  // sample index
    int spNum;  // 1-4 reflection channel number
    int prn_code;  // 1-32 GPS PRN
    double inc_angle;  // specular point incidence angle
    double fresnel_coeff2;  // square of reflectivity
    double sp_rx_gain;

    int ddm_ant;  // antenna type 2 = starboard, 3 = port
    double ant_temperature_cels;
    double gps_eirp_watt;
    double noise_figure;

    double ddm_sp_delay_row;
    double ddm_sp_dopp_col;

    double rx_position_ecef_m[3];
    double rx_velocity_ecef_ms[3];
    double tx_position_ecef_m[3];
    double tx_velocity_ecef_ms[3];
    double sp_position_ecef_m[3];
    double sp_lat, sp_lon;
    double sc_att_rad[3];

    double DDM_power[17][11];
    //double DDM_brcs[17][11];
};

void readL1data(char L1dataFilename[], int sampleIndex, int ddm_index, struct CYGNSSL1 *l1data);

int readnc_int_1d(int ncid, char varName[], int index);
int readnc_int_2d(int ncid, char varName[], int index, int ddm_index);
long int readnc_long_2d(int ncid, char varName[], int index, int ddm_index);
float readnc_float_1d(int ncid, char varName[], int index);
float readnc_float_2d(int ncid, char varName[], int index, int ddm_index);

void DDMobs_saveToFile(struct CYGNSSL1 l1data, int index, int pathType, char saveDir[1000]);

#endif //CFORWARDMODEL_CYGNSS_H
