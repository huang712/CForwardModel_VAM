//
// Created by Feixiong Huang on 1/22/18.
//

//strutures to store CYGNSS measurement data and related functions
#ifndef CFORWARDMODEL_CYGNSS_H
#define CFORWARDMODEL_CYGNSS_H

//#define DIM 85633 //number of sample in L1 data, hard code this to avoid dynamic array

struct CYGNSSL1
{
    int quality_flags;
    int utc_sec;
    int sc_num;  //1-8
    int index; //index of measurement of the day
    int spNum; //1-4 reflection number
    int prn_code; //GPS PRN 1-32

    int ddm_ant; // antenna type 2 = starboard, 3 = port
    double ant_temperature_cels;
    double gps_eirp_watt;
    double noise_figue;

    int ddm_peak_delay_row;
    int ddm_peak_dopp_col;
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
    double DDM_brcs[17][11];
};

void readL1data(char L1dataFilename[], int sampleIndex, int ddm_index, struct CYGNSSL1 *l1data);

int readnc_int_1d(int ncid, char varName[], int index);
int readnc_int_2d(int ncid, char varName[], int index, int spNum);
float readnc_float_1d(int ncid, char varName[], int index);
float readnc_float_2d(int ncid, char varName[], int index, int ddm_index);

void DDMobs_saveToFile(struct CYGNSSL1 l1data, int index, int pathType);

#endif //CFORWARDMODEL_CYGNSS_H
