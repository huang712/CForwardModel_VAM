//
// Created by Feixiong Huang on 1/22/18.
//
#include <stdio.h>
#include "cygnss.h"
#include <netcdf.h>

#define ERR(e) {printf("Error: %s\n", nc_strerror(e));}

void readL1data(char L1dataFilename[], int sampleIndex, int ddm_index, struct CYGNSSL1 *l1data) {
    printf("read CYGNSS L1 data\n");

    int ncid;
    int sample_id, sample_num, sc_num_id;
    int retval;  //Error handling

    // Open the file
    if ((retval = nc_open(L1dataFilename, NC_NOWRITE, &ncid))) ERR(retval);

    //read and check number of samples of the data
    if ((retval = nc_inq_dimid(ncid, "sample", &sample_id))) ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, sample_id, &sample_num))) ERR(retval);
    if (sampleIndex > sample_num-1) printf("sample index is larger than length of data\n");

    //read CYGNSS spacecraft ID 1-8
    if ((retval = nc_inq_varid(ncid, "spacecraft_num", &sc_num_id))) ERR(retval);
    if ((retval = nc_get_var_int(ncid, sc_num_id, &l1data->sc_num))) ERR(retval);

    //read DDM power
    size_t start[4]={sampleIndex,ddm_index,0,0};
    size_t count[4]={1,1,17,11};
    int power_analog_id;
    if ((retval = nc_inq_varid(ncid, "power_analog", &power_analog_id))) ERR(retval);
    if ((retval = nc_get_vara_double(ncid, power_analog_id, start, count, &l1data->DDM_power[0][0]))) ERR(retval);

    //read variables
    l1data->rx_position_ecef_m[0] = readnc_int_1d(ncid, "sc_pos_x", sampleIndex);
    l1data->rx_position_ecef_m[1] = readnc_int_1d(ncid, "sc_pos_y", sampleIndex);
    l1data->rx_position_ecef_m[2] = readnc_int_1d(ncid, "sc_pos_z", sampleIndex);
    l1data->rx_velocity_ecef_ms[0] = readnc_int_1d(ncid, "sc_vel_x", sampleIndex);
    l1data->rx_velocity_ecef_ms[1] = readnc_int_1d(ncid, "sc_vel_y", sampleIndex);
    l1data->rx_velocity_ecef_ms[2] = readnc_int_1d(ncid, "sc_vel_z", sampleIndex);
    l1data->tx_position_ecef_m[0] = readnc_int_2d(ncid, "tx_pos_x", sampleIndex, ddm_index);
    l1data->tx_position_ecef_m[1] = readnc_int_2d(ncid, "tx_pos_y", sampleIndex, ddm_index);
    l1data->tx_position_ecef_m[2] = readnc_int_2d(ncid, "tx_pos_z", sampleIndex, ddm_index);
    l1data->tx_velocity_ecef_ms[0] = readnc_int_2d(ncid, "tx_vel_x", sampleIndex, ddm_index);
    l1data->tx_velocity_ecef_ms[1] = readnc_int_2d(ncid, "tx_vel_y", sampleIndex, ddm_index);
    l1data->tx_velocity_ecef_ms[2] = readnc_int_2d(ncid, "tx_vel_z", sampleIndex, ddm_index);
    l1data->sp_position_ecef_m[0] = readnc_int_2d(ncid, "sp_pos_x", sampleIndex, ddm_index);
    l1data->sp_position_ecef_m[1] = readnc_int_2d(ncid, "sp_pos_y", sampleIndex, ddm_index);
    l1data->sp_position_ecef_m[2] = readnc_int_2d(ncid, "sp_pos_z", sampleIndex, ddm_index);
    l1data->sp_lat = readnc_float_2d(ncid, "sp_lat", sampleIndex, ddm_index);
    l1data->sp_lon = readnc_float_2d(ncid, "sp_lon", sampleIndex, ddm_index);
    l1data->sc_att_rad[0] = readnc_float_1d(ncid, "sc_pitch", sampleIndex);
    l1data->sc_att_rad[1] = readnc_float_1d(ncid, "sc_roll", sampleIndex);
    l1data->sc_att_rad[2] = readnc_float_1d(ncid, "sc_yaw", sampleIndex);

    l1data->utc_sec = readnc_int_1d(ncid, "ddm_timestamp_utc", sampleIndex);
    l1data->quality_flags = readnc_int_2d(ncid, "quality_flags", sampleIndex, ddm_index);
    l1data->ddm_peak_delay_row = readnc_int_2d(ncid, "brcs_ddm_peak_bin_delay_row", sampleIndex, ddm_index);
    l1data->ddm_peak_dopp_col = readnc_int_2d(ncid, "brcs_ddm_peak_bin_dopp_col", sampleIndex, ddm_index);
    l1data->ddm_sp_delay_row = readnc_float_2d(ncid, "brcs_ddm_sp_bin_delay_row", sampleIndex, ddm_index);
    l1data->ddm_sp_dopp_col = readnc_float_2d(ncid, "brcs_ddm_sp_bin_dopp_col", sampleIndex, ddm_index);
    l1data->index = sampleIndex;
    l1data->spNum = ddm_index;
    l1data->prn_code = readnc_int_2d(ncid, "prn_code", sampleIndex, ddm_index);
    l1data->gps_eirp_watt = readnc_float_2d(ncid, "gps_eirp", sampleIndex, ddm_index);
    l1data->ddm_ant = readnc_int_2d(ncid, "ddm_ant", sampleIndex, ddm_index);
    l1data->noise_figue = readnc_float_2d(ncid, "lna_noise_figure", sampleIndex, ddm_index);

    if (l1data->ddm_ant == 2){
        l1data->ant_temperature_cels = readnc_float_1d(ncid, "lna_temp_nadir_starboard", sampleIndex);
    }
    else if (l1data->ddm_ant == 3){
        l1data->ant_temperature_cels = readnc_float_1d(ncid, "lna_temp_nadir_port", sampleIndex);
    }
    else {printf("no antenna temperature\n");}

    if ((retval = nc_close(ncid))) ERR(retval);
}

int readnc_int_1d(int ncid, char varName[], int index){
    int retval;
    int var_id;
    int var;
    size_t start[1]={index};
    size_t count[1]={1};
    if ((retval = nc_inq_varid(ncid, varName, &var_id))) ERR(retval);
    if ((retval = nc_get_vara_int(ncid, var_id, start, count, &var))) ERR(retval);
    return var;
}

int readnc_int_2d(int ncid, char varName[], int index, int ddm_index){
    int retval;
    int var_id;
    int var;
    size_t start[2]={index,ddm_index};
    size_t count[2]={1,1};
    if ((retval = nc_inq_varid(ncid, varName, &var_id))) ERR(retval);
    if ((retval = nc_get_vara_int(ncid, var_id, start, count, &var))) ERR(retval);
    return var;
}

float readnc_float_1d(int ncid, char varName[], int index){
    int retval;
    int var_id;
    float var;
    size_t start[1]={index};
    size_t count[1]={1};
    if ((retval = nc_inq_varid(ncid, varName, &var_id))) ERR(retval);
    if ((retval = nc_get_vara_float(ncid, var_id, start, count, &var))) ERR(retval);
    return var;
}

float readnc_float_2d(int ncid, char varName[], int index, int ddm_index){
    int retval;
    int var_id;
    float var;
    size_t start[2]={index,ddm_index};
    size_t count[2]={1,1};
    if ((retval = nc_inq_varid(ncid, varName, &var_id))) ERR(retval);
    if ((retval = nc_get_vara_float(ncid, var_id, start, count, &var))) ERR(retval);
    return var;
}

void DDMobs_saveToFile(struct CYGNSSL1 l1data, int index, int pathType) {

    double val;
    char filename[50];

    FILE *outp;
    switch(pathType){
        case 0:
            outp = fopen("DDMobs.dat","wb");
            break;
        case 1:
            sprintf(filename, "DDMobs/DDMobs%d.dat", index);
            outp = fopen(filename, "wb");
            break;
    }

    for (int i = 0; i < 11 ; i++) {
        for (int j = 0; j < 17 ; j++){
            val = l1data.DDM_power[j][i];
            fwrite(&val, 1, sizeof(double), outp);
        }
    }
    fclose(outp);
    printf("save CYGNSS DDM into file\n");
}