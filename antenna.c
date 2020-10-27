//---------------------------------------------------------------------------
//
// This file had antenna related functions.
// Currently, we just support a fixed Tx EIRP and a file-based Rx pattern.
//
//****************************************************************************/

#include "gnssr.h"
#include "forwardmodel.h"

double RXG_DB, TXG_DB;
int convertAzElToIndex( double azimuth_rad, double elevation_rad);
double get_GPS_satAnt_gainPattern_dB( double elevation_rad );
int antennaInputType;

void antenna_initialize(struct powerParm pp){
    // This is to initilize the receiver antenna pattern
    // For GPS transmitter, we have pp->Tx_eirp_watt

    RXG_DB = pp.Rx_antennaGain_dB;  // excess gain, zero by default
    TXG_DB = pp.Tx_antennaGain_dB;  // not used, zero by default
    antennaInputType = pp.AntennaType;

    //double TXP_DB   = pp.Tx_Power_dB;
    //double ATTEN_DB = pp.AtmosphericLoss_dB;

    int numEl = pp.Rx_numEl;
    int numAz = pp.Rx_numAz;

    // Allocate memory
    antenna.data = (antennaDataPixel *)calloc(numEl * numAz, sizeof(antennaDataPixel));
    if(antenna.data == NULL) {
        printf("Error: Not enough memory to load antenna data\n");
        exit(1);
    }

    // Read from data file into memory
    int N = numEl * numAz;
    int M = 3;  // num fields

    for(int i=0;i<N;i++){
        antenna.data[i].gain_dB = pp.data[i].gain_dB;
        antenna.data[i].el_deg  = pp.data[i].el_deg;
        antenna.data[i].az_deg  = pp.data[i].az_deg;
    }
    antenna.numAz = numAz;
    antenna.numEl = numEl;

}

double antenna_getGain_abs(  antennaType antType, polarizationType pol, double angles_rad[2] ){
    // Not working for GPS transmitter

    double azimuth_rad   = angles_rad[0];
    double elevation_rad = angles_rad[1];

    switch (antType) {

        case CYGNSS_NADIR_ANT:
            switch (pol) {
                case RHCP: return -INFINITY;
                case LHCP:
                    if( antennaInputType == 0 )
                        return pow(10,( ( RXG_DB + antenna.data[convertAzElToIndex( azimuth_rad, elevation_rad )].gain_dB )/10));
                    if( antennaInputType == 1 )
                        return pow(10,( ( RXG_DB )/10));
            }
        case CYGNSS_ZENITH_ANT:
            switch (pol) {
                case RHCP: return -INFINITY;
                case LHCP: return -INFINITY;
            }
        case GPS_SAT_ANT: // not working
            switch (pol) {
                case RHCP: return pow(10,(TXG_DB/10));  // get_GPS_satAnt_gainPattern_dB( elevation_rad );
                case LHCP: return -INFINITY;
            }
        default:
            printf("Error: bad antenna type");
            exit(0);
    }
}

int convertAzElToIndex( double azimuth_rad, double elevation_rad){
    // In Matlab, we upsample the antenna pattern by this amount
    // to get a nice smooth antenna pattern over the surface grid.

    if( azimuth_rad < 0 )     { azimuth_rad += 2*pi; }
    if( azimuth_rad >=  2*pi ){ azimuth_rad -= 2*pi; }
    int az = (int)floor(azimuth_rad * R2D * upsamplingFactor);
    int el = (int)floor(elevation_rad * R2D * upsamplingFactor) +  90 * upsamplingFactor;

    int numEl = 180 * upsamplingFactor + 1;
    int numAz = 360 * upsamplingFactor + 1;

    return az*numEl + el;
}


