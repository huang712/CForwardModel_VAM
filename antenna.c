//
// Created by Feixiong Huang on 10/23/17.
//

//***************************************************************************
// This file is part of the CYGNSS E2ES.
//
// CYGNSS E2ES is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// CYGNSS E2ES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CYGNSS E2ES.  If not, see <http://www.gnu.org/licenses/>.
//
//---------------------------------------------------------------------------
//
//  This file had antenna related functions.  Currently, we just support a
//  fixed gain Tx pattern and a file-based Rx pattern.  TODO: implement a
//  few other antenna options.
//
//****************************************************************************/


#include "gnssr.h"
#include "forwardmodel.h"

double RXG_DB, TXG_DB;
int convertAzElToIndex( double azimuth_rad, double elevation_rad);
double get_GPS_satAnt_gainPattern_dB( double elevation_rad );
int antennaInputType;

void antenna_initialize(struct powerParm pp){

    // note that Rx_antennaGain_dB in the config file is now the excess gain!!!
    // TODO: make this more clear

    RXG_DB = pp.Rx_antennaGain_dB; //excess gain
    TXG_DB = pp.Tx_antennaGain_dB;
    antennaInputType = pp.AntennaType;

    double TXP_DB   = pp.Tx_Power_dB;
    double ATTEN_DB = pp.AtmosphericLoss_dB;

    //int plotAntennaData   = 0;
    int numEl = pp.Rx_numEl;
    int numAz = pp.Rx_numAz;

    // allocate memory
    antenna.data = (antennaDataPixel *)calloc(numEl * numAz, sizeof(antennaDataPixel));
    if(antenna.data == NULL) {
        printf("Error: Not enough memory to load antenna data\n");
        exit(1);
    }

    // read from data file into memory
    int N = numEl * numAz;
    int M = 3; // num fields

    for(int i=0;i<N;i++){
        antenna.data[i].gain_dB = pp.data[i].gain_dB;
        antenna.data[i].el_deg  = pp.data[i].el_deg;
        antenna.data[i].az_deg  = pp.data[i].az_deg;
    }
    antenna.numAz = numAz;
    antenna.numEl = numEl;

}

double antenna_getGain_abs(  antennaType antType, polarizationType pol, double angles_rad[2] ){
    // this is pretty simple now, but will be expanded ...

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
        case GPS_SAT_ANT:
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
    // in Matlab, we upsample the antenna pattern by this amount
    // to get a nice smooth antenna pattern over the surface grid.
    // we should really move this upsampling into the initialization
    // code above.  TODO.

    int upsamplingPower = 3;

    //int upsamplingFactor = pp.upsamplingFactor;
    //int upsamplingFactor = pow(2,upsamplingPower);

    //if( azimuth_rad < -pi ){ azimuth_rad += 2*pi; }
    //if( azimuth_rad >=  pi ){ azimuth_rad -= 2*pi; }
    //int az = floor(azimuth_rad * R2D * upsamplingFactor)   + 180 * upsamplingFactor;

    if( azimuth_rad < 0 )     { azimuth_rad += 2*pi; }
    if( azimuth_rad >=  2*pi ){ azimuth_rad -= 2*pi; }
    int az = (int)floor(azimuth_rad * R2D * upsamplingFactor);
    int el = (int)floor(elevation_rad * R2D * upsamplingFactor) +  90 * upsamplingFactor;

    int numEl = 180 * upsamplingFactor + 1;
    int numAz = 360 * upsamplingFactor + 1;

    // correct for fact that I should have gone 0-360 in the file instead of 0-359
    // TODO: fix this! (fixed now)
    //if( az >= numAz ){ az = numAz - 1; }

    return az*numEl + el;
}


/****************************************************************************/
//  Load antenna data from binary file
/****************************************************************************/

/*
void antenna_loadfile(const char *filename ) {

    FILE *file;

    // open file
    file = fopen(filename,"rb");
    if (!file) {
        printf("Error: Couldn't open file %s in antenna_loadfile\n", filename);
        exit(0);
    }

    // read windfield dimensions
    double a,b;
    fread(&a,sizeof(double),1,file);
    fread(&b,sizeof(double),1,file);
    int numEl = floor(a);
    int numAz = floor(b);

    // allocate memory
    antenna.data = (antennaDataPixel *)calloc(numEl * numAz, sizeof(antennaDataPixel));
    if(antenna.data == NULL) {
        printf("Error: Not enough memory to load antenna data\n");
        exit(1);
    }

    // read from data file into memory
    int N = numEl * numAz;
    int M = 3; // num fields
    double *tempBuffer = (double *)calloc(N*M, sizeof(double));
    int numBytesReadFromFile = (int)fread(tempBuffer,sizeof(double),N*M,file);
    if (numBytesReadFromFile != (M*N)) {
        printf("Error: Reading error in antenna data: read %d, expected %d \n", numBytesReadFromFile, N*M);
        free(antenna.data);
        exit(3);
    }
    for(int i=0;i<N;i++){
        antenna.data[i].gain_dB = tempBuffer[0*N+i];
        antenna.data[i].el_deg  = tempBuffer[1*N+i];
        antenna.data[i].az_deg  = tempBuffer[2*N+i];
    }
    antenna.numAz = numAz;
    antenna.numEl = numEl;
    free(tempBuffer);
    fclose(file);

}
*/

/****************************************************************************/
//  Plot antenna data to .png image (for debugging)
/****************************************************************************/

/*
void antenna_save2PNG( void ) {
    unsigned height = antenna.numAz;
    unsigned width  = antenna.numEl;
    unsigned char* image;
    unsigned x, y, idx;
    double val, min, max;
    char filename[100];

    double* vals  = calloc(width * height, sizeof(double));

    fprintf(outputPtr,"\n Plotting Antenna Data --------------------------------\n");

    for(int type = 1; type <= 8; type++){

        switch (type){
            case 1:  min = 0; max = 0; strcpy(filename, "ant_gain_dB.png"); break;
            case 2:  min = 0; max = 0; strcpy(filename, "ant_el_deg.png"); break;
            case 3:  min = 0; max = 0; strcpy(filename, "ant_az_deg.png"); break;
            case 4:  min = 0; max = 0; strcpy(filename, "ant_idx_x.png"); break;
            case 5:  min = 0; max = 0; strcpy(filename, "ant_idx_y.png"); break;
            case 6:  min = 0; max = 0; strcpy(filename, "ant_idx.png"); break;
            case 7:  min = 0; max = 0; strcpy(filename, "ant_idxInt.png"); break;
            case 8:  min = 0; max = 0; strcpy(filename, "ant_idxIntGain.png"); break;
        }

        for(y = 0; y < height; y++){
            for(x = 0; x < width; x++) {
                idx = x*height + y;
                switch (type){
                    case 1:  val = antenna.data[idx].gain_dB; break;
                    case 2:  val = antenna.data[idx].el_deg;  break;
                    case 3:  val = antenna.data[idx].az_deg;  break;
                    case 4:  val = x;  break;
                    case 5:  val = y;  break;
                    case 6:  val = idx;  break;
                    case 7:  val = convertAzElToIndex( antenna.data[idx].az_deg * D2R, antenna.data[idx].el_deg * D2R ); break;
                    case 8:  val = antenna.data[convertAzElToIndex( antenna.data[idx].az_deg * D2R, antenna.data[idx].el_deg * D2R )].gain_dB; break;
                }
                vals[idx] = val;
            }
        }

        fprintf(outputPtr,"   ... (%d) plotting %s ", type, filename);
        fflush(outputPtr);
        image_createImageFromDouble(&image, vals, width, height, min, max);
        image_flipud(width, height, image);
        image_plotFigure(filename, image, width, height);
        free(image);

    }
    fprintf(outputPtr,"\n");
}
*/

//******************************************************************************/
// GPS satellite antenna pattern (nadir angle deg, gain dBi)

const double GPS_satAnt_gainPattern[91][2] = {
        { -30.000000, -80.422230 },
        { -29.333333, -71.773414 },
        { -28.666667, -63.670095 },
        { -28.000000, -56.094762 },
        { -27.333333, -49.029905 },
        { -26.666667, -42.458014 },
        { -26.000000, -36.361578 },
        { -25.333333, -30.723089 },
        { -24.666667, -25.525035 },
        { -24.000000, -20.749906 },
        { -23.333333, -16.380193 },
        { -22.666667, -12.398385 },
        { -22.000000, -8.786971 },
        { -21.333333, -5.528443 },
        { -20.666667, -2.605289 },
        { -20.000000, 0.000000 },
        { -19.333333, 2.304935 },
        { -18.666667, 4.327025 },
        { -18.000000, 6.083782 },
        { -17.333333, 7.592714 },
        { -16.666667, 8.871333 },
        { -16.000000, 9.937148 },
        { -15.333333, 10.807669 },
        { -14.666667, 11.500407 },
        { -14.000000, 12.032872 },
        { -13.333333, 12.422573 },
        { -12.666667, 12.687022 },
        { -12.000000, 12.843727 },
        { -11.333333, 12.910200 },
        { -10.666667, 12.903951 },
        { -10.000000, 12.841219 },
        { -9.333333, 12.731376 },
        { -8.666667, 12.581484 },
        { -8.000000, 12.398600 },
        { -7.333333, 12.189784 },
        { -6.666667, 11.962093 },
        { -6.000000, 11.722587 },
        { -5.333333, 11.478323 },
        { -4.666667, 11.236362 },
        { -4.000000, 11.003760 },
        { -3.333333, 10.787578 },
        { -2.666667, 10.594873 },
        { -2.000000, 10.432704 },
        { -1.333333, 10.308130 },
        { -0.666667, 10.228209 },
        { 0.000000, 10.200000 },
        { 0.666667, 10.228209 },
        { 1.333333, 10.308130 },
        { 2.000000, 10.432704 },
        { 2.666667, 10.594873 },
        { 3.333333, 10.787578 },
        { 4.000000, 11.003760 },
        { 4.666667, 11.236362 },
        { 5.333333, 11.478323 },
        { 6.000000, 11.722587 },
        { 6.666667, 11.962093 },
        { 7.333333, 12.189784 },
        { 8.000000, 12.398600 },
        { 8.666667, 12.581484 },
        { 9.333333, 12.731376 },
        { 10.000000, 12.841219 },
        { 10.666667, 12.903951 },
        { 11.333333, 12.910200 },
        { 12.000000, 12.843727 },
        { 12.666667, 12.687022 },
        { 13.333333, 12.422573 },
        { 14.000000, 12.032872 },
        { 14.666667, 11.500407 },
        { 15.333333, 10.807669 },
        { 16.000000, 9.937148 },
        { 16.666667, 8.871333 },
        { 17.333333, 7.592714 },
        { 18.000000, 6.083782 },
        { 18.666667, 4.327025 },
        { 19.333333, 2.304935 },
        { 20.000000, -0.000000 },
        { 20.666667, -2.605289 },
        { 21.333333, -5.528443 },
        { 22.000000, -8.786971 },
        { 22.666667, -12.398385 },
        { 23.333333, -16.380193 },
        { 24.000000, -20.749906 },
        { 24.666667, -25.525035 },
        { 25.333333, -30.723089 },
        { 26.000000, -36.361578 },
        { 26.666667, -42.458014 },
        { 27.333333, -49.029905 },
        { 28.000000, -56.094762 },
        { 28.666667, -63.670095 },
        { 29.333333, -71.773414 },
        { 30.000000, -80.422230 }
};

double get_GPS_satAnt_gainPattern_dB( double elevation_rad ){

    // convert to angle of ascention
    double asc = 90 - elevation_rad * R2D;

    // this should never happen
    if( fabs(asc) >= 30 ){
        printf("Error: transmit pattern angle unexpected range\n");
        exit(0);
    }

    // linearly interpolate data
    for(int i=0; i < 91; i++){
        if( GPS_satAnt_gainPattern[i][0] >= asc ){
            double d = ( GPS_satAnt_gainPattern[i][0] - GPS_satAnt_gainPattern[i-1][0] );
            double f = (GPS_satAnt_gainPattern[i][0] - asc) / d;
            return ( GPS_satAnt_gainPattern[i-1][0] * f + GPS_satAnt_gainPattern[i][0] * (1 - f) );
        }
    }

    return -80;
}

