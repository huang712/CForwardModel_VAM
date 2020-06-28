//---------------------------------------------------------------------------
//
// Code for loading the wind field data and converting to mss.  Surface.c has
// code which requests the wind field at different locations.
//
//****************************************************************************/


#include "gnssr.h"
#include "forwardmodel.h"

/*
void ddmaLUT_initialize(void){
    //added by Feixiong
    //load DDMA LUT from CYGNSS
    //ws = 0.05:0.1:69.95
    //theta = 1:90
    FILE *file;
    char *filename = "../../Data/ddmaLUTv11.bin";
    file = fopen(filename,"rb");
    if (file == NULL){
        printf("fail to open ddmaLUT file\n");
        exit(1);
    }
    fread(ddmaLUT,sizeof(double),63000,file);

    fclose(file);
}
*/

void wind_interpolate(windField *wf,struct Geometry geom, struct inputWindField iwf, double grid_resolution){

    //interpolate from iwf.data[] to wf.data[]
    printf("Interpolate wind field into surface frame\n");

    double r_sp; //earth radius at specular point
    double d; //grid resolution
    double dphi, dtheta, phi, theta;  //all in unit of rad
    int numX, numY, numPts, ind;

    //printf ("geom sp = %f %f %f\n",geom.sp_position_ecef_m[0],geom.sp_position_ecef_m[1],geom.sp_position_ecef_m[2]);
    r_sp = vector_norm(geom.sp_position_ecef_m);
    d = grid_resolution;
    dphi = d/r_sp; //in rad
    dtheta = d/r_sp;

    numX = wf->numGridPtsX;
    numY = wf->numGridPtsX;
    numPts = wf->numGridPts;

    double *PUTx, *PUTy, *PUTz; //coordinates of all patches on specular frame
    PUTx = (double *)calloc(numPts,sizeof(double));
    PUTy = (double *)calloc(numPts,sizeof(double));
    PUTz = (double *)calloc(numPts,sizeof(double));


    //First rorate along Y for phi, then along X for theta
    for (int i = 0; i < numY; i++){
        for(int j = 0;j < numX;j++){
            phi = (numY/2-i) * dphi;
            theta= (numX/2-j) * dtheta;
            ind = i * numY + j;
            PUTx[ind] = r_sp * sin(phi);
            PUTy[ind] = -r_sp * cos(phi) * sin(theta);
            PUTz[ind] = r_sp * cos(phi) * cos(theta);
        }
    }

    //SPEC TO ECEF
    double *PUTx_ECEF, *PUTy_ECEF, *PUTz_ECEF;
    PUTx_ECEF = (double *)calloc(numPts,sizeof(double));
    PUTy_ECEF = (double *)calloc(numPts,sizeof(double));
    PUTz_ECEF = (double *)calloc(numPts,sizeof(double));

    double R_ECEFtoSPEC[9],R_SPEC2ECEF[9];

    getECEF2SpecularFrameXfrm(geom.rx_position_ecef_m,geom.tx_position_ecef_m,geom.sp_position_ecef_m,R_ECEFtoSPEC);
    matrix_transpose(3,3,R_ECEFtoSPEC,R_SPEC2ECEF);
    for (int i = 0; i < numPts; i++){
        PUTx_ECEF[i] = PUTx[i] * R_SPEC2ECEF[0] + PUTy[i] * R_SPEC2ECEF[1] + PUTz[i] * R_SPEC2ECEF[2];
        PUTy_ECEF[i] = PUTx[i] * R_SPEC2ECEF[3] + PUTy[i] * R_SPEC2ECEF[4] + PUTz[i] * R_SPEC2ECEF[5];
        PUTz_ECEF[i] = PUTx[i] * R_SPEC2ECEF[6] + PUTy[i] * R_SPEC2ECEF[7] + PUTz[i] * R_SPEC2ECEF[8];
    }

    //ECEF TO LLH
    double *PUT_LAT, *PUT_LON, *PUT_H;
    PUT_LAT = (double *)calloc(numPts,sizeof(double));
    PUT_LON = (double *)calloc(numPts,sizeof(double));
    PUT_H = (double *)calloc(numPts,sizeof(double));

    for (int i = 0; i < numPts; i++){
        double pos_ecef[3] = {PUTx_ECEF[i],PUTy_ECEF[i],PUTz_ECEF[i]};
        double pos_lla[3];
        wgsxyz2lla(pos_ecef,pos_lla);
        PUT_LAT[i] = pos_lla[0];
        PUT_LON[i] = pos_lla[1];
        if (PUT_LON[i] < 0) PUT_LON[i]=PUT_LON[i]+360; //convert from -180 -180 to 0-360
        PUT_H[i] = pos_lla[2];

    }

    int savePUT;
    savePUT = 0;
    if (savePUT==1){
        FILE *outp = fopen("PUT.dat", "wb");
        for (int i = 0;i< numPts;i++) {
            fwrite(&PUT_LAT[i], sizeof(double), 1, outp);
            fwrite(&PUT_LON[i], sizeof(double), 1, outp);
        }
        fclose(outp);
    }

    //Interpolate iwf.data[] to wf.data[]
    int *positions;
    positions = (int *)calloc(numPts,sizeof(int));

    double *lon_vec, *lat_vec;
    lon_vec = (double *)calloc(iwf.numPtsLon,sizeof(double));
    lat_vec = (double *)calloc(iwf.numPtsLat,sizeof(double));

    //PUT_LON is in -180 to 180; iwf.lon_min_deg is in 0-360

    for (int i = 0; i < iwf.numPtsLon; i++){
            lon_vec[i] = iwf.lon_min_deg + i*iwf.resolution_lon_deg;
    }
    for (int i = 0; i < iwf.numPtsLat; i++){
        lat_vec[i] = iwf.lat_min_deg + i*iwf.resolution_lat_deg;
    }

    //bilinear interpolation
    int bi_index[4];
    double bi_weight[4];
    for (int i = 0; i<numPts; i++){
        bilinear_interp(lon_vec, lat_vec, iwf.numPtsLon, iwf.numPtsLat,
                               PUT_LON[i], PUT_LAT[i], bi_index, bi_weight, fabs(iwf.resolution_lat_deg));

        wf->data[i].windSpeed_U10_ms = bi_weight[0]*iwf.data[bi_index[0]].windSpeed_ms
                                       +bi_weight[1]*iwf.data[bi_index[1]].windSpeed_ms
                                       +bi_weight[2]*iwf.data[bi_index[2]].windSpeed_ms
                                       +bi_weight[3]*iwf.data[bi_index[3]].windSpeed_ms;
        wf->data[i].windSpeed_V10_ms = 0;
        wf->data[i].freezingHeight_m = bi_weight[0]*iwf.data[bi_index[0]].freezingHeight_m
                                       +bi_weight[1]*iwf.data[bi_index[1]].freezingHeight_m
                                       +bi_weight[2]*iwf.data[bi_index[2]].freezingHeight_m
                                       +bi_weight[3]*iwf.data[bi_index[3]].freezingHeight_m;
        wf->data[i].rainRate_mmhr = bi_weight[0]*iwf.data[bi_index[0]].rainRate_mmhr
                                    +bi_weight[1]*iwf.data[bi_index[1]].rainRate_mmhr
                                    +bi_weight[2]*iwf.data[bi_index[2]].rainRate_mmhr
                                    +bi_weight[3]*iwf.data[bi_index[3]].rainRate_mmhr;
        for (int j = 0 ; j<4; j++){
            bi_index0[i][j]=bi_index[j]; //store into global variable
            bi_weight0[i][j]=bi_weight[j];
        }
    }

    /* nearest interpoation
    int ind_lat, ind_lon;
    for (int i = 0; i < numPts; i++){
        ind_lat = find_nearest(lat_vec, iwf.numPtsLat, PUT_LAT[i]);
        ind_lon = find_nearest(lon_vec, iwf.numPtsLon, PUT_LON[i]);
        positions[i] = ind_lon * iwf.numPtsLat + ind_lat;
        wf->data[i].windSpeed_U10_ms = iwf.data[positions[i]].windSpeed_U10_ms;
        wf->data[i].windSpeed_V10_ms = iwf.data[positions[i]].windSpeed_V10_ms;
        wf->data[i].freezingHeight_m = iwf.data[positions[i]].freezingHeight_m;
        wf->data[i].rainRate_mmhr = iwf.data[positions[i]].rainRate_mmhr;
    }
     */
    //all checked, correct!

    free(PUTx);free(PUTy);free(PUTz);
    free(PUTx_ECEF);free(PUTy_ECEF);free(PUTz_ECEF);
    free(PUT_LAT);free(PUT_LON);free(PUT_H);
    free(lon_vec);free(lat_vec);free(positions);

    double mss[5];
    double mag_abs, dir_rad;
    for(int i = 0; i < numPts; i++){   //1:8100
        wind_convertWindXY2MagDir( wf->data[i].windSpeed_U10_ms, wf->data[i].windSpeed_V10_ms, &mag_abs, &dir_rad);  //from U10, V10 to wind magnitude and direction
        wf->data[i].windSpeed_ms = mag_abs;
        wf->data[i].windDir_rad  = dir_rad;
        wind_converWindToMSS( wf->data[i].windSpeed_ms, wf->data[i].windDir_rad * R2D, mss );   //from wind speed to surface slope variances and correlation Katzberg Model
        wf->data[i].mss_perp     = mss[0];
        wf->data[i].mss_para     = mss[1];
        wf->data[i].mss_x        = mss[2];
        wf->data[i].mss_y        = mss[3];
        wf->data[i].mss_b        = mss[4];
    }


}

void wind_initialize(windField *wf, struct metadata meta, struct Geometry geom, struct inputWindField iwf){
    wf->numGridPtsX   = meta.numGridPoints[0];	//120
    wf->numGridPtsY   = meta.numGridPoints[1];   //120
    wf->resolutionX_m = meta.grid_resolution_m; //1000
    wf->resolutionY_m = meta.grid_resolution_m;
    wf->numGridPts    = wf->numGridPtsX * wf->numGridPtsY;

    //wf->isLoaded = 1;
    int N = wf->numGridPts;//8100
    wf->data = (windFieldPixel *)calloc(N, sizeof(windFieldPixel)); //initiallize this array with length of N

    //initialize the index
    for(int y = 0; y < wf->numGridPtsY; y++){
        for(int x = 0; x < wf->numGridPtsX; x++) {
            wf->data[y * wf->numGridPtsX + x].x = x;
            wf->data[y * wf->numGridPtsX + x].y = y;
        }
    }

    wf->type = 2;//non-uniform wind
}

void wind_getWindFieldAtXY( windField *wf, double x_m, double y_m, windFieldPixel *value ){  //copy wf to value
    int x_idx, y_idx, idx;

    switch( wf->type ){  //type=2
        case 1: // uniform wind field
            memcpy( value, &(wf->data[wf->locCurrentPt]), sizeof(windFieldPixel) );
            break;

        case 2: // non-uniform wind field
            x_idx = floor(  (x_m - 0) / wf->resolutionX_m );
            y_idx = floor(  (y_m - 0) / wf->resolutionY_m );
            //printf("XX: %d %d\n", x_idx, y_idx );
            //idx   = y_idx * wf->numGridPtsX + x_idx;
            idx = x_idx * wf->numGridPtsY + y_idx;

            //printf("idx= %d \n", idx);

            if( (idx < 0) || (idx >= wf->numGridPts) ){
                //printf("Error: windfield out of range\n");
                //exit(0);
                idx = 0;
            }
            memcpy( value, &(wf->data[idx]), sizeof(windFieldPixel) );
            break;

        default:
            printf("Error: bad wind field type in wind_getWindFieldAtXY\n");
            exit(0);
    }
}

/****************************************************************************/
//  Load binary wind field file (Nature Run)
/****************************************************************************/

void wind_loadWindField(const char *filename, windField *wf ) {

    FILE *file;
    size_t numBytesReadFromFile;
    double wind_Xcomp, wind_Ycomp, mag_abs, dir_rad, mss[5];

    /*
    if( wf->isLoaded == 1 ) {
        free(wf->data);
        wf->isLoaded = 0;
    }*/

    // open file
    fprintf(outputPtr,"   Loading windfield file: %s\n", filename);
    file = fopen(filename,"rb");
    if (!file) {
        printf("Error: Couldn't open file in wind_loadWindField\n");
        return;
    }

    // read windfield dimensions
    double a,b;
    fread(&a,sizeof(double),1,file);
    fread(&b,sizeof(double),1,file);
    wf->numGridPtsX   = floor(a);
    wf->numGridPtsY   = floor(b);
    wf->resolutionX_m = 1000;
    wf->resolutionY_m = 1000;
    wf->numGridPts    = wf->numGridPtsX * wf->numGridPtsY;
    fprintf(outputPtr,"   Headers says its %d x %d (X x Y) \n", wf->numGridPtsX, wf->numGridPtsY);

    // allocate memory
    fprintf(outputPtr,"   Allocating %4.2f MB for windfield data ... \n", sizeof(windFieldPixel) * wf->numGridPts / 1048576.0 );
    wf->data = (windFieldPixel *)calloc(wf->numGridPts, sizeof(windFieldPixel));
    if(wf->data == NULL) {
        fputs("Error: Not enough memory to load file in wind_loadWindField\n",stderr);
        exit(1);
    }

    // read from data file into memory
    int N = wf->numGridPts;
    int M = 8; // num fields
    double *tempBuffer = (double *)calloc(N*M, sizeof(double));
    numBytesReadFromFile = fread(tempBuffer,sizeof(double),N*M,file);
    if (numBytesReadFromFile != (M*N)) {
        printf("Error: Reading error in wind_loadWindField: read %lu, expected %d \n", numBytesReadFromFile, N*M);
        free(wf->data);
        exit(3);
    }
    for(int i=0;i<N;i++){
        // wind and rain field data
        wf->data[i].windSpeed_U10_ms = tempBuffer[0*N+i];
        wf->data[i].windSpeed_V10_ms = tempBuffer[1*N+i];
        wf->data[i].rainRate_mmhr    = tempBuffer[2*N+i];
        wf->data[i].freezingHeight_m = tempBuffer[3*N+i];

        // used for debugging purposes
        wf->data[i].lat_deg          = tempBuffer[4*N+i];
        wf->data[i].lon_deg          = tempBuffer[5*N+i];
        wf->data[i].rowIdx           = tempBuffer[6*N+i];
        wf->data[i].colIdx           = tempBuffer[7*N+i];
    }
    free(tempBuffer);
    fclose(file);
    wf->isLoaded = 1;

    // set x,y surface coordinates (for debugging purposes)
    for(int y = 0; y < wf->numGridPtsY; y++){
        for(int x = 0; x < wf->numGridPtsX; x++) {
            wf->data[y * wf->numGridPtsX + x].x = x;
            wf->data[y * wf->numGridPtsX + x].y = y;
        }
    }

    // get magnitude and direction
    for(int i = 0; i < wf->numGridPts; i++){
        wind_Xcomp = wf->data[i].windSpeed_U10_ms;
        wind_Ycomp = wf->data[i].windSpeed_V10_ms;
        wind_convertWindXY2MagDir( wind_Xcomp, wind_Ycomp, &mag_abs, &dir_rad );
        wf->data[i].windSpeed_ms = mag_abs;
        wf->data[i].windDir_rad  = dir_rad;
    }

    // limit minimum wind
    if( wf->minimumWindSpeed_ms > 0){
        for(int i = 0; i < wf->numGridPts; i++){
            if( wf->data[i].windSpeed_ms < wf->minimumWindSpeed_ms){
                double c = wf->minimumWindSpeed_ms / wf->data[i].windSpeed_ms;
                wf->data[i].windSpeed_ms     *= c;
                wf->data[i].windSpeed_U10_ms *= c;
                wf->data[i].windSpeed_V10_ms *= c;
            }
        }
    }

    // solve for mss
    for(int i = 0; i < wf->numGridPts; i++){
        wind_converWindToMSS( wf->data[i].windSpeed_ms, wf->data[i].windDir_rad * R2D, mss );
        wf->data[i].mss_perp     = mss[0];
        wf->data[i].mss_para     = mss[1];
        wf->data[i].mss_x        = mss[2];
        wf->data[i].mss_y        = mss[3];
        wf->data[i].mss_b        = mss[4];
    }
}

/****************************************************************************/
//  Katzberg Model - Wind Speed to MSS
/****************************************************************************/

void wind_converWindToMSS( double windSpeedMag_ms, double windDirectionAngle_deg, double mss[5] ) {	//wind to MSS (62)-(67)

    double f, sigma2_sx0, sigma2_sy0, sigma2_sx, sigma2_sy, sxsy, b_xy, phi0_rad;

    //*********************************************************************************
    // Katzberg Model

    if(windSpeedMag_ms <= 3.49)
        f = windSpeedMag_ms;
    else
    if( (windSpeedMag_ms > 3.49 ) & (windSpeedMag_ms <= 46) )
        f = 6*log(windSpeedMag_ms) - 4.0;
    else
        f = ((1.855e-4) * windSpeedMag_ms+0.0185) / (3.16e-3) / 0.45;

    sigma2_sx0 = 0.45*(0.003 + (1.92e-3)*f);  // perp component
    sigma2_sy0 = 0.45*((3.16e-3)*f);          // parallel component

    //*********************************************************************************
    // Rotate mss into X-Y coordinate system using equations (44),(45),(46) and (42)
    // from Zavorotny & Voronovich 2000.

    phi0_rad   = windDirectionAngle_deg * D2R;

    sigma2_sx  = sigma2_sx0 * pow(cos(phi0_rad),2) + sigma2_sy0 * pow(sin(phi0_rad),2);
    sigma2_sy  = sigma2_sy0 * pow(cos(phi0_rad),2) + sigma2_sx0 * pow(sin(phi0_rad),2);
    sxsy       = (sigma2_sy0 - sigma2_sx0)*cos(phi0_rad)*sin(phi0_rad);
    b_xy       = sxsy / sqrt( sigma2_sx*sigma2_sy );

    mss[0] = sigma2_sx0; //mss_perp
    mss[1] = sigma2_sy0; //mss_para
    mss[2] = sigma2_sx;  //mss_x
    mss[3] = sigma2_sy;  //mss_y
    mss[4] = b_xy;       //mss_b
}

void wind_convertWindXY2MagDir( double x, double y, double *mag, double *dir_rad ){
    // wind vector w = [ x y ], y-axis vector y = [ 0 1 ];
    // evalutates  angle = sign(w(1)) * acos( dot(w,y) / norm(w) )
    double sign_x = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
    *mag = sqrt(pow(x,2) + pow(y,2));
    *dir_rad = sign_x * acos( y / *mag );
}
