//---------------------------------------------------------------------------
//
// In the E2ES, all of the satellite and specular geometry is stored in a single geometry
// structure (a "orbitGeometryStruct").  The Rx\Tx\Sx positions and velocities in ECEF coordinates
// are loaded from the configuration file using the geom_initialize function, which then
// precalculates a lot of secondary geometric information.  The E2ES eventually sends this
// geometry information to the surface_calcGeom function, which uses it to calculate
// geometry parameters over the entire gridded surface.
//
//***************************************************************************

#include "gnssr.h"

void geom_initialize(geometryData *gd, struct Geometry geom){
    // Read geometry information(ECEF) to gd

    gd->geomStartIdx    = 0;
    gd->geomEndIdx      = 0;
    gd->numGeometries   = 1;

    gd->g = (orbitGeometryStruct *)calloc(gd->numGeometries, sizeof(orbitGeometryStruct));

    for(int i=0;i<gd->numGeometries;i++){
        //printf("read geometry data\n");

        gd->g[i].rx_pos_ecef[0] = geom.rx_position_ecef_m[0];
        gd->g[i].rx_pos_ecef[1] = geom.rx_position_ecef_m[1];
        gd->g[i].rx_pos_ecef[2] = geom.rx_position_ecef_m[2];

        gd->g[i].tx_pos_ecef[0] = geom.tx_position_ecef_m[0];
        gd->g[i].tx_pos_ecef[1] = geom.tx_position_ecef_m[1];
        gd->g[i].tx_pos_ecef[2] = geom.tx_position_ecef_m[2];

        gd->g[i].rx_vel_ecef[0] = geom.rx_velocity_ecef_ms[0];
        gd->g[i].rx_vel_ecef[1] = geom.rx_velocity_ecef_ms[1];
        gd->g[i].rx_vel_ecef[2] = geom.rx_velocity_ecef_ms[2];

        gd->g[i].tx_vel_ecef[0] = geom.tx_velocity_ecef_ms[0];
        gd->g[i].tx_vel_ecef[1] = geom.tx_velocity_ecef_ms[1];
        gd->g[i].tx_vel_ecef[2] = geom.tx_velocity_ecef_ms[2];

        gd->g[i].sx_pos_ecef[0] = geom.sp_position_ecef_m[0];
        gd->g[i].sx_pos_ecef[1] = geom.sp_position_ecef_m[1];
        gd->g[i].sx_pos_ecef[2] = geom.sp_position_ecef_m[2];

        gd->g[i].sx_vel_ecef[0] = geom.sp_velocity_ecef_m[0];
        gd->g[i].sx_vel_ecef[1] = geom.sp_velocity_ecef_m[1];
        gd->g[i].sx_vel_ecef[2] = geom.sp_velocity_ecef_m[2];

        gd->g[i].sc_att_rad[0] = geom.sc_att_rad[0];
        gd->g[i].sc_att_rad[1] = geom.sc_att_rad[1];
        gd->g[i].sc_att_rad[2] = geom.sc_att_rad[2];

        // Convert ECEF to specular frame
        // Compute sx_pos, RX, TX, Dopper, Delay
        geom_calculateSecondaryGeometry(&(gd->g[i]));
    }

}

orbitGeometryStruct *geom_getOrbitData(geometryData *gd, int geomIdx ){
    if( (geomIdx < gd->geomStartIdx) || (geomIdx > gd->geomEndIdx) )
    fatalError("geomIdx out of range, %d not between %d %d", geomIdx, gd->geomStartIdx, gd->geomEndIdx);
    return &(gd->g[geomIdx - gd->geomStartIdx]);
}


/****************************************************************************/
// Given Rx\Ts position & velocity, calculate all other geometry parameters
/****************************************************************************/

void geom_calculateSecondaryGeometry( orbitGeometryStruct *g ){

    // Solve for specular point (now use sp_pos from L1data)
    // SolveSpecularPtPosition(g->rx_pos_ecef, g->tx_pos_ecef, g->sx_pos_ecef, 0, 100);  //sx_pos_ecef---specular position in ECEF

    // Calculate Lat Lon Height (altitude above WGS84 geoid)
    wgsxyz2lla( g->rx_pos_ecef, g->rx_pos_llh );
    wgsxyz2lla( g->tx_pos_ecef, g->tx_pos_llh );
    wgsxyz2lla( g->sx_pos_ecef, g->sx_pos_llh );

    // Propagate satellites forward in time to solve for specular point velocity
    double tempRx[3],tempTx[3],tempSx[3],tempSxLLH[3]; // position at next time
    double timeInc_s = 1;
    for(int i=0;i<3;i++){ //Tx and Rx pos at next time
        tempRx[i] = g->rx_pos_ecef[i] + g->rx_vel_ecef[i] * timeInc_s;
        tempTx[i] = g->tx_pos_ecef[i] + g->tx_vel_ecef[i] * timeInc_s;
        tempSx[i] = g->sx_pos_ecef[i] + g->sx_vel_ecef[i] * timeInc_s;
    }
    //solveSpecularPtPosition(tempRx, tempTx, tempSx, 0, 100); // it may not be very accurate; sometimes fail angle check

    wgsxyz2lla( tempSx, tempSxLLH );
    //for(int i=0;i<3;i++) g->sx_vel_ecef[i] = (tempSx[i]    - g->sx_pos_ecef[i]) / timeInc_s;
    for(int i=0;i<3;i++) g->sx_vel_llh[i]  = (tempSxLLH[i] - g->sx_pos_llh[i] ) / timeInc_s;

    // Convert ECEF to specular frame
    getECEF2SpecularFrameXfrm(g->rx_pos_ecef, g->tx_pos_ecef, g->sx_pos_ecef, g->ECEF_TO_SPEC_FRAME );  // get Transformation matrix
    matrixVectorMult3x3(g->ECEF_TO_SPEC_FRAME, g->rx_pos_ecef, g->rx_pos );
    matrixVectorMult3x3(g->ECEF_TO_SPEC_FRAME, g->tx_pos_ecef, g->tx_pos );
    matrixVectorMult3x3(g->ECEF_TO_SPEC_FRAME, g->sx_pos_ecef, g->sx_pos );  // sx_pos: specular position at specular frame
    matrixVectorMult3x3(g->ECEF_TO_SPEC_FRAME, g->rx_vel_ecef, g->rx_vel );
    matrixVectorMult3x3(g->ECEF_TO_SPEC_FRAME, g->tx_vel_ecef, g->tx_vel );
    matrixVectorMult3x3(g->ECEF_TO_SPEC_FRAME, g->sx_vel_ecef, g->sx_vel );

    // Specular-to-ECEF frame
    matrix_invert_3x3(g->ECEF_TO_SPEC_FRAME, g->SPEC_TO_ECEF_FRAME);

    // Calculate Rx\Tx\Sx speeds
    g->rx_speed_ms = vector_norm(g->rx_vel);
    g->tx_speed_ms = vector_norm(g->tx_vel);
    g->sx_speed_ms = vector_norm(g->sx_vel);

    // Calculate Rx\Tx range and Doppler contributions in Specular Frame
    double vector_sx_to_rx[9],vector_sx_to_tx[9],RSx_unit[3],TSx_unit[3];
    vector_subtract(g->rx_pos, g->sx_pos, vector_sx_to_rx);  // vector from Sx to Rx
    vector_subtract(g->tx_pos, g->sx_pos, vector_sx_to_tx);  // vector from Sx to Tx
    vector_unit( vector_sx_to_rx, RSx_unit);  // unit vector from Sx to Rx
    vector_unit( vector_sx_to_tx, TSx_unit);  // unit vector from Sx to Tx
    g->rx_range_m         = vector_norm(vector_sx_to_rx);
    g->tx_range_m         = vector_norm(vector_sx_to_tx);
    g->rx_LOS_velocity    = vector_dot_product(g->rx_vel,RSx_unit);
    g->tx_LOS_velocity    = vector_dot_product(g->tx_vel,TSx_unit);

    // Doppler contributions of Rx and Tx to reflected signal's Doppler
    g->dopplerRx_Hz  = -1 * g->rx_LOS_velocity * (L1)/(speedlight);
    g->dopplerTx_Hz  = -1 * g->tx_LOS_velocity * (L1)/(speedlight);

    // Calculate total path delay and Doppler at specular point
    g->specularDistance_m = g->rx_range_m   + g->tx_range_m;
    g->specularDoppler_Hz = g->dopplerTx_Hz + g->dopplerRx_Hz;

    // Calculate elevation angles (hopefully these should all be the same...)
    g->rx_angle_rad = asin(RSx_unit[2]);
    g->tx_angle_rad = asin(TSx_unit[2]);

    // Incidence angle at specular point
    g->sx_angle_rad = acos(vector_dot_product(TSx_unit,RSx_unit))/2;

    // Calculate Rx and Tx orbit frames transforms (as defined by doc148-0336-1)
    // Updated by Feixiong
    // ECEF to Orbit
    getECEF2OrbitFrameXfrm( g->rx_pos_ecef, g->rx_vel_ecef, g->ECEF_TO_RX_ORB_FRAME );  //receiver
    getECEF2OrbitFrameXfrm( g->tx_pos_ecef, g->tx_vel_ecef, g->ECEF_TO_TX_ORB_FRAME );  //transmitter

    // SPEC to Orbit
    matrix_multiply(3,3,3,g->ECEF_TO_RX_ORB_FRAME,g->SPEC_TO_ECEF_FRAME,g->SPEC_TO_RX_ORB_FRAME);
    matrix_multiply(3,3,3,g->ECEF_TO_TX_ORB_FRAME,g->SPEC_TO_ECEF_FRAME,g->SPEC_TO_TX_ORB_FRAME);

    // Orbit to Body (only for RX)
    // Fot Tx, Body = Orbit
    getOrbit2BodyFrameXfrm (g->sc_att_rad[0],g->sc_att_rad[1],g->sc_att_rad[2],g->ORB_TO_RX_BODY_FRAME ); //by Feixiong

    // SPEC to Body (only for RX)
    matrix_multiply(3,3,3,g->ORB_TO_RX_BODY_FRAME,g->SPEC_TO_RX_ORB_FRAME,g->SPEC_TO_RX_BODY_FRAME);

    // Get distance, Doppler, and relative angles between Rx and Tx
    double vector_tx_to_rx[9],vector_rx_to_tx[9],RTx_unit[3],TRx_unit[3];
    vector_subtract(g->rx_pos, g->tx_pos, vector_tx_to_rx);  // vector from Tx to Rx
    vector_subtract(g->tx_pos, g->rx_pos, vector_rx_to_tx);  // vector from Rx to Tx
    vector_unit( vector_tx_to_rx, TRx_unit);  // unit vector from Tx to Rx
    vector_unit( vector_rx_to_tx, RTx_unit);  // unit vector from Rx to Tx
    g->directPathDistance_m = vector_norm(vector_tx_to_rx);
    g->directPathDoppler_Hz = +1 * ( vector_dot_product(g->rx_vel,RTx_unit) +
                                     vector_dot_product(g->tx_vel,TRx_unit) ) * (L1)/(speedlight);
    g->pathDifference_chips = (g->specularDistance_m - g->directPathDistance_m) / speedlight * chipRate_cs;
    geom_getRelativeAngleInFrame(g->rx_pos, g->tx_pos, g->SPEC_TO_RX_ORB_FRAME, g->angleTxFromRx_rad );
    geom_getRelativeAngleInFrame(g->tx_pos, g->rx_pos, g->SPEC_TO_TX_ORB_FRAME, g->angleRxFromTx_rad );

    // Get relative angle to specular point as seen from Rx or Tx
    geom_getRelativeAngleInFrame(g->rx_pos, g->sx_pos, g->SPEC_TO_RX_BODY_FRAME, g->angleSxFromRx_rad ); //Body frame
    geom_getRelativeAngleInFrame(g->tx_pos, g->sx_pos, g->SPEC_TO_TX_ORB_FRAME, g->angleSxFromTx_rad );

    // Rx & Tx antenna gains for both polarizations (not all are used currently)
    // Specular point gain of Tx and Rx
    g->antennaRxGainAtSx_RHCP_dB = 10*log10(antenna_getGain_abs( CYGNSS_NADIR_ANT,  RHCP, g->angleSxFromRx_rad ) );
    g->antennaRxGainAtSx_LHCP_dB = 10*log10(antenna_getGain_abs( CYGNSS_NADIR_ANT,  LHCP, g->angleSxFromRx_rad ) );
    g->antennaTxGainAtSx_RHCP_dB = 10*log10(antenna_getGain_abs( GPS_SAT_ANT,       RHCP, g->angleSxFromTx_rad ) );
    g->antennaTxGainAtSx_LHCP_dB = 10*log10(antenna_getGain_abs( GPS_SAT_ANT,       LHCP, g->angleSxFromTx_rad ) );
    g->antennaRxGainAtTx_RHCP_dB = 10*log10(antenna_getGain_abs( CYGNSS_ZENITH_ANT, RHCP, g->angleTxFromRx_rad ) );
    g->antennaRxGainAtTx_LHCP_dB = 10*log10(antenna_getGain_abs( CYGNSS_ZENITH_ANT, LHCP, g->angleTxFromRx_rad ) );
    g->antennaTxGainAtRx_RHCP_dB = 10*log10(antenna_getGain_abs( GPS_SAT_ANT,       RHCP, g->angleRxFromTx_rad ) );
    g->antennaTxGainAtRx_LHCP_dB = 10*log10(antenna_getGain_abs( GPS_SAT_ANT,       LHCP, g->angleRxFromTx_rad ) );

    printf("Calculated Rx gain = %f\n",g->antennaRxGainAtSx_LHCP_dB);

    // Path loss (UNDER CONSTRUCTION, not used now)
    g->pathloss_Rx2Sx_dB = 20*log10( g->rx_range_m ) + 20*log10( L1 ) - 147.55;
    g->pathloss_Tx2Sx_dB = 20*log10( g->tx_range_m ) + 20*log10( L1 ) - 147.55;
    g->pathloss_Rx2Tx_dB = 20*log10( g->directPathDistance_m ) + 20*log10( L1 ) - 147.55;
    g->pathloss_Sx_dB    = 20*log10( g->specularDistance_m ) + 20*log10( L1 ) - 147.55;

    // Calculate motion of specular bin in DDM
    vector_subtract(tempRx, tempSx, vector_sx_to_rx);  // vector from Sx to Rx
    vector_subtract(tempTx, tempSx, vector_sx_to_tx);  // vector from Sx to Tx
    vector_unit( vector_sx_to_rx, RSx_unit);  // unit vector from Sx to Rx
    vector_unit( vector_sx_to_tx, TSx_unit);  // unit vector from Sx to Tx
    double deltaDoppler_Hz = (vector_dot_product(g->rx_vel_ecef,RSx_unit) +
                              vector_dot_product(g->tx_vel_ecef,TSx_unit)) * -1 * (L1)/(speedlight)
                             - g->specularDoppler_Hz;
    double deltaDelay_s    = ( (vector_norm(vector_sx_to_rx) + vector_norm(vector_sx_to_tx))
                               - g->specularDistance_m ) / speedlight ;
    g->ddm_delayRate_chips_per_s = deltaDelay_s * 1.023e6 / timeInc_s;
    g->ddm_dopplerRate_kH_per_s  = deltaDoppler_Hz / 1000 / timeInc_s;

}

void geom_getRelativeAngleInFrame( double origin[3], double pos[3] , double M[9], double angles_rad[2] ){
    // Used primarily for finding the angle to the specular pt from the perspective of the satellites
    // origin, pos are in specular frame
    // M is the rotation matrix from SPEC to orbit

    double temp0[3], temp[3], tempSph[3];
    vector_subtract(pos, origin, temp0);   // temp0 = pos-origin vector from origin to pos (coordinates in SPEC frame)
    matrixVectorMult3x3( M, temp0, temp ); // temp = M*temp0 put it in orbit frame (coordinates in orbit frame)
    cart2sph(temp,tempSph);
    angles_rad[0] = tempSph[1];  // phi (azimuth)
    angles_rad[1] = tempSph[2];  // theta (elevation)
}


/****************************************************************************/
//  Transforms between coordinate frames
/****************************************************************************/

void getECEF2SpecularFrameXfrm( double rx_pos_ecef[3], double tx_pos_ecef[3],
                                double sx_pos_ecef[3], double xfrmMatrix[9] ){
    // This is Scott's specular frame definition
    // z-hat points normal to surface of Earth at specular point
    // x-hat points from Tx to Rx (the component orthog to z)
    // y-hat follows from right hand system

    double tempx[3],tempy[3],tempz[3];
    vector_unit(sx_pos_ecef, tempz);
    vector_subtract(rx_pos_ecef, tx_pos_ecef, tempx);
    vector_orthoNorm(tempz, tempx);
    vector_cross_product(tempz, tempx, tempy);
    matrix_form3x3(tempx, tempy, tempz, xfrmMatrix); //form a matrix by combining three row vectors
}

void getECEF2OrbitFrameXfrm( double sat_pos_ecef[3], double sat_vel_ecef[3],
                                           double xfrmMatrix[9] ){
    // Matrix for column vector [x1;y1;z1]=M*[x0;y0;z0];
    // This transform defines the orbit frame of a satellite (Rx or Tx)
    // Writen by Feixiong; accounting for Earth rotation according to the TDS1 doc
    // y-hat is negative orbit normal
    // z-hat is towards to the center of earth
    // x-hat is along-track by right hand system

    double temp[3],tempx[3],tempy[3],tempz[3];
    double sat_vel_inertial[3];
    double we[3]={0,0,omega0};  //earth rotation vector

    vector_cross_product(we,sat_pos_ecef,temp);  // temp = we x sat_pos
    vector_add(sat_vel_ecef,temp,sat_vel_inertial);  // sat_vel_inertial = sat_vel+temp

    vector_cross_product(sat_pos_ecef,sat_vel_inertial,temp);  // temp = sat_pos x sat_vel_inertial
    vector_scale(temp,tempy,-1/vector_norm(temp));  // tempy = temp * -1/vector_norm(temp)

    vector_scale(sat_pos_ecef,tempz,-1/vector_norm(sat_pos_ecef));  // tempz
    vector_cross_product(tempy, tempz, tempx);

    matrix_form3x3(tempx, tempy, tempz, xfrmMatrix);
}

void getOrbit2BodyFrameXfrm (double pitch, double roll, double yaw,double xfrmMatrix[9] ){
    // Added by Feixiong
    // from Orbit frame to Body frame using attitude dynamics

    double theta = pitch;  // pitch-Y
    double phi = roll;  // roll-X
    double psi = yaw;  // yaw-Z

    double R_X[9],R_Y[9],R_Z[9],R_ZX[9];
    RotationMatrix_X (phi, R_X);
    RotationMatrix_Y (theta, R_Y);
    RotationMatrix_Z (psi, R_Z);
    matrix_multiply(3,3,3,R_Z,R_X,R_ZX);  // 2-1-3 sequence
    matrix_multiply(3,3,3,R_ZX,R_Y,xfrmMatrix);
}

/****************************************************************************/
//  Write geometry summary table to a text file
/****************************************************************************/

void geom_printToLog(FILE *outputPtr, int index, orbitGeometryStruct *g){

    fprintf(outputPtr,"\n  *******************\n");
    fprintf(outputPtr,"  *  Geometry #%d   *\n", index);
    fprintf(outputPtr,"  *******************\n\n");
    fprintf(outputPtr,"   Rx Position (ECEF Frame)            (%f,%f,%f) (m)\n",  g->rx_pos_ecef[0],g->rx_pos_ecef[1],g->rx_pos_ecef[2]);
    fprintf(outputPtr,"   Tx Position (ECEF Frame)            (%f,%f,%f) (m)\n",  g->tx_pos_ecef[0],g->tx_pos_ecef[1],g->tx_pos_ecef[2]);
    fprintf(outputPtr,"   Sx Position (ECEF Frame)            (%f,%f,%f) (m)\n",  g->sx_pos_ecef[0],g->sx_pos_ecef[1],g->sx_pos_ecef[2]);
    fprintf(outputPtr,"   Rx Velocity (ECEF Frame)            (%f,%f,%f) (m/s)\n",g->rx_vel_ecef[0],g->rx_vel_ecef[1],g->rx_vel_ecef[2]);
    fprintf(outputPtr,"   Tx Velocity (ECEF Frame)            (%f,%f,%f) (m/s)\n",g->tx_vel_ecef[0],g->tx_vel_ecef[1],g->tx_vel_ecef[2]);
    fprintf(outputPtr,"   Sx Velocity (ECEF Frame)            (%f,%f,%f) (m/s)\n\n",g->sx_vel_ecef[0],g->sx_vel_ecef[1],g->sx_vel_ecef[2]);
    fprintf(outputPtr,"   Rx Position (Spec Frame)            (%f,%f,%f) (m)\n",  g->rx_pos[0],g->rx_pos[1],g->rx_pos[2]);
    fprintf(outputPtr,"   Tx Position (Spec Frame)            (%f,%f,%f) (m)\n",  g->tx_pos[0],g->tx_pos[1],g->tx_pos[2]);
    fprintf(outputPtr,"   Sx Position (Spec Frame)            (%f,%f,%f) (m)\n",  g->sx_pos[0],g->sx_pos[1],g->sx_pos[2]);
    fprintf(outputPtr,"   Rx Velocity (Spec Frame)            (%f,%f,%f) (m/s)\n",g->rx_vel[0],g->rx_vel[1],g->rx_vel[2]);
    fprintf(outputPtr,"   Tx Velocity (Spec Frame)            (%f,%f,%f) (m/s)\n",g->tx_vel[0],g->tx_vel[1],g->tx_vel[2]);
    fprintf(outputPtr,"   Sx Velocity (Spec Frame)            (%f,%f,%f) (m/s)\n\n",g->sx_vel[0],g->sx_vel[1],g->sx_vel[2]);
    fprintf(outputPtr,"   Rx Position (Lat,Lon,Alt)           (%f,%f,%f) (deg,deg,km)\n",g->rx_pos_llh[0],g->rx_pos_llh[1],g->rx_pos_llh[2]/1000);
    fprintf(outputPtr,"   Tx Position (Lat,Lon,Alt)           (%f,%f,%f) (deg,deg,km)\n",g->tx_pos_llh[0],g->tx_pos_llh[1],g->tx_pos_llh[2]/1000);
    fprintf(outputPtr,"   Sx Position (Lat,Lon,Alt)           (%f,%f,%f) (deg,deg,km)\n",g->sx_pos_llh[0],g->sx_pos_llh[1],g->sx_pos_llh[2]/1000);
    fprintf(outputPtr,"   Sx Velocity (Lat,Lon,Alt)           (%f,%f,%f) (deg/s,deg/s,km/s)\n\n",g->sx_vel_llh[0],g->sx_vel_llh[1],g->sx_vel_llh[2]/1000);
    fprintf(outputPtr,"   Distance from Sx to Rx              %f (km)\n",  g->rx_range_m / 1000);
    fprintf(outputPtr,"   Distance from Sx to Tx              %f (km)\n",  g->tx_range_m / 1000);
    fprintf(outputPtr,"   Distance along complete path        %f (km)\n\n",  g->specularDistance_m / 1000);
    fprintf(outputPtr,"   Rx Speed                            %f (km/s)\n",g->rx_speed_ms / 1000);
    fprintf(outputPtr,"   Tx Speed                            %f (km/s)\n",g->tx_speed_ms / 1000);
    fprintf(outputPtr,"   Sx Speed                            %f (km/s)\n\n",g->sx_speed_ms / 1000);
    fprintf(outputPtr,"   Rx LOS Vel.                         %f (m/s)\n", g->rx_LOS_velocity );
    fprintf(outputPtr,"   Tx LOS Vel.                         %f (m/s)\n", g->tx_LOS_velocity );
    fprintf(outputPtr,"   Doppler of Rx at Sx                 %f (kHz)\n", g->dopplerRx_Hz / 1000);
    fprintf(outputPtr,"   Doppler of Tx at Sx                 %f (kHz)\n", g->dopplerTx_Hz / 1000);
    fprintf(outputPtr,"   Doppler of reflected signal         %f (kHz)\n\n", g->specularDoppler_Hz / 1000);
    fprintf(outputPtr,"   Elevation Angle of Rx               %f (deg)\n", g->rx_angle_rad * R2D);
    fprintf(outputPtr,"   Elevation Angle of Tx               %f (deg)\n", g->tx_angle_rad * R2D);
    fprintf(outputPtr,"   Specular Angle (incidencen)          %f (deg)\n", g->sx_angle_rad * R2D);
    fprintf(outputPtr,"   Specular Angle (elevation)          %f (deg)\n\n", 90 - g->sx_angle_rad * R2D);
    fprintf(outputPtr,"   Angle to Sx from Rx Orb Frame (az)  %f (deg)\n",  g->angleSxFromRx_rad[0] * R2D);
    fprintf(outputPtr,"   Angle to Sx from Rx Orb Frame (el)  %f (deg)\n\n",g->angleSxFromRx_rad[1] * R2D);
    fprintf(outputPtr,"   Angle to Sx from Tx Orb Frame (az)  %f (deg)\n",g->angleSxFromTx_rad[0] * R2D);
    fprintf(outputPtr,"   Angle to Sx from Tx Orb Frame (el)  %f (deg)\n",g->angleSxFromTx_rad[1] * R2D);
    fprintf(outputPtr,"   Angle to Sx from Tx Orb Frame (asc) %f (deg)\n\n",90 - g->angleSxFromTx_rad[1] * R2D);
    fprintf(outputPtr,"   Direct Path Distance                %f (km)\n",g->directPathDistance_m / 1000);
    fprintf(outputPtr,"   Direct Path Doppler                 %f (kHz)\n",g->directPathDoppler_Hz / 1000);
    fprintf(outputPtr,"   Direct/Reflected Path Delay Diff    %f, %f (chips, chips mod 1023)\n",g->pathDifference_chips, fmod(g->pathDifference_chips,1023));
    fprintf(outputPtr,"   Direct/Reflected Path Doppler Diff  %f (kHz)\n\n",(g->specularDoppler_Hz - g->directPathDoppler_Hz)/1000);
    fprintf(outputPtr,"   Angle to Tx from Rx Orb Frame (az)  %f (deg)\n",g->angleTxFromRx_rad[0] * R2D);
    fprintf(outputPtr,"   Angle to Tx from Rx Orb Frame (el)  %f (deg)\n",g->angleTxFromRx_rad[1] * R2D);
    fprintf(outputPtr,"   Angle to Rx from Tx Orb Frame (az)  %f (deg)\n",g->angleRxFromTx_rad[0] * R2D);
    fprintf(outputPtr,"   Angle to Rx from Tx Orb Frame (el)  %f (deg)\n\n",g->angleRxFromTx_rad[1] * R2D);
    fprintf(outputPtr,"   Rx Antenna Gain at Sx (LHCP)        %f (dB)\n", g->antennaRxGainAtSx_LHCP_dB );
    fprintf(outputPtr,"   Rx Antenna Gain at Sx (RHCP)        %f (dB)\n", g->antennaRxGainAtSx_RHCP_dB );
    fprintf(outputPtr,"   Tx Antenna Gain at Sx (RHCP)        %f (dB)\n", g->antennaTxGainAtSx_RHCP_dB );
    fprintf(outputPtr,"   Rx Antenna Gain at Tx (RHCP)        %f (dB)\n",   g->antennaRxGainAtTx_RHCP_dB );
    fprintf(outputPtr,"   Tx Antenna Gain at Rx (RHCP)        %f (dB)\n\n", g->antennaTxGainAtRx_RHCP_dB );
    fprintf(outputPtr,"   Specular Bin Motion in DDM          (%f,%f) (chips/s,kHz/s)\n\n", g->ddm_delayRate_chips_per_s, g->ddm_dopplerRate_kH_per_s );
    fflush(outputPtr);
}
