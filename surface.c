//---------------------------------------------------------------------------
//
//  This portion of the E2ES is responsible for creating a gridded surface
//  and evaluating a number of parameters over that surface
//
//****************************************************************************/

#include "forwardmodel.h"
#include "gnssr.h"

double dmdx(double ws);

void getSurfaceFieldDataName( int type, char *filename );
void getSurfaceFieldData( int type, double *vals, char *filename, double minmax[2]);


double dmdx (double ws){
    //compute derivative of MSS respect to wind speed in Katzberg model
    if(ws>=0 && ws<=3.49) {return 1.143e-3;}
    else if(ws>3.49 && ws<=46) {return 6.858e-3/ws;}
    else if(ws>46) {return 4.69773e-4;}
    else {
        printf("Negative wind speed\n"); exit(0);
    }
}


void surface_initialize(struct metadata meta){
    surface.numGridPtsX = meta.numGridPoints[0];
    surface.numGridPtsY = meta.numGridPoints[1];
    surface.resolution_m = meta.grid_resolution_m;
    surface.surfaceCurvatureType = meta.surfaceCurvatureType;
    surface.rainOnOff = 0; //turn rain off
    surface.width_m = surface.resolution_m * surface.numGridPtsX;
    surface.height_m = surface.resolution_m * surface.numGridPtsY;
    surface.numGridPts = surface.numGridPtsX  * surface.numGridPtsY;
    surface.specularLoactionX_m = floor( surface.numGridPtsX / 2.0 ) * surface.resolution_m;
    surface.specularLoactionY_m = floor( surface.numGridPtsY / 2.0 ) * surface.resolution_m;

    // allocate surface buffers
    surface.data     = (surfacePixel *) calloc( surface.numGridPts, sizeof(surfacePixel) );
    surface.windData = (windFieldPixel *) calloc( surface.numGridPts, sizeof(windFieldPixel) );
    surface_resetToZero();  //initialize surface.Data

}

void surface_cleanup(void){
    free(surface.data);
}

void surface_calcGeomOverSurface(orbitGeometryStruct *geometry, int surfType, struct powerParm pp){   //use surfType=0    3.1.1   use geometry data (ECEF)

    double temp_vec[9];
    double PUT[3],RSx_unit[3],TSx_unit[3], q_vec[3], n_vec[3];
    double R1,R2,RxG,TxG,TxP,lambda,Ti, Area_dS,powerFactor,extra,normal,pathloss;
    double dopplerRx_Hz,dopplerTx_Hz,doppler_Hz,sxangle_rad;
    double angleSxFromRx_rad[2], angleSxFromTx_rad[2];

    double *rx_pos = geometry->rx_pos;	//specular frame
    double *tx_pos = geometry->tx_pos;
    double *rx_vel = geometry->rx_vel;
    double *tx_vel = geometry->tx_vel;
    double *sx_pos = geometry->sx_pos;  //specular position at specular frame

    double TXP_DB = pp.Tx_Power_dB;
    double ATTEN_DB = pp.AtmosphericLoss_dB;

    surface.specularGridPt_x_idx = (int)floor(surface.numGridPtsX / 2.0 );	//N_theta/2  (48)
    surface.specularGridPt_y_idx = (int)floor(surface.numGridPtsY / 2.0 );	//N_phi/2	(49)

    // we have a quick mode that only does four corners
    int i_inc, j_inc, i0, i1, j0, j1;
    switch(surfType){  //surfType=0
        case 1:
            i0    = 0;
            j0    = 0;
            i1    = surface.numGridPtsX;
            j1    = surface.numGridPtsY;
            i_inc = surface.numGridPtsX - 1;
            j_inc = surface.numGridPtsY - 1;
            break;

        case 2:
            i0    = surface.specularGridPt_x_idx;
            j0    = surface.specularGridPt_y_idx;
            i1    = i0 + 1;
            j1    = j0 + 1;
            i_inc = 1;
            j_inc = 1;
            break;

        default:  //case 0
            i0    = 0;
            j0    = 0;
            i1    = surface.numGridPtsX;
            j1    = surface.numGridPtsY;
            i_inc = 1;
            j_inc = 1;
            break;
    }

    for (int i = i0; i < i1; i+=i_inc) {
        for (int j = i0; j < j1; j+=j_inc) {

            //printf("i = %d, j = %d \n",i,j);

            // get the coordinates of the surface grid pt, PUT
            switch (surface.surfaceCurvatureType) { //Type=1=spherical
                case 1: grid_getGridPt_sphericalEarth(i, j, sx_pos, PUT, vector_norm(geometry->sx_pos)); break; //sx_pos(specular position), PUT(specular frame coordinates of the patch i,j; norm(earth radius))
                case 2: grid_getGridPt_flatEarth(i, j, sx_pos, PUT); break;
                default: fprintf(errPtr,"Error: Bad surfaceCurvatureType in surface_calcGeom"); exit(0);
            }

            // correction for coordinate system (should technically correct direction too)
            //PUT[0] = -1 * PUT[0];

            vector_unit(PUT, n_vec);

            vector_subtract(rx_pos,PUT,temp_vec);
            R2 = vector_norm(temp_vec);         // distance between Rx and PUT
            vector_unit(temp_vec, RSx_unit);    // scattering point to receiver unit vector

            vector_subtract(tx_pos,PUT,temp_vec);
            R1 = vector_norm(temp_vec);         // distance between Tx and PUT
            vector_unit(temp_vec, TSx_unit);    // scattering point to transmitter unit vector

            // Doppler components
            dopplerRx_Hz = -1*vector_dot_product(rx_vel,RSx_unit)*(L1)/(speedlight);
            dopplerTx_Hz = -1*vector_dot_product(tx_vel,TSx_unit)*(L1)/(speedlight);
            doppler_Hz   = dopplerTx_Hz + dopplerRx_Hz;

            // scattering vector and incidence angle
            surface_getScatteringVector(TSx_unit, RSx_unit, PUT, q_vec);// get q_vec
            sxangle_rad = acos(vector_dot_product(TSx_unit,RSx_unit))/2; //incidence angle, corrected by Feixiong

            int useOldAntennaAngles = 0;

            if( useOldAntennaAngles ){
                // get relative angle to grid point as seen from Rx or Tx
                geom_getRelativeAngleInFrame(geometry->rx_pos, PUT, geometry->SPEC_TO_RX_ORB_FRAME_AARON, angleSxFromRx_rad );
                geom_getRelativeAngleInFrame(geometry->tx_pos, PUT, geometry->SPEC_TO_TX_ORB_FRAME_AARON, angleSxFromTx_rad );

                // apply correction to Aaron's angles
                angleSxFromRx_rad[0] -= 180 * D2R;
                if( angleSxFromRx_rad[0] < -1*pi )
                    angleSxFromRx_rad[0] += 360*D2R;
                if( angleSxFromRx_rad[0] > pi )
                    angleSxFromRx_rad[0] -= 360*D2R;
            }else{
                // get relative angle to grid point as seen from Rx or Tx
                geom_getRelativeAngleInFrame(geometry->rx_pos, PUT, geometry->SPEC_TO_RX_ORB_FRAME, angleSxFromRx_rad );
                geom_getRelativeAngleInFrame(geometry->tx_pos, PUT, geometry->SPEC_TO_TX_ORB_FRAME, angleSxFromTx_rad );
            }


            // get Rx & Tx antenna gains for specific angles
            RxG = antenna_getGain_abs(CYGNSS_NADIR_ANT,   LHCP, angleSxFromRx_rad );
            TxG = antenna_getGain_abs(GPS_SAT_ANT,        RHCP, angleSxFromTx_rad );

            // get relative angle to grid point as seen from Rx or Tx (angles for rain ...)
            // TODO: don't reuse variables
            geom_getRelativeAngleInFrame(geometry->rx_pos, PUT, geometry->SPEC_TO_RX_ORB_FRAME, angleSxFromRx_rad );
            geom_getRelativeAngleInFrame(geometry->tx_pos, PUT, geometry->SPEC_TO_TX_ORB_FRAME, angleSxFromTx_rad );

            // calculate the geometry-dependent, windfield-independent power factor for
            // the scattering equation.
            Area_dS  = pow(surface.resolution_m,2);
            normal   = 1 / n_vec[2]; // Account for change in surface area due to Earth curvature
            TxP      = pow(10,(TXP_DB/10)); // Tx Power
            lambda   = L1_WAVELENGTH;
            Ti       = ddm.cohIntegrationTime_s;
            extra    = pow(10,((-ATTEN_DB)/10));
            pathloss = 1 / ( pow(R1,2) * pow(R2,2) );
            // removed pow(Ti,2) factor
            powerFactor    = TxP * pow(lambda,2) * TxG * RxG * Area_dS * normal * extra * pathloss / pow(4*pi,3); //PowerFactor

            // powerfactor associated with 25km footprint
            if( surfType == 2 )
                powerFactor = TxP * pow(lambda,2) * TxG * RxG * (pow((30e3) / 2,2)*pi) * normal * extra * pathloss / pow(4*pi,3);

            // save surface data of eahc grid to buffer
            int idx = SURFINDEX(i, j);   //=i * surface.numGridPtsY + j   i*90+j
            surface.data[idx].delay_s       = ((R1+R2)- geometry->specularDistance_m) / speedlight;
            //surface.data[idx].doppler_Hz    = geometry->specularDoppler_Hz - doppler_Hz;
            surface.data[idx].doppler_Hz    =  doppler_Hz - geometry->specularDoppler_Hz;
            surface.data[idx].q[0]          = q_vec[0];
            surface.data[idx].q[1]          = q_vec[1];
            surface.data[idx].q[2]          = q_vec[2];
            surface.data[idx].sx_angle_rad  = sxangle_rad;
            surface.data[idx].powerFactor   = powerFactor;

            // to solve for rain atten. we'll need the elevation angles across surface
            surface.data[idx].rx_elevationAngle_rad = geometry->rx_angle_rad;// asin(RSx_unit[2]); //geometry->rx_angle_rad;
            surface.data[idx].tx_elevationAngle_rad = geometry->tx_angle_rad;// asin(TSx_unit[2]); //geometry->tx_angle_rad;

            // these parameters aren't used except for debugging purposes
            surface.data[idx].i             = i;
            surface.data[idx].j             = j;
            surface.data[idx].position[0]   = PUT[0];
            surface.data[idx].position[1]   = PUT[1];
            surface.data[idx].position[2]   = PUT[2];
            surface.data[idx].antennaGainRx_abs = RxG;
            surface.data[idx].antennaGainTx_abs = TxG;
            surface.data[idx].pathloss      = pathloss;

            // antenna pattern angles over surface
            surface.data[idx].angleSxFromRx_theta_rad = angleSxFromRx_rad[0];
            surface.data[idx].angleSxFromRx_phi_rad   = angleSxFromRx_rad[1];
            surface.data[idx].angleSxFromTx_theta_rad = angleSxFromTx_rad[0];
            surface.data[idx].angleSxFromTx_phi_rad   = angleSxFromTx_rad[1];

            // the x,y location of each surface grid for purposes of
            // looking up the wind field values
            surface.data[idx].windFieldLocation_x_m      = (i - surface.specularGridPt_x_idx) * surface.resolution_m;	//discrete set of angles (48)  (i-45)*1000
            surface.data[idx].windFieldLocation_y_m      = (j - surface.specularGridPt_y_idx) * surface.resolution_m;	//(49)
            surface.data[idx].windFieldLocationRange_km  = sqrt( pow(surface.data[idx].windFieldLocation_x_m,2) +
                                                                 pow(surface.data[idx].windFieldLocation_y_m,2) ) / 1000;

            // convert location from specular frame to ECEF and LLH
            double pos_ecef[3], pos_spec[3], pos_llh[3];
            memcpy( pos_spec, surface.data[idx].position, 3 * sizeof(double) );
            matrixVectorMult3x3(geometry->SPEC_TO_ECEF_FRAME, pos_spec, pos_ecef );
            wgsxyz2lla( pos_ecef, pos_llh );
            memcpy( surface.data[idx].pos_spec, pos_spec, 3 * sizeof(double) );
            memcpy( surface.data[idx].pos_ecef, pos_ecef, 3 * sizeof(double) );
            memcpy( surface.data[idx].pos_llh,  pos_llh , 3 * sizeof(double) );
        }
    }

    if( surfType == 0 ){

        // determine the valid DDM region for this surface size and geometry
        surface_createSurfaceMask();

        // speckle requires geometry, so initialize it here
        surface_initSpeckle();

    }
}




void surface_getScatteringVector(double TSx_unit[3], double RSx_unit[3], double PUT[3], double q_vec_new[3]){
    // SG's old scattering vector code:
    //      double temp_vec[3],Q;
    //      vector_add(TSx_unit,RSx_unit,temp_vec);
    //      Q = ((2*pi*L1)/speedlight);
    //      q_vec[0]=Q*temp_vec[0];
    //      q_vec[1]=Q*temp_vec[1];
    //      q_vec[2]=Q*temp_vec[2];
    // The implementation below if from JTJ based on VZ for curved Earth

    double n_vec[3], q_vec[3], qt_vec[3], p_vec[3], t_vec[3];
    double qn, qtp, qtd, K;

    vector_add(TSx_unit,RSx_unit,q_vec);

    // n_vec is unit vector normal to surface at grid point PUT
    vector_unit(PUT, n_vec);

    // qt is tangential component of q
    qn=vector_dot_product(n_vec,q_vec);
    qt_vec[0]=q_vec[0]-n_vec[0]*qn;
    qt_vec[1]=q_vec[1]-n_vec[1]*qn;
    qt_vec[2]=q_vec[2]-n_vec[2]*qn;

    // phi vector tangent to sphere; nominally in yhat direction
    p_vec[0]=0.0;
    p_vec[1]= n_vec[2];
    p_vec[2]=-n_vec[1];
    vector_unit(p_vec, p_vec);

    // theta vector tangent to sphere; nominally in xhat direction
    vector_cross_product(p_vec, n_vec, t_vec);

    qtp = vector_dot_product(qt_vec,p_vec);
    qtd = vector_dot_product(qt_vec,t_vec);

    K = ((2*pi*L1)/speedlight);

    q_vec_new[0] = K*qtd;
    q_vec_new[1] = K*qtp;
    q_vec_new[2] = K*qn;
}



/****************************************************************************/
//  Evaluate Sigma0 Over Surface
/****************************************************************************/

void surface_calcSigma0OnSurface(int windModelType){	//compute RCS at surface (68)
    // this function assumes that the geometric and windfield properties of the
    // surface data have already been filled in.  It calculates the sigma0
    // using a bivariate Gaussian slope pdf

    double ws,mss_x,mss_y,mss_b,sxangle,q_vec[3],sigma0,sigma0_dP,x,y,P,Q4,R2,dP;

    for (int idx = 0; idx < surface.numGridPts; idx++) {

        // evaluate slope pdf (Eqn 40 [ZV 2000])
        ws = surface.windData[idx].windSpeed_ms;
        mss_x    = surface.windData[idx].mss_x;
        mss_y    = surface.windData[idx].mss_y;
        mss_b    = surface.windData[idx].mss_b;
        sxangle  = surface.data[idx].sx_angle_rad;  //incidence angle
        q_vec[0] = surface.data[idx].q[0];
        q_vec[1] = surface.data[idx].q[1];
        q_vec[2] = surface.data[idx].q[2];
        x        = -q_vec[0]/q_vec[2];
        y        = -q_vec[1]/q_vec[2];
        double mss_iso=(mss_x+mss_y)/2;

        // evaluate sigma0 (Eqn 34 [ZV 2000])
        R2     = pow(cabs(reflectionCoef(sxangle)),2);   //reflection coefficient
        Q4     = pow(vector_norm(q_vec),4) / pow(q_vec[2],4);

        if (windModelType == 0){    //isotropic model
            mss_b=0;
            P=1/(2*pi*mss_iso) * exp(-(pow(x,2)+pow(y,2))/(2*mss_iso));

            //P = 1/(2*pi*sqrt(mss_iso*mss_iso)*sqrt(1-pow(mss_b,2))) *
            //    exp(  -1/(2*(1-pow(mss_b,2))) * ( pow(x,2) / mss_iso -
            //                                      2*mss_b* (x*y)/sqrt(mss_iso*mss_iso) + pow(y,2) / mss_iso )); // PDF (61)
            sigma0 = pi * R2 * Q4 * P;    //RCS  (68)
        }
        if (windModelType == 1){    ////anisotropic model
            P = 1/(2*pi*sqrt(mss_x*mss_y)*sqrt(1-pow(mss_b,2))) *
                exp(  -1/(2*(1-pow(mss_b,2))) * ( pow(x,2) / mss_x -
                                                  2*mss_b* (x*y)/sqrt(mss_x*mss_y) + pow(y,2) / mss_y )); // PDF (61)
            sigma0 = pi * R2 * Q4 * P;    //RCS  (68)
        }

        dP       = -1/(2*pi) * ( 1/pow(mss_iso,2) +                   //dP for isotropic mss
                                 (-1.0/2 * (pow(x,2) + pow(y,2))/pow(mss_iso,3) ) ) *
                   exp(  -1.0/2 * ( ( pow(x,2) + pow(y,2) )/ mss_iso )) * dmdx(ws);

        sigma0_dP = pi * R2 * Q4 * dP;  //derivative of sigma0 (used for H matrix)

        //printf("s2 = %f\n",sqrt(pow(x,2) + pow(y,2)));
        //printf("test = %f\n",-1/2.0);
        // if use DDMA LUT, need to first ddmaLUT_initialize();
        /*
        double sigma0_GMF;
        int index_ws, index_theta;
        double ws = surface.windData[idx].windSpeed_ms;
        index_ws = (int) round((ws-0.05)/0.1)+1; //index=1,2,3...
        index_theta = (int) round(sxangle * 180/pi);
        sigma0_GMF = ddmaLUT[(index_ws-1)*90+index_theta-1];
        */

        /*
        if (idx==0 ||idx == 7260 || idx==14399){
            printf("idx = %d\n",idx);
            printf("q_vec = %f %f %f\n",q_vec[0],q_vec[1],q_vec[2]);
            printf("mss_iso = %f \n",mss_iso);
            printf("R2 = %f, Q4 = %f, P = %f\n",R2, Q4,P);
            printf("wind = %f theta_deg = %f \n",ws, sxangle*180/pi);
            printf("sigma0 = %f, sigma_GMF = %f\n", sigma0, sigma0_GMF);
            printf("\n");
        }
        */

        // store sigma parameters
        surface.data[idx].sigma0    = sigma0;
        surface.data[idx].sigma0_R2 = R2;
        surface.data[idx].sigma0_P  = P;
        surface.data[idx].sigma0_dP  = sigma0_dP;
        surface.data[idx].sigma0_Q4 = Q4;
    }
    surface_calcRainAttenOnSurface();
}

complex double reflectionCoef(double Sxangle) {
    // Based on Valery powsp.for polarization reflection coefficient
    // (see Eqs.(A11-A13) in [1]):
    double st,ct,tet;
    double complex eps,esq,ev1,ev2,b1,b2;
    tet = Sxangle;
    eps = 74.62 + I*51.92;
    st  = sin(tet);
    ct  = cos(tet);
    esq = csqrt(eps-pow(st,2));
    ev1 = (eps*ct-esq)/(eps*ct+esq); //corrected by Feixiong
    ev2 = (ct-esq)/(ct+esq);        //corrected by Feixiong
    b1  = (ev1-ev2)/2; // Right-to-left circular polarization (LHCP)
    b2  = (ev1+ev2)/2; // Right-to-right circular polarization (RHCP)
    return(b1);
}

/****************************************************************************/
//  Put all the pieces together to find the Total Scattered Power
/****************************************************************************/

void surface_composeTotalScatPowrOnSurface(int type){  //type=1
    // in some cases, we want sigma0 and in others, we want sqrt(sigma0) with phase.
    // I don't like this function, but here is how we handle both cases currently.
    switch (type) {
        case 1: // no speckle.  no mask
            for (int idx = 0; idx < surface.numGridPts; idx++){
                //surface.data[idx].total = surface.data[idx].powerFactor * surface.data[idx].sigma0
                //* surface.data[idx].mask * surface.data[idx].rainAtten_abs ;
                //surface.data[idx].total_dP = surface.data[idx].powerFactor * surface.data[idx].sigma0_dP
                //* surface.data[idx].mask * surface.data[idx].rainAtten_abs ;
                surface.data[idx].total = surface.data[idx].powerFactor * surface.data[idx].sigma0
                                          * surface.data[idx].rainAtten_abs ;
                surface.data[idx].total_dP = surface.data[idx].powerFactor * surface.data[idx].sigma0_dP
                                             * surface.data[idx].rainAtten_abs ;
            }
            break;
        case 2: // speckle, so sqrt of scattered power and phase factor - this method has problems
            for (int idx = 0; idx < surface.numGridPts; idx++){
//                surface.data[idx].total =
//                        sqrt(surface.data[idx].powerFactor * surface.data[idx].sigma0 * surface.data[idx].rainAtten_abs)
//                        * surface.data[idx].phaseShiftFactor0 * surface.data[idx].phaseShiftFactor1 * surface.data[idx].mask;
//                surface.data[idx].total_dP =
//                        sqrt(surface.data[idx].powerFactor * surface.data[idx].sigma0_dP * surface.data[idx].rainAtten_abs)
//                        * surface.data[idx].phaseShiftFactor0 * surface.data[idx].phaseShiftFactor1 * surface.data[idx].mask;
                surface.data[idx].total =
                        sqrt(surface.data[idx].powerFactor * surface.data[idx].sigma0 * surface.data[idx].rainAtten_abs)
                        * surface.data[idx].phaseShiftFactor0 * surface.data[idx].phaseShiftFactor1;
                surface.data[idx].total_dP =
                        sqrt(surface.data[idx].powerFactor * surface.data[idx].sigma0_dP * surface.data[idx].rainAtten_abs)
                        * surface.data[idx].phaseShiftFactor0 * surface.data[idx].phaseShiftFactor1;
            }
            break;
        default:
            fprintf(errPtr,"Error: Bad type in surface_composeTotalPower \n");
            exit(0);
            break;
    }
}

/****************************************************************************/
//  Surface Wind Field
/****************************************************************************/

void surface_loadSurfWindfield(windField *wf, int wfNum){	//load wind to surface frame   //wfNum=0;
    double x_m, y_m;

    if( wfNum > wf->locNumPts ){ //locNumPts=1;
        printf("Error: bad wind field loc index in surface_loadSurfWindfield"); exit(0);
    }

    wf->locCurrentPt = wfNum;  //0

    if( wf->type == 1 )
        fprintf(outputPtr,"Uniform wind case (%f m/s) ...\n\n", wf->data[wf->locCurrentPt].windSpeed_ms );

    //printf("%lf \n", surface.specularLoactionX_m);  //45000
    //printf("%lf \n", surface.specularLoactionY_m);  //45000

    //surface.specularLoactionX_m = wf->loc_rowIdx[wfNum] * wf->resolutionX_m;  //250*1000
    //surface.specularLoactionY_m = wf->loc_colIdx[wfNum] * wf->resolutionY_m;  //250000

    //printf("%lf \n", surface.specularLoactionX_m);
    //printf("%lf \n", surface.specularLoactionY_m);

    for(int i=0; i<surface.numGridPts; i++){ //0:8100
        //surface.data[i].windFieldLocation_x_m=(ix-45)*1000
        x_m = surface.specularLoactionX_m + surface.data[i].windFieldLocation_x_m;
        y_m = surface.specularLoactionY_m + surface.data[i].windFieldLocation_y_m;

        wind_getWindFieldAtXY( wf, x_m, y_m, &(surface.windData[i]) );  //important: copy wf to surface.windData[i]


    }
}

/****************************************************************************/
//  Surface Rain Field
/****************************************************************************/

void surface_calcRainAttenOnSurface(void){

    double theta1_rad, theta2_rad, rainRate_mmhr, freezeht_km, rainAtten_abs;

    if( surface.rainOnOff == 1 ){
        for (int idx = 0; idx < surface.numGridPts; idx++) {
            theta1_rad     = surface.data[idx].rx_elevationAngle_rad;
            theta2_rad     = surface.data[idx].tx_elevationAngle_rad;
            rainRate_mmhr  = surface.windData[idx].rainRate_mmhr;
            freezeht_km    = surface.windData[idx].freezingHeight_m / 1000;
            rainAtten_abs  = getRainAtten_abs( theta1_rad, theta2_rad, freezeht_km, rainRate_mmhr);
            surface.data[idx].rainAtten_abs = rainAtten_abs;
        }
    }
    else
        for (int idx = 0; idx < surface.numGridPts; idx++)
            surface.data[idx].rainAtten_abs = 1;
}

double getRainAtten_abs( double theta1_rad, double theta2_rad, double h_km, double R_mmhr){
    // theta1_rad: elevation angle to Rx in Sx Frame
    // theta2_rad: elevation angle to Tx in Sx Frame
    // h_km:       freezing height,
    // R_mmhr:     rain rate (mm/hr)
    // coef's:
    //   a = 24.312e-5; b = 0.9567;   ITU R838-3  http://www.itu.int/rec/R-REC-P.838-3-200503-I/en
    //   a = 22.326e-5; b = 1.15;     email/paper
    //   (ITU predicts lower attenuation, but is probably better, so we'll use it)

    double a             = 24.312e-5;
    double b             = 0.9567;
    double alpha         = a * pow(R_mmhr,b); // (Np/km)
    double rainAtten_abs = exp(-1 * alpha * h_km * ( csc(theta1_rad) + csc(theta2_rad) ));
    return rainAtten_abs;
}

/****************************************************************************/
//  Speckle
/****************************************************************************/

void surface_initSpeckle(void){

    double updatePeriod_s = ddm.cohIntegrationTime_s;

    surface.speckleType = 2; // hardcoded for now

    switch (surface.speckleType) {
        case 0: // speckle off, do nothing
            break;
        case 1: // speckle is just random process on surface (only used for testing)
            for(int i=0; i < surface.numGridPts; i++){
                surface.data[i].phaseShiftFactor0 = cexp(I * uniformRandf() * 2 * pi);
                surface.data[i].phaseShiftFactor1 = cexp(I * uniformRandf() * 2 * pi);
                surface.data[i].phaseShiftFactor2 = 1;
            }
            break;
        case 2: // speckle is based on geometry
            for(int i=0; i < surface.numGridPts; i++){
                surface.data[i].phaseShiftFactor0 = cexp(I * uniformRandf() * 2 * pi);
                surface.data[i].phaseShiftFactor1 = cexp( I * 2 * pi * surface.data[i].doppler_Hz * updatePeriod_s );
                surface.data[i].phaseShiftFactor2 = cexp( I * 2 * pi * surface.data[i].doppler_Hz * updatePeriod_s );
            }
            break;
        default:
            fprintf(errPtr,"Error: bad speckleType in surface_initSpeckle\n");
            exit(0);
    }
}

void surface_updateSpeckle(void){

    switch (surface.speckleType) {
        case 0: // speckle off, do nothing
            break;
        case 1: // random phase each time  (only used for testing)
            for(int i=0; i < surface.numGridPts; i++)
                surface.data[i].phaseShiftFactor1 = cexp(I * uniformRandf() * 2 * pi);
            break;
        case 2: // update phase based on geometry changes
            for(int i=0; i < surface.numGridPts; i++){
                surface.data[i].phaseShiftFactor1 *= surface.data[i].phaseShiftFactor2;
            }
            break;
        default:
            fprintf(errPtr,"Error: bad speckleType in surface_initSpeckle\n");
            exit(0);
    }
}


/****************************************************************************/
//  Misc
/****************************************************************************/

void surface_createSurfaceMask(void){
    int i, j=0;//add =0
    double val, maxDelay_s,maxDoppler_Hz,minDoppler_Hz;
    double maxDelayEllipse_s;
    double maxDopplerEllipse_Hz, minDopplerEllipse_Hz;

    // find max delay & Doppler on surface
    maxDelay_s    = surface.data[SURFINDEX(0, j)].delay_s;     //just initial value
    minDoppler_Hz = surface.data[SURFINDEX(0, j)].doppler_Hz;  //just initial value
    maxDoppler_Hz = surface.data[SURFINDEX(0, j)].doppler_Hz;  //just initial value

    for (i = 0; i < surface.numGridPts; i++) {
        val = surface.data[i].delay_s;
        if( val > maxDelay_s ) maxDelay_s = val;
        val = surface.data[i].doppler_Hz;
        if( val > maxDoppler_Hz ) maxDoppler_Hz = val;
        if( val < minDoppler_Hz ) minDoppler_Hz = val;
    }

    // determine DDM working range (i.e the max complete
    // delay "ellipse" on the surface by looking at the edges of the surface
    int cx = (int)floor( (1.0 * surface.numGridPtsX ) / 2 );//45
    int cy = (int)floor( (1.0 * surface.numGridPtsY ) / 2 );//45
    double val1  = surface.data[SURFINDEX(0, cy)].delay_s;//[0 45]
    double val2  = surface.data[SURFINDEX(surface.numGridPtsX - 1, cy)].delay_s;//[89 45]
    double val3  = surface.data[SURFINDEX(cx,0)].delay_s;//[45 0]
    double val4  = surface.data[SURFINDEX(cx, surface.numGridPtsY - 1)].delay_s;//[45 89]
    maxDelayEllipse_s = val1;
    if( val2 < maxDelayEllipse_s ) maxDelayEllipse_s = val2;//change to max delay >
    if( val3 < maxDelayEllipse_s ) maxDelayEllipse_s = val3;
    if( val4 < maxDelayEllipse_s ) maxDelayEllipse_s = val4;

    // now we need to move at least a chip shorter from the edge
    //maxDelayEllipse_s = maxDelayEllipse_s - 2.0/chipRate_cs;

    // find the range of Doppler that is within that max delay ellipse
    maxDopplerEllipse_Hz = minDoppler_Hz;//just initilize
    minDopplerEllipse_Hz = minDoppler_Hz;
    for (i = 0; i < surface.numGridPtsX; i++) {
        for (j = 0; j < surface.numGridPtsY; j++) {
            if( surface.data[SURFINDEX(i, j)].delay_s <= maxDelayEllipse_s ){
                val = surface.data[SURFINDEX(i, j)].doppler_Hz;
                if( val > maxDopplerEllipse_Hz )  maxDopplerEllipse_Hz = val;
                if( val < minDopplerEllipse_Hz )  minDopplerEllipse_Hz = val;
            }
        }
    }

    // create mask on the surface so that only points within
    // the max valid delay are integrated

    for (i = 0; i < surface.numGridPts; i++) {
        //if( surface.data[i].delay_s <= (maxDelayEllipse_s + 2.0/chipRate_cs) )
        if( surface.data[i].delay_s <= maxDelayEllipse_s )
            surface.data[i].mask = 1;
        else
            surface.data[i].mask = 0;
    }

    /*
    int saveSurfaceLogFile = getParamFromConfigFileWDefault_int("run.saveSurfLogFile", 0);
    if (saveSurfaceLogFile == 1) {
        FILE *outputPtr=fopen("info.txt","w+");//output file
        fprintf(outputPtr,"Surface & DDM Resolution Check -------------------------------------------- \n");
        fprintf(outputPtr,"  Maximum delay on surface:                              %4.2f (chips)\n", maxDelay_s * 1.023e6);
        fprintf(outputPtr,"  Doppler range on surface:                              %4.2f to %4.2f (kHz)\n", minDoppler_Hz / 1000, maxDoppler_Hz / 1000);
        fprintf(outputPtr,"  Largest complete delay ellipse:                        %4.2f (chips)\n", maxDelayEllipse_s * 1.023e6);
        fprintf(outputPtr,"  Doppler range within delay ellipse:                    %4.2f to %4.2f (kHz) \n", minDopplerEllipse_Hz/1000, maxDopplerEllipse_Hz/1000);
        fprintf(outputPtr,"  Minimum working DDM size (no masking):                 %4.2f chips by %f kHz \n",
                ceil( 2 + maxDelay_s * 1.023e6 ), ceil(2 + maxDoppler_Hz/1000));
        fprintf(outputPtr,"  Minimum working DDM size (masking):                    %4.2f chips by %f kHz \n",
                ceil( maxDelayEllipse_s * 1.023e6 ), ceil(2 + (maxDopplerEllipse_Hz-minDopplerEllipse_Hz)/1000));
        fprintf(outputPtr,"  Current working DDM extent is:                         %4.2f chips by %4.2f kHz \n\n",
                ddm.delayRes_chips*ddm.numDelayBins, ddm.numDoppBins*ddm.dopplerRes_Hz/1000);
        fclose(outputPtr);
    }
    */

}

void surface_resetToZero(void){

    surfacePixel zeroPixel;
    zeroPixel.mask             = 1;
    zeroPixel.delay_s          = 0;
    zeroPixel.doppler_Hz       = 0;
    zeroPixel.sx_angle_rad     = 0;
    zeroPixel.q[0]             = 0;
    zeroPixel.q[1]             = 0;
    zeroPixel.q[2]             = 0;
    zeroPixel.powerFactor      = 0;
    zeroPixel.phaseShiftFactor0= 0;
    zeroPixel.phaseShiftFactor1= 0;
    zeroPixel.phaseShiftFactor2= 0;
    zeroPixel.total            = 0;
    zeroPixel.bin_index        = 0;
    zeroPixel.i                = 0;
    zeroPixel.j                = 0;
    zeroPixel.sigma0           = 0;
    zeroPixel.sigma0_R2        = 0;
    zeroPixel.sigma0_P         = 0;
    zeroPixel.sigma0_Q4        = 0;
    zeroPixel.position[0]      = 0;
    zeroPixel.position[1]      = 0;
    zeroPixel.position[2]      = 0;
    zeroPixel.antennaGainRx_abs = 0;
    zeroPixel.antennaGainTx_abs = 0;

    for(int i=0; i < surface.numGridPts; i++){
        surface.data[i] = zeroPixel;
    }
}


void surface_saveWindToFile(void) {  //added by Feixiong
    FILE *outp = fopen("surfaceWind.dat", "wb");
    for (int i = 0;i<surface.numGridPts;i++) {
        fwrite(&surface.windData[i].windSpeed_ms, sizeof(double), 1, outp);
    }
    fclose(outp);
}

void surface_saveDopplerToFile(void){  //added by Feixiong
    FILE *outp = fopen("surfaceDoppler.dat","wb");
    for(int i=0;i<surface.numGridPts;i++){
        fwrite(&surface.data[i].doppler_Hz, sizeof(double),1,outp);
    }
    fclose(outp);
}

void surface_saveDelayToFile(void) {  //added by Feixiong
    FILE *outp = fopen("surfaceDelay.dat", "wb");
    for (int i = 0;i < surface.numGridPts;i++) {
        fwrite(&surface.data[i].delay_s, sizeof(double), 1, outp);
    }
    fclose(outp);
}
