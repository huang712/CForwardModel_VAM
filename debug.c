//
// Created by Feixiong Huang on 11/14/17.
//

#include "gnssr.h"

void printfGeometry(geometryData orbitGeometry){
    printf("tx_pos = %f   ",orbitGeometry.g[0].tx_pos_ecef[0]);
    printf("%f   ",orbitGeometry.g[0].tx_pos_ecef[1]);
    printf("%f   ",orbitGeometry.g[0].tx_pos_ecef[2]);
    printf("\n");
    printf("rx_pos = %f   ",orbitGeometry.g[0].rx_pos_ecef[0]);
    printf("%f   ",orbitGeometry.g[0].rx_pos_ecef[1]);
    printf("%f   ",orbitGeometry.g[0].rx_pos_ecef[2]);
    printf("\n");
    printf("sx_pos = %f   ",orbitGeometry.g[0].sx_pos_ecef[0]);
    printf("%f   ",orbitGeometry.g[0].sx_pos_ecef[1]);
    printf("%f   ",orbitGeometry.g[0].sx_pos_ecef[2]);
    printf("\n");
    //printf("sx_vel = %f   ",orbitGeometry.g[0].sx_vel_ecef[0]);
    //printf("%f   ",orbitGeometry.g[0].sx_vel_ecef[1]);
    //printf("%f   ",orbitGeometry.g[0].sx_vel_ecef[2]);
    //printf("\n");
}

/*
printf("delay = %f\n",surface.data[0].delay_s);
printf("sigma0 = %f\n",surface.data[0].sigma0);

printf("total = %e\n", surface.data[1000].total);
printf("total = %e\n", surface.data[3000].total);

printf("DDM = %e %e %e \n",DDM[1],DDM[100],DDM[1000]);

printf("delay = %e\n",jacob->data[0].value);
*/


/*
void ddm_Hmatrix_old(struct metadata meta, struct inputWindField iwf, struct Jacobian *jacob){


    unsigned height = (unsigned) surface.numGridPts;  //14400
    unsigned width  = (unsigned) ddm.numBins;	//160000
    int startDelay_bin    = meta.resample_startBin[0];
    int startDoppelr_bin  = meta.resample_startBin[1];
    int resDelay_bins     = meta.resample_resolution_bins[0];
    int resDoppelr_bins   = meta.resample_resolution_bins[1];
    int numDelayBins      = meta.resample_numBins[0];
    int numDoppelrBins    = meta.resample_numBins[1];


    double temp = 1.0 * numDelayBins * numDoppelrBins;

    _ddm_zero( H, width );

    int i, j, surface_index;
    int numBins = numDelayBins * numDoppelrBins; //187
    int numSurfacePt = meta.numGridPoints[0] * meta.numGridPoints[1] / 100; //144
    double **H0; //2D array numBins * numSurfacePt 187x144
    H0 = (double**)malloc(sizeof(double*) * numBins);//
    for (i = 0; i < numBins; i++){//
        H0[i] = (double*)malloc(sizeof(double)*numSurfacePt);
    }

    double *H0_lat_vec = (double *)calloc(numSurfacePt, sizeof(double)); //144
    double *H0_lon_vec = (double *)calloc(numSurfacePt, sizeof(double)); //144

    int i0=0; //index of ddmbin
    int j0=0; //index of surfacePts
    //int lat_index, lon_index;
    for(int m = 0 ; m < surface.numGridPtsX/10 ; m = m + 1)
    {
        for(int n = 0 ; n < surface.numGridPtsY/10 ; n = n + 1)
        {
            surface_index = (m * 10 + 4) * surface.numGridPtsY + n * 10 + 4; //resolution from 1km to 10km
            H0_lat_vec[j0] = surface.data[surface_index].pos_llh[0];
            H0_lon_vec[j0] = surface.data[surface_index].pos_llh[1];

            //printf("bin_index = %d\n",surface.data[surface_index].bin_index);
            if (surface.data[surface_index].bin_index < 0){ //if bin_index=-1, the surface point is outside the glistening zone, dHdm=0
                for (i0=0;i0<numBins;i0++){
                    H0[i0][j0]=0;
                }
                j0 = j0+1;

            }
            if (surface.data[surface_index].bin_index >= 0){
                H[surface.data[surface_index].bin_index] = surface.data[surface_index].total_dP;
                ddm_convolveH_FFT(2); //convolve H with ambiguity function, save into H

                i0=0;
                // Resample: from high resolution DDM to 17x11 DDM
                for (int k=startDoppelr_bin; k < (numDoppelrBins*resDoppelr_bins + startDoppelr_bin); k+=resDoppelr_bins) {
                    for (int l=startDelay_bin; l < (numDelayBins*resDelay_bins + startDelay_bin); l+=resDelay_bins) {
                        H0[i0][j0] = creal(H[DDMINDEX(k,l)]) * 100;
                        //jacob->data[i].value = creal(H[DDMINDEX(k,l)]) * 100;//H is derivative respect to pixel in 10km resolution
                        i0 = i0+1;
                    }
                }
                _ddm_zero( H, width ); // length of H = 187
                j0 = j0+1;
            }
        }
    }

    //************** Compute H matrix for points on lat/lon ******************
    //create lat/lon vector
    double *lon_vec, *lat_vec;
    lon_vec = (double *)calloc(iwf.numPtsLon,sizeof(double));
    lat_vec = (double *)calloc(iwf.numPtsLat,sizeof(double));
    for (i = 0; i < iwf.numPtsLon; i++){
        lon_vec[i] = iwf.lon_min_deg + i*iwf.resolution_lon_deg-360;   //-180 ~ 180
    }
    for (i = 0; i < iwf.numPtsLat; i++){
        lat_vec[i] = iwf.lat_min_deg + i*iwf.resolution_lat_deg;
    }

    //initialize matrix bi_weight[numSurfacePt][4]  bi_index[numSurfacePt][4]  in the bilinear interpolation
    double **bi_weight;
    int **bi_index;
    bi_weight = (double**)malloc(sizeof(double*) * numSurfacePt);
    for (i = 0; i < numSurfacePt; i++){
        bi_weight[i] = (double*)malloc(sizeof(double)*4);
    }
    bi_index = (int**)malloc(sizeof(int*) * numSurfacePt);
    for (i = 0; i < numSurfacePt; i++){
        bi_index[i] = (int*)malloc(sizeof(int)*4);
    }

    //bilinear interpolation and get bi_index, bi_weight
    for (i=0; i< numSurfacePt; i++){
        bilinear_interp(lat_vec, lon_vec, iwf.numPtsLat, iwf.numPtsLon, H0_lat_vec[i], H0_lon_vec[i], bi_index[i], bi_weight[i], fabs(iwf.resolution_lat_deg));
    }

    int k = 0;
    int *bi_index1 = (int *)calloc(numSurfacePt*4,sizeof(int));  //reshape bi_index to a 1D array
    for (i = 0; i< numSurfacePt; i++){
        for (j=0;j<4;j++){
            bi_index1[i*4+j]=bi_index[i][j];
        }
    }
    bubble(bi_index1,numSurfacePt*4); // put in order

    // throw repeated index and find the length N
    for (i = 1; i < numSurfacePt*4; i++)
    {
        if (bi_index1[k] != bi_index1[i])
        {
            bi_index1[k + 1] = bi_index1[i];
            k++;  //index of the non-repeated array
        }
    }
    int numPt_LL=k+1; //length of M (num of points on lat/lon) = N = numPt_LL

    int *indexLL = (int *)calloc(numPt_LL,sizeof(int));
    for (i=0; i<numPt_LL; i++){
        indexLL[i]=bi_index1[i];
    }

    //construct T matrix T: 144 x N
    double **T; // T[144][N]interpolation transformation matrix: fill it with bi_weight[144][4]
    T = (double**)malloc(sizeof(double*) * numSurfacePt);
    for (i = 0; i < numSurfacePt; i++){
        T[i] = (double*)malloc(sizeof(double)*numPt_LL);
    }

    //initialize T
    for (i = 0; i<numSurfacePt; i++){
        for(j = 0; j<numPt_LL; j++){
            T[i][j]=0;
        }
    }
    //fill T matrix
    for (i = 0; i<numSurfacePt; i++){
        for (j = 0; j<numPt_LL; j++){
            for (k=0; k<4; k++){
                if(bi_index[i][k]==indexLL[j]){
                    T[i][j]=bi_weight[i][k];
                    //printf("T= %f, i= %d, j=%d\n",T[i][j],i,j);
                }
            }
        }
    }


    int saveH0,saveT;
    saveH0 = 1;
    saveT = 1;
    if (saveH0 == 1){
        FILE *outp = fopen("H0.dat", "wb");
        for (j = 0;j< 144;j++) {
            for (i = 0; i < 187; i++){
                if(fabs(H0[i][j])>1e-10){
                    printf("something wrong in H0 second!! i= %d, j= %d H=%e\n",i,j,H0[i][j]);
                }
                fwrite(&H0[i][j], sizeof(double), 1, outp);
            }
        }
        fclose(outp);
    }

    if (saveT == 1){
        FILE *outp1 = fopen("T.dat", "wb");
        for (j = 0;j < numPt_LL;j++) {
            for (i = 0; i < numSurfacePt; i++){
                fwrite(&T[i][j], sizeof(double), 1, outp1);
            }
        }
        fclose(outp1);
    }


    //matrix multiplication H[187][110] = H0[187][144] * T[144][110]
    double **H_LL; //H matrix respect to lat/lon 187x110
    H_LL = (double**)malloc(sizeof(double*) * numBins);//
    for (i = 0; i < numBins; i++){//
        H_LL[i] = (double*)malloc(sizeof(double)*numPt_LL);
    }
    for (i = 0; i < numBins; i++){
        for (j = 0; j<numPt_LL; j++){
            H_LL[i][j]=0;
            for (k=0; k<numSurfacePt; k++){
                H_LL[i][j] += H0[i][k] * T[k][j];
            }
        }
    }

    //save H_LL to jabob.data[187x144]
    for (i = 0; i < numBins; i++){
        for (j = 0; j<numPt_LL; j++){
            jacob->data[j*numBins+i].value = H_LL[i][j];
            jacob->data[j*numBins+i].lat_deg = iwf.data[indexLL[j]].lat_deg;
            jacob->data[j*numBins+i].lon_deg = iwf.data[indexLL[j]].lon_deg;
        }
    }
    jacob->numDDMbins = numBins;
    jacob->numPts_LL = numPt_LL;

    for(i = 0; i< numPt_LL; i++){
        jacob->Pts_lat_vec[i] = iwf.data[indexLL[i]].lat_deg;
        jacob->Pts_lon_vec[i] = iwf.data[indexLL[i]].lon_deg;
        jacob->Pts_ind_vec[i] = indexLL[i];
    }

    //free menmory
    for (i = 0; i < numBins; ++i){
        free(H0[i]); free(H_LL[i]);
    }
    free(H0); free(H_LL);
    free(H0_lat_vec); free(H0_lon_vec); free(lat_vec); free(lon_vec);
    free(bi_index1); free(indexLL);
    for (i = 0; i < numSurfacePt; ++i){
        free(bi_index[i]); free(bi_weight[i]); free(T[i]);
    }
    free(bi_index); free(bi_weight); free(T);
}
*/