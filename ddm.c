//
// Created by Feixiong Huang on 10/22/17.
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
// The DDM parameters are stored in a global structure "ddm", while buffers for
// the actual DDM, the ambiguity function, etc are in separate globals, "DDM",
// "DDM_avg", "DDM_amb", etc ... The parameters are read in using the
// "ddm_initializeFromConfigFile" function, which also allocates memory for the
// DDM and precalculates the ambiguity function.  Then one calls various functions
// to assemble the DDM.  See ddm_* calls from main.c
//
//****************************************************************************/

#include "gnssr.h"
#include "forwardmodel.h"

//prototypes only used within wind.c
//int find_nearest(double *vec, int size, double value);

void ddm_initialize(struct metadata meta) {

    //  Initialization of DDM:  Read in DDM params from config file and alloc
    //  and initialize DDM buffers

    // reads DDM parameters from the configuration file into global struct
    ddm.numDelayBins         = meta.numDelaybins;   //400
    ddm.numDoppBins          = meta.numDopplerbins; //400
    ddm.delayOffset_bins     = meta.specular_delayBinIdx; //6.18 * 5
    ddm.dopplerOffset_bins   = meta.specular_dopplerBinIdx; //5.02 * 20
    ddm.dopplerRes_Hz        = meta.dopplerRes_Hz;
    ddm.delayRes_chips       = meta.delayRez_chips;
    ddm.cohIntegrationTime_s = 0.001;

    ddm.numBins              = ddm.numDelayBins * ddm.numDoppBins;//160000
    ddm.chipsPerSec          = 1.023e6;
    ddm.refPower_dB          = -200; // when plotting log, this is what is used for log10(0)

    // allocate & init DDM buffers
    _ddm_alloc(&DDM,      ddm.numBins);
    _ddm_alloc(&H,        ddm.numBins);
    _ddm_alloc(&DDM_temp, ddm.numBins);
    _ddm_alloc(&DDM_avg,  ddm.numBins);
    _ddm_alloc(&DDM_amb,  ddm.numBins);
    _ddm_alloc(&DDM_amb1, ddm.numBins);
    _ddm_alloc(&DDM_amb2, ddm.numBins);
    _ddm_alloc(&DDM_store,ddm.numBins);

    _ddm_zero( DDM,      ddm.numBins );
    _ddm_zero( H,        ddm.numBins );
    _ddm_zero( DDM_temp, ddm.numBins );
    _ddm_zero( DDM_avg,  ddm.numBins );
    _ddm_zero( DDM_amb,  ddm.numBins );
    _ddm_zero( DDM_amb1, ddm.numBins );
    _ddm_zero( DDM_amb2, ddm.numBins );
    _ddm_zero( DDM_store,ddm.numBins );

    // optimize FFT based on requested DDM size
    ddm.FFTWPLAN  = fftw_plan_dft_2d(ddm.numDoppBins, ddm.numDelayBins, (fftw_complex*) DDM, (fftw_complex*) DDM, FFTW_FORWARD, FFTW_MEASURE);
    ddm.IFFTWPLAN = fftw_plan_dft_2d(ddm.numDoppBins, ddm.numDelayBins, (fftw_complex*) DDM, (fftw_complex*) DDM, FFTW_BACKWARD, FFTW_MEASURE);
    h.FFTWPLAN  = fftw_plan_dft_2d(ddm.numDoppBins, ddm.numDelayBins, (fftw_complex*) H, (fftw_complex*) H, FFTW_FORWARD, FFTW_MEASURE);
    h.IFFTWPLAN = fftw_plan_dft_2d(ddm.numDoppBins, ddm.numDelayBins, (fftw_complex*) H, (fftw_complex*) H, FFTW_BACKWARD, FFTW_MEASURE);


    // precalculate and store the ambiguity function
    ddm_initACF();
    ddm_initAmbFuncBuffers(meta.prn_code);

    // calculate the thermal noise power
    ddm_initThermalNoise(meta);

}

void ddm_cleanup(void){

    // free DDM buffers
    _ddm_free(DDM);
    _ddm_free(DDM_temp);
    _ddm_free(DDM_avg);
    _ddm_free(DDM_amb);
    _ddm_free(DDM_amb1);
    _ddm_free(DDM_amb2);
    _ddm_free(DDM_store);

    fftw_destroy_plan(ddm.FFTWPLAN);
    fftw_destroy_plan(ddm.IFFTWPLAN);
}

/****************************************************************************/
//  Surface Binning and Mapping to DDM
/****************************************************************************/

void ddm_binSurface(void) {
    // loop over surface and use the delay and Doppler value to determine what
    // DDM bin that surface patch maps to
    int delayBin,dopplerBin;
    for(int i=0; i<surface.numGridPts; i++) {
        delayBin   = (int)(round(surface.data[i].delay_s * ddm.chipsPerSec / ddm.delayRes_chips)) + ddm.delayOffset_bins;
        dopplerBin = (int)(round(surface.data[i].doppler_Hz/ddm.dopplerRes_Hz)) + ddm.dopplerOffset_bins;

        if((dopplerBin < ddm.numDoppBins) && (dopplerBin >= 0) &&
           (delayBin < ddm.numDelayBins) && (delayBin >= 0 )) {
            surface.data[i].bin_index = DDMINDEX(dopplerBin, delayBin);
        }
        else
            surface.data[i].bin_index = -1;
    }
}

void ddm_Hmatrix(struct metadata meta, struct inputWindField iwf, struct Jacobian *jacob){

    printf("compute Jacobian matrix respect to points on lat/lon\n");
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
    int numSurfacePt1 = meta.numGridPoints[0] * meta.numGridPoints[1]; //14400  num point in 1 km resolution
    int numSurfacePt10 = (int)(numSurfacePt1/100); //144  num point in 10 km resolution

    //H0: 2D array numBins * numSurfacePt 187x144
    double **H0 = (double **)calloc(numBins, sizeof(double *));
    for (i = 0; i < numBins; i++){//
        H0[i] = (double *)calloc(numSurfacePt10,sizeof(double));
    }

    double *H0_lat_vec = (double *)calloc(numSurfacePt10, sizeof(double)); //144
    double *H0_lon_vec = (double *)calloc(numSurfacePt10, sizeof(double)); //144

    int i0=0; //index of ddmbin
    int j0=0; //index of surfacePts
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
    //x'=Mx  M[14400][K]
    //X=Tx   T[144][K]
    //First create M[14400][K] matrix from bi_index0[14400][4] and bi_weight0[14400][4]

    int *bi_index1 = (int *)calloc(numSurfacePt1*4,sizeof(int));  //reshape bi_index to a 1D array
    for (i = 0; i< numSurfacePt1; i++){
        for (j=0;j<4;j++){
            bi_index1[i*4+j]=bi_index0[i][j];
        }
    }

    //throw repeated index in bi_index1
    int *a = (int *)calloc(iwf.numPts,sizeof(int)); //array to store occurence times of each index;
    memset(a,0,sizeof(a));
    //int a[2000]={0};  //array to store occurence times of each index; length must be larger than all bi_index1[
    for (i=0;i<numSurfacePt1*4;i++){
        a[bi_index1[i]]=a[bi_index1[i]]+1;
    }

    //find the lenth numPt_LL and store index into the indexLL
    int numPt_LL = 0;
    for (i=0;i<iwf.numPts;i++){
        if(a[i]>0) numPt_LL++;
    }

    int *indexLL = (int *)calloc(numPt_LL,sizeof(int));
    j=0;
    for (i=0;i<iwf.numPts;i++){
        if(a[i]>0){
            indexLL[j]=i;
            j++;
        }
    }
    free(a);

//    int save_indexLL=1;
//    if(save_indexLL==1){
//        FILE *outp = fopen("/users/fax/CYGNSS/VAM/VAM_MATLAB/indexLL.dat", "wb");
//        fwrite(indexLL, sizeof(int), numPt_LL, outp);
//        fclose(outp);
//    }

    //construct M matrix M: 14400 x K
    double **M; // T[14400][K] interpolation transformation matrix: fill it with bi_weight[14400][4]
    M = (double**)calloc(numSurfacePt1, sizeof(double*));
    for (i = 0; i < numSurfacePt1; i++){
        M[i] = (double*)calloc(numPt_LL,sizeof(double));
    }

    //initialize M
    for (i = 0; i<numSurfacePt1; i++){
        for(j = 0; j<numPt_LL; j++){
            M[i][j]=0;
        }
    }

    int k;
    //fill M matrix
    for (i = 0; i<numSurfacePt1; i++){
        for (j = 0; j<numPt_LL; j++){
            for (k=0; k<4; k++){
                if(bi_index0[i][k]==indexLL[j]){
                    M[i][j]=bi_weight0[i][k];
                    //printf("M= %f, i= %d, j=%d\n",T[i][j],i,j);
                }
            }
        }
    }

    int saveM;
    saveM=0;
    if(saveM==1){
        FILE *outp = fopen("M.dat", "wb");
        for (j = 0;j< numPt_LL;j++) {
            for (i = 0; i < 14400; i++){
                fwrite(&M[i][j], sizeof(double), 1, outp);
            }
        }
        fclose(outp);
    }

    double **T; // T[144][K]
    T = (double**)calloc(numSurfacePt10,sizeof(double*));
    for (i = 0; i < numSurfacePt10; i++){
        T[i] = (double*)calloc(numPt_LL,sizeof(double));
    }

    //initialize T
    for (i = 0; i<numSurfacePt10; i++){
        for(j = 0; j<numPt_LL; j++){
            T[i][j]=0;
        }
    }

    //from M to T : 14400 to 144  average
    //M[14400][numPt_LL]  T[144][numPt_LL]
    int index;
    i=0;
    for (int m = 0 ; m < surface.numGridPtsX/10 ; m = m + 1) {
        for (int n = 0; n < surface.numGridPtsY / 10; n = n + 1) {

            for (int i1=0; i1<10; i1++){
                for (int j1=0; j1<10; j1++){
                    index = (10*m+i1)*surface.numGridPtsY +10*n+j1;
                    for (int k=0;k<numPt_LL;k++){
                        T[i][k]+=M[index][k];
                        if(i1==9&&j1==9) {T[i][k]=T[i][k]/100;}
                    }
                }
            }
            i++;  //i=0:143
        }
    }

    //printf("%d %d %d %d\n", bi_index[0][0],bi_index[0][1],bi_index[0][2],bi_index[0][3]);
    //printf("%f %f %f %f\n", bi_weight[0][0],bi_weight[0][1],bi_weight[0][2],bi_weight[0][3]);
    //printf("%d %d %d %d\n", indexLL[0],indexLL[1],indexLL[2],indexLL[3]);
    //printf("%f %f %f %f\n", T[0][0],T[0][1],T[0][2],T[0][3]);


    int saveH0,saveT;
    saveH0 = 0;
    saveT = 0;
    if (saveH0 == 1){
        FILE *outp = fopen("H0.dat", "wb");
        for (j = 0;j< numSurfacePt10;j++) {
            for (i = 0; i < numBins; i++){
                fwrite(&H0[i][j], sizeof(double), 1, outp);
            }
        }
        fclose(outp);
    }

    if (saveT == 1){
        FILE *outp1 = fopen("T.dat", "wb");
        for (j = 0;j < numPt_LL;j++) {
            for (i = 0; i < numSurfacePt10; i++){
                fwrite(&T[i][j], sizeof(double), 1, outp1);
            }
        }
        fclose(outp1);
    }

    //matrix multiplication H[187][110] = H0[187][144] * T[144][110]
    double **H_LL; //H matrix respect to lat/lon 187x110
    H_LL = (double**)calloc(numBins,sizeof(double*));//
    for (i = 0; i < numBins; i++){//
        H_LL[i] = (double*)calloc(numPt_LL,sizeof(double));
    }
    for (i = 0; i < numBins; i++){
        for (j = 0; j<numPt_LL; j++){
            H_LL[i][j]=0;
            for (k=0; k<numSurfacePt10; k++){
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

    //save indexLL into jacob
    for(i = 0; i< numPt_LL; i++){
        jacob->Pts_lat_vec[i] = iwf.data[indexLL[i]].lat_deg;
        jacob->Pts_lon_vec[i] = iwf.data[indexLL[i]].lon_deg;
        jacob->Pts_ind_vec[i] = indexLL[i];
    }

    //free menmory
    for (i = 0; i < numBins; i++){
        free(H0[i]);
    }
    free(H0);

    for (i = 0; i < numBins; i++){
        free(H_LL[i]);
    }
    free(H_LL);

    free(H0_lat_vec); free(H0_lon_vec);
    free(bi_index1); free(indexLL);

    for (i = 0; i < numSurfacePt1; i++){
        free(M[i]);
    }
    free(M);

    for (i = 0; i < numSurfacePt10; i++){
        free(T[i]);
    }
    free(T);
}



void ddm_Hmatrix_old(struct metadata meta, struct inputWindField iwf, struct Jacobian *jacob){

    printf("compute Jacobian matrix respect to points on lat/lon\n");
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
    int numSurfacePt1 = meta.numGridPoints[0] * meta.numGridPoints[1]; //14400  num point in 1 km resolution
    int numSurfacePt10 = (int)(numSurfacePt1/100); //144  num point in 10 km resolution

    //double **H0; //2D array numBins * numSurfacePt 187x144
    //H0 = (double**)malloc(sizeof(double*) * numBins);//
    //for (i = 0; i < numBins; i++){//
    //    H0[i] = (double*)malloc(sizeof(double)*numSurfacePt10);
    //}
    //memset(H0,0, sizeof(H0));
    double **H0 = (double **)calloc(numBins, sizeof(double *));
    for (i = 0; i < numBins; i++){//
        H0[i] = (double *)calloc(numSurfacePt10,sizeof(double));
    }

    double *H0_lat_vec = (double *)calloc(numSurfacePt10, sizeof(double)); //144
    double *H0_lon_vec = (double *)calloc(numSurfacePt10, sizeof(double)); //144

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

    printf("coordinate transformation\n");
    //************** Compute H matrix for points on lat/lon ******************
    //x'=Mx  M[14400][K]
    //X=Tx   T[144][K]
    //First create M[14400][K] matrix from bi_index0[14400][4] and bi_weight0[14400][4]

    int *bi_index1 = (int *)calloc(numSurfacePt1*4,sizeof(int));  //reshape bi_index to a 1D array
    for (i = 0; i< numSurfacePt1; i++){
        for (j=0;j<4;j++){
            bi_index1[i*4+j]=bi_index0[i][j];
        }
    }

    double start, end;
    start = clock();

    bubble(bi_index1,numSurfacePt1*4); // put in order
    end =clock();
    printf("bubble time: %f seconds\n", (end-start)/CLOCKS_PER_SEC);

    start = clock();
    // throw repeated index and find the length K
    int k = 0;
    for (i = 1; i < numSurfacePt1*4; i++)
    {
        if (bi_index1[k] != bi_index1[i])
        {
            bi_index1[k + 1] = bi_index1[i];
            k++;  //index of the non-repeated array
        }
    }
    int numPt_LL=k+1; //length of M (num of points on lat/lon) = K = numPt_LL

    int *indexLL = (int *)calloc(numPt_LL,sizeof(int));
    for (i=0; i<numPt_LL; i++){
        indexLL[i]=bi_index1[i];
    }

    int save_indexLL=1;
    if(save_indexLL==1){
        FILE *outp = fopen("/users/fax/CYGNSS/VAM/MATLAB/indexLL.dat", "wb");
        fwrite(indexLL, sizeof(int), numPt_LL, outp);
        fclose(outp);
    }

    //construct M matrix M: 14400 x K
    double **M; // T[14400][K] interpolation transformation matrix: fill it with bi_weight[14400][4]
    M = (double**)calloc(numSurfacePt1, sizeof(double*));
    for (i = 0; i < numSurfacePt1; i++){
        M[i] = (double*)calloc(numPt_LL,sizeof(double));
    }

    //initialize M
    for (i = 0; i<numSurfacePt1; i++){
        for(j = 0; j<numPt_LL; j++){
            M[i][j]=0;
        }
    }


    //fill M matrix
    for (i = 0; i<numSurfacePt1; i++){
        for (j = 0; j<numPt_LL; j++){
            for (k=0; k<4; k++){
                if(bi_index0[i][k]==indexLL[j]){
                    M[i][j]=bi_weight0[i][k];
                    //printf("M= %f, i= %d, j=%d\n",T[i][j],i,j);
                }
            }
        }
    }

    int saveM;
    saveM=0;
    if(saveM==1){
        FILE *outp = fopen("M.dat", "wb");
        for (j = 0;j< numPt_LL;j++) {
            for (i = 0; i < 14400; i++){
                fwrite(&M[i][j], sizeof(double), 1, outp);
            }
        }
        fclose(outp);
    }

    double **T; // T[144][K]
    T = (double**)calloc(numSurfacePt10,sizeof(double*));
    for (i = 0; i < numSurfacePt10; i++){
        T[i] = (double*)calloc(numPt_LL,sizeof(double));
    }

    //initialize T
    for (i = 0; i<numSurfacePt10; i++){
        for(j = 0; j<numPt_LL; j++){
            T[i][j]=0;
        }
    }

    //from M to T : 14400 to 144  average
    //M[14400][numPt_LL]  T[144][numPt_LL]
    int index;
    i=0;
    for (int m = 0 ; m < surface.numGridPtsX/10 ; m = m + 1) {
        for (int n = 0; n < surface.numGridPtsY / 10; n = n + 1) {

            for (int i1=0; i1<10; i1++){
                for (int j1=0; j1<10; j1++){
                    index = (10*m+i1)*surface.numGridPtsY +10*n+j1;
                    for (int k=0;k<numPt_LL;k++){
                        T[i][k]+=M[index][k];
                        if(i1==9&&j1==9) {T[i][k]=T[i][k]/100;}
                    }
                }
            }
            i++;  //i=0:143
        }
    }



    //printf("%d %d %d %d\n", bi_index[0][0],bi_index[0][1],bi_index[0][2],bi_index[0][3]);
    //printf("%f %f %f %f\n", bi_weight[0][0],bi_weight[0][1],bi_weight[0][2],bi_weight[0][3]);
    //printf("%d %d %d %d\n", indexLL[0],indexLL[1],indexLL[2],indexLL[3]);
    //printf("%f %f %f %f\n", T[0][0],T[0][1],T[0][2],T[0][3]);


    int saveH0,saveT;
    saveH0 = 0;
    saveT = 0;
    if (saveH0 == 1){
        FILE *outp = fopen("H0.dat", "wb");
        for (j = 0;j< numSurfacePt10;j++) {
            for (i = 0; i < numBins; i++){
                fwrite(&H0[i][j], sizeof(double), 1, outp);
            }
        }
        fclose(outp);
    }

    if (saveT == 1){
        FILE *outp1 = fopen("T.dat", "wb");
        for (j = 0;j < numPt_LL;j++) {
            for (i = 0; i < numSurfacePt10; i++){
                fwrite(&T[i][j], sizeof(double), 1, outp1);
            }
        }
        fclose(outp1);
    }

    //matrix multiplication H[187][110] = H0[187][144] * T[144][110]
    double **H_LL; //H matrix respect to lat/lon 187x110
    H_LL = (double**)calloc(numBins,sizeof(double*));//
    for (i = 0; i < numBins; i++){//
        H_LL[i] = (double*)calloc(numPt_LL,sizeof(double));
    }
    for (i = 0; i < numBins; i++){
        for (j = 0; j<numPt_LL; j++){
            H_LL[i][j]=0;
            for (k=0; k<numSurfacePt10; k++){
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
    for (i = 0; i < numBins; i++){
        free(H0[i]);
    }
    free(H0);

    for (i = 0; i < numBins; i++){
        free(H_LL[i]);
    }
    free(H_LL);

    free(H0_lat_vec); free(H0_lon_vec);
    free(bi_index1); free(indexLL);

    for (i = 0; i < numSurfacePt1; i++){
        free(M[i]);
    }
    free(M);

    for (i = 0; i < numSurfacePt10; i++){
        free(T[i]);
    }
    free(T);
    end =clock();
    printf("transformation time: %f seconds\n", (end-start)/CLOCKS_PER_SEC);
}


void ddm_mapSurfaceToDDM(void) {
    // creates a DDM by looping over the surface and adding the scattered
    // power/complex value to a particular DDM bin
    _ddm_zero( DDM, ddm.numBins );
    for(int i=0; i<surface.numGridPts; i++){
        if (surface.data[i].bin_index >= 0){
            DDM[surface.data[i].bin_index] += surface.data[i].total;
        }
    }
}

void ddm_mapRegionToDDM(void) {
    // Averages markers to produce DDM
    _ddm_zero( DDM,      ddm.numBins );
    _ddm_zero( DDM_temp, ddm.numBins );
    for(int i=0; i<surface.numGridPts; i++){
        if (surface.data[i].bin_index >= 0){
            //DDM_temp[surface.data[i].bin_index] += 1;
            DDM[surface.data[i].bin_index]  = surface.data[i].regionMarker;
        }
    }
    //for (int i = 0; i < ddm.numBins; i++)
    //    if( cabs(DDM_temp[surface.data[i].bin_index]) > 0 )
    //        DDM[surface.data[i].bin_index] /= DDM_temp[surface.data[i].bin_index];
}

void ddm_setSingleBin(int dopplerBin, int delayBin) {
    // test case for debugging.  Just sets single bin of DDM to one
    // and all the others zero  (good way to plot the amb func)

    _ddm_zero( DDM, ddm.numBins );
    if((dopplerBin < ddm.numDoppBins) && (dopplerBin >= 0) &&
       (delayBin < ddm.numDelayBins) && (delayBin >= 0 ))
        DDM[DDMINDEX(dopplerBin, delayBin)] = 1;
    else{
        fprintf(errPtr,"Error: requested bin out of range in ddm_setSingleBin\n");
        exit(1);
    }
}

void ddm_setBox(int centerDopplerBin, int centerDelayBin, int dopplerHalfWidth, int delayHalfWidth) {
    // test case for debugging.  Just sets a rect region of DDM to one
    // and all the others zero  (for use with ddm_mapDDMToSurface)

    for (int l=0; l < ddm.numDelayBins; l++) {
        for (int k=0; k < ddm.numDoppBins; k++) {
            if( ( abs(k - centerDopplerBin) < dopplerHalfWidth ) && ( abs(l - centerDelayBin) < delayHalfWidth )  )
                DDM[DDMINDEX(k, l)] = 1;
        }
    }
}

void ddm_mapDDMToSurface(void) {
    // performs the "inverse" of ddm_mapSurfaceToDDM, for testing the spatial
    // extent of certain DDM bins.  Use ddm_setSingleBin or ddm_setBox
    // before calling this.

    for(int i=0; i<surface.numGridPts; i++){
        if (surface.data[i].bin_index >= 0){
            surface.data[i].total = DDM[surface.data[i].bin_index];
        }else
            surface.data[i].total = 0;
    }

}


/****************************************************************************/
//  FFTs for convolution (uses FFTW)
/****************************************************************************/

void ddm_convolveFFT(int ambFuncType){   //type=2
    // convolves the working DDM with the ambiguity function.  Different
    // versions of the amb func are used for different purposes
    ddm_fft();
    switch (ambFuncType) {
        case 0:  for(int i = 0; i<ddm.numBins; i++) { DDM[i] = DDM_amb[i]  * DDM[i]; } break; //FFT DDM_amb version
        case 1:  for(int i = 0; i<ddm.numBins; i++) { DDM[i] = DDM_amb1[i] * DDM[i]; } break; //sqrt DDM_amb version
        case 2:  for(int i = 0; i<ddm.numBins; i++) { DDM[i] = DDM_amb2[i] * DDM[i]; } break; //mag squared DDM_amb version
        default: fprintf(errPtr,"Error: bad ambFuncType in ddm_convolveFFT"); break;
    }
    ddm_ifft();
}

void ddm_convolveH_FFT(int ambFuncType){
    // convolves the working DDM with the ambiguity function.  Different
    // versions of the amb func are used for different purposes
    ddm_h_fft();
    switch (ambFuncType) {
        case 0:  for(int i = 0; i<ddm.numBins; i++) { H[i] = DDM_amb[i]  * H[i]; } break;
        case 1:  for(int i = 0; i<ddm.numBins; i++) { H[i] = DDM_amb1[i] * H[i]; } break;
        case 2:  for(int i = 0; i<ddm.numBins; i++) { H[i] = DDM_amb2[i] * H[i]; } break;
        default: fprintf(errPtr,"Error: bad ambFuncType in ddm_convolveHFFT"); break;
    }
    ddm_h_ifft();
}


void ddm_fft(void){  fftw_execute(ddm.FFTWPLAN); }

void ddm_ifft(void){ fftw_execute(ddm.IFFTWPLAN); ddm_scale( (1.0/ddm.numBins)); }

void ddm_h_fft(void){  fftw_execute(h.FFTWPLAN); }

void ddm_h_ifft(void){ fftw_execute(h.IFFTWPLAN); ddm_h_scale( (1.0/ddm.numBins)); }



void ddm_fftshift(void){
    circshift(DDM_temp, DDM , ddm.numDoppBins,  ddm.numDelayBins, (ddm.numDoppBins/2),  (ddm.numDelayBins/2));
    for(int i = 0; i<ddm.numBins; i++) { DDM[i] = DDM_temp[i];}
}

void circshift(DDMtype *out, const DDMtype *in, int xdim, int ydim, int xshift, int yshift){
    int ii,jj;
    for (int i =0; i < xdim; i++) {
        ii = (i + xshift) % xdim;
        for (int j = 0; j < ydim; j++) {
            jj = (j + yshift) % ydim;
            out[ii * ydim + jj] = in[i * ydim + j];
        }
    }
}

/****************************************************************************/
//  Ambiguity Function (analytic, unfiltered)
//  based on eqns (20) and (21) [Zavorotny & Voronovich 2000]
/****************************************************************************/

void ddm_initACF(void){
    //added by Feixiong
    // load PRN ACF matrix
    FILE *file;
    char *ACF_filename = PRN_ACF_FILE;
    file = fopen(ACF_filename,"rb");
    if (file == NULL){
        printf("fail to open PRN ACF file\n");
        exit(1);
    }
    fread(PRN_ACF,sizeof(double),32736,file);
    fclose(file);
}

void ddm_initAmbFuncBuffers(int prn_code) {
    // called from ddm_initialize

    // generate amb func and align it to first bin, store a copy in temp buffer
    ddm_genAmbFunc(prn_code);  //computer ambiguity function values for each delay/Doppler bin and store in DDM[]
    ddm_fftshift();   //circular shift DDM[] in delay and Doppler dimension (why?)
    for(int i = 0; i<ddm.numBins; i++) {  DDM_temp[i] = DDM[i]; }

    // save FFT'd version in DDM_amb buffer
    ddm_fft();
    for(int i = 0; i<ddm.numBins; i++) {  DDM_amb[i] = DDM[i]; }

    // save a normalized, FFT'd version in DDM_amb1 buffer
    for(int i = 0; i<ddm.numBins; i++) {  DDM_amb1[i] = sqrt(cabs(DDM_amb[i])); }

    // take mag squared, FFT it, and save it in DDM_amb2 buffer
    for(int i = 0; i<ddm.numBins; i++) {  DDM[i] = DDM_temp[i]; }
    ddm_magSqr();
    ddm_fft(); // why?
    for(int i = 0; i<ddm.numBins; i++) {  DDM_amb2[i] = DDM[i]; }  //use this one

    // reset DDM buffer to zero when we are done
    _ddm_zero( DDM, ddm.numBins );
}

void ddm_genAmbFunc(int prn_code){
    // generate (non-squared) amb func. centered in DDM buffer
    double dtau_s, dfreq_Hz;

    double cohIntTime_s  = ddm.cohIntegrationTime_s;  // 0.001
    double tauChip_s     = 1 / ddm.chipsPerSec; // 1/1.023e6
    int centerDelayBin   = (int) floor(ddm.numDelayBins/2); //200
    int centerDopplerBin = (int) floor(ddm.numDoppBins/2);  //200

    for (int l=0; l < ddm.numDelayBins; l++) {  //l=0:399
        dtau_s =  (l - centerDelayBin) * ddm.delayRes_chips * tauChip_s; //(l-200)*0.05*tauChip_s
        for (int k=0; k < ddm.numDoppBins; k++) {
            dfreq_Hz = (k - centerDopplerBin) * ddm.dopplerRes_Hz;
            //DDM[DDMINDEX(k,l)] = lambda(dtau_s,tauChip_s,cohIntTime_s) * S(dfreq_Hz,cohIntTime_s);
            DDM[DDMINDEX(k,l)] = lambda_prn(prn_code, dtau_s,tauChip_s,cohIntTime_s) * S(dfreq_Hz,cohIntTime_s);
        }
    }
}

double lambda( double dtau_s, double tauChip_s, double cohIntTime_s ) {
    //perfect trianglular on ZV paper
    return (fabs(dtau_s) <= (tauChip_s*(1+tauChip_s/cohIntTime_s))) ?
           (1 - fabs(dtau_s)/tauChip_s) : -tauChip_s/cohIntTime_s;
}

double lambda_prn(int prn_code, double dtau_s, double tauChip_s, double cohIntTime_s ){
    //ACF for each PRN
    double bin;
    int bin0, bin1;
    double ACF, ACF0, ACF1;

    bin = dtau_s/tauChip_s; //-10chips to 9.95chips
    bin0 = (int)floor(bin);
    bin1 = (int)floor(bin)+1;
    ACF0 = PRN_ACF[1023*(prn_code-1)+511+bin0];
    ACF1 = PRN_ACF[1023*(prn_code-1)+511+bin1];
    ACF = (bin1-bin)*ACF0+(bin-bin0)*ACF1; //linear interpolation

    return ACF;
    //printf("%f %f %f \n",ACF,ACF0, ACF1);

}

complex double S(double dfreq_Hz, double cohIntTime_s) {
    double x       = dfreq_Hz*pi*cohIntTime_s;
    double ang_rad = -1*pi*dfreq_Hz*cohIntTime_s;
    return (x==0) ? 1 : (sin(x)/x) * (cos(ang_rad) + I*sin(ang_rad));
}


/****************************************************************************/
//  Thermal Noise
/****************************************************************************/

void ddm_initThermalNoise(struct metadata meta){
    // called from ddm_initialize

    // get thermal noise params from config file
    double k           = -228.599167840;	// Boltzman, dBW/K/Hz
    double temp_K      = meta.temp_K;
    double noiseFig_dB = meta.noiseFigure_dB;
    double Ti          = 0.001;
    double BW_Hz       = 1/Ti;

    double noisePower_abs =  pow(10,(k/10)) * temp_K * BW_Hz * pow(10,(noiseFig_dB/10)); //
    double noisePower_dBW = 10*log10(noisePower_abs);

    ddm.thermalNoisePwr_abs = sqrt(noisePower_abs);  //sigma of noise
    ddm.thermalNoisePwr_abs = ddm.thermalNoisePwr_abs/1000; // because of incoherent summation
}

void ddm_addGaussianNoise(void){
    // Add complex, Gaussian noise to DDM (colored with ambiguity function)
    // Currently, the convolution is done separately from the DDM/speckle convolution
    // due to the normalization.
    printf("Add thermal noise\n");
    for(int i = 0; i<ddm.numBins; i++) {  DDM_temp[i] = DDM[i]; DDM[i] = 0;  }
    ddm_addWhiteGaussianNoise(ddm.thermalNoisePwr_abs);
    //ddm_convolveFFT(1);

    //printf("thermalNoisePwr_abs = %e\n",ddm.thermalNoisePwr_abs);
    //printf("DDM noise = %e\n",DDM[100]);

    for(int i = 0; i<ddm.numBins; i++) { DDM[i] += DDM_temp[i]; }
}

void ddm_addWhiteGaussianNoise(double sigma){
    // Add complex, *white* Gaussian noise to DDM (mean = 0 and variance = sigma^2)
    #define RAND_GEN_TYPE 1
    double U1,U2;
    double R,T1,T2;

    #if (RAND_GEN_TYPE == 1) // standard Box-Muller Method
        for(int i = 0; i<ddm.numBins; i++) {
            U1 = uniformRandf(); //uniform random number (0,1]
            U2 = uniformRandf();
            R  = sigma*sqrt(-2.0*log(U1))*(1/sqrt(2));
            T1 = R*cos(2*pi*U2);
            T2 = R*sin(2*pi*U2);
            //DDM[i] += T1 + I*T2;
            DDM[i] += T1*T1 + T2*T2;  //changed by Feixiong
        }
    #endif

    #if (RAND_GEN_TYPE == 2) // polar form (supposedly faster ...)
        double S;
        for(int i = 0; i<ddm.numBins; i++) {
            do {
                U1=2 * uniformRandf() - 1; /* U1=[-1,1] */
                U2=2 * uniformRandf() - 1; /* U2=[-1,1] */
                S=U1 * U1 + U2 * U2;
            } while(S >= 1);
            R  = sigma*sqrt(-2.0*log(S) / S)*(1/sqrt(2));
            T1 = R*U1;
            T2 = R*U2;
            DDM[i] += T1 + I*T2;
        }
    #endif
}

// uniform random number (0,1]. never = 0 so that we can safely take log of it
double uniformRandf( void ) { return ((double)rand() + 1.0)/((double)RAND_MAX + 1.0); }

void ddm_addRandomPhase(void){
    // adds a uniform random phase to the DDM
    // (for testing a poor man's speckle noise)
    for(int i = 0; i<ddm.numBins; i++) { DDM[i] *= cexp(I * uniformRandf() * 2 * pi); }
}

void testNoisePowerLevels(void){
    // used internally for debugging noise powers
    printf("Noise Power 1 is:  %f (dBW)\n", 20*log10(ddm.thermalNoisePwr_abs) );

    _ddm_zero( DDM,      ddm.numBins );
    ddm_addGaussianNoise();
    printf("Noise Power 2 is:  %f (dBW)\n", 20*log10(ddm_getRMS()) );

    ddm_resetRunningAvg();
    for( int i=0; i<100; i++){
        _ddm_zero( DDM,      ddm.numBins );
        ddm_addGaussianNoise();
        ddm_magSqr();
        ddm_addToRunningAvg();
    }
    ddm_getRunningAvg();
    printf("Noise Power 3 is:  %f (dBW)\n", 10*log10(ddm_getRMS()) );
}


/****************************************************************************/
//  Operations on DDMs:  These are various operations performed on the DDM.
//  Note that there are several DDM buffers.  "DDM" is the working DDM buffer
//  where most of the work happens, but there is "DDM_avg" for running averages
//  (i.e. non-coherent integration) & "DDM_store" for storing a DDM for later
//  use.
/****************************************************************************/

void ddm_normalize(void){
    double max = ddm_getMax();
    double min = ddm_getMin();
    if(max == 0) return;
    for (int i = 0; i < ddm.numBins; i++){
        DDM[i] = (cabs(DDM[i])-min)/(max-min);
    }
}

double ddm_getMax(void){
    double val, max = creal(DDM[0]);
    for (int i = 0; i < ddm.numBins; i++){
        val = creal(DDM[i]);
        if (val > max){ max = val;}
    }
    return max;
}

double ddm_getMin(void){
    double val, min = creal(DDM[0]);
    for (int i = 0; i < ddm.numBins; i++){
        val = creal(DDM[i]);
        if (val < min){ min = val;}
    }
    return min;
}

double ddm_getRMS(void)       { double val=0; for (int i = 0; i < ddm.numBins; i++){val += pow(cabs(DDM[i]),2);} return sqrt(val/ddm.numBins); }
double ddm_integrate(void)    { double val=0; for (int i = 0; i < ddm.numBins; i++){val += cabs(DDM[i]);} return val; }
void ddm_mag(void)            { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = cabs(DDM[i]); } }
void ddm_angle(void)          { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = carg(DDM[i]); } }
void ddm_real(void)           { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = creal(DDM[i]); } }
void ddm_imag(void)           { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = cimag(DDM[i]); } }
void ddm_magSqr(void)         { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = pow(cabs(DDM[i]),2);} }
void ddm_sqrt(void)           { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = csqrt(DDM[i]);} }
void ddm_convertTodB(void)    { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = (DDM[i] == 0) ? ddm.refPower_dB : 10*log10(cabs(DDM[i])); } }
void ddm_scale(double c)      { for(int i = 0; i<ddm.numBins; i++) { DDM[i] *= c; } }
void ddm_h_scale(double c)      { for(int i = 0; i<ddm.numBins; i++) { H[i] *= c; } }
void ddm_store(void)          { for(int i = 0; i<ddm.numBins; i++) { DDM_store[i] = DDM[i]; } }
void ddm_restore(void)        { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = DDM_store[i]; } }
void ddm_addToRunningAvg(void){ for(int i = 0; i<ddm.numBins; i++) { DDM_avg[i] = DDM_avg[i] + DDM[i]; } DDM_avgCount++; }
void ddm_getRunningAvg(void)  { for(int i = 0; i<ddm.numBins; i++) { DDM[i] = DDM_avg[i] / DDM_avgCount; } }
void ddm_resetRunningAvg(void){ DDM_avgCount = 0; _ddm_zero( DDM_avg,  ddm.numBins ); }
int  ddm_checkNAN(void)       { for(int i = 0; i<ddm.numBins; i++) { if( isnan(DDM[i]) ) return 1; } return 0; }

/****************************************************************************/
//  Save DDM to binary file or .png image
/****************************************************************************/

void ddm_save(struct metadata meta, struct DDMfm *ddm_fm, int realOrComplex){
    int startDelay_bin    = meta.resample_startBin[0];//0
    int startDoppelr_bin  = meta.resample_startBin[1];//100
    int resDelay_bins     = meta.resample_resolution_bins[0];//5
    int resDoppelr_bins   = meta.resample_resolution_bins[1];//20
    int numDelayBins      = meta.resample_numBins[0];//17
    int numDoppelrBins    = meta.resample_numBins[1];//11
    double val;
    complex double valc;

    int i = 0;
    switch(realOrComplex){  //case=1
        case 1: // real
            //for(int i = 0; i<ddm.numBins; i++) {
            //    val = creal(DDM[i]);
            //    fwrite(&val, 1, sizeof(double), outp);
            //}

            for (int k=startDoppelr_bin; k < (numDoppelrBins*resDoppelr_bins + startDoppelr_bin); k+=resDoppelr_bins) {
                for (int l=startDelay_bin; l < (numDelayBins*resDelay_bins + startDelay_bin); l+=resDelay_bins) {
                    ddm_fm->data[i].power = creal(DDM[DDMINDEX(k,l)]);
                    ddm_fm->data[i].delay = (l - meta.specular_delayBinIdx) * meta.delayRez_chips;
                    ddm_fm->data[i].Doppler = (k - meta.specular_dopplerBinIdx) * meta.dopplerRes_Hz;
                    i=i+1;
                }
            }
            break;

        case 2: // complex under construction
            for(int i = 0; i<ddm.numBins; i++) {
                valc = DDM[i];
            }
            break;
        default:
            fprintf(errPtr,"Error: bad realOrComplex in ddm_saveToFile");
            exit(0);
    }
}


double getImagePixelValue( int x, int y ){
    return creal(DDM[DDMINDEX(x, y)]);
}
/*
void ddm_save2PNG(const char* filename, int plotType, int limType, double min, double max ) {
    unsigned height = ddm.numDelayBins;
    unsigned width  = ddm.numDoppBins;
    unsigned char* image = malloc(width * height * 4);
    unsigned x, y;
    unsigned char color[3];
    double m1,m2;

    // save a copy in temp buffer
    for (int i = 0; i < ddm.numBins; i++){  DDM_temp[i] = DDM[i]; }

    // select component of DDM to plot
    switch (plotType){
        case 1: ddm_mag(); break;
        case 2: ddm_angle(); break;
        case 3: ddm_real(); break;
        case 4: ddm_imag(); break;
        case 5: ddm_convertTodB(); break;
        case 6:
            for(y = 0; y < height; y++)
                for(x = 0; x < width; x++)
                    DDM[DDMINDEX(x, y)] = x;
            break;
        case 7:
            for(y = 0; y < height; y++)
                for(x = 0; x < width; x++)
                    DDM[DDMINDEX(x, y)] = y;
            break;
        case 8: _ddm_zero( DDM, ddm.numBins ); break;
        default:
            fprintf(errPtr,"Error: Bad plotType in ddm_save2PNG");
            exit(0);
    }

    // select color limits to use (specified limits or auto limits)
    switch (limType) {
        case 1:  m1 = min; m2 = max; break;
        case 2:  m1 = ddm_getMin(); m2 = ddm_getMax(); break;
        default:  fprintf(errPtr,"Error: Bad limType in ddm_save2PNG"); exit(0);
    }

    // plot to png
    for(y = 0; y < height; y++)
        for(x = 0; x < width; x++) {
            image_mapFloat2RGB( m1, m2, getImagePixelValue(x,y), color );
            image[4 * width * y + 4 * x + 0] = color[0];
            image[4 * width * y + 4 * x + 1] = color[1];
            image[4 * width * y + 4 * x + 2] = color[2];
            image[4 * width * y + 4 * x + 3] = 255;
        }
    image_encodePNG(filename, image, width, height);
    free(image);

    // restore original DDM value from temp buffer
    for (int i = 0; i < ddm.numBins; i++){  DDM[i] = DDM_temp[i]; }
}
*/

/****************************************************************************/
// internal functions for alloc/dealloc DDM buffers
/****************************************************************************/

void _ddm_alloc( DDMtype **ptr, int numBins ) {
    *ptr = (DDMtype *) calloc( numBins, sizeof(DDMtype) );
    if (*ptr==NULL) {
        fprintf(errPtr,"Error allocating memory for DDM \n");
        exit (1);
    }
}

void _ddm_free( DDMtype *ptr ) {
    if (ptr!=NULL) {
        free(ptr);
    }
}

void _ddm_zero( DDMtype *ptr, int numBins  ) {
    if (ptr==NULL) {
        printf ("Error: Null pointer to _ddm_zero \n ");
        exit (1);
    }
    for(int i = 0; i<numBins; i++) {
        ptr[i] = 0 + I*0;
    }
}


