// be careful about SURF index (details in my Thesis)
// 1. wind_initialize
// 2. wind_interpolated
// 3. surface_calcGeomOverSurface
// 4. gnssr.h SURFINDEX = i*numGridPtsY + j

// be careful about sxangle or sxangle_rad or ax_angle (incidence angle)
// 1. geom.c
// 2. surface.c
// 3. GMF.c

#ifndef CFORWARDMODEL_GNSSR_H
#define CFORWARDMODEL_GNSSR_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <sys/time.h>
#include "forwardmodel.h"

#define speedlight 2.99792458e8 // speed of light (m/s)
#define pi 3.141592653589793
#define L1 1575.42e6
#define L1_WAVELENGTH 0.1904
#define GPS_CA_CHIP_LENGTH 293.05   // meters
#define R2D 57.295779513082320876798154814105  // radians to degrees
#define D2R 1.0/R2D // degrees to radians
#define chipRate_cs 1.023e6 // chips per second
#define omega0 7.2921158553e-5 //earth rotation in rad/s

#define upsamplingFactor 10  // upsamplingFactor of receiver antenna pattern

FILE *outputPtr, *errPtr, *consolePtr;


//***** Global variable *******************************************************************/
//grid transformation
int bi_index0[14400][4];  // for 1km in SURF
double bi_weight0[14400][4]; // for 1km in SURF
int sp_index; //index of specular point on the 120 X 120 grid (global variable)
double cyg_R2; // CYGNSS square of reflectivity at SP (very slight change on incidence angle)

//******************************************************************************/
//GMF.c
double GMF_converWindToMSS( double windSpeedMag_ms, double R2);
void GMF_init(double sp_sxangle);
double GMF_getCoef(char [], double inc_angle_degree);
double get_gfds(double u);
double get_gyslf(double u);
double get_g(double u);
double get_dgfds(double u);
double get_dgyslf(double u);
double get_dg(double u);
double get_dmdx_GMF(double R2, double u);

struct
{
    double a0,a1,a2,b0,b1,b2,c0,c1,c2,d0,d1;
    double trans_u_fds, trans_u_yslf;
} GMF;

// Geometry (geom.c)

// global structure for holding Rx\Tx\Sx geometry
typedef struct {

    // Rx\Tx\SX positions & velocities in specular frame
    double rx_pos[3], tx_pos[3], sx_pos[3];
    double rx_vel[3], tx_vel[3], sx_vel[3];

    // Rx\Tx\SX positions & velocities in ECEF frame
    double rx_pos_ecef[3], tx_pos_ecef[3], sx_pos_ecef[3];
    double rx_vel_ecef[3], tx_vel_ecef[3], sx_vel_ecef[3];

    // Rx\Tx\SX positions & velocities in LLH coordinates
    double rx_pos_llh[3], tx_pos_llh[3], sx_pos_llh[3];
    double rx_vel_llh[3], tx_vel_llh[3], sx_vel_llh[3];

    // attitude dynamic
    double sc_att_rad[3];

    // Rx\Tx to Sx parameters
    double rx_range_m, tx_range_m;
    double rx_speed_ms, tx_speed_ms, sx_speed_ms;
    double rx_angle_rad, tx_angle_rad, sx_angle_rad;
    double rx_LOS_velocity, tx_LOS_velocity;
    double dopplerRx_Hz, dopplerTx_Hz;
    double specularDistance_m, specularDoppler_Hz;

    // transforms between frames
    double ECEF_TO_SPEC_FRAME[9], SPEC_TO_ECEF_FRAME[9];
    double SPEC_TO_RX_ORB_FRAME[9], SPEC_TO_TX_ORB_FRAME[9];
    double SPEC_TO_RX_ORB_FRAME_AARON[9], SPEC_TO_TX_ORB_FRAME_AARON[9]; //not used
    double ECEF_TO_RX_ORB_FRAME[9], ECEF_TO_TX_ORB_FRAME[9];
    double ORB_TO_RX_BODY_FRAME[9], SPEC_TO_RX_BODY_FRAME[9];

    // various angles between things
    double angleSxFromRx_rad[2], angleSxFromTx_rad[2];
    double angleTxFromRx_rad[2], angleRxFromTx_rad[2];
    double distance_rx2tx_m;
    double doppler_rx2tx_Hz;

    // antenna gains
    double antennaRxGainAtSx_RHCP_dB, antennaRxGainAtSx_LHCP_dB;
    double antennaTxGainAtSx_RHCP_dB, antennaTxGainAtSx_LHCP_dB;
    double antennaRxGainAtTx_RHCP_dB, antennaRxGainAtTx_LHCP_dB;
    double antennaTxGainAtRx_RHCP_dB, antennaTxGainAtRx_LHCP_dB;

    // path loss
    double pathloss_Rx2Sx_dB, pathloss_Tx2Sx_dB, pathloss_Rx2Tx_dB, pathloss_Sx_dB;

    // quantify motion of specular bin on DDM
    double ddm_delayRate_chips_per_s, ddm_dopplerRate_kH_per_s;
    double directPathDistance_m,directPathDoppler_Hz,pathDifference_chips;

    // transmit power (this doesnt really belong here ...)
    double txTransmitPower_dB;

} orbitGeometryStruct;

// geometryData contains a list of geometries to be run
typedef struct {
    orbitGeometryStruct *g;
    int numGeometries, geomStartIdx, geomEndIdx;
} geometryData;

void geom_printToLog(FILE *outputPtr, int index, orbitGeometryStruct *g); // keep it as it may be useful

void geom_calculateSecondaryGeometry( orbitGeometryStruct *g );
void geom_initialize(geometryData *gd, struct Geometry geom);
void geom_getRelativeAngleInFrame( double origin[3], double pos[3] , double M[9], double angles_rad[2] );
orbitGeometryStruct *geom_getOrbitData(geometryData *gd, int geomIdx );

void getECEF2SpecularFrameXfrm( double rx_pos_ecef[3], double tx_pos_ecef[3], double sx_pos_ecef[3], double xfrmMatrix[9]);
void getECEF2OrbitFrameXfrm( double sat_pos_ecef[3], double sat_vel_ecef[3], double xfrmMatrix[9] ); // by Feixiong
void getOrbit2BodyFrameXfrm (double pitch, double roll, double yaw, double xfrmMatrix[9]); // by Feixiong
//******************************************************************************/
// specular.c

int solveSpecularPtPosition(double rx_pos_ecef[3], double tx_pos_ecef[3],
                            double sx_pos_ecef[3], double r_init_m, int maxIterations);

void solveSpecularPt(double rx_pos_ecef[3], double tx_pos_ecef[3], double sx_pos_ecef[3],
                     double rx_vel_ecef[3], double tx_vel_ecef[3], double sx_vel_ecef[3], double *sx_valid, double r_init_m );

//******************************************************************************/
//coord.c
void cart2sph( double x[3], double y[3] );
void wgslla2xyz(double x_llh[3], double x_ecef[3]);
void wgsxyz2enu(double p_e[3] ,double x_llh[3], double x_enu[3]);
void wgsxyz2lla( double *x_ecef, double *x_lla );

//******************************************************************************/
// wind.c

typedef struct{
    double windSpeed_U10_ms, windSpeed_V10_ms;
    double rainRate_mmhr, freezingHeight_m;
    double windSpeed_ms, windDir_rad;
    double mss_perp, mss_para, mss_x, mss_y, mss_b;
    int x,y;
    double lat_deg,lon_deg;
    int rowIdx, colIdx;
    int domainIdx;
    double weight;
} windFieldPixel;

typedef struct{
    int  numGridPtsX, numGridPtsY, numGridPts;  // num pixels
    double resolutionX_m, resolutionY_m;

    windFieldPixel *data;   //structure in structrure
    int type; // 1 = uniform, 2 = file

    // for specular locations
    int locType, locLoaded, locNumPts, locCurrentPt, locStartIdx, locEndIdx;

} windField;

void wind_interpolate(windField *wf,struct Geometry geom, struct inputWindField iwf, double grid_resolution);
void wind_initialize(windField *wf, struct metadata meta, struct Geometry geom, struct inputWindField iwf);
void wind_converWindToMSS( double windSpeedMag_ms, double windDirectionAngle_deg, double mss[5] );
void wind_convertWindXY2MagDir( double x, double y, double *mag, double *dir_rad );
void wind_getWindFieldAtXY( windField *wf, double x_m, double y_m, windFieldPixel *value );

//******************************************************************************/
// Surface

// each grid point (i.e. pixel) on the surface has the following parameters:
typedef struct{

    // surface gridding params
    double pos_spec[3], pos_ecef[3], pos_llh[3], normal_spec[3], area_m2, distance_m;
    double position[3];
    int i,j;

    // these params are filled in using surface_calcGeom
    double delay_s, doppler_Hz;
    double q[3];
    double sx_angle_rad;
    double powerFactor, pathloss;

    // antenna params
    double antennaGainRx_abs, antennaGainTx_abs;
    double angleSxFromRx_theta_rad, angleSxFromRx_phi_rad, angleSxFromTx_theta_rad, angleSxFromTx_phi_rad;

    // bin_index is the index into the DDM buffer.  It is calculated using ddm_binSurface
    int bin_index;

    complex double phaseShiftFactor0, phaseShiftFactor1, phaseShiftFactor2;
    complex double total, total_dP;

    // wind field X,Y locations relative to specular on surface
    double windFieldLocation_x_m, windFieldLocation_y_m, windFieldLocationRange_km;

    // solving for rain attenuations needs to know elevation angles across surface
    double rx_elevationAngle_rad, tx_elevationAngle_rad, rainAtten_abs;

    // sigma0 and components of it for debugging
    double sigma0, sigma0_R2, sigma0_P, sigma0_dP, sigma0_Q4;

    // mask portion of surface outside complete annuli
    double mask;

} surfacePixel;

// global that holds surface parameters and data
struct {
    int numGridPtsX, numGridPtsY, numGridPts;
    double width_m, height_m, resolution_m;
    int surfaceCurvatureType;
    int speckleType;

    // X,Y coordinates of specular point in meteres in the master wind field
    double specularLoactionX_m, specularLoactionY_m;
    int specularGridPt_x_idx, specularGridPt_y_idx;

    surfacePixel *data;
    windFieldPixel *windData;

    int rainOnOff;

} surface;

#define SURFINDEX(XINDEX,YINDEX) XINDEX * surface.numGridPtsY + YINDEX

// surface.c prototypes

void surface_initialize(struct metadata meta);
void surface_cleanup(void);

void surface_calcGeomOverSurface(orbitGeometryStruct *geometry, int fourCornerTest, struct powerParm pp);
void surface_getScatteringVector(double TSx_unit[3], double RSx_unit[3], double PUT[3], double q_vec_new[3]);
void surface_createSurfaceMask(void);

void surface_calcSigma0OnSurface(int windModelType);
void surface_calcRainAttenOnSurface(void);
void surface_composeTotalScatPowrOnSurface(int type);
complex double reflectionCoef(double Sxangle);

void surface_loadSurfWindfield(windField *wf, int wfNum);

void surface_initSpeckle(void);
void surface_updateSpeckle(void);

void surface_resetToZero(void);

double getRainAtten_abs( double theta1_rad, double theta2_rad, double h_km, double R_mmhr);

void surface_saveWindToFile(void);
void surface_saveDopplerToFile(void);
void surface_saveDelayToFile(void);
//******************************************************************************/
// grid.c

void grid_getGridPt_sphericalEarth(int i, int j, double *sx_pos, double *grid_pos, double magSx);
void grid_getGridPt_flatEarth(int i, int j, double *sx_pos, double *grid_pos);

//******************************************************************************/
// DDM

// global for holding DDM data
struct
{
    int numDelayBins, numDoppBins, numBins;
    double dopplerRes_Hz, delayRes_chips;
    int delayOffset_bins, dopplerOffset_bins;
    double chipsPerSec;
    double refPower_dB;
    fftw_plan FFTWPLAN,IFFTWPLAN;
    double cohIntegrationTime_s;
    double thermalNoisePwr_abs;
} ddm;

struct
{
    int numDelayBins, numDoppBins, numBins;
    double dopplerRes_Hz, delayRes_chips;
    int delayOffset_bins, dopplerOffset_bins;
    double chipsPerSec;
    double refPower_dB;
    fftw_plan FFTWPLAN,IFFTWPLAN;
    double cohIntegrationTime_s;
    double thermalNoisePwr_abs;
} h;

// DDM buffers
typedef complex double DDMtype;
DDMtype *DDM,*DDM_temp,*DDM_avg,*DDM_amb,*DDM_amb1,*DDM_amb2,*DDM_store, *H;
int DDM_avgCount;

#define DDMINDEX(DOPPLERBINIDX,DELAYBINIDX) DOPPLERBINIDX * ddm.numDelayBins + DELAYBINIDX

double PRN_ACF[32736]; // PRN aotocorrelation matrix

// ddm.c prototypes

void ddm_initialize(struct metadata meta);
void ddm_cleanup(void);

void ddm_Hmatrix(struct metadata meta, struct inputWindField iwf, struct Jacobian *jacob);
void ddm_mapSurfaceToDDM(void);
void ddm_binSurface(void);

void ddm_convolveFFT(int ambFuncType);
void ddm_convolveH_FFT(int ambFuncType);
void ddm_fft(void);
void ddm_ifft(void);
void ddm_h_fft(void);
void ddm_h_ifft(void);
void ddm_fftshift(void);
void circshift(DDMtype *out, const DDMtype *in, int xdim, int ydim, int xshift, int yshift);

void ddm_initACF(void);
void ddm_genAmbFunc(int prn_code);
void ddm_initAmbFuncBuffers(int prn_code);
double lambda( double tau_s, double tauChip_s, double cohIntTime_s );
double lambda_prn(int prn_code, double tau_s, double tauChip_s, double cohIntTime_s );
complex double S(double dfreq_Hz, double cohIntTime_s);

void ddm_save(struct metadata meta, struct DDMfm *ddm_fm, int realOrComplex);

double ddm_getMax(void);
double ddm_getMin(void);

void ddm_mag(void);
void ddm_angle(void);
void ddm_real(void);
void ddm_imag(void);
void ddm_magSqr(void);
void ddm_addToRunningAvg(void);
void ddm_resetRunningAvg(void);
void ddm_convertTodB(void);
void ddm_getRunningAvg(void);
void ddm_normalize(void);
void ddm_sqrt(void);
void ddm_store(void);
void ddm_restore(void);
void ddm_scale(double c);
void ddm_h_scale(double c);
double ddm_getRMS(void);
double ddm_integrate(void);

void ddm_initThermalNoise(struct metadata meta);
double uniformRandf( void );
void ddm_addWhiteGaussianNoise(double sigma);
void ddm_addGaussianNoise(void);
void ddm_addRandomPhase(void);
void testNoisePowerLevels(void);

void ddm_resetToZero(void);
void _ddm_alloc( DDMtype **, int );
void _ddm_free( DDMtype * );
void _ddm_zero( DDMtype *, int );

//******************************************************************************/
// Antennas (antenna.c)

typedef struct{
    double az_deg, el_deg, gain_dB;
} antennaDataPixel;

typedef enum { CYGNSS_ZENITH_ANT, CYGNSS_NADIR_ANT, GPS_SAT_ANT } antennaType;
typedef enum { RHCP, LHCP } polarizationType;

struct{
    antennaDataPixel *data;
    antennaType type;
    double gain, excessGain;
    double numAz, numEl;
} antenna;

void antenna_initialize(struct powerParm pp);
double antenna_getGain_abs( antennaType at, polarizationType pt, double angles_rad[2] );

//******************************************************************************/
// math.c prototypes
void bubble(int *a,int n);
void bilinear_interp(double *x_vec, double *y_vec, int size_x, int size_y, double x, double y, int *bi_index, double *bi_weight, double resolution);
int find_nearest(double *vec, int size, double value);
double linear_interp( double a, double b, int direction, double time_01 );
void cubic_interpolation( double f0, double f1, double df0, double df1,
                          double t, double *ft, double *dft, double timeInterval_s);
void cubic_interpolation_3vector( double f0[3], double f1[3], double df0[3], double df1[3],
                                  double t, double ft[3], double dft[3], double timeInterval_s);

double cot(double z);
double sec(double z);
double csc(double z);

void RotationMatrix_X (double theta, double xfrmMatrix[9]);
void RotationMatrix_Y (double theta, double xfrmMatrix[9]);
void RotationMatrix_Z (double theta, double xfrmMatrix[9]);

void vector_cross_product(double a[3], double b[3],double c[3]);
void vector_unit(double *v, double *u);
double vector_dot_product (double *a, double *b);
double vector_norm(double v[3]);
void vector_subtract(double a[3], double b[3],double c[3]);
void vector_add(double a[3], double b[3],double c[3]);
void vector_scale( double a[3], double b[3], double s );
void vector_constrainToPlane( double a[3], double b[3], double c[3] );
void vector_orthoNorm( double a[3], double b[3] );

void matrix_multiply(unsigned rowsA,unsigned columnsA,unsigned columnsB,double *A,double *B,double *C);
void matrix_transpose(unsigned rows,unsigned columns,double *A,double *B);
void matrix_form3x3( double row1[3], double row2[3], double row3[3], double *M );
void matrixVectorMult3x3( double M[9], double x[3], double y[3] );
double matrix_det_3x3(double m[3][3]);
void matrix_scaleAdjoint_3x3(double a[3][3],double s,double m[3][3]);
void matrix_invert_3x3(double A[9], double invA[9]);

//******************************************************************************/

#define fatalError(...) { char str[1000]; sprintf(str, __VA_ARGS__); printf( "%s (%s in %s, line %d)\n", str, __func__, __FILE__, __LINE__); exit(EXIT_FAILURE); }

//******************************************************************************/
//debug.c
void printfGeometry(geometryData orbitGeometry);
void printGMF();

#endif //CFORWARDMODEL_GNSSR_H


