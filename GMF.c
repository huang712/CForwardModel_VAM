// Functions for modified wind-MSS model from CYGNSS GMF
// GMF is a global structure containing coefficients

#include "gnssr.h"

double GMF_converWindToMSS( double windSpeedMag_ms, double sp_sxangle, double sxangle){
    //Modified model from CYGNSS GMF
    //sxangle: elevation angle at current point
    //sp_sxangle: elevation angle at specular point
    //mss: mss_iso = mss_x = mss_y

    double u,R2,mss;
    double g;

    R2 = pow(cabs(reflectionCoef(sxangle)),2);
    u = windSpeedMag_ms;
    g = get_g(u);  // CYGNSS sigma0
    mss=R2/(2*g); // compute MSS

//    printf("R2 = %f\n",R2);
//    printf("ws = %f\n",u);
//    printf("sp_inc = %f, inc = %f\n",90-R2D*sp_sxangle, 90-R2D*sxangle);
//    printf("mss = %f\n",mss);

    return mss;
}


void GMF_init(double sp_sxangle){
    // read GMF coefficients from data
    // sp_sxangle: elevation angle in rad

    printf("Read CYGNSS GMF model\n");
    double inc_angle_degree = 90-R2D*sp_sxangle; // incidence angle
    GMF.a0=GMF_getCoef("a0",inc_angle_degree);
    GMF.a1=GMF_getCoef("a1",inc_angle_degree);
    GMF.a2=GMF_getCoef("a2",inc_angle_degree);
    GMF.b0=GMF_getCoef("b0",inc_angle_degree);
    GMF.b1=GMF_getCoef("b1",inc_angle_degree);
    GMF.b2=GMF_getCoef("b2",inc_angle_degree);
    GMF.c0=GMF_getCoef("c0",inc_angle_degree);
    GMF.c1=GMF_getCoef("c1",inc_angle_degree);
    GMF.c2=GMF_getCoef("c2",inc_angle_degree);
    GMF.d0=GMF_getCoef("d0",inc_angle_degree);
    GMF.d1=GMF_getCoef("d1",inc_angle_degree);
    GMF.trans_u_fds=GMF_getCoef("trans_ws_fds",inc_angle_degree);
    GMF.trans_u_yslf=GMF_getCoef("trans_ws_yslf",inc_angle_degree);

    //printf("inc angle = %f\n",inc_angle_degree);
    //printGMF();

}

double GMF_getCoef(char name[],double inc_angle_degree){
    FILE *file;
    double val_inc[90],val; // value for each incidence angle
    int index;
    char filename[1000]=GMF_PATH;
    strcat(filename,name);
    strcat(filename,".dat");
    file = fopen(filename,"rb");
    if (file == NULL){
        printf("fail to open GMF file");
        exit(1);
    }
    fread(val_inc,sizeof(double),90,file);
    index = (int) round(inc_angle_degree)-1;
    if (index<0){
        index=0;
    }
    val = val_inc[index];
    return val;

}

double get_g(double u){
    // compute CYGNSS DDMA from wind speed
    double g_fds, g_yslf, g;

    // compute g_fds and g_yslf
    g_fds = get_gfds(u); //FDS
    g_yslf = get_gyslf(u); //YSLF

    // compute combined g
    if (u <= 10){
        g=g_fds;
    }
    else if ((u > 10) && (u < 20)){
        g=(20.0-u)/10.0*g_fds+(u-10.0)/10.0*g_yslf;
    }
    else{
        g=g_yslf;
    }

    return g;
}
double get_gfds(double u){
    double a0,a1,a2,b0,b1,b2,trans_u_fds;
    double g_fds;
    a0 = GMF.a0; a1 = GMF.a1; a2 = GMF.a2;
    b0 = GMF.b0; b1 = GMF.b1; b2 = GMF.b2;
    trans_u_fds = GMF.trans_u_fds;
    if (u <= trans_u_fds){
        g_fds=a0+a1/u+a2/pow(u,2);
    }
    else{
        g_fds=b0+b1*u+b2*pow(u,2);
    }

    return g_fds;
}

double get_gyslf(double u){
    double c0,c1,c2,d0,d1,trans_u_yslf;
    double g_yslf;
    c0 = GMF.c0; c1 = GMF.c1; c2 = GMF.c2;
    d0 = GMF.d0; d1 = GMF.d1;
    trans_u_yslf = GMF.trans_u_yslf;

    if (u <= trans_u_yslf){
        g_yslf=c0+c1/u+c2/pow(u,2);
    }
    else{
        g_yslf=d0+d1*u;
    }

    return g_yslf;
}

double get_dmdx_GMF(double R2, double u){
    double g, dg, dmdx;
    g = get_g(u);
    dg = get_dg(u);
    dmdx = -R2/2 * 1/pow(g,2) *dg;
    return dmdx;
}

double get_dg(double u){
    double g_fds, g_yslf, dg_fds, dg_yslf, dg;
    g_fds = get_gfds(u);
    g_yslf = get_gyslf(u);
    dg_fds = get_dgfds(u);
    dg_yslf = get_dgyslf(u);

    if (u <= 10){
        dg=dg_fds;
    }
    else if ((u > 10) && (u < 20)){
        dg = -0.1*g_fds + (20-u)/10*dg_fds + 0.1*g_yslf + (u-10)/10*dg_yslf;
    }
    else {
        dg=dg_yslf;
    }
    return dg;
}

double get_dgfds(double u){
    // derivative of g_fds respect to wind speed
    double a0,a1,a2,b0,b1,b2,trans_u_fds;
    double dg_fds;
    a0 = GMF.a0; a1 = GMF.a1; a2 = GMF.a2;
    b0 = GMF.b0; b1 = GMF.b1; b2 = GMF.b2;
    trans_u_fds = GMF.trans_u_fds;

    if (u <= trans_u_fds){
        dg_fds = -a1/pow(u,2) - 2*a2/pow(u,3);
    }
    else{
        dg_fds = b1 + 2*b2*u;
    }
    return dg_fds;
}

double get_dgyslf(double u){
    // derivative of g_yslf respect to wind speed
    double c0,c1,c2,d0,d1,trans_u_yslf;
    double dg_yslf;
    c0 = GMF.c0; c1 = GMF.c1; c2 = GMF.c2;
    d0 = GMF.d0; d1 = GMF.d1;
    trans_u_yslf = GMF.trans_u_yslf;

    if (u <= trans_u_yslf){
        dg_yslf = -c1/pow(u,2) - 2*c2/pow(u,3);
    }
    else{
        dg_yslf = d1;
    }
    return dg_yslf;
}