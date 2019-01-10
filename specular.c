//
// Created by Feixiong Huang on 11/9/17.
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
//  This file is dedicated to solving for the specular point location given
//  the Rx and Tx positions in ECEF coordinates.  This file has a single point
//  of entry (solveSpecularPt).  Note that it accurately accounts for an
//  WGS-84 ellipsoidal Earth.
//
/****************************************************************************/

#include "gnssr.h"

// prototypes
int solveSpecularPt_sphericalEarth(double rx_pos_ecef[3], double tx_pos_ecef[3],
                                   double sx_pos_ecef_unit[3], double r);
double minPositiveRoot( double polynomCoefs[5] );
double get_earthRadius_m( double x_ecef[3] );
int determineLOSSatelliteVis( double rx_pos_ecef[3], double tx_pos_ecef[3] );


/****************************************************************************/
//  Solve for Specular Point position and velocity in ECEF on WGS-84 Earth.

void solveSpecularPt(double rx_pos_ecef[3], double tx_pos_ecef[3], double sx_pos_ecef[3],
                     double rx_vel_ecef[3], double tx_vel_ecef[3], double sx_vel_ecef[3],
                     double *sx_valid, double r_init_m )
{
    int result = solveSpecularPtPosition(rx_pos_ecef, tx_pos_ecef, sx_pos_ecef,
                                         r_init_m, 100);

    if( result <= 0 ){
        // specular point should exist, but solver failed
        // or specular point does not exist
        *sx_valid = 0;
        sx_pos_ecef[0] = 0;
        sx_pos_ecef[1] = 0;
        sx_pos_ecef[2] = 0;
        sx_vel_ecef[0] = 0;
        sx_vel_ecef[1] = 0;
        sx_vel_ecef[2] = 0;
        return;
    }

    //-------------------------------------------------------------------------
    // solve for specular point velocity by propagating orbits along velocity
    // vectors for a small amount of time and solving for new specular point

    double rx_pos2_ecef[3],tx_pos2_ecef[3], sx_pos2_ecef[3],  dtime;
    dtime = 1;
    for(int i=0; i < 3; i++){
        rx_pos2_ecef[i] = rx_pos_ecef[i] + dtime * rx_vel_ecef[i];
        tx_pos2_ecef[i] = tx_pos_ecef[i] + dtime * tx_vel_ecef[i];
    }

    result = solveSpecularPtPosition(rx_pos2_ecef, tx_pos2_ecef, sx_pos2_ecef,
                                     vector_norm(sx_pos_ecef), 100);

    if( result <= 0 ){
        // specular point should exist, but solver failed
        // or specular point does not exist
        *sx_valid = -1;
        sx_pos_ecef[0] = 0;
        sx_pos_ecef[1] = 0;
        sx_pos_ecef[2] = 0;
        sx_vel_ecef[0] = 0;
        sx_vel_ecef[1] = 0;
        sx_vel_ecef[2] = 0;
        return;
    }

    for(int i=0; i < 3; i++)
        sx_vel_ecef[i] = ( sx_pos2_ecef[i] - sx_pos_ecef[i] ) / dtime;

    *sx_valid = 1;
}

/****************************************************************************/
//  Solve for Specular Point on WGS-84 Ellipsoid
//
// Solve for specular point on WGS84 elllipsoidal Earth by making iterative
// calls to Joel's spherical Earth spcular point solver.
// r_init_m is the initial guess at the Earth radius at the specular point.
// If it is <= 0, then this funciton uses the Earth radius at the Rx as the
// initial guess.  The better this initial guess, the fewer iterations will be
// needed.  Set maxIterations = 1 if you are really confident in r_init_m.
// If there is no specular point, then it returns a zero vector, so be sure
// to check for it.
/****************************************************************************/

int solveSpecularPtPosition(double rx_pos_ecef[3], double tx_pos_ecef[3],
                            double sx_pos_ecef[3], double r_init_m, int maxIterations){

    // first, check LOS visibility between Rx and Tx
    if( ! determineLOSSatelliteVis( rx_pos_ecef, tx_pos_ecef ) ){
        // no specular point for this geometry, so return a zero vector
        sx_pos_ecef[0] = 0;
        sx_pos_ecef[1] = 0;
        sx_pos_ecef[2] = 0;
        return(-1);
    }

    double r;
    double S[3], Stemp[3], temp[3];
    double correction;       // difference between iterations, meters
    double tol        = 1;   // convergence tolerance, meters
    double iterations;

    if( r_init_m <= 0 ){
        // first iteration: use Earth radius at receiver
        r = get_earthRadius_m(rx_pos_ecef);
    }else{
        r = r_init_m;
    }

    // do first iteration
    if( ! solveSpecularPt_sphericalEarth( rx_pos_ecef, tx_pos_ecef, S, r) )
        return(0);

    iterations = 1;

    correction = 100;
    while( (correction > tol) && (iterations < maxIterations)){

        // sucessive iterations: use Earth radius at estimated specular point
        r = get_earthRadius_m(S);
        if( ! solveSpecularPt_sphericalEarth(rx_pos_ecef, tx_pos_ecef, Stemp, r) )
            return(0);

        vector_subtract(S, Stemp, temp);
        correction = r * abs(vector_norm(temp));

        S[0] = Stemp[0];
        S[1] = Stemp[1];
        S[2] = Stemp[2];

        if(iterations++ > 100)
        fatalError("Error: In solveSpecularPt, couldnt converge.");
    }

    // constrain to Earth Surface
    r = get_earthRadius_m(S);
    vector_scale(S, sx_pos_ecef, r);
    return(1);
}

double get_earthRadius_m( double x_ecef[3] ){
    double WGS84_a = 6378137;
    double WGS84_f = 1/298.257223563;
    //double WGS84_b = WGS84_a*(1-WGS84_f);
    double WGS84_e = sqrt(2*WGS84_f - pow(WGS84_f,2));
    double theta   = asin(x_ecef[2]/vector_norm(x_ecef));  // latitude angle

    // ellipse radius at that latitude
    return WGS84_a*sqrt( (1 - pow(WGS84_e,2)) / (1 - (pow(WGS84_e,2)*pow(cos(theta),2))) );
}

/****************************************************************************/
//  Solve for Specular Point on Spherical Earth
/****************************************************************************/

int solveSpecularPt_sphericalEarth(double rx_pos_ecef[3], double tx_pos_ecef[3],
                                   double sx_pos_ecef_unit[3], double r){
    // Joel's analytic solution of specular point assuming spherical Earth
    // (some approximation due to truncation of the polynomial, on the 100m scale)

    double Rmag, Tmag, Rhat[3], That[3], gam, r1, r2, A, B, C, D, E, alp;

    Rmag = vector_norm(rx_pos_ecef);
    Tmag = vector_norm(tx_pos_ecef);
    vector_unit(rx_pos_ecef, Rhat);
    vector_unit(tx_pos_ecef, That);

    gam = acos(vector_dot_product(Rhat, That));

    r1  = r/Rmag;
    r2  = r/Tmag;

    A = (16-r1)/24*sin(gam);
    B = ((8-r1)*cos(gam)-r2)/6;
    C = (r1-4)*sin(gam)/2;
    D = cos(gam)*(r1-2)+r2;
    E = (1-r1)*sin(gam);

    double polynomCoefs[5] = { A, B, C, D, E };
    alp = minPositiveRoot(polynomCoefs);

    double Stemp[3];
    Stemp[0] = (Rhat[0] * sin(gam-alp) + That[0] * sin(alp)) / sin(gam);
    Stemp[1] = (Rhat[1] * sin(gam-alp) + That[1] * sin(alp)) / sin(gam);
    Stemp[2] = (Rhat[2] * sin(gam-alp) + That[2] * sin(alp)) / sin(gam);
    vector_unit(Stemp,sx_pos_ecef_unit);

    if(0){
        printf(" rx pos: %f %f %f\n", rx_pos_ecef[0], rx_pos_ecef[1], rx_pos_ecef[2] );
        printf(" tx pos: %f %f %f\n", tx_pos_ecef[0], tx_pos_ecef[1], tx_pos_ecef[2] );
        printf(" sx pos: %f %f %f\n", sx_pos_ecef_unit[0], sx_pos_ecef_unit[1], sx_pos_ecef_unit[2] );
    }

#define MIN_SPECULAR_ANGLE_TOL_DEG 0.01
#define MIN_SPECULAR_OFF_PLANE     0.001

    // check/constrain to Rx\Tx plane
    vector_constrainToPlane( rx_pos_ecef, tx_pos_ecef, sx_pos_ecef_unit );
    if( ( 1 - vector_norm(sx_pos_ecef_unit) ) > MIN_SPECULAR_OFF_PLANE ){
        printf("Error: Specular point failed in-plane check.\n");
        return(0);
    }
    vector_unit(sx_pos_ecef_unit, sx_pos_ecef_unit);

    // do verification on angles
    double rs_unit[3], ts_unit[3];
    vector_subtract(rx_pos_ecef, sx_pos_ecef_unit, rs_unit); vector_unit(rs_unit, rs_unit);
    vector_subtract(tx_pos_ecef, sx_pos_ecef_unit, ts_unit); vector_unit(ts_unit, ts_unit);
    double angleDiff_deg =  abs( acos(vector_dot_product(rs_unit,sx_pos_ecef_unit))
                                 - acos(vector_dot_product(ts_unit,sx_pos_ecef_unit)) ) * R2D;
    if( angleDiff_deg > MIN_SPECULAR_ANGLE_TOL_DEG ){
        printf("Error: Specular point failed angle check. (difference = %f deg)\n",angleDiff_deg);
        return(0);
    }

    return(1);
}


/****************************************************************************/
// find the minimum positive root of a 4th order polynomial

typedef struct FCOMPLEX {float r,i;} fcomplex;
void zroots(fcomplex a[], int m, fcomplex roots[5], int polish);
void laguer(fcomplex a[], int m, fcomplex *x, int *its);

double minPositiveRoot( double polynomCoefs[5] ){
    fcomplex roots[5], a[5];

    // Matlab uses descending coefs, but zroots expects ascending
    for(int i = 0; i< 5; i++){
        a[i].r = polynomCoefs[4-i];
        a[i].i = 0.0;
        roots[i].r = 0;
        roots[i].i = 0;
    }

    zroots(a, 4, roots, 1);

    double r, minPosRoot = DBL_MAX;
    for(int i = 0; i<=4; i++ ){
        r = roots[i].r;
        if( (r > 0) && (r < minPosRoot) )
            minPosRoot = r;
    }

    return minPosRoot;
}

/****************************************************************************/
//  Check LOS visibility between Rx and Tx

int determineLOSSatelliteVis( double rx_pos_ecef[3], double tx_pos_ecef[3] ){

    // define WGS84 Ellipsoid ---------------
    double a = 6378137.0;      // semi-major axis, meters  (x-axis)
    double b = a;
    double c = 6356752.314245;  // semi-minor axis, meters (z-axis)

    // define ray between Rx and Tx
    double v[3],rx[3],tx[3];

    rx[0] = (1/a) * rx_pos_ecef[0];
    rx[1] = (1/b) * rx_pos_ecef[1];
    rx[2] = (1/c) * rx_pos_ecef[2];
    tx[0] = (1/a) * tx_pos_ecef[0];
    tx[1] = (1/b) * tx_pos_ecef[1];
    tx[2] = (1/c) * tx_pos_ecef[2];

    vector_subtract(tx,rx,v); // from rx to tx
    vector_unit(v, v);

    // coefs of quadratic formula
    double A = pow(vector_norm(v),2);
    double B = 2*vector_dot_product(rx, v);
    double C = pow(vector_norm(rx),2) - 1;

    // if all roots are imaginary, then no intersection between ray and
    // ellipsoid, so LOS exists
    if( (pow(B,2) - 4*A*C) < 0 )
        return 1;
    else{
        // if both solutions are negative, then nothing between rx and tx
        double t1 = (-B + sqrt(pow(B,2) - 4*A*C)) / (2*A);
        double t2 = (-B - sqrt(pow(B,2) - 4*A*C)) / (2*A);
        if( (t1 < 0) && (t2 < 0 ) )
            return 1;
        else
            return 0;
    }
}

/****************************************************************************/
//  Roots of Polynomials, from Numerical Recipes C9-5
//  TODO: Check licensing issues
/****************************************************************************/

#define EPSS 1.0e-5 // was e-7, but I changed it to prevent "TOO MANY ITERATIONS" error
#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define EPS 2.0e-6
#define MAXM 100

fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(float re, float im);
fcomplex Conjg(fcomplex z);
fcomplex Cdiv(fcomplex a, fcomplex b);
float Cabs(fcomplex z);
fcomplex Csqrt(fcomplex z);
fcomplex RCmul(float x, fcomplex a);

void zroots(fcomplex a[], int m, fcomplex roots[5], int polish){
    int i,its,j,jj;
    fcomplex x,b,c,ad[MAXM];
    for (j=0;j<=m;j++) ad[j]=a[j];
    for (j=m;j>=1;j--) {
        x = Complex(0.0,0.0);
        laguer(ad,j,&x,&its);
        if (fabs(x.i) <= 2.0*EPS*fabs(x.r)) x.i = 0.0;
        roots[j]=x;
        b=ad[j];
        for (jj=j-1;jj>=0;jj--) {
            c=ad[jj];
            ad[jj]=b; b = Cadd(Cmul(x,b),c);
        }
    }
    if (polish)
        for (j=1;j<=m;j++)
            laguer(a,m,&roots[j],&its);
    for (j=2;j<=m;j++) {
        x=roots[j];
        for (i=j-1;i>=1;i--) {
            if (roots[i].r <= x.r) break;
            roots[i+1]=roots[i];
        }
        roots[i+1]=x;
    }
}

void laguer(fcomplex a[], int m, fcomplex *x, int *its){
    int iter,j;
    float abx,abp,abm,err;
    fcomplex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
    static float frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

    for(iter = 1; iter <= MAXIT; iter++){
        *its=iter;
        b=a[m];
        err=Cabs(b); d=f=Complex(0.0,0.0); abx=Cabs(*x);
        for (j=m-1;j>=0;j--) {
            f=Cadd(Cmul(*x,f),d);
            d=Cadd(Cmul(*x,d),b);
            b=Cadd(Cmul(*x,b),a[j]);
            err=Cabs(b)+abx*err;
        }
        err *= EPSS;
        if (Cabs(b) <= err) return;
        g=Cdiv(d,b); //The generic case: use Laguerre’s formula.
        g2=Cmul(g,g);
        h=Csub(g2,RCmul(2.0,Cdiv(f,b)));
        sq=Csqrt(RCmul((float) (m-1),Csub(RCmul((float) m,h),g2)));
        gp=Cadd(g,sq);
        gm=Csub(g,sq);
        abp=Cabs(gp);
        abm=Cabs(gm);
        if (abp < abm) gp=gm;
        dx=((fmax(abp,abm) > 0.0 ? Cdiv(Complex((float) m,0.0),gp)
                                 : RCmul(1+abx,Complex(cos((float)iter),sin((float)iter)))));
        x1=Csub(*x,dx);
        if (x->r == x1.r && x->i == x1.i) return; //Converged.
        if (iter % MT) *x=x1;
        else *x=Csub(*x,RCmul(frac[iter/MT],dx));
    }
    printf("roots: too many iterations in laguer\n");
    return;
}

/****************************************************************************/
// From Numerical Recipes Complex.c.  These are just used in the root finder
// above.  No where else in the simulator.  Probably best to replace with
// standard complex doubles eventually
/****************************************************************************/

fcomplex Cadd(fcomplex a, fcomplex b)
{
    fcomplex c;
    c.r=a.r+b.r;
    c.i=a.i+b.i;
    return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
{
    fcomplex c;
    c.r=a.r-b.r;
    c.i=a.i-b.i;
    return c;
}


fcomplex Cmul(fcomplex a, fcomplex b)
{
    fcomplex c;
    c.r=a.r*b.r-a.i*b.i;
    c.i=a.i*b.r+a.r*b.i;
    return c;
}

fcomplex Complex(float re, float im)
{
    fcomplex c;
    c.r=re;
    c.i=im;
    return c;
}

fcomplex Conjg(fcomplex z)
{
    fcomplex c;
    c.r=z.r;
    c.i = -z.i;
    return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b)
{
    fcomplex c;
    float r,den;
    if (fabs(b.r) >= fabs(b.i)) {
        r=b.i/b.r;
        den=b.r+r*b.i;
        c.r=(a.r+r*a.i)/den;
        c.i=(a.i-r*a.r)/den;
    } else {
        r=b.r/b.i;
        den=b.i+r*b.r;
        c.r=(a.r*r+a.i)/den;
        c.i=(a.i*r-a.r)/den;
    }
    return c;
}

float Cabs(fcomplex z)
{
    float x,y,ans,temp;
    x=fabs(z.r);
    y=fabs(z.i);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
    } else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
    }
    return ans;
}

fcomplex Csqrt(fcomplex z)
{
    fcomplex c;
    float x,y,w,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
        c.r=0.0;
        c.i=0.0;
        return c;
    } else {
        x=fabs(z.r);
        y=fabs(z.i);
        if (x >= y) {
            r=y/x;
            w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        } else {
            r=x/y;
            w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.r >= 0.0) {
            c.r=w;
            c.i=z.i/(2.0*w);
        } else {
            c.i=(z.i >= 0) ? w : -w;
            c.r=z.i/(2.0*c.i);
        }
        return c;
    }
}

fcomplex RCmul(float x, fcomplex a)
{
    fcomplex c;
    c.r=x*a.r;
    c.i=x*a.i;
    return c;
}
