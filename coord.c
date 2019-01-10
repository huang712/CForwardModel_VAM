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
// Spatial & Temportal Coordinate System Transforms
//
//***************************************************************************


#include "gnssr.h"

/****************************************************************************/
//

void cart2sph( double x[3], double y[3] ){
    // Matlab / mathematicians definition cartesian to spherical coords
    y[0] = sqrt (pow(x[0],2) + pow(x[1],2) + pow(x[2],2)); //R
    y[1] = atan2(x[1], x[0]); // Theta  (azimuth)
    y[2] = atan2(x[2], sqrt(pow(x[0],2) + pow(x[1],2))); //Phi (elevation)
}

/****************************************************************************/
//

void wgslla2xyz(double x_llh[3], double x_ecef[3]){
    // LLH to ECEF

    const double f = 1/298.257223563;        //   WGS-84 Flattening.
    const double e = sqrt(f*(2 - f));        //   Eccentricity.
    const double R_0 = 6378137;              //   WGS-84 equatorial radius (m).

    //   Compute East-West Radius of curvature at current position
    const double R_E = R_0 / sqrt(1 - pow(e*sin(x_llh[0]),2));

    //   Compute ECEF coordinates
    x_ecef[0] = (R_E + x_llh[2])*cos(x_llh[0])*cos(x_llh[1]);
    x_ecef[1] = (R_E + x_llh[2])*cos(x_llh[0])*sin(x_llh[1]);
    x_ecef[2] = ((1 - pow(e,2))*R_E + x_llh[2])*sin(x_llh[0]);
}

void wgsxyz2enu(double p_ecef[3] ,double x_llh[3], double pe_enu[3]){
    // ECEF to ENU

    double delta_xyz[3],x_ecef[3];

    // Calculate the relative position vector
    wgslla2xyz(x_llh, x_ecef);
    vector_subtract(p_ecef, x_ecef, delta_xyz);

    // Calculate ENU coordinates
    pe_enu[2] = cos(x_llh[0])*cos(x_llh[1])*delta_xyz[0] +
                cos(x_llh[0])*sin(x_llh[1])*delta_xyz[1] +
                sin(x_llh[0])*delta_xyz[2];
    pe_enu[0] = -sin(x_llh[1])*delta_xyz[0] + cos(x_llh[1])*delta_xyz[1];
    pe_enu[1] = -sin(x_llh[0])*cos(x_llh[1])*delta_xyz[0] -
                sin(x_llh[0])*sin(x_llh[1])*delta_xyz[1] +
                cos(x_llh[0])*delta_xyz[2];
}

void wgsxyz2lla( double *x_ecef, double *x_lla ) {
    // ECEF to LLH

    double f,e,omega_ie,R_0,R_P,mu_E;
    double x,y,z,lon;
    double p,E,F,G,c,s,P,Q,k_1,k_2,k_3,k_4,r_0,k_5,U,V;
    double alt,z_0,e_p,lat;

    // paramters describing the WGS-84 reference ellipsoid
    f = 1/298.257223563;        //   WGS-84 Flattening.
    e = sqrt(f*(2 - f));        //   Eccentricity.
    omega_ie = 7.292115e5;      //   WGS-84 Earth rate (rad/s).
    R_0 = 6378137;              //   WGS-84 equatorial radius (m).
    R_P = R_0*(1 - f);          //   Polar radius (m).
    mu_E = 3.986004418e14;      //   WGS-84 Earth's gravitational

    //  Place holder for output and temporary variables
    x = x_ecef[0];
    y = x_ecef[1];
    z = x_ecef[2];

    //  Calculate longitude
    lon = atan2(y,x)*R2D;

    p = sqrt(pow(x,2) + pow(y,2));
    E = sqrt(pow(R_0,2) - pow(R_P,2));
    F = 54*pow((R_P*z),2);
    G = pow(p,2) + (1 - pow(e,2))*pow(z,2) - pow((e*E),2);
    c = pow(e,4)*F*pow(p,2)/pow(G,3);
    s = pow((1 + c + sqrt(pow(c,2) + 2*c)),(1/3));
    P = (F/(3*pow(G,2)))/(pow((s + (1/s) + 1),2));
    Q = sqrt(1 + 2*pow(e,4)*P);
    k_1 = -P*pow(e,2)*p/(1 + Q);
    k_2 = 0.5*pow(R_0,2)*(1 + 1/Q);
    k_3 = -P*(1 - pow(e,2))*pow(z,2)/(Q*(1 + Q));
    k_4 = -0.5*P*pow(p,2);
    r_0 = k_1 + sqrt(k_2 + k_3 + k_4);
    k_5 = (p - pow(e,2)*r_0);
    U = sqrt(pow(k_5,2) + pow(z,2));
    V = sqrt(pow(k_5,2) + (1 - pow(e,2))*pow(z,2));
    alt = U*(1 - (pow(R_P,2)/(R_0*V)));

    //   Compute additional values required for computing
    //   latitude
    z_0 = (pow(R_P,2)*z)/(R_0*V);
    e_p = (R_0/R_P)*e;
    lat = atan((z + z_0*pow(e_p,2))/p) * R2D;

    x_lla[0] = lat;
    x_lla[1] = lon;
    x_lla[2] = alt;
}
