//---------------------------------------------------------------------------
//
//  This portion of the E2ES is responsible for creating a gridded surface
//  over the Earth centered at the specular point.  UNDER CONSTRUCTION
//
//****************************************************************************/

#include "gnssr.h"

typedef enum { FLAT_EARTH, SPHERICAL_EARTH, WGS84_EARTH } gridType;

/****************************************************************************/
//

void grid_getGridPt_sphericalEarth(int i, int j, double *sx_pos, double *grid_pos, double magSx ){
    // curved Earth using Scott's rectangular approach
    // given an i,j index and the specular position sx_pos, determine the
    // location of the surface grid position in the specular frame
    double T_II2x_A[9],T_II2x_B[9];
    double temp_vec[9];

    // magSx =  meters per radian (i.e. radius)  = earth radius at specular point
    double theta_rad = ((surface.width_m/magSx)/2.0) - j*(surface.resolution_m/magSx);   //(X/2-j)*dtheta
    double phi_rad   = ((surface.height_m/magSx)/2.0) - i*(surface.resolution_m/magSx);

    // find the coordinates of the scattering point, using rotation matrices
    // transformation with respect X axis. width angle
    T_II2x_B[0] = 1;                // A(1,1)
    T_II2x_B[1] = 0;                // A(1,2)
    T_II2x_B[2] = 0;                // A(1,3)
    T_II2x_B[3] = 0;                // A(2,1)
    T_II2x_B[4] = cos(theta_rad);   // A(2,2)
    T_II2x_B[5] = sin(theta_rad);  // A(2,3)
    T_II2x_B[6] = 0;                // A(3,1)
    T_II2x_B[7] = -1*sin(theta_rad);   // A(3,2)
    T_II2x_B[8] = cos(theta_rad);   // A(3,3)

    // transformation with respect Y axis. length angle
    T_II2x_A[0] = cos(phi_rad);  // A(1,1)
    T_II2x_A[1] = 0;                // A(1,2)
    T_II2x_A[2] = -1*sin(phi_rad);     // A(1,3)
    T_II2x_A[3] = 0;                // A(2,1)
    T_II2x_A[4] = 1;                // A(2,2)
    T_II2x_A[5] = 0;                // A(2,3)
    T_II2x_A[6] = sin(phi_rad);     // A(3,1)
    T_II2x_A[7] = 0;                // A(3,2)
    T_II2x_A[8] = cos(phi_rad);     // A(3,3)

    // coordinates of the scattering point (grid element), PUT = Sx*T_II2x_A*T_II2x_B;
    matrix_multiply(3,3,3,T_II2x_A,T_II2x_B,temp_vec);  //temp_vec=II2x_A * TT2x_B
    matrix_multiply(1,3,3,sx_pos,temp_vec,grid_pos);  //grid_pos=PUT  PUT: 1x3
}

void grid_getGridPt_flatEarth(int i, int j, double *sx_pos, double *grid_pos){
    grid_pos[0] = sx_pos[0] -(surface.width_m/2)  + j * surface.resolution_m;
    grid_pos[1] = sx_pos[1] -(surface.height_m/2) + i * surface.resolution_m;
    grid_pos[2] = sx_pos[2];
}