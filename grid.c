//
// Created by Feixiong Huang on 11/13/17.
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
//  This portion of the E2ES is responsible for creating a gridded surface
//  over the Earth centered at the specular point.  UNDER CONSTRUCTION
//
//****************************************************************************/


#include "gnssr.h"


typedef enum { FLAT_EARTH, SPHERICAL_EARTH, WGS84_EARTH } gridType;

void grid_analyzeSurface( geometryData *gd ){
    // for each geometry, determine how the surface should be
    // gridded for the desired DDM size


}

void grid_determineSurfaceSize( orbitGeometryStruct *geometry, double maxDelay, double minDoppler_hz, double maxDoppler_Hz){

}

void grid_construct(orbitGeometryStruct *geometry){

    double grid_pos_specular[3], grid_pos_ecef[3], grid_pos_llh[3], normal[3], area, distance;

    double *sx_pos = geometry->sx_pos;

    int Ni   = surface.numGridPtsX;
    int Nj   = surface.numGridPtsY;
    gridType type = surface.surfaceCurvatureType;
    double R = vector_norm(geometry->sx_pos);
    double Area_dS = pow(surface.resolution_m,2);


    for (int i = 0; i < Ni; i++) {
        for (int j = 0; j < Nj; j++) {

            // get the coordinates of the surface grid point
            switch (type) {

                case WGS84_EARTH:
                    // spare no expense
                    grid_getGridPt_sphericalEarth(i, j, sx_pos, grid_pos_specular, R);
                    matrixVectorMult3x3(geometry->SPEC_TO_ECEF_FRAME, grid_pos_specular, grid_pos_ecef );
                    wgsxyz2lla( grid_pos_ecef, grid_pos_llh );
                    grid_pos_llh[2] = 0;
                    wgslla2xyz(grid_pos_llh, grid_pos_ecef);
                    matrixVectorMult3x3(geometry->ECEF_TO_SPEC_FRAME, grid_pos_ecef, grid_pos_specular );

                    // approx normal
                    vector_unit(grid_pos_specular, normal);

                    // approx.  Account for change in surface area due to Earth curvature
                    area = Area_dS * (1 / normal[2]);
                    break;

                case SPHERICAL_EARTH:
                    grid_getGridPt_sphericalEarth(i, j, sx_pos, grid_pos_specular, R);
                    matrixVectorMult3x3(geometry->SPEC_TO_ECEF_FRAME, grid_pos_specular, grid_pos_ecef );
                    wgsxyz2lla( grid_pos_ecef, grid_pos_llh );
                    vector_unit(grid_pos_specular, normal);

                    // Account for change in surface area due to Earth curvature
                    area = Area_dS * (1 / normal[2]);
                    break;

                case FLAT_EARTH:
                    grid_getGridPt_flatEarth(i, j, sx_pos, grid_pos_specular);
                    matrixVectorMult3x3(geometry->SPEC_TO_ECEF_FRAME, grid_pos_specular, grid_pos_ecef );
                    wgsxyz2lla( grid_pos_ecef, grid_pos_llh );

                    // z-directed normal
                    normal[0] = 0;
                    normal[1] = 0;
                    normal[2] = 1;

                    area = Area_dS;
                    break;
            }


            int idx = SURFINDEX(i, j);
            surface.data[idx].i              = i;
            surface.data[idx].j              = j;
            surface.data[idx].pos_spec[0]    = grid_pos_specular[0];
            surface.data[idx].pos_spec[1]    = grid_pos_specular[1];
            surface.data[idx].pos_spec[2]    = grid_pos_specular[2];
            surface.data[idx].pos_ecef[0]    = grid_pos_ecef[0];
            surface.data[idx].pos_ecef[1]    = grid_pos_ecef[1];
            surface.data[idx].pos_ecef[2]    = grid_pos_ecef[2];
            surface.data[idx].pos_llh[0]     = grid_pos_llh[0];
            surface.data[idx].pos_llh[1]     = grid_pos_llh[1];
            surface.data[idx].pos_llh[2]     = grid_pos_llh[2];
            surface.data[idx].normal_spec[0] = normal[0];
            surface.data[idx].normal_spec[1] = normal[1];
            surface.data[idx].normal_spec[2] = normal[2];
            surface.data[idx].area_m2        = area;
            surface.data[idx].distance_m     = distance;
        }
    }
}


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