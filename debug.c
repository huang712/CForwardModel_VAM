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
