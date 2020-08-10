//
// Created by Feixiong Huang on 11/14/17.
//

#include "gnssr.h"
void printGMF(){
    printf("a0 = %f  ",GMF.a0);
    printf("a1 = %f  ",GMF.a1);
    printf("a2 = %f  ",GMF.a2);
    printf("\n");
    printf("b0 = %f  ",GMF.b0);
    printf("b1 = %f  ",GMF.b1);
    printf("b2 = %f  ",GMF.b2);
    printf("\n");
    printf("c0 = %f  ",GMF.c0);
    printf("c1 = %f  ",GMF.c1);
    printf("c2 = %f  ",GMF.c2);
    printf("\n");
    printf("d0 = %f  ",GMF.d0);
    printf("d1 = %f  ",GMF.d1);
    printf("\n");
    printf("tras_u_fds = %f  ",GMF.trans_u_fds);
    printf("trans_u_yslf = %f  ",GMF.trans_u_yslf);
    printf("\n");
}

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
    printf("rx_vel = %f   ",orbitGeometry.g[0].rx_vel_ecef[0]);
    printf("%f   ",orbitGeometry.g[0].rx_vel_ecef[1]);
    printf("%f   ",orbitGeometry.g[0].rx_vel_ecef[2]);
    printf("\n");
    printf("tx_vel = %f   ",orbitGeometry.g[0].tx_vel_ecef[0]);
    printf("%f   ",orbitGeometry.g[0].tx_vel_ecef[1]);
    printf("%f   ",orbitGeometry.g[0].tx_vel_ecef[2]);
    printf("\n");
    printf("sx_vel = %f   ",orbitGeometry.g[0].sx_vel_ecef[0]);
    printf("%f   ",orbitGeometry.g[0].sx_vel_ecef[1]);
    printf("%f   ",orbitGeometry.g[0].sx_vel_ecef[2]);
    printf("\n");
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
