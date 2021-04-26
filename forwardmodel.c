//---------------------------------------------------------------------------
//
// Main function of the forward model
// Created by Feixiong Huang on 10/21/17
//
//***************************************************************************

#include "forwardmodel.h"
#include "gnssr.h"

void forwardModel(struct metadata meta, struct powerParm pp, struct inputWindField iwf,
                  struct Geometry geom, struct DDMfm *ddm_fm, struct Jacobian *jacob, struct option opt)
{
    // opt.JacobOnOff = 0, compute only DDM;
    // opt.JacobOnOff = 1, compute DDM + Jacobian

    printf("Running forward model...\n");
    surface_initialize(meta);
    ddm_initialize(meta);
    antenna_initialize(pp);

    // Load wind field into structure wf
    windField wf;
    wind_initialize(&wf,meta, geom, iwf);
    wind_interpolate(&wf,geom,iwf,meta.grid_resolution_m);

    // Initialize geometry & specular point
    geometryData orbitGeometry;
    geom_initialize(&orbitGeometry, geom);
    //printfGeometry(orbitGeometry);

    int wfNum = 0;
    int geomInd = 0;  // just process one geometry
    int windModelType = 0;  // 0 for isotropic, 1 for anisotropic

    surface_calcGeomOverSurface(geom_getOrbitData(&orbitGeometry,geomInd),0,pp);
    surface_loadSurfWindfield(&wf,wfNum);
    surface_calcSigma0OnSurface(windModelType);

    // Calculate the DDM forward model
    ddm_binSurface();
    surface_composeTotalScatPowrOnSurface(1);  // 1 for no speckle and no mask; speckle not working now
    //surface_effArea();

    ddm_mapSurfaceToDDM();
    ddm_convolveFFT(2);
    ddm_mag();  //
    if (opt.thermalNoiseOnOff == 1) {ddm_addGaussianNoise();}  // add thermal noise

    ddm_save(meta,ddm_fm,1);  // resample DDM to 17x11 and save to structure ddm_fm
    if (opt.JacobOnOff == 1) {ddm_Hmatrix(meta, iwf, jacob);}  // compute and save to structure jacob

    // free memory
    free(wf.data);
    free(antenna.data);
    free(orbitGeometry.g);
    surface_cleanup(); // free surface.data
    ddm_cleanup(); // free DDM, H, ...

//    surface_saveWindToFile();
//    surface_saveDopplerToFile();
//    surface_saveDelayToFile();
}

