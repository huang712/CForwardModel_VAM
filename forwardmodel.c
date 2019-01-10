//
// Created by Feixiong Huang on 10/21/17.
//
#include "forwardmodel.h"
#include "gnssr.h"

void forwardModel(struct metadata meta, struct powerParm pp, struct inputWindField iwf,
                  struct Geometry geom, struct DDMfm *ddm_fm, struct Jacobian *jacob, int option)
{
    //option = 0, only compute DDM; option = 1, compute DDM+Jacobian
    printf("Running forward model...\n");
    surface_initialize(meta);
    ddm_initialize(meta);
    antenna_initialize(pp);
    //ddmaLUT_initialize();

    //load wind field into structure wf
    windField wf;
    wind_initialize(&wf,meta, geom, iwf);
    wind_interpolate(&wf,geom,iwf,meta.grid_resolution_m);

    //initialize geometry & specular point
    geometryData orbitGeometry;
    geom_initialize(&orbitGeometry, geom);

    //printfGeometry(orbitGeometry);

    int wfNum = 0;
    int geomInd = 0; //just process one geometry
    int windModelType = 0; //0 for isotropic, 1 for anisotropic

    surface_calcGeomOverSurface(geom_getOrbitData(&orbitGeometry,geomInd),0,pp);
    surface_loadSurfWindfield(&wf,wfNum);
    surface_calcSigma0OnSurface(windModelType);

    //calculate the DDM forward model
    ddm_binSurface();
    surface_composeTotalScatPowrOnSurface(1); //1 for no speckle.  no mask

    ddm_mapSurfaceToDDM();
    ddm_convolveFFT(2);
    ddm_mag(); //save to DDM[i];

    ddm_save(meta,ddm_fm,1);  //save to structure ddm_fm

    if (option == 1){
        ddm_Hmatrix(meta, iwf, jacob);  //compute and save to structure jacob
    }


    //surface_saveWindToFile();
    //surface_saveDopplerToFile();
    //surface_saveDelayToFile();
}

