# CForwardModel_VAM
This is the GNSS-R delay-Doppler map (DDM) forward model developed by Feixiong Huang from Purdue University. (contact email: huang712@purdue.edu)

Reading a gridded wind field and CYGNSS Level 1 data, the forward model can produce a simulated DDM and a Jacobian matrix.

## source code explaination: 

main.c : main functions of the program for data processing 

cygnss.h : define the structure of CYGNSS Level 1 data (observed DDM & parameters used by forward model)

cygnss.c : include the function readL1data() to read CYGNSS Level 1 data

forwardmodel.h : define all input/output structures of the forward model

initialization.c : initialization of input/output structures of forward model, including functions to read CYGNSSL1 structure, HWRF data and antenna file.

forwardmodel.c : include main function of forward model with the help of gnssr.h

gnssr.h : define all global values, structures and functions for use of forward model (from E2ES)

wind.c surface.c specular.c math.c grid.c geom.c ddm.c cood.c antenna.c debug.c : functions defined by gnssr.h
