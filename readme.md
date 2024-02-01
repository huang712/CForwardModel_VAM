# CForwardModel_VAM
This is the GNSS-R delay-Doppler map (DDM) forward model developed by Feixiong Huang from Purdue University. (contact email: fxhunag55@foxmail.com)

Reading a gridded wind field, a CYGNSS Level 1 product file and CYGNSS receiver antenna patterns, the forward model can produce a simulated DDM and a Jacobian matrix. For more information, please refer to the paper:

*Huang, Feixiong, et al. "A Forward Model for Data Assimilation of GNSS Ocean Reflectometry Delay-Doppler Maps." IEEE Transactions on Geoscience and Remote Sensing (2020). https://doi.org/10.1109/TGRS.2020.3002801*

Example demonstrations and a user guide document are provided on the Code Ocean:

*Feixiong Huang, Andrew O’Brien, Nereida Rodriguez-Alvarez, James Garrison (2020) GNSS-R DDM Forward Model and Variational Analysis Method for Data Assimilation [Source Code]. https://doi.org/10.24433/CO.5369859.v2*

The code is under the terms of GNU General Public License, Version 3. Please give approriate reference when you use the code for any publications.

## source code explaination: 

main.c : main functions of the program for data processing 

cygnss.h : define the structure of CYGNSS Level 1 data (observed DDM & parameters used by forward model)

cygnss.c : include the function readL1data() to read CYGNSS Level 1 data

forwardmodel.h : define all input/output structures of the forward model

initialization.c : initialization of input/output structures of forward model, including functions to read CYGNSS L1 structure, wind field data and antenna file.

forwardmodel.c : include the main function of forward model 

gnssr.h : define all global values, structures and functions for use of forward model (some from E2ES)

wind.c surface.c specular.c math.c grid.c geom.c ddm.c cood.c antenna.c debug.c GMF.c: functions defined by gnssr.h
