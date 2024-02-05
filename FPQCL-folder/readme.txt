
General notes on the use of the Matlab code. Also red carefully the several comments presents in the code itself.
The units of all the constants and variables are specified in the code.
(Note: The code can be further optimised in terms of CPU time)

Written by Carlo Silvestri and Lorenzo Columbo. DECEMBER 2022

This code allows to integrate the Effective Semiconductor Maxwell Bloch Equations for the case of a mid-IR QCL as reported in the paper C. Silvestri et al., Opt. Express. 28, 23850 (2020) (ref. [1] in the following]).

1. Code structure

The complete code consists in three Matlab scripts named:
-'Main_TDTW.m'
-'TDTW_rungeFP.m'
-'DeviceData_OptExpr_FP.m'

The first one is the main script to be used to run the simulations.
The second one contains the algorithm to integrate the dynamical PDEs.
The third one contains the value of the main physical and geometrical parameters of the considered device
(NOTE: In particular the file 'DeviceData_OptExpr_FP.m' allows to reproduce the data in Fig. 2 of ref. [1])

2. Code execution

The (command window) instruction to run a single simulation is of the following type:

ResTD=Main_TDTW(pwd,600,0,300,'')

where in this case a bias current I=600mA and an integration interval between between 0 and 300ns have been chosen. 

In this example the simulation produces as output a structure 'ResTD' stored in a file named 'Res_600mAFrom0nsTo300ns.mat' and a file named 'Status_600mA_300ns.mat_600mAAt300ns'. In the latter the values of field, polarization and carriers variables as a functions of z at the end of the simulation are stored (in this case this corresponds to t=300ns).
To run a new simulation that lasts for example 400ns starting from the final state of the previous simulation you have to write in the command window the following (command window) instruction:

ResTD=Main_TDTW(pwd,600,0,400,'Status_600mA_300ns.mat_300mAAt300ns.mat')

By specifying the 'status' file with the final state to start with as initial condition.

3. Code visualisation

The results of the simulations are visualised by a) loading the file with the stored structure 'ReTD' for example with the (command window) instruction:

load('Res_600mAFrom0nsTo300ns.mat')

and b) running the post processing script named 'post_proc_simple.m' with the (command window) instruction:

post_proc_simple


This procedure generates a figure with the following subplots:
Power vs time
Optical Spectrum
Power Spectrum
Carriers at order 0 vs time
Carriers grating vs time
Instantaneous frequency vs time
 
and a separated plot with power and instantaneous frequency superimposed.

Note that in 'post_proc_simple.m' you can select the temporal trace that you want to analyse by properly setting the values of the start time 'T0' and the end time 'TEnd' to be considered.

The Matlab igure obtained by post processing the results of a simulation run with the current version of the code and device parameters is named 'OptExpress_Fig2_Comb_600mA.fig' and shows the results partially reported in Fig. 2 of ref. [1].

