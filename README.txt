README

This document is a brief  user guide to run the RIConduit 5.4 version of the code (13 october 2024)


%%%%%%%%%%
%1. Files%
%%%%%%%%%%

viscosity.m: This file contains a function that calculate the viscosity of melt using chemical composition, temperature and disolved water content. It is based on the paper of Giordano et al., (2008)

density.m: This file contains a function that calculate the density of melt, using chemical composition, temperature, disolved water content and pressure. It is based on the code of Tom D Perling, University of Sheffield,
which is based on the paper of Bottinga and Weill (1970).

fvrel.m: This file contains a function that calculates the relative viscosity due to crystals. The user can chose between 5 models: Einstein and Roscoe, 1952; Costa et al., (2009); Vona et al., (2011); Costa 2009, but with parameters of Vona et al., (2011)
and effective medium with 2 populations of crystals (phenocrysts and microlites)

radrub4: This file constains a function that calculates overpressure and dike width using the equations of Rubin (1995)

critara3: This file contains a function that iterates to find the conduit radius and solution that satisfies the criteria of Aravena et al. (2017) in order to find a estable conduit

histograph: This file contaiins a function that generate a frequency graph with the range of solutions found by the program.

RIconduitef5_4c...: These files contain the main part of the code and solve the differential equations of the conduit model to find Pressure, exsolved gas content, velocity of melt and gas, crystal content and bubble number density.
These files try to finf the effusive solution. The differences are: RIconduitef5_4c: Solves for the mainconduit5 file. In this case, only a single run is performed with the input parameters from the input script. RIconduit5_4crange
find solutions to a range of initial input parameters, specified in the calcrange file. RIconduitef5_4csmrange: It is the same than RIconduitef5_4crange, but with a different shooting method that applies when bubble coalescence is 
relevant (e.g. villarrica2015 case)

RIconduitex5_4c...: Same than the RIconduitef5_4c files, but here the program tries to find the explosive solution.

mainconduit5: Is the main file that run the code. The user should run this file to get a solution. In this file the program uses a conduit radius and overpressure specified by the user.

mainconduit5b: Similar to mainconduit5, but here, the program calculate the conduit dimensions and overpressure at the inlet of the conduit.

mainconduit5figp: Similar to mainconduit5b, but with different graphs. This program was used to generate te figures listed in the paper

calcrange: Similar to mainconduit5b, but here the program tries to find solutions to a range of initial conditions for crystal content, temperature and water content.



%%%%%%%%%%%%%%%%%
%2. Input script%
%%%%%%%%%%%%%%%%%

Input scripts are named following the volcano name and year of the eruption (e.g. calbuco2015d.m). Here, the user should specify the initial parameters. The main inputs are: Chemical composition of the melt (glass) in order to 
calculate viscosity and density; magma temperature, initial crystal content, depth of conduit inlet and water content. Others parameters not tested with respect to their sensitivity yet include characteristic time of crystallization,
aspect ratio of crystals, model for relative viscosity, relative size of bubbles, etc. This file should be in the same folder as the rest of files.


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. How to use the program%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

1. Fill the input script with the adequate parameters and save it with an specific name.

2. Open one of the following files: mainconduit5, mainconduit5b, mainconduit5figp or calcrange

3. Replace the name_parameters value with the input script name (line 5)

4. Run the code

Notes:

In the mainconduit5... codes, the output is a graph with the solution. The value of each single point is found in the output matrices soluef and soluex, where:

soluef(:,1)=magma pressures
soluef(:,2)=exsolved gas volume fraction
soluef(:,3)=bubble number density
soluef(:,4)=crystal content
soluef(:,5)=melt velocity
soluef(:,6)=exsolved gas velocity
zsolef=vertical coordinate points

The eruption rates for the explosive and effusive solution (m3/s) are saved in the Qex and Qef variables

If the user use the calcrange file, only the eruption rates (m3/s) are saved as outputs in the matrices resultsex and resultsef (m3/s).

In order to use the histograph file, the user should save the resultsex and resultsef matrices. For example, using the command: "save('calbuco2015test01.mat','resultsex')" and then put in line 1: A=load('calbuco2015test01.mat');
and in line 3: "B=A.resultsex" or "B=A.resultsef"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Written by Angelo Castruccio
Universidad de Chile
Santiago, Chile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






