# CO2BLOCK: A tool for estimates of storage resource capacity in saline aquifers

CO2BLOCK provides the automatic estimate of the storage capacity of an aquifer, based on the geomechanical limits to pressure build-up. 

Users only need to open and modify the script "CO2BLOCK.m". There is no need to open or modify the other scripts, which is instead not recommended.

Please refer to the user guide for more detailed instructions. 
If you want to learn more about the theory behind CO2BLOCK, please refer to the following papers:
- De Simone and Krevor (2021) A tool for first order estimates and optimisation of dynamic storage resource capacity in saline aquifers”.  International Journal of Greenhouse Gas Control, 106, 103258. 
- De Simone S., Jackson S. J. & Krevor S. (2019). The error in using superposition to estimate pressure during multi-site subsurface CO2 storage. Geophysical Research Letters 46: 6525-6533 

Should you have any question/comment/suggestion, please send an email to Silvia De Simone at silviadesi@gmail.com

CO2BLOCK is a free licence software.
If you use it for academic purposes, please cite the reference paper:

- De Simone and Krevor (2021) A tool for first order estimates and optimisation of dynamic storage resource capacity in saline aquifers”.  International Journal of Greenhouse Gas Control, 106, 103258. https://doi.org/10.1016/j.ijggc.2021.103258

###########

THIS VERSION ("CO2BLOCK_mod" release from 2026) CONTAINS SOME MODIFICATIONS WITH RESPECT TO THE ORIGINAL VERSION ("1st release 1.1" from 2021).

The UserGuide pdf does not contemplate these modifications, and the use of the original version is suggested in case of doubts.

Modifications:
1) the pressure build-up is calculated through the full Nordbotten solution (three-zone), not through the two-zone simplified formulation - differences are minor
2) there is the possibility of activating the Srivastava & Guzman correction, which improves accuracy close to the pressure front 
3) Nordbotten_solution now includes the case in which the CO2 plume extends outside of the pressurized region (R<csi)


