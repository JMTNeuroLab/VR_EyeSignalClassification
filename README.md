# HcTask_EyeSignalProcessing
Matlab scripts for processing eye signals from NHP collected when using dynamic stimuli which can generate smooth pursuits, such as curing virtual navigation. 

Can be used on data where blinks and non-physiological signals have been removed, and data is in degrees visual angle, collected at 500Hz, but might be possible to use at lower sampling rates with minor modifications.

This work is described in:
Benjamin W. Corrigan, Roberto A. Gulli, Guillaume Doucet, Julio C. Martinez-Trujillo; Characterizing eye movement behaviors and kinematics of non-human primates during virtual navigation tasks. Journal of Vision 2017;17(12):15. doi: https://doi.org/10.1167/17.12.15.

Includes meanangle.m from Jeff Dunne and circ_rtest.m from Philip Berens' Circular statistics toolbox for matlab: https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics licenses for which are included in the separate txt.

