# Radar-Project
Project for ananlysing simulation and experimnetal data from FMCW radar for range detection.


Current Version 1.2

Features:

Takes Data and Creates Basic Range Angle Map
Create Basic Geometric Simulation and get correct range information
correct angle info too, but need to reverse angle axes for simulated data, not sure why yet

Now Has All targets at correct positions

Now 4 ray reflections off floor

Used TE Fresnels equations to calculate power loss from each reflection

Consider different paths seperatley now taking into account their number of reflections/ distances travelled

Added interpolation, but currently limited to linear

Added Experiments file to load different experiment data

Added zeros to sim data and now it almost perfectly matches the sim data and the interpolation is working a dreams

Still Need:

Noise?
RCS? from numerical stuff by propagation group
Finish Profiles of objects in rooms for different data sets E.g. RCS for big reflector

Add emtpy array before and after real data to get the same sinc functions

