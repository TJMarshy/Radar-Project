# Radar-Project
Project for ananlysing simulation and experimnetal data from FMCW radar for range detection.


Current Version 1.1

Features:

Takes Data and Creates Basic Range Angle Map
Create Basic Geometric Simulation and get correct range information
correct angle info too, but need to reverse angle axes for simulated data, not sure why yet

Now Has All targets at correct positions

Now 4 ray reflections off floor

Used TE Fresnels equations to calculate power loss from each reflection

Consider different paths seperatley now taking into account their number of reflections/ distances travelled

Added interpolation, but currently limited to linear so may use interp2d to get other interpolation types

Still Need:

Noise?
RCS? from numerical stuff by propagation group
Creat Profiles of objects in rooms for different data sets

Notice angles are wrong for further out targets due to lower precision, so gonna test with more antenna just for curiosity
