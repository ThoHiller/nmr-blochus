## BLOCHUS

**BLOCH** **U**niversal **S**imulator

### About

**BLOCHUS** is a set of MATLAB<sup>TM</sup> tools, that allow some basic simulations of (S)NMR spin dynamics based on the Bloch equations. The Bloch equations are solved in the laboratory frame of reference with MATLABs built-in `ode45` solver. Because it was developed within the scope of a near surface SNMR project, its main features are the simulation of (1) pre-polarization switch-off ramps and (2) excitation pulses. The main front-end to the underlying simulation tools is a graphical user interface (GUI) that allows playing around with the different features and helps to understand the basic concepts of (S)NMR spin dynamics.

#### Basic features:
1. Choose between different protons (e.g. *Hydrogen*, *Helium*, *Flourine*, etc.)
2. Choose between different pre-polarization switch-off ramp shapes (e.g. *exponential*, *linear*, *half cosine*,  etc.) with arbitrary ramp time 
3. Choose between different excitation pulses (*pi/2*, *pi*, *adiabatic half passage*, etc.) with arbitrary off-resonant frequency, or in case of the adiabatic pulses, arbitrary frequency and current modulation.

- - -

### Requirements

In order to work properly you need to meet the following requirements:

1. The [Mathworks](https://www.mathworks.com) MATLAB<sup>TM</sup> software development environment (tested with R2017b and newer)
2. The GUI Layout Toolbox (get it from [FEX](https://de.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)) (<span style="color:red">required</span>)
3. `findjobj` (get it from [FEX](https://de.mathworks.com/matlabcentral/fileexchange/14317-findjobj-find-java-handles-of-matlab-graphic-objects)) (<span style="color:red">required</span>)

#### Operating System

I tested it successfully under Windows 7 (64bit) and 10 (64bit) with Matlab R2018a.

**NOTE:** So far I did not test anything on Linux or a Mac. If you get it to work on either of the two systems (which it basically should I guess) please let me know.

- - -

### Installation

1. It is recommended to install the GUI Layout Toolbox directly into MATLAB<sup>TM</sup> via the mltbx-file (but it should also work via the old-school way of adding the toolbox folders to the MATLAB<sup>TM</sup> path)
2. To use **BLOCHUS** you just need to place the `blochus` folder from  the git repository on your hard drive and use the start script `startBLOCHUS` (within this script all necessary **BLOCHUS** folders are added to the MATLAB<sup>TM</sup> path)

- - -

### Usage

1. By executing the start scripts (see above)
2. Simply type `BLOCHUS` on the MATLAB<sup>TM</sup> prompt (make sure the `blochus` folder is on the MATLAB<sup>TM</sup> path)
3. Check the example scripts for the usage of the core functions without the GUI (inside the `scripts` folder)

- - -

### Documentation

A basic documentation to **BLOCHUS** can be found in the `blochus\doc` folder. Just open the `index.html` in the web browser of your choice. The documentation was created with [m2html](https://www.artefact.tk/software/matlab/m2html/) by Guillaume Flandin.

- - -

### TODO

In no particular order and without guarantee that it will ever happen :-) :

1. Adapt the core functionality in a Python module (this is on top of my agenda)
2. A real Manual
3. ...

- - -

### References

1. Hiller, T., Dlugosch, R. and Müller-Petke, M., "Utilizing pre-polarization to enhance SNMR signals - effect of imperfect switch-off", Geophysical Journal International **0**(0), p.1-15, 2020, [DOI](https://doi.org/10.1093/gji/ggaa216)

- - -
<p style="text-align: center;"> MATLAB is a registered trademark of The Mathworks, Inc. </p>