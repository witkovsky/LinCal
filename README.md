# LinCal
LinCal: Linear Calibration

MATLAB files and supplementary material to the working paper:

Palenčár J., Palenčár R., Chytil M., Witkovský V., Wimmer G. Jr., and Wimmer G.: On Linear Comparative Calibration and Uncertainty Analysis of the Measurement Result Obtained With the Calibrated Instrument. Working Paper 2022.

For current status of the MATLAB toolbox see the CharFunTool development available at

- https://github.com/witkovsky/LinCal

About
=====

The Linear Calibration Toolbox (LinCal) consists of a set of algorithms (LinearCalibration.m) and example data (Simulations.m) for linear comparative calibration as presented in the working paper Palenčár et al. 2022. Moreover, the PolyCal is an algorithm suggested for fitting a polynomial calibration and for computing the exact coverage intervals (if used with cf2DistGP algorithm of the package CharFunTool).

Installation and requirements
=============================

To install, you can either clone the directory with Git or download a .zip file. 

## Option 1: Download .zip file

Download a .zip of LinCal from

- https://github.com/witkovsky/LinCal/archive/master.zip

After unzipping, you will need to add LinCal to the MATLAB path. You can do this either (a) by typing
```
addpath(LinCalRoot), savepath
```
where `LinCalRoot` is the path to the unzipped directory, (b) by selecting the `LinCal` directory with the `pathtool` command, or (c) though the File > Set Path... dialog from the MATLAB menubar.

## Option 2: Clone with Git

To clone the LinCal repository, first navigate in a terminal to where you want the repository cloned, then type
```
git clone https://github.com/witkovsky/LinCal.git
```
To use LinCal in MATLAB, you will need to add the `LinCal` directory to the MATLAB path as above.

Getting started
===============

We recommend taking a look at the Simulations file and the detailed helps of the included functions. 

* To get a taste of what computing with LinCal is like, try simple example,
```
   x   = [1.2 1.9 2.9 4.0 4.7 5.9]';
   ux  = [0.2 0.2 0.2 0.2 0.2 0.2]';
   y   = [3.4 4.4 7.2 8.5 10.8 13.5]'; 
   uy  = [0.2 0.2 0.2 0.4 0.4 0.4]';
   y0  = [9 11];
   uy0 = [0.4 0.4];
   [ab,Uab,x0a,ux0a] = LinearCalibration(x,y,ux,uy,[],y0,uy0);
```

License
=======

See `LICENSE.txt` for CharFunTool licensing information.
