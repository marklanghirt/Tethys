function [vao] = tethysparams()
% [jsonpath]=tethysparams()
% Description: this routine provides a template/example for writing the 
%  parameter .json-file read in by the "Tethys" energy flux model wrapper
%
% Copyright: ©2025, Mark Adam Langhirt
% License: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT

arguments (Output, Repeating); vao; end
mpth = string(mfilename("fullpath")); 
jsnpth = split(mpth,filesep); jsnpth = jsnpth(1:end-1);
name = "wej1zslc_dbg"; jsnpth = join([jsnpth;name+".json"],filesep);
optcell = { "clear", 1, %(bool) clears intermediary variables during execution to save memory
            "debug", 0, %(bool) pauses executation at the end of the main computation
            "parType", "for", %(string) type of parallelization: "vec","for","background","process","thread","gpu"; see "chspool.m"
            "omitCore", 1, %(bool) omit one core when using thread/process parallelization
            "calcRays", 0, %(bool) simulate ray trajectories and include in output data
            "ntRays", 10, %(int) number of theta (elevation) angles for rays
            "npRays", 10, %(int) number of phi (azimuth) angles for rays
            "botLoss", 1, %(bool) calculate accumulated bottom loss along horizontal trajectories
            "ntSamp", 47, %(int) number of theta sample points for interpolations
            "nxSamp", 48, %(int) number of x-range sample points for interpolations
            "nySamp", 49, %(int) number of y-range sample points for interpolations
            "nIsoLim", 50, %(int) uniform max array length for iso-modenumber k_xy-profiles
            "dirKern", 1, %(bool) calculate Lloyd Mirror directivity kernels for src/rcv near srf/bot
            "maxLamb", 20, %(num) max wavelength distance from boundary for lloyd mirror directivity kernel
            "convType", 1, %(int) type of convergence factor calculation: (0)off, (1)periodicGaussian, (2)dirichletKernel/cosineSummation 
            "yConv", 1, %(bool) compute transverse convergence factor
            "ymmd", 16, %(num) transverse max modenumber difference for convergence factor
            "zConv", 0, %(bool) compute vertical convergence factor
            "zmmd", 16, %(num) vertical max modenumber difference for convergence factor
            "wkbLim", 1, %(bool) limit WKB modal amplitudes using rough approx. of Airy function peak
            "finFreq", 1, %(bool) limit |theta|_min according to finite freq. approx: (0)off, (1)vac-rgd, (2)median, (3)vac-vac
            "splitTh", 1, %(bool) split theta (elevation) angle domain in final integration
            "splitPh", 1 %(bool) split phi (azimuth) angle domain in final integration
          }.';
geocell = { "z", [nan,nan], %(num) list of sediment layer thicknesses, first (watercolumn) and last (basement) are always NaN
            "cp", [1500,1700], %(num) compressive sound speed of layers in [m/s]
            "ap", [0,0], %(num) compressive attenuation of layers in [dB/kmHz]
            "rho", [1,50000] %(num) medium density of layers in [g/cc]
          }.';
envcell = { "c0", 1500, %(num) watercolumn sound speed in [m/s]
            "y", [-4000, 0, 10000], %(num) transverse y-range nodes for bathymetry profile
            "H", [0, 199.999634236445,699.998719827557], %(num) transverse depth values for bathymetry profile
            "geo", struct(geocell{:})
          }.';
srccell = { "thEps", 0.000001, %(num) epsilon (small number) theta (elevation) buffer to avoid integrating over 0deg and 90deg
            "phEps", 0.000001, %(num) epsilon (small number) phi (azimuth) buffer to avoid integrating over 0deg and 90deg
            "monoAmp", 12.566370614359172, %(num) monopole amplitude, usually 4*pi to normalize pressure to 1Pa at 1m distance from source
            "f", 25, %(num) temporal frequency of acoustic source signal
            "x", 0, %(num) x-range position of source in [m], not used...
            "y", 0, %(num) y-range position of source in [m]
            "z", 100, %(num) z-depth position of source in [m]
            "pthN", 104, %(int) number of absolute theta (elevation) angles for angular integration domain
            "pthDegMin", 0, %(num) absolute theta (elevation) angle minimum for angular integration domain
            "pthDegMax", 89.9, %(num) absolute theta (elevation) angle maximum for angular integration domain
            "pphN", 105, %(int) number of absolute phi (azimuth) angles for angular integration domain
            "pphDegMin", 0, %(num) absolute phi (azimuth) angle minimum for angular integration domain
            "pphDegMax", 89.9 %(num) absolute phi (azimuth) angle maximum for angular integation domain
          }.';
rcvcell = { "xN", 101, %(int) number of longitudinal x-range positions for receiver grid
            "xMin", 100, %(num) minimum longitudinal x-range position in [m] for receiver grid
            "xMax", 25000, %(num) maximum longitudinal x-range position in [m] for receiver grid
            "yN", 102, %(int) number of transverse y-range positions for receiver grid
            "yMin", -3999.9, %(num) minimum transverse y-range position in [m] for receiver grid
            "yMax", 5999.9, %(num) maximum transverse y-range position in [m] for receiver grid
            "zN", 1, %(int) number of vertical z-depth positions for receiver grid
            "zMin", 0, %(num) minimum vertical z-depth position in [m] for receiver grid
            "zMax", 30 %(num) maximum vertical z-depth position in [m] for receiver grid
          }.';
s.opt = struct(optcell{:});
s.env = struct(envcell{:});
s.src = struct(srccell{:});
s.rcv = struct(rcvcell{:});
writestruct(s,jsnpth);
if nargout == 0; fprintf("Tethys params written to:  "+jsnpth+"\n"); 
else vao{1} = jsnpth; end
end
