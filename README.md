# Tethys: A three-dimensional solid-angle energy flux ocean acoustic propagation model

This code is for research purposes and is distributed as open-source software, with absolutely no warranty or guarantee of any kind.

COPYRIGHT: ©2025, Mark Adam Langhirt \
LICENSE: This research code is distributed under the MIT license, a copy of which is included in the code's repository or alternatively can be viewed at: [https://opensource.org/license/MIT](https://opensource.org/license/MIT).

---

This model is the culmination of a Ph.D. dissertation project.  The full derivation of this model and a detailed explanation of the physical and mathematical concepts is written up in a pre-print draft manuscript and hosted on arXiv at the following link: [http://arxiv.org/abs/2506.03325](http://arxiv.org/abs/2506.03325).  References used in the model's derivation are provided in the manuscript.

---

## Overview

![Example horizontal slice of transmission loss field at 30m depth in the canonical ideal wedge problem](/pic/wej1zslcdbg_Cy_fld2d.png)
*Example horizontal slice of transmission loss field at 30m depth in the canonical ideal wedge problem*

Tethys is an ocean acoustic propagation model built upon energy flux methods analytically extended to handle three-dimensional propagation problems with asymmetric environmental range-dependence.  The model's primary output is the semi-coherent transmission loss field sampled on an arbitrary receiver grid in Cartesian coordinates (x,y,z).  A demonstration of simulated 3D adiabatic ray trajectories for the canonical wedge problem using the cycle-tracking method is also included.

The model's derivation may be briefly summarized as follows. A solution to the 3D range-dependent Helmholtz equation is expanded as an adiabatic double Wentzel-Kramers-Brillouin (WKB) mode sum.  The discrete mode sum is transformed into an integration across a continuum of modes, and considerations of the model's differential topology permits deriving a Jacobian to transform differential modenumbers to differential propagation angles.
 The continuum of vertical and transverse modes is mapped bijectively to the source propagation angles in elevation and azimuth, and the adiabatic approximation enables smooth mappings of propagation angles via the ray invariant as the waves travel throughout the waveguide to the receiver positions.  Inter-modal interference from the double mode sum cross product is partially retained and approximated with Taylor series expansions to form a convergence factor integration kernel, which acts as an eigenray filter to focus the modal energy when a particular ray/mode completes the appropriate number of transverse and vertical cycles to arrive at a receiver.  Additional integration kernels are also incorporated to: simulate mode stripping as waves pass through the transparent boundaries of the transverse computational domain, accumulate bottom attenuation along the horizontal trajectory of the transverse wave cycles, and apply Lloyd mirror type directivity patterns to the source and receivers due to boundary reflection with the surface and seafloor.  The end result is an integration of modal amplitudes across source propagation solid-angles that accounts for geometric spreading, horizontal refraction with convergence zones and shadow zones, and bottom loss that is accumulated along each unique transverse cycle trajectory.

This model is implemented in Matlab, and the computations are primarily performed as vectorized operations on multidimensional arrays.  Sections of the computation that are particularly intensive have been isolated for parallelization, but the user has the option to perform these calculations using fully vectorized operations on a single process, for loops on a single process, or parfor loops on a parallel pool of processes or threads.  There are several types of computational grids used by the routine: the transverse environmental grid used to specify the bathymetry profile and transverse boundaries of the computational domain, the source angular grid of propagation angles over which the main calculation is integrated, the Cartesian grid of receiver coordinates at which to calculate the transmission loss field, and some additional temporary sampling grids for interpolation of wavenumber profiles and partial cycle integrals.  All interpolations and mappings are piecewise linear, and the final integration over solid-angles is computed with a composite Simpson's rule quadrature scheme.

There are two wrapper functions included, and these are the main functions that a user should call.  These routines expect one argument called "name", which is the filename of the JSON or MAT file to load without its extension. \
  (1) "*runtethys.m*" for loading parameters from JSON files, calling the main Tethys calculation routine, and saving the results \
  (2) "*plottethys.m*" for plotting 2D transmission loss field slices

Example invocation: \
  runtethys("myTestCase") \
  plottethys("myTestCase")

There is a separate included routine called "*isoWedge_quadRay_cycleTest.m*", that demonstrates the cycle-tracking method to emulate ray trajectories using the adiabatic wave cycles.  This script uses simplified analytic expressions that are particular to the ideal wedge environment, and thus is not intended as a general purpose routine for generating wave trajectories.  Trajectories for generalized environments can be extracted from the main Tethys calculation routine's cycle-tracking algorithm, but this feature is still in development (there may be bugs).

![Demonstration of adiabatic cycle trajectories in the ideal wedge environment generated by isoWedge_quadRay_cycleTest.m](/pic/quadRayWedgeDemo.png)
*Demonstration of adiabatic cycle trajectories in the ideal wedge environment with a rigid wall at y=10km, generated by isoWedge_quadRay_cycleTest.m*

---

## Capabilities and Limitations

The model implemented here does not include all of the generality of what the analytic derivation can capture, and the typical limitations of energy flux modeling also apply.

![Example of environmental geometry and coordinates](/pic/environmentCoords.png)
*Example of environmental geometry and coordinates*

- Capabilities:
  - This energy flux model computes the transmission loss field at arbitrary receiver locations and captures horizontal refraction due to transverse range-dependence.  It is an approximate frequency-domain solution to a 3D range-dependent Helmholtz equation and may be interpreted as the locally averaged harmonic acoustic intensity subject to horizontal refraction.
  - Theoretical capability for environmental variability from then analytic derivation:
    - The environment is assumed to be longitudinally range-independent, i.e. translational symmetry along x, the forward propagating direction
    - Vertical depth-dependence and transverse range-dependence for water-column profile: k(z|y)
    - Transverse range-dependence for the bathymetry profile, H(y), and bottom reflection coefficient, R(y)
  - Implementation's capability for environmental variability for demonstration purposes:
    - Water-column is assumed to be isovelocity
    - Transverse range-dependence for bathymetry profile only: H(y)
    - Sediment properties are assumed to be homogeneous across the environment
  - A WKB amplitude limiter is included to cap WKB modal amplitudes based roughly on the peak amplitude of the Airy function solution in the vicinity of the turning points.
  - A convergence factor integration kernel acts as an eigenray filter to focus modal energy when the adequate number of wave cycles are traversed for a particular receiver location.
  - Transverse boundaries of the computational domain are transparent, and the model has a mechanism for stripping modal energy from the solution after it has escaped the computational domain.  The combined effect of the transverse mode stripping kernel and the convergence factor kernel is the predominant mechanism simulating geometric spreading.
  - A bottom attenuation kernel accumulates losses by a product integral along the horizontal transverse cycle trajectory, y(x).
  - A Lloyd's mirror directivity kernel adds boundary reflection interference at the source and receiver positions.
- General Limitations:
  - Energy flux models compute an incoherent or semi-coherent transmission loss field by design for efficiency, and thus are incapable of resolving fine-scale interference structure or fully coherent phase fronts. Some of the inter-modal coherence on the order of the cycle distances is incorporated in the convergence factor in order to capture refractive convergence zones and shadow zones, hence the model is considered semi-coherent.
  - All bulk, reduced, and effective wavenumber profiles are assumed to be convex. This is to avoid dealing with bimodal potential wells that would represent bifurcating ray trajectory families or modefunctions that exhibit wave tunneling across a potential barrier.  This also simplifies checking for wave arrivals at the boundaries of the computational domain.
  - Since the environment is assumed to be range-independent in the forward longitudinal x-direction, there is no mechanism in this model that can account for total refraction in the backward direction, i.e. no back propagation.
  - Additional physical acoustic propagation phenomena that are not incorporated in this model derivation include: diffraction, wave tunneling, volume attenuation, volume and boundary scattering, reverberation, and acousto-elastic wave mechanics.

---

## Input Parameters

The model ingests a specific set of input parameters and options that are packaged into structs loaded from a JSON file.  There is a template of the JSON file included called "*template.json*", a selection of JSON files for the canonical ideal wedge environment and lossy penetrable wedge environment, and an example of how to generate these JSON files from the structs is provided in "*tethysparams.m*".  The following is a list of the model's input parameters with descriptions.

- "opt" (struct):
  - "clear" (bool):  clears intermediary variables during execution to save memory
  - "debug" (bool):  pauses executation at the end of the main computation
  - "parType" (string):  type of parallelization: "vec","for","background","process","thread","gpu"; see "chspool.m"
  - "omitCore" (bool):  omit one core when using thread/process parallelization
  - "calcRays" (bool):  simulate ray trajectories and include in output data
  - "ntRays" (int):  number of theta (elevation) angles for rays
  - "npRays" (int):  number of phi (azimuth) angles for rays
  - "botLoss" (bool):  calculate accumulated bottom loss along horizontal trajectories
  - "ntSamp" (int):  number of theta sample points for interpolations
  - "nxSamp" (int):  number of x-range sample points for interpolations
  - "nySamp" (int):  number of y-range sample points for interpolations
  - "nIsoLim" (int):  uniform max array length for iso-modenumber k_xy-profiles
  - "dirKern" (bool):  calculate Lloyd Mirror directivity kernels for src/rcv near srf/bot
  - "maxLamb" (num):  max wavelength distance from boundary for lloyd mirror directivity kernel
  - "convType" (int):  type of convergence factor calculation: (0)off, (1)periodicGaussian, (2)dirichletKernel/cosineSummation
  - "yConv" (bool):  compute transverse convergence factor
  - "ymmd" (num):  transverse max modenumber difference for convergence factor
  - "zConv" (bool):  compute vertical convergence factor
  - "zmmd" (num):  vertical max modenumber difference for convergence factor
  - "wkbLim" (bool):  limit WKB modal amplitudes using rough approx. of Airy function peak
  - "finFreq" (int):  limit |theta|_min according to finite freq. approx: (0)off, (1)vac-rgd, (2)median, (3)vac-vac
  - "splitTh" (bool):  split theta (elevation) angle domain in final integration
  - "splitPh" (bool):  split phi (azimuth) angle domain in final integration
- "env" (struct):
  - "c0" (num):  watercolumn sound speed in [m/s]
  - "y" [(num),...,(num)]:  transverse y-range nodes for bathymetry profile
  - "H" [(num),...,(num)]:  transverse depth values for bathymetry profile
  - "geo" (struct):
    - "z" [NaN,(num),...,(num),NaN]: list of sediment layer thicknesses, first (watercolumn) and last (basement) are always NaN
    - "cp" [(num),(num),...,(num),(num)]:  compressive sound speed of layers in [m/s]
    - "ap" [(num),(num),...,(num),(num)]:  compressive attenuation of layers in [dB/kmHz]
    - "rho" [(num),(num),...,(num),(num)]:  medium density of layers in [g/cc]
- "src" (struct):
  - "thEps" (num):  epsilon (small number) theta (elevation) buffer to avoid integrating over 0deg and 90deg
  - "phEps" (num):  epsilon (small number) phi (azimuth) buffer to avoid integrating over 0deg and 90deg
  - "monoAmp" (num):  monopole amplitude, usually 4*pi to normalize pressure to 1Pa at 1m distance from source
  - "f" (num):  temporal frequency of acoustic source signal
  - "x" (num):  x-range position of source in [m], not used...
  - "y" (num):  y-range position of source in [m]
  - "z" (num):  z-depth position of source in [m]
  - "pthN" (int):  number of absolute theta (elevation) angles for angular integration domain
  - "pthDegMin" (num):  absolute theta (elevation) angle minimum for angular integration domain
  - "pthDegMax" (num):  absolute theta (elevation) angle maximum for angular integration domain
  - "pphN" (int):  number of absolute phi (azimuth) angles for angular integration domain
  - "pphDegMin" (num):  absolute phi (azimuth) angle minimum for angular integration domain
  - "pphDegMax" (num):  absolute phi (azimuth) angle maximum for angular integation domain
- "rcv" (struct):
  - "xN" (int):  number of longitudinal x-range positions for receiver grid
  - "xMin" (num):  minimum longitudinal x-range position in [m] for receiver grid
  - "xMax" (num):  maximum longitudinal x-range position in [m] for receiver grid
  - "yN" (int):  number of transverse y-range positions for receiver grid
  - "yMin" (num):  minimum transverse y-range position in [m] for receiver grid
  - "yMax" (num):  maximum transverse y-range position in [m] for receiver grid
  - "zN" (int):  number of vertical z-depth positions for receiver grid
  - "zMin" (num):  minimum vertical z-depth position in [m] for receiver grid
  - "zMax" (num):  maximum vertical z-depth position in [m] for receiver grid
