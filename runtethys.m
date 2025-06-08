function [] = runtethys(name)
% []=runtethys(name)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: a wrapper routine to load parameters for and execute "Tethys",
% a three-dimensional solid-angle energy flux ocean acoustic propagation model

%% directories and paths
mflnm = mfilename("fullpath"); 
tetdir = strrep(mflnm,"runtethys","");
prmdir = tetdir + "prm" + filesep;  
datdir = tetdir + "dat" + filesep;
prmpth = prmdir + name + ".json";
datpth = datdir + name + ".mat";

%% load input parameters from json file
opt = readstruct(prmpth,"StructSelector","/opt");
env = readstruct(prmpth,"StructSelector","/env");
src = readstruct(prmpth,"StructSelector","/src");
rcv = readstruct(prmpth,"StructSelector","/rcv");

%% OPT SETUP
opt.name = name;
opt.yFF = 1 + 2*opt.ymmd; % gaussian focusing for y-coherence
opt.zFF = 1 + 2*opt.zmmd; % gaussian focusing for z-coherence

%% ENV SETUP
% put geo matrix back together
tmpgeo = env.geo; env = rmfield(env,"geo");
env.geo = [tmpgeo.z(:),tmpgeo.cp(:),tmpgeo.ap(:),tmpgeo.rho(:)];

%% SRC SETUP
% H: water depth at source
src.H = linppval(env.y(:),env.H(:),1,src.y);
% THETA: elevation angle array
if src.pthDegMin < src.thEps; src.pthDegMin = src.thEps; end
if src.pthDegMax > (90-src.thEps); src.pthDegMax = 90-src.thEps; end
% finfreq
switch opt.finFreq;
case 1; % srf=vac & bot=rgd:  1/4 * lam_z0 = H0
  ffElvMin = asind(env.c0/4/src.f/src.H);
case 2; % halfway wavelength:  3/8 * lam_z0 = H0
  ffElvMin = asind(3*env.c0/8/src.f/src.H);
case 3; % srf=vac & bot=vac:  1/2 * lam_z0 = H0
  ffElvMin = asind(env.c0/2/src.f/src.H);
otherwise; % no finite frequency limiting for min(theta_0)
  ffElvMin = pthDegMin; end
[src.thN,src.thd,~] = arrUnfold(2.*src.pthN,max(src.pthDegMin,ffElvMin),src.pthDegMax);
src.pthDeg = src.thd(src.pthN+1:end); src.pthDegDif = src.pthDeg(2) - src.pthDeg(1);
src.pth = src.pthDeg/180*pi; src.th = src.thd/180*pi; src.pthDif = src.pthDegDif/180*pi;
% PHI: azimuthal angle array
if src.pphDegMin < src.phEps; src.pphDegMin = src.phEps; end
if src.pphDegMax > (90-src.phEps); src.pphDegMax = 90-src.phEps; end
src.pphDegDif = (src.pphDegMax - src.pphDegMin) ./ (src.pphN-1); % [Deg]
src.pphDeg = src.pphDegMin:src.pphDegDif:src.pphDegMax; % [Deg]
src.pph=src.pphDeg./180.*pi; src.pphDif=src.pphDegDif./180.*pi; % [Rad]
src.ph = [-flip(src.pph),src.pph]; src.phd=src.ph.*180./pi;
%% RCV SETUP
% X: longitudinal position array
rcv.x = linspace(rcv.xMin,rcv.xMax,rcv.xN);
rcv.xDif = (rcv.x(end)-rcv.x(1))/(numel(rcv.x)-1);
% Y: transverse position array
rcv.y = linspace(rcv.yMin,rcv.yMax,rcv.yN);
rcv.yDif = (rcv.y(end)-rcv.y(1))/(numel(rcv.y)-1);
% Z: vertical position array
rcv.z = linspace(rcv.zMin,rcv.zMax,rcv.zN); 
rcv.zDif = (rcv.z(end)-rcv.z(1))/(numel(rcv.z)-1);

%% Run main calculation script
tic; [fld,ray] = tethys(env,src,rcv,opt); runtime.secs = toc;
[runtime.dhms(1),runtime.dhms(2),runtime.dhms(3),runtime.dhms(4),runtime.sstr] = S2DHMS(runtime.secs);
runtime.lstr = sprintf('%d days   %02d:%02d:%04.1f (HH:MM:SS.s)',runtime.dhms(1),runtime.dhms(2),runtime.dhms(3),runtime.dhms(4));
fprintf('Elapsed Time:   %s\n',runtime.lstr);
%% Save results from main calculation
fprintf("Saving data: %s \n", datpth);
save(datpth,'rcv','src','env','opt','fld','ray','runtime');
