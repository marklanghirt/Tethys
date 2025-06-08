function [intCz]=intCz_pwLinC(ppz,ppc,ppdim,cr)
%[intCz]=intCz_pwLinC(ppz,ppc,ppdim,cr)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: Integrates c_z profiles from piecewise linear c(z) profiles.  
%  This function will integrate over the entire profiles and sum the real parts omitting nans.
%
% INPUTS:
%   ppz: [real nd-array] depths/nodes/breakpoints of piecewise polynomials (must be linear, 1st-order)
%   ppc: [real nd-array] speeds/values/samples of piecewise polynomials (must be linear, 1st-order)
%   ppdim: [integer scalar] the array dimension along which the piecewise polynomials extend
%   cr: [real nd-array] the constant speeds/values to subtract from c, i.e. 1/cz = sqrt((1/c)^2-(1/cr)^2);
%         Cr must be singleton in ppdim (constant for the profile); but it can be different for all 
%         other dimensions/elements of the pp-arrays, and will be duplicated to the sizes of ppz/ppc
%         when any dimension besides ppdim is also singleton, ie. binary array expansion.
% OUTPUTS:
%   intCz: results from integrating along ppdim, note: w is omega = 2 pi freq
%          Int(w/kz)dz = Int(cz)dz = Int( 1/sqrt( (1/c)^2-(1/cr)^2 ) )dz


%% Section 1: Binary array expansion and validation of input parameters
z = ppz; c = ppc; szz = size(z); szc = size(c); szcr = size(cr);
d = ppdim; ndz = numel(szz); ndc = numel(szcr); md = max(ndz,ndc);
if sum(szz~=szc); error('ppz and ppc must have same size'); end
szz = [szz,ones([1,md-ndz])]; szcr = [szcr,ones([1,md-ndc])];
if szcr(d)~=1; error('cr must be singleton in ppdim'); end
zrepsize = ceil(szcr./szz); crepsize = ceil(szz./szcr); crepsize(d) = 1;
try z = repmat(z,zrepsize); c = repmat(c,zrepsize); cr = repmat(cr,crepsize);
catch; error('incompatible size(ppz) and size(cr)'); end

%% Section 2: Prepare arrays for the calculation (reshapes and predefines intermediaries)
%  ppdim becomes dim1, all other dimensions are extruded into dim2, so each column is a separate profile
z = permshift(z,md,d-1); z = z(:,:); % pp depths/nodes
nz = size(z,1); ni = nz-1; %  nz: num nodes;   ni: num intervals
c = permshift(c,md,d-1); c = c(:,:); % pp speeds/values
cr = permshift(cr,md,d-1); szc = size(cr); % intermediary reshape size before returning calculations
cr = cr(:,:); nc = size(cr,2); % constant speeds/values orthogonal to cz

%  intermediary variables for simplifying formulas
z1 = z(1:ni,:); z2 = z(2:nz,:); dz = z2-z1; 
c1 = c(1:ni,:); c2 = c(2:nz,:); dc = c2-c1; dcdz = dc./dz;
g1 = sqrt(cr.^2-c1.^2); g2 = sqrt(cr.^2-c2.^2); % these are convenience variables, from algebraic manipulations

%% Section 3: Calculation routine
intCzArr = zeros([ni,nc]); % pre-populate calculation array
iIso = dcdz==0; % logical mask for isovelocity intervals, a limiting case with a separate formula
iLin = ~iIso; % logical mask for linear gradient intervals
iInf = isinf(cr); % logical mask for infinite orthogonal speed component, another limiting case with a separate formula
iLI = iLin & iInf; iLF = iLin & ~iInf; % masks for linear gradient with infinite/finite cr
cr = repmat(cr,[ni,ones(1,md-1)]); % duplicate cr along profiles for direct array arithmetic
% formulas for integrating Int(cz)dz in each segment of profiles
intCzArr(iIso) = cr(iIso).*c1(iIso).*dz(iIso)./g1(iIso);    % when c is isovelocity
intCzArr(iLI) = (dcdz(iLI)./2.*dz(iLI) + c1(iLI)).*dz(iLI); % when cr is infinite
intCzArr(iLF) = cr(iLF)./dcdz(iLF).*(g1(iLF) - g2(iLF));    % when cr is finite

%% Section 4: Reshape calculation results to match input sizes before returning
% Sums along dim1 only the real-values and omits NaNs, then reshape to intermediary size szc
intCz = reshape(sum(real(intCzArr),1,'omitnan'),szc);
intCz = permshift(intCz,md,md-(d-1)); % finally reshape to match original input array sizes
end % intCz_pwLinC()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [B] = permshift(A,ndim,shift)
% PERMSHIFT: a version of shiftdim() using permute()
%   This treats array dimensions in a more consistent manner to:
%     -treat 1d-vectors with the same manipulation rules as nd-arrays
%     -avoid collapsing leading singleton dimensions
% INPUTS:
%   A: [anyType nd-array] the array to rotate
%   ndim: [integer scalar] the max number of intended dimensions in the array
%   shift: [integer scalar] the number of dimensions by which to rotate the array
% OUTPUTS:
%   B: [anyType nd-array] the rotated array
B = permute(A,circshift(1:ndim,-shift,2));
end % permshift()