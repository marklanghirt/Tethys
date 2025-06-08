function [Ims, Imb] = lloydMirrorDP(f,z,c,H,th,Rs,Rb,nLam)
%[Ims,Imb]=lloydMirrorDP(f,z,c,H,th,Rs,Rb,[nLam])
%
% ORIGINAL AUTHOR: Charles W. Holland
% EDITED AND ADAPTED BY: Mark Adam Langhirt
% COPYRIGHT: ©2025, Charles Holland, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: Computes Lloyd mirror pattern for sea surface And bottom boundaries 
%  Assumes isoveloity, i.e. no refraction 
%  sea surface reflection interference AND 
%  bottom reflection (downgoing paths) relative to source 
%  this is equivalent to having two integrals one 0->pi/2 and one -pi/2->0 that are symmetric except for the Lloyd mirror effects
%  must take out factor of two for symmetry (is already taken care of in sin(...) expression
%
% INPUTS:
%  freq [Hz] (scalar)
%  c local sound speed [m/s] (vector dimN2)
%  z source or receiver depth [m] (vector dimN3)
%  H water depth [m] (vector dimN2)
%  th grazing angle [rad] corresponding to Rs, Rb, (dimN1 and dimN2)
%  Rs surface pressure reflection coefficient (complex) (vector dimN2)
%  Rb bottom pressure reflection coefficient (complex) (vector dimN2)
%
% OUTPUTS: (apply to intensity or pressure-squared)
%  Im the product of sea surface and seabed boundary intereference
%  Ims the sea surface intereference (Lloyd mirror)
%  Imb the seabed interference


arguments
  f {mustBeNumeric} % dimN1
  z {mustBeNumeric} % dimN2
  c {mustBeNumeric} % dimN3
  H {mustBeNumeric} % dimN3
  th {mustBeNumeric} % dimN3 and dimN4
  Rs {mustBeNumeric} % dimN3 and dimN4
  Rb {mustBeNumeric} % dimN3 and dimN4
nLam {mustBeNumeric} = inf; end
%% Validation
% make sure singleton dimensions are replicated between arrays that share dims
% make sure only (c,H) and (th,Rs,Rb) can share only 1 dim
fsz=size(f); zsz=size(z); Hsz=size(H); Rssz=size(Rs); Rbsz=size(Rb);
c=c.*ones(Hsz); csz=size(c); H=H.*ones(csz); Hsz=size(H);
th=th.*ones(Rssz).*ones(Rbsz); tsz=size(th); 
Rs=Rs.*ones(tsz).*ones(Rbsz); Rssz=size(Rs); 
Rb=Rb.*ones(tsz).*ones(Rssz); Rbsz=size(Rb);
ndf=numel(fsz); ndz=numel(zsz); ndc=numel(csz); ndH=numel(Hsz); ndt=numel(tsz); 
ndRs=numel(Rssz); ndRb=numel(Rbsz); md = max([ndf,ndz,ndc,ndH,ndt,ndRs,ndRb]); 
fsz=[fsz,ones(1,md-ndf)]; zsz=[zsz,ones(1,md-ndz)]; csz=[csz,ones(1,md-ndc)]; 
Hsz=[Hsz,ones(1,md-ndH)]; tsz=[tsz,ones(1,md-ndt)];
Rssz = [Rssz,ones(1,md-ndRs)]; Rbsz = [Rbsz,ones(1,md-ndRb)]; 
sz = [fsz;zsz;csz;Hsz;tsz;Rssz;Rbsz]; ld = sz>1; id = 1:md; vd = cell(7,1);
for ii=1:7; vd{ii,1}=id(ld(1,:)); end; nd = sum(ld,2);
nm = ["f","z","c","H","th","Rs","Rb"]; ndlim = [1,1,1,1,2,2,2];
for ii = 1:7; if nd(ii)>ndlim(ii)
error("%s can have numel()>1 for %.0f dim(s)",nm(ii),ndlim(ii)); end; end
if any(ld(1,:).*(ld(2,:)+ld(3,:)+ld(5,:))); error("f cannot share dim with others"); end
if any(ld(2,:).*(ld(3,:)+ld(5,:))); error("z cannot share dim with others"); end
if nd(5)==2; if sum(ld(4,:).*ld(5,:))~=1; error("2D (th,Rs,Rb) must share a dim with (c,H)"); end
if nd(4)==1; if sz(4,vd{4})~=sz(5,vd{4}); error("(thd,Rs,Rb) and (c,H) must have same size in shared dim"); 
end; end; end

%% Calculation
% avoid calculating when k*z > kdLim, or numWavelengths > kdLim/2pi
H = H.*ones(zsz); zUp = z.*ones(Hsz); % must replicate for difference (H-z)
zDn = H-zUp; blwMsk = zDn<0; kdLim = 2.*pi.*nLam; k = 2.*pi.*f./c; 
Ims = 1./2.*abs(1+Rs.*exp(1i.*2.*k.*zUp.*sin(th))).^2; szIs = size(Ims);
Imb = 1./2.*abs(1+Rb.*exp(1i.*2.*k.*zDn.*sin(th))).^2; szIb = size(Imb); %corrected 12 April 2010 phase consistent with Refl Coeff. phase
upMsk = logical((k.*zUp > kdLim).*ones(szIs));
dnMsk = logical((k.*zDn > kdLim).*ones(szIb));
blwMsk = logical(blwMsk.*ones(szIs)); if szIs ~= szIb; error("szIs ~= szIb"); end
Ims(upMsk)=1; Ims(blwMsk)=0; Imb(dnMsk)=1; Imb(blwMsk)=0;