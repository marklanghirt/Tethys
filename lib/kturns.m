function [zmin,kozmin,zmax,kozmax] = kturns(ppz,ppk,ppdim,turnk)
%[zmin,cozmin,zmax,cozmax]=kturns(ppz,ppk,ppdim,turnk)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: finds vertexing/turning points of a piecewise-polynomial k-profile
%
% INPUTS:
%  ppz: array of depth values
%  ppk: array of wavenumber values
%  turnk: array of vertexing wavenumbers
% OUTPUTS: 
%  zmin: minimum z-values where k=t
%  kozmin: k of min z-values where k=t
%  zmax: maximum z-values where k=t
%  kozmax: k of max z-values where k=t


z = ppz; k=ppk; t = turnk; szz=size(ppz); szk=size(ppk); szt=size(t);
d = ppdim; ndz=numel(szz); ndt=numel(szt); md=max(ndz,ndt);
if sum(szz~=szk); error('ppz and ppk must have same size'); end
szz=[szz,ones([1,md-ndz])]; szt=[szt,ones([1,md-ndt])];
if szt(d) ~= 1; error('turnk must be singleton in ppdim'); end
zrepsize = ceil(szt./szz); trepsize=ceil(szz./szt); trepsize(d)=1; 
try z=repmat(z,zrepsize); k=repmat(k,zrepsize); t=repmat(t,trepsize);
catch; error('incompatible size(ppz) and size(turnk)'); end 
z=permshift(z,md,d-1); nz=szz(d); %znn=~isnan(z);
np=nz-1; %iz = (1:nz).'; ip = (1:np).';
k=permshift(k,md,d-1); t=permshift(t,md,d-1); 
szz=size(z); szt=[1,szz(2:end)]; %n = prod(szz(2:end));
z=z(:,:); z1 = z(1:np,:); z2 = z(2:nz,:); 
k=k(:,:); k1 = k(1:np,:); k2 = k(2:nz,:); t=t(:,:);
[ifz,fz] = dimFind(z,1,1,'first','num');
[ilz,lz] = dimFind(z,1,1,'last','num');
fk = dimIndex(k,1,ifz); lk = dimIndex(k,1,ilz);
dk = diff(k,1,1); pdk = dk>=0; ndk = dk<0;
tgek = t>=k1; tlek = t<=k2; tink = tgek&tlek;
ifp = dimFind(tink&pdk,1,1,'first'); 
iln = dimFind(tink&ndk,1,1,'last');
ifpin = isnan(ifp); ifpwon = ifp; ifpwon(ifpin) = ifz(ifpin);
ilnin = isnan(iln); ilnwon = iln; ilnwon(ilnin) = ilz(ilnin)-1;
%htz1 = dimIndex(z1,1,ifpwon); htz2 = dimIndex(z2,1,ifpwon);
%ltz1 = dimIndex(z1,1,ilnwon); ltz2 = dimIndex(z2,1,ilnwon);
%htk1 = dimIndex(k1,1,ifpwon); htk2 = dimIndex(k2,1,ifpwon);
%ltk1 = dimIndex(k1,1,ilnwon); ltk2 = dimIndex(k2,1,ilnwon);
ltz1 = dimIndex(z1,1,ifpwon); ltz2 = dimIndex(z2,1,ifpwon);
htz1 = dimIndex(z1,1,ilnwon); htz2 = dimIndex(z2,1,ilnwon);
ltk1 = dimIndex(k1,1,ifpwon); ltk2 = dimIndex(k2,1,ifpwon);
htk1 = dimIndex(k1,1,ilnwon); htk2 = dimIndex(k2,1,ilnwon);
htz = (htz2-htz1)./(htk2-htk1).*(t-htk1) + htz1;
ltz = (ltz2-ltz1)./(ltk2-ltk1).*(t-ltk1) + ltz1;
htz(ilnin) = lz(ilnin); htk = t; htk(ilnin) = lk(ilnin);
ltz(ifpin) = fz(ifpin); ltk = t; ltk(ifpin) = fk(ifpin);
htz = permshift(reshape(htz,szt),md,md-d+1); zmax = htz;
htk = permshift(reshape(htk,szt),md,md-d+1); kozmax = htk;
ltz = permshift(reshape(ltz,szt),md,md-d+1); zmin = ltz; 
ltk = permshift(reshape(ltk,szt),md,md-d+1); kozmin = ltk;