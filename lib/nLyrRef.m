function ref = nLyrRef(f,thb,geo)
%[ref]=nLyrRef(f,thb,geo) 
%
% ORIGINAL AUTHOR: Charles W. Holland - June 2003
% EDITED AND ADAPTED BY: Mark Adam Langhirt
% COPYRIGHT: ©2025, Charles Holland, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: compute the plane wave reflection coefficient for arbitrary number of isospeed layers
% -> from Brekovskikh/Godin Waves in Layered Media I Eq 2.5.1-2.5.4  
% Eg: [ref] = nLyrRef([200 1000 4000],1:90,Nx4geo)
%
% INPUTS:
%  "f" is frequencies
%  "thb" is grazing angle in radians
%  "geo" contains environmental data, e.g., 
%    thickness    speed    attenuation    density
%      (m)        (m/s)    (dB/m/kHz)      (g/cc)
%  geo=[NaN        1511          0          1.03
%      1.7         1500        0.015        1.45
%      2.8         1555         0.1         1.7 
%      0.6         1600         0.1         2.0
%      NaN         1565         0.01        1.9];
%
% NB: ref is the (complex) plane wave pressure reflection coefficient (BL =-20 log10 (ref))
% Aug. 16, 2016: eliminated do loop on frequency by adding another dimension 
%  to several matrices, minor code cleanup


nf=length(f); nc=length(geo(:,2)); nt=length(thb);
dB2nep=2*pi*20*1000/log(10);    
d=geo(2:end-1,1); c=geo(:,2); a=geo(:,3); r=geo(:,4);
nlyr=nc-1; v=1./(1./c + a*1i/dB2nep); % force radiation condition to be satisfied

zz=zeros(nc,nt); th=zz; %th(1,:)=th1 .ne. 0 as is; but never used at th(1,)
zz(1,:)=r(1)*v(1)./sin(thb);  %changed to this on 22 Dec 2006 from %zz(1,:)=r(1)*c(1)./sin(th1); 
for m=2:length(c); th(m,:)=acos( cos(thb)*v(m)/v(1) ); %changed to this on 22 Dec 2006 from %th(m,:)=acos( cos(th1)*v(m)/c(1) ); 
zz(m,:)=r(m)*v(m)./sin(th(m,:)); end %could do this in vector math e.g., thm=acos(v/v(1)*cos(th1));which would give correct values at m=1

z=repmat(zz,1,1,nf);
k=(2*pi*repmat(f,nlyr+1,1)./repmat(v,1,nf)).';

if nlyr==1 % sediment halfspace 
  rfh=(z(2,:,:) - z(1,:,:) )./(z(2,:,:) + z(1,:,:));
  if nf==1;ref=rfh.'; else; ref=squeeze(rfh); end

else %COMPUTE input impedance
  tankd= tan( repmat(k(:,nlyr),1,nt).*d(nlyr-1).*repmat(sin(th(nlyr,:)),nf,1 ) );
  znlay=reshape(z(nlyr,:,:),nt,nf);  znlay1=reshape(z(nlyr+1,:,:),nt,nf);
  zin = znlay.*(znlay1-1i*znlay.*tankd.')./(znlay - 1i*znlay1.*tankd.');
  for m=nlyr:-1:3
    tankd = tan(repmat(k(:,m-1),1,nt).*d(m-2).*repmat(sin(th(m-1,:)),nf,1));
    zm1 = reshape(z(m-1,:,:),nt,nf);
    zin = zm1.*(zin-1i*zm1.*tankd.')./(zm1-1i*zin.*tankd.'); end
  z1=reshape(z(1,:,:),nt,nf);
  ref=(zin-z1)./(zin+z1); end
ref=ref.'; ref(isnan(ref))=-1; %at 0 degrees, reflection coefficeint = -1