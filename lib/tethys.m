function [fld,rays] = tethys(env,src,rcv,opt)
%[fld,rays]=tethys(env,src,rcv,opt)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
% 
% DESCRIPTION: Main calculation routine for the Tethys three-dimensional 
%  solid-angle energy flux ocean acoustic propagation model
%
% 20230625: 3D solid-angle Tethys split into vectorized/parallel sections
%  (I)  VECtorized pre-calculations evaluated to/from 3D-arrays or smaller.
%  (II) (PAR)FOR loop over receiver coordinates, where calculations that evaluate to
%       or are dependent on 4D-arrays or larger are computed.
%  NB:  For the final integration, theta and phi dependence must be included in
%       the arrays stored in memory, but Rcvr locations are independent.
% COMPUTATION DIMS:
%   dim 1: theta (src)
%   dim 2: phi (src)
%   dim 3: x (rcv)
%   dim 4: y (rcv)
%   dim 5: z (rcv)
%   dim 6: extra (interpolation)
%   dim 7: extra (interpolation)

% Initialize Options and Parameters
mainTic = tic; 
fprintf('%s:  Initializing parameters, constants, and arrays\n',S2DHMS(toc(mainTic)));
dbg = opt.debug; clr = opt.clear; ray = opt.calcRays; mlam = opt.maxLamb;
prtp = opt.parType; omcr = opt.omitCore;
Np = opt.ymmd; Nm = opt.zmmd; Fp = opt.yFF; Fm = opt.zFF; 
Copt = opt.convType; CyOn = opt.yConv; CzOn = opt.zConv;
botloss = opt.botLoss; dirKern = opt.dirKern; 
wkblim = opt.wkbLim; splitth = opt.splitTh; splitph = opt.splitPh;
qtN = opt.ntSamp; % sampling/query theta length for iso-Iz mapping and bottom loss
qyN = opt.nySamp; % sampling/query y-rng length for iso-Iz mapping and transverse y-cycle profiles (qLyL, etc.)
qxN = opt.nxSamp; % sampling/query x-rng length for cycle-tracking calculation
isoN = opt.nIsoLim; % uniform array length for iso-Iz (rayInv/modenum) kxy(y)-profiles

% Initialize Arrays and Enforce Dimensions
st = src.th(:);                  % dim 1
std = src.thd(:);                %#ok<NASGU>
sigt = sign(st);                 % dim 1
negt = sigt < 0;                 % dim 1
post = ~negt;                    % dim 1
sp = permute(src.ph(:),[2,1]);   % dim 2
spd = permute(src.phd(:),[2,1]); %#ok<NASGU>
sigp = sign(sp);                 % dim 2
negp = sigp < 0;                 % dim 2
posp = ~negp;                    % dim 2
rx = permute(rcv.x(:),3:-1:1);   % dim 3
ry = permute(rcv.y(:),4:-1:1);   % dim 4
rz = permute(rcv.z(:),5:-1:1);   % dim 5
ey = permute(env.y(:),6:-1:1);   % dim 6
eH = permute(env.H(:),6:-1:1);   % dim 6
% Initialize Input Variables and Independent Constants
tN = numel(st); xN = numel(rx); yN = numel(ry); pi2 = 2.*pi; pN = numel(sp); 
zN = numel(rz); sy = src.y; sz = src.z; dth = src.pthDif; dph = src.pphDif;
ec = env.c0; egeo = env.geo; f = src.f; w = pi2.*f;
k = w./ec; kk = k.^2; sS = src.monoAmp; sS2 = sS.^2;
FFp_2pi = Fp.^2./pi2; FFm_2pi = Fm.^2./pi2;
% Start pool if necessary
[pp,parType,pfW] = chspool('type',prtp,'omitcore',omcr);


%% Phase I Calculation (vectorized) >>>
%% Specify kxy0(|m) to calculate kz(z|y,m) and kxy(y|m) Profiles
fprintf('%s:  Calculating source horiz. k_xy and sampling ray invariant\n',S2DHMS(toc(mainTic)));
sH = linppval(ey(:),eH(:),1,sy);
if sz >= sH; error('z_s >= H_s'); end
rH = linppval(ey,eH,6,ry); bMsk = rz>rH;
skxy = k.*cos(st); skz = k.*sin(st); %skkxy = skxy.^2;
qy = linspace(ey(1),ey(end),qyN).'; % dim1 (independent) for mkisopp
qkxy = k-k.*logspace(-6,0,qtN); % dim2 (dependent) for mkisopp
qH = linppval(ey,eH,6,qy);
qIz = sqrt(kk-qkxy.^2).*qH;
fprintf('%s:  Mapping horiz. k_xy-profiles using ray invariant\n',S2DHMS(toc(mainTic)));
[ym,kxym] = mkisopp(qy,qkxy,qIz,sH.*skz((tN/2+1):tN));                                              if clr; clear qy qkxy qIz; end
if isoN; fIdx = dimFind(ym,2,1,'first','num'); lIdx = dimFind(ym,2,1,'last','num');
  newym = downsamp(ym,2,isoN,fIdx,lIdx);                                                            if clr; clear fIdx lIdx; end
  newkxym = dimSwap(linppval(ym,kxym,2,dimSwap(newym,2,3)),2,3);
  ym = dimSwap(newym,2,6); kxym = dimSwap(newkxym,2,6);                                             if clr; clear newym newkxym; end
else; ym = dimSwap(ym,2,6); kxym = dimSwap(kxym,2,6); end
ym = [flip(ym,1);ym]; kxym = [flip(kxym,1);kxym];
dkxym = dimDiff(kxym,6)./dimDiff(ym,6); rdkxydy = linppval(ym,dkxym,6,ry);                          if clr; clear dkxym; end
rkxy = linppval(ym,kxym,6,ry); rkkxy = rkxy.^2; rt = acos(rkxy./k);
rkz = sqrt(kk-rkkxy); rkz(imag(rkz)~=0) = nan;
%% Specify kx(|m,p) to calculate kxy(y|m,p) and ky(y|m,p) Profiles
fprintf('%s:  Calculating transverse k_y-profiles\n',S2DHMS(toc(mainTic)));
kx = skxy.*cos(sp); kkx = kx.^2; cx = w./kx;
rky = sqrt(rkkxy-kkx); rky(imag(rky)~=0) = nan;                                                     if clr; clear rkkxy; end
[ympL,kxympL,ympH,kxympH] = kturns(ym,kxym,6,kx);
ympL(ympL>sy) = sy; ympH(ympH<sy) = sy;
ymp = ympL + (ympH - ympL).*dimExtrude(linspace(0,1,qyN).',7);
ymp(:,:,1) = ympL; ymp(:,:,end) = ympH;                                                             if clr; clear ympL ympH; end
Hmp = linppval(ey,eH,6,ymp);
kxymp = linppval(ym,kxym,6,ymp); kkxymp = kxymp.^2;                                                 if clr; clear ym kxym; end
%kxympL = kxymp(:,:,1); kxympH = kxymp(:,:,end);
kxymp(:,:,1) = kxympL; kxymp(:,:,end) = kxympH;
kzmp = sqrt(kk - kkxymp);
%% Z-Cycle Integrals
fprintf('%s:  Calculating vertical z-cycle integrals\n',S2DHMS(toc(mainTic)));
sLz = abs(sH./skz); sLzL = sz./skz;
Lzmp = abs(Hmp./kzmp); Lzmp2 = Lzmp.^2;                                                             if clr; clear kzmp Lzmp Hmp; end
rLz = abs(rH./rkz); rLzL = rz./rkz;                                                                 if clr; clear skz; end
kxyEff = sqrt( (kkxymp-kkx).*Lzmp2 + kkx );                                                         if clr; clear kkx Lzmp2 kkxymp; end

%% Phase II Calculation (parfor over y-dim for y-cycle calculations) >>>
switch parType
case "vec"; fprintf('%s:  Vectorized y-queries for y-cycle calculations\n',S2DHMS(toc(mainTic))); lpNm="VEC";
case "off"; fprintf('%s:  For-loop over y-queries for y-cycle calculations\n',S2DHMS(toc(mainTic))); lpNm="FOR";
case "gpu"; fprintf('%s:  ERROR - GPU parallelization not yet implemented\n',S2DHMS(toc(mainTic))); lpNm="GPU";
case {"bkg" "prc" "thr"}; lpNm = "PARFOR";
  fprintf('%s:  Parfor-loop over y-queries for y-cycle calculations\n',S2DHMS(toc(mainTic))); end
qLyL = zeros(tN,pN,1,1,1,1,qyN); qLyH = zeros(tN,pN,1,1,1,1,qyN);
qLyzL = zeros(tN,pN,1,1,1,1,qyN); qLyzH = zeros(tN,pN,1,1,1,1,qyN);
Q = parallel.pool.DataQueue; afterEach(Q,@(~)checkprog);
checkprog('reset',1,'N',qyN,'label',lpNm+" progress:");
%% Y-Cycle Integrals
switch parType
case "vec"; % vectorized version
	[qympL,qkxympL,qympH,qkxympH] = slcLpp(dimSwap(ymp,6,7),dimSwap(kxymp,6,7),6,ymp); % for calculating Ly
	[ ~, qkxyEffmpL, ~, qkxyEffmpH] = slcLpp(dimSwap(ymp,6,7),dimSwap(kxyEff,6,7),6,ymp); % for calculating Lyz
  qLyL = intCz_pwLinC(qympL,w./qkxympL,6,cx)./w;
	qLyH = intCz_pwLinC(qympH,w./qkxympH,6,cx)./w;
	qLyzL = intCz_pwLinC(qympL,w./qkxyEffmpL,6,cx)./w;
	qLyzH = intCz_pwLinC(qympH,w./qkxyEffmpH,6,cx)./w;
case "off"; % for loop as if using parfor
  for ii = 1:qyN; iymp = ymp(:,:,ii);
    [qLyL(:,:,1,1,1,1,ii), qLyH(:,:,1,1,1,1,ii),	qLyzL(:,:,1,1,1,1,ii), qLyzH(:,:,1,1,1,1,ii)] ...
			= parforfun1(ymp,kxymp,kxyEff,iymp,w,cx); send(Q,[]); end
case "gpu"; error('GPU parallelization not yet implemented');
case {"bkg" "prc" "thr"} % parfor loop
	parfor (ii = 1:qyN,pfW); iymp = ymp(:,:,ii);
	[qLyL(:,:,1,1,1,1,ii), qLyH(:,:,1,1,1,1,ii),	qLyzL(:,:,1,1,1,1,ii), qLyzH(:,:,1,1,1,1,ii)] ...
			= parforfun1(ymp,kxymp,kxyEff,iymp,w,cx); send(Q,[]); end; %#ok<PFBNS>
end;                                                                                                

%% Phase III Calculations (vectorized) >>>
%% Source Cycle Integrals
fprintf('%s:  Interpolating source transverse cycle integrals\n',S2DHMS(toc(mainTic)))
%sIyH = +sigp.*linppval(ymp,qIyH,7,sy);
sLyL = +sigp.*linppval(ymp,qLyL,7,sy);
sLyH = -sigp.*linppval(ymp,qLyH,7,sy);
sLyzL = +sigp.*linppval(ymp,qLyzL,7,sy);
sLyzH = -sigp.*linppval(ymp,qLyzH,7,sy);                                                            if clr; clear qLyzH; end
switch 1; % TEST SWITCH
  case 0; % original
    Ly = abs(sLyL-sLyH); Ly2 = 2.*Ly; Lyz = abs(sLyzL-sLyzH); 
  case 1; % Test 1
    Ly = intCz_pwLinC(ymp,w./kxymp,7,cx)./w; Ly2 = 2.*Ly;
    Lyz = intCz_pwLinC(ymp,w./kxyEff,7,cx)./w; end;                                                 if clr; clear cx sLyzH kxyEff; end
if  sum(Ly<0,'all')>0; error('negative Ly'); end
if sum(Lyz<0,'all')>0; error('negative Lyz'); end                                                  
%% WKB Modal Amplitude-Squared per Solid Angle:  Psi(y|m,n)
fprintf('%s:  Calculating WKB limiter / modal amplitude\n',S2DHMS(toc(mainTic)))
if ~wkblim; Psi = skxy./(kx.*rky.*rkz.*Ly.*rLz);
else; rkyLim = max(rky,(3.*pi.*rkxy.*abs(rdkxydy)).^(1./3)); 
rkzLim = rkz; Psi = skxy./(kx.*rkyLim.*rkzLim.*Ly.*rLz); end;                                       if clr; clear skxy rky rkz rkxy rdkxydy rkyLim rkzLim; end
%% Transverse Y-Mode Stripping Kernel:  Ry(x|m,p,sgn(ph))
fprintf('%s:  Calculating transverse mode-stripping kernel\n',S2DHMS(toc(mainTic)))
rhoL = kxympL<=kx; rhoH = kxympH<=kx;                                                               if clr; clear kxympL kxympH; end
rho1st = dimIndex(cat(3,rhoL,rhoH),3,(sigp>=0)+1); 
sLy1st = dimIndex(cat(3,abs(sLyL),abs(sLyH)),3,(sigp>=0)+1);                                        if clr; clear sigp sLyH; end
dx2bdry = kx.*(sLy1st + rho1st.*Ly)./~(rhoL&rhoH);                                                  if clr; clear sLy1st rhoL rhoH rho1st; end
Ry = rx < dx2bdry;
%% Rcvr Cycle Integrals
fprintf('%s:  Interpolating receiver transverse cycle integrals\n',S2DHMS(toc(mainTic)))
rLyL = linppval(ymp,qLyL,7,ry); rLyL(rLyL<0) = 0;                                                   if clr; clear qLyH; end
rLyzL = linppval(ymp,qLyzL,7,ry); rLyzL(rLyzL<0) = 0;                                               if clr; clear qLyzL; end
%% Transverse Y-Cycle tracking:  y(x|m,p,sgn(ph))
if ~(botloss || ray);                                                                                       if clr; clear Ly2 qLyL ymp kxymp; end
else; fprintf('%s:  Transverse cycle tracking\n',S2DHMS(toc(mainTic)))
  Dy = kx.*Ly2; cqx = dimExtrude(0:qxN-1,6).*Dy./(qxN-1);
  qdx = cqx(:,:,:,:,:,2)-cqx(:,:,:,:,:,1); dcqp = cqx./Dy;
  cqp = dcqp + sLyL./Ly2; pcqp = mod(cqp,1);                                                        if clr; clear cqp; end
	mcqp = .5-abs(pcqp-.5); scqp = .5-(pcqp>=.5);                                                     if clr; clear pcqp; end
	yofcqx = linppval(abs(qLyL),ymp,7,Ly2.*mcqp);                                                     if clr; clear Ly2 qLyL mcqp; end
  kxyofcqx = linppval(ymp,kxymp,7,yofcqx); thofcqx = acos(kxyofcqx./k);                             if clr; clear ymp kxymp; end
  dcp = rx./Dy;  ncyc = floor(dcp); rcyc = dcp-ncyc;                                                if clr; clear Dy dcp; end
	%y_x = linppval(dcqp,yofcqx,6,rcyc);
	sigy_x = linppval(dcqp,scqp,6,rcyc); 
	%sigy_x(sigy_x>=0) = +1; sigy_x(sigy_x<0) = -1;
	%LyH_x = +sigy_x.*linppval(ymp,qLyH,7,y_x);
	%LyzH_x = +sigy_x.*linppval(ymp,qLyzH,7,y_x);
 end
%% Pseudo Y-Cycle Ray Tracing
rays = struct("x",[],"y",[],"xi",[],"up",[],"th",[],"ph",[],"bnd",[]); 
if ray; xiofcqx = kx./kxyofcqx; upofcqx = scqp.*2.*sqrt(1-xiofcqx.^2);                              if clr; clear scqp; end
  rays.x = permshift(cqx,6,-1); rays.y = permshift(yofcqx,6,-1);
  rays.xi = permshift(xiofcqx,6,-1); rays.up = permshift(upofcqx,6,-1);                             if clr; clear xiofcqx; end
  rays.th = permshift(thofcqx,6,-1); rays.ph = permshift(asin(upofcqx),6,-1);                       if clr; clear upofcqx; end
	rays.bnd = permshift(cqx<dx2bdry,6,-1); % logical for ray within domain bounds
end;                                                                                                if clr; clear scqp pcqp cqx dx2bdry; end
%% Convergence Factor Eigenray Filter Kernel (precalcs):  C(x,y,z|m,p)
% Xys, Xzs, and Xyzs have to be calculated even if not being used, because they could
% potentially be used in the parfor loop so they are broadcast variables
Xys = + pi .* (rx./kx./Ly + sLyL./Ly);                                                 
Xzs = + pi .* (rx.*Lyz./kx./Ly + sLzL./sLz);
%Xzs2 = + (rLyzH_x./Lyz - rLyH_x./Ly);
%Xzs3 = - (sLyzL./Lyz - sLyL./Ly);                                                                   if clr; clear sLyL sLyzL; end
Xyzs = + (Lyz.*sLyL./Ly - sLyzL);                                                                   if clr; clear sLyL sLyzL; end
%Xzs2 = + sLzL./sLz;                                                                                 if clr; clear sLz; end
%Xzs = pi .* (Lyz .* (Xzs1 + Xzs2 + Xzs3) + Xzs4);
%Xzs = pi .* (Lyz .* (Xzs1 + Xzs3) + Xzs4);                                                          if clr; clear Xzs1 Xzs3 Xzs4; end

%% Shared resource for botLoss and dirKern
if botloss || dirKern; qthBot = dimExtrude(linspace(0,pi./2,qtN),7);
qRofthBot = dimExtrude(nLyrRef(f,qthBot,egeo),7); end
%% Cumulative bottom intensity RefCoeff along y(x):  Rz(x|m,p,sgn(ph))
if ~botloss; Rz = ones(1,1,xN);                                                                     
else; fprintf('%s:  Calculating bottom loss along transverse cycles\n',S2DHMS(toc(mainTic)))
  Rofcqx = linppval(qthBot,abs(qRofthBot),7,thofcqx);                                               
  Hofcqx = linppval(dimExtrude(ey,7),dimExtrude(eH,7),7,yofcqx);                                    if clr; clear yofcqx; end
  kzofcqx = sqrt(kk-kxyofcqx.^2);                                                                   
  Lzofcqx = Hofcqx./kzofcqx;                                                                        if clr; clear Hofcqx kzofcqx; end
  Rzofcqx = exp(dimCumTrapz(log(Rofcqx)./Lzofcqx,6,[],qdx)./kx);                                    if clr; clear Rofcqx Lzofcqx; end
  Rz = Rzofcqx(:,:,:,:,:,end).^ncyc .* linppval(dcqp,Rzofcqx,6,rcyc); end;                          if clr; clear ncyc rcyc dcqp Rzofcqx kx qdx yofcqx kxyofcqx thofcqx; end
%% Lloyd-Mirror Directivity Pattern Kernel:  D(y,z,z0|m)
if ~dirKern; DP = 1;                                                                                if clr; clear qthBot qRofthBot rH rt; end
else; fprintf('%s:  Calculating Lloyd mirror directivity kernel\n',S2DHMS(toc(mainTic)))
  Rst = linppval(qthBot,qRofthBot,7,abs(st));
  Rrt = linppval(qthBot,qRofthBot,7,rt);
  [sDPs,sDPb] = lloydMirrorDP(f,sz,ec,sH,abs(st),-1,Rst,mlam);                                      if clr; clear Rst; end
  [rDPs,rDPb] = lloydMirrorDP(f,rz,ec,rH,rt,-1,Rrt,mlam);                                           if clr; clear rH rt Rrt; end
  DP = sDPs.*sDPb.*rDPs.*rDPb; end;                                                                 if clr; clear sDPs sDPb rDPs rDPb; end


%% Phase IV Calculation: calc CF and integrate EF (vectorized or parallel) >>>
% s = [-1,1]; sg = dimSwap(s([mod(0:3,2);(0:3)>=2]+1),2,6); 
% ysg = sg(2,:,:,:,:,:); zsg = sg(1,:,:,:,:,:);
ysg = cat(6,-1,-1,1,1);
zsg = cat(6,-1,1,-1,1);
switch parType
case "vec"; lpNm = "VEC";
  fprintf('%s:  Vectorized convergence factor calculation over receivers and final integrations\n',S2DHMS(toc(mainTic)));
case "off"; lpNm = "FOR";
  fprintf('%s:  For-loop over receivers for convergence factor calculation and final integration\n',S2DHMS(toc(mainTic)));
case "gpu"; fprintf('%s:  ERROR - GPU parallelization not yet implemented\n',S2DHMS(toc(mainTic))); lpNm="GPU";
case {"bkg" "prc" "thr"}; lpNm = "PARFOR";
  fprintf('%s:  Parfor-loop over receivers for convergence factor calculation and final integration\n',S2DHMS(toc(mainTic))); end
Q = parallel.pool.DataQueue; afterEach(Q,@(~)checkprog);
% prepare arrays and indices for parallel sections
ixy = allsubs([xN,yN]); xyN = size(ixy,1); 
PP = nan(xyN,zN); vix = ixy(:,1); viy = ixy(:,2);                                                   if clr; clear ixy; end
checkprog('reset',1,'N',xyN,'label',lpNm+" progress:");
switch parType
case "vec"; 
	if Copt; if CyOn; Xy = Xys - ysg.*pi.*rLyL./Ly; end
		%if CzOn; Xz = Xzsp - zsg.*pi.*rLzHp./rLzp; end
  	%if CzOn; Xz = Xzs + pi.*( ysg.*(rLyzL - Lyz.*rLyL./Ly) - zsg.*rLzL./rLz ); end
  	%if CzOn; Xz = Xzs + pi.*( ysg.*sigy_x.*(rLyzL - Lyz.*rLyL./Ly) - zsg.*rLzL./rLz ); end
  	if CzOn; Xz = Xzs + pi.*( Xyzs - ysg.*sigy_x.*(rLyL./Ly.*Lyz - rLyzL) - zsg.*rLzL./rLz ); end
    % Convergence Factor (Periodic Gaussian)
  	if Copt == 1; CFy = ones([1,1,1,1,1,4]); CFz = CFy;
    	if CyOn; CFy = Fp.*exp(FFp_2pi.*(cos(Xy)-1)); end
    	if CzOn; CFz = Fm.*exp(FFm_2pi.*(cos(Xz)-1)); end
    	CF = sum(CFy.*CFz,6)./16;
  	% Convergence Factor (Cosine Summation)
  	elseif Copt==2; CF = ones([1,1,1,1,1,4]);
			dp = 1:Np; dm = 1:Nm; [ddp,ddm] = ndgrid(dp,dm); dp=dimExtrude(dp,7); 
			dm=dimExtrude(dm,7); ddp=dimExtrude(ddp,7); ddm=dimExtrude(ddm,7);
    	if CyOn; CFy = 2.*sum(cos(dp.*Xy),7); end
			if CzOn; CFz = 2.*sum(cos(dm.*Xz),7); end
			if CyOn&&CzOn; CFyz = 4.*sum(cos(ddp.*Xy+ddm.*Xz),7); end 
    	CF = sum(CF+CFy+CFz+CFyz,6)./16; end
	else; CF = 1./4; end; 
	phInt = Psi .* DP .* Ry .* Rz .* CF;
	% Final Integration over Solid Angle:  TL(x,y,z) = Int Int [] dPh dTh
	if splitph; thInt = (compSimpQuad(phInt(:,negp,:,:,:),2) ...
                   	+ compSimpQuad(phInt(:,posp,:,:,:),2)).*dph;
	else; thInt = compSimpQuad(phInt,2).*dph; end
	if splitth; PP = (compSimpQuad(thInt(negt,:,:,:,:),1) ...
                 	+ compSimpQuad(thInt(post,:,:,:,:),1)).*dth;
	else; PP = compSimpQuad(thInt,1).*dth; end
case "gpu"; error('gpu processing not yet implemented');
case "off" % for loop as if using parfor
  clear parforfun2 ixy
	for ixy = 1:xyN
  	ix = vix(ixy); iy = viy(ixy); % init
  	if bMsk(:,:,:,iy,:); send(Q,[]); continue; end  % short-circuit bail-out
  	% Slice variables by receiver x,y,z
  	Xysp = Xys(:,:,ix); Xzsp = Xzs(:,:,ix); Xyzsp = Xyzs(:,:); rLyLp = rLyL(:,:,:,iy); sigy_xp = sigy_x(:,:,ix);
    rLyzLp = rLyzL(:,:,:,iy,:); rLzLp = rLzL(:,:,:,iy,:); rLzp = rLz(:,:,:,iy);
  	phInt = Psi(:,:,:,iy) .* DP(:,:,:,iy,:) .* Ry(:,:,ix) .* Rz(:,:,ix);
		PP(ixy,:) = parforfun2(Ly,rLyLp,rLzp,rLzLp,Lyz,rLyzLp,Xysp,Xzsp,Xyzsp,phInt,ysg,zsg,sigy_xp, ...
			Np,Nm,Fp,Fm,FFp_2pi,FFm_2pi,Copt,CyOn,CzOn,negp,posp,negt,post,dph,dth,splitph,splitth);
	  send(Q,[]); end % end for
case {"bkg" "prc" "thr"} % parfor version
	clear parforfun2 ixy
	parfor (ixy = 1:xyN,pfW)
		ix = vix(ixy); iy = viy(ixy); % init
  	if bMsk(:,:,:,iy,:); send(Q,[]); continue; end %#ok<PFBNS> % short-circuit bail-out
  	% Slice variables by receiver x,y,z
  	Xysp = Xys(:,:,ix); Xzsp = Xzs(:,:,ix); Xyzsp = Xyzs(:,:,:,iy); sigy_xp = sigy_x(:,:,ix,iy); %#ok<PFBNS>
  	rLyLp = rLyL(:,:,:,iy); rLyzLp = rLyzL(:,:,:,iy,:); %#ok<PFBNS>
		rLzLp = rLzL(:,:,:,iy,:); rLzp = rLz(:,:,:,iy) %#ok<PFBNS>
  	phInt = Psi(:,:,:,iy) .* DP(:,:,:,iy,:) .* Ry(:,:,ix) .* Rz(:,:,ix); %#ok<PFBNS>
		PP(ixy,:) = parforfun2(Ly,rLyLp,rLzp,rLzLp,Lyz,rLyzLp,Xysp,Xzsp,Xyzsp,phInt,ysg,zsg,sigy_xp, ...
			Np,Nm,Fp,Fm,FFp_2pi,FFm_2pi,Copt,CyOn,CzOn,negp,posp,negt,post,dph,dth,splitph,splitth);
	  send(Q,[]); end % end parfor
end % switch partype
if exist('pp','var'); delete(pp); end; % fprintf('\n');                                             
if clr; clear rLz Ly Lyz Xys Xzs rLyH rDzLamH rIyH Psi DP Ry Rz bMsk pi_Ly pi_sDzLam rLyzL rLyL rLzL; end  
%TL = -10 * ( log10( reshape(PP,[xN,yN,zN]) ) + log10( 4*sS2/(4*pi).^2 ) );                          if clr; clear PP; end
pp = 4*sS2/(4*pi).^2 * reshape(PP,[xN,yN,zN]); tl = -10 * log10(pp); 
iitl = isinf(tl); tl(iitl) = 1e9; intl = isnan(tl); tl(intl) = 1e9;
fld.tl = tl; fld.pp = pp;                                         if clr; clear PP pp tl iitl; end
if dbg; keyboard; end %#ok<KEYBOARDFUN>
end % end function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% parforfun1
function [iqLyL,iqLyH,iqLyzL,iqLyzH] = parforfun1(ymp,kxymp,kxyEff,iymp,w,cx)
[qympL,qkxympL,qympH,qkxympH] = slcLpp(dimSwap(ymp,6,7),dimSwap(kxymp,6,7),6,iymp); % for calculating Ly
	[ ~, qkxyEffmpL, ~, qkxyEffmpH] = slcLpp(dimSwap(ymp,6,7),dimSwap(kxyEff,6,7),6,iymp); % for calculating Lyz
  iqLyL = intCz_pwLinC(qympL,w./qkxympL,6,cx)./w;
	iqLyH = intCz_pwLinC(qympH,w./qkxympH,6,cx)./w;
	iqLyzL = intCz_pwLinC(qympL,w./qkxyEffmpL,6,cx)./w;
	iqLyzH = intCz_pwLinC(qympH,w./qkxyEffmpH,6,cx)./w;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% parforfun2
function [PP] = parforfun2(Ly,rLyLp,rLzp,rLzLp,Lyz,rLyzLp,Xysp,Xzsp,Xyzsp,phInt,ysg,zsg,sigy_xp, ...
		Np,Nm,Fp,Fm,FFp_2pi,FFm_2pi,Copt,CyOn,CzOn,negp,posp,negt,post,dph,dth,splitph,splitth)
persistent ipfn2; if isempty(ipfn2); ipfn2 = 0; end
% Convergence Factor Eigenray Filter Kernel (part II):  C(x,y,z|m,p)
if Copt; if CyOn; Xy = Xysp - ysg.*pi.*rLyLp./Ly; end %#ok<>
  %if CzOn; Xz = Xzsp - zsg.*pi.*rLzHp./rLzp; end
  %if CzOn; Xz = Xzsp + pi.*( ysg.*(rLyzLp - Lyz.*rLyLp./Ly) - zsg.*rLzLp./rLzp ); end
  %if CzOn; Xz = Xzsp + pi.*( ysg.*sigy_xp.*(rLyzLp - Lyz.*rLyLp./Ly) - zsg.*rLzLp./rLzp ); end
  if CzOn; Xz = Xzsp + pi.*( Xyzsp - ysg.*sigy_xp.*(rLyLp./Ly.*Lyz - rLyzLp) - zsg.*rLzLp./rLzp ); end
  % Convergence Factor (Periodic Gaussian)
  if Copt == 1; CFy = ones([1,1,1,1,1,4]); CFz = CFy;
    if CyOn; CFy = Fp.*exp(FFp_2pi.*(cos(Xy)-1)); end
    if CzOn; CFz = Fm.*exp(FFm_2pi.*(cos(Xz)-1)); end
    CF = sum(CFy.*CFz,6)./16;
  % Convergence Factor (Cosine Summation)
  elseif Copt==2; CFy = ones([1,1,1,1,1,4])./2; CFz = CFy;
    if CyOn; for dp = 1:Np; addCFy = cos(dp.*Xy);
      CFy = CFy + addCFy; end; CFy = 2.*CFy; end
    if CzOn; for dm = 1:Nm; addCFz = cos(dm.*Xz);
      CFz = CFz + addCFz; end; CFz = 2.*CFz; end
    CF = sum(CFy.*CFz,6)./16; end
else; CF = 1./4; end; phInt = phInt.*CF;
%% Final Integration over Solid Angle:  TL(x,y,z) = Int Int [] dPh dTh
if splitph; thInt = (compSimpQuad(phInt(:,negp,:,:,:),2) ...
                   + compSimpQuad(phInt(:,posp,:,:,:),2)).*dph;
else; thInt = compSimpQuad(phInt,2).*dph; end
if splitth; PP_ = (compSimpQuad(thInt(negt,:,:,:,:),1) ...
                 + compSimpQuad(thInt(post,:,:,:,:),1)).*dth;
else; PP_ = compSimpQuad(thInt,1).*dth; end
PP = PP_(:).'; 
% ipfn2 = ipfn2 + 1; if ipfn2 == 500; keyboard; end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% sub S2DHMS
function [dhmsStr] = S2DHMS(SS)
%[D,H,M,S,dhmsStr] = S2DHMS(SS)
%   Converts total input seconds to days hours minutes and seconds
D = floor(SS./86400); SS = SS-D.*86400;
H = floor(SS./3600); SS = SS-H.*3600;
M = floor(SS./60); S = SS-M.*60;
dhmsStr = sprintf('(%02.0fd:%02.0fh:%02.0fm:%02.0fs)',D,H,M,S); 
end