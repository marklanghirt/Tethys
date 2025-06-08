function [ppLb,ppLv,ppHb,ppHv] = slcLpp(ppb,ppv,ppdim,slc)
%[ppL,ppH]=slcLpp(ppb,ppv,ppdim,slc)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: splits nd-array of profiles at slice points
%  each profile can have a different slice point to split at
%
% INPUTS:
%  ppb: nd-array of piecewise linear polynomial breaks (nodes)
%  ppv: nd-array of piecewise linear polynomial values
%  ppdim: dimension of nd-arrays along which profiles extend
%  slc: ND-array of slice points slc(i) to slice pp(i) in two
% OUTPUTS:
%  ppLb: lower pp breaks before slice points
%  ppLv; lower pp values before slice points
%  ppHb: upper pp breaks after slice points
%  ppHv: upper pp values after slice points

b=ppb; v=ppv; s=slc; szb=size(b); szv=size(v); szs=size(s);
d=ppdim; ndb=numel(szb); nds=numel(szs); md=max(ndb,nds);
if sum(szb~=szv); error('ppb and ppv must have same size'); end
szb=[szb,ones([1,md-ndb])]; szs=[szs,ones([1,md-nds])];
if szs(d) ~= 1; error('slc must be singleton in ppdim'); end
brepsize = ceil(szs./szb); srepsize=ceil(szb./szs); srepsize(d)=1;
try b=repmat(b,brepsize); v=repmat(v,brepsize); s=repmat(s,srepsize);
catch; error('incompatible size(ppb) and size(slc)'); end
b=permshift(b,md,d-1); v=permshift(v,md,d-1); s=permshift(s,md,d-1);
vofsl = linppval(b,v,1,s); vofsl(isnan(vofsl)) = 0; vofsl=vofsl(:,:);
[ib1st,b1st] = dimFind(b,1,1,'first','num'); ib1st=ib1st(:,:); b1st=b1st(:,:);
[ibEnd,bEnd] = dimFind(b,1,1,'last','num'); ibEnd=ibEnd(:,:); bEnd=bEnd(:,:);
nb=size(b,1); ni=nb-1; szs=size(s); s=s(:,:); b=b(:,:); v=v(:,:);
ns=size(s,2); ibArr = (1:nb).'; szb=[nb,szs(2:end)]; % icArr = (1:nc).'; 
b1 = b(1:ni,:,:); b2 = b(2:end,:,:);
v1st = linppval(b,v,1,b1st); 
vEnd = linppval(b,v,1,bEnd); 
sgeb1 = s >= b1; sltb2 = s < b2; seqbEnd = s == bEnd; 
seqb = sum(s == b,1,'omitnan')>0; 
sgeb1altb2oeqbE = (sgeb1&sltb2) | ((2:nb).'==ibEnd&seqbEnd);
intIdx = dimFind(sgeb1altb2oeqbE,1,1,'first','log'); 
below = s <= b1st; above = s >= bEnd; ext = below|above;
intIdxNoNan = intIdx; intIdxNoNan(isnan(intIdx)) = 0;
ppL_lastbIdx = (~ext).*(intIdxNoNan + ~seqb) + below.*ib1st + above.*ibEnd; 
ppH_firstbIdx = (~ext).*intIdxNoNan + below.*ib1st + above.*ibEnd; 
ppL_bNanIdx = ppL_lastbIdx + 1; ppH_bNanIdx = ppH_firstbIdx-1;
%ppL_lastcIdx = (~ext).*(intIdxNoNan - seqb) + below.*ib1st + above.*(ibEnd-1); 
%ppH_firstcIdx = (~ext).*intIdxNoNan + below.*ib1st + above.*(ibEnd-1); 
%ppL_cNanIdx = ppL_lastcIdx + 1 - below; ppH_cNanIdx = ppH_firstcIdx - 1 + above;
ppLb = b; ppHb = b; ppLv = v; ppHv = v; 
%ppL_co = co; ppH_co = co;
ppL_lastbLog = ibArr==ppL_lastbIdx; ppH_firstbLog = ibArr==ppH_firstbIdx;
%ppL_lastcLog = icArr==ppL_lastcIdx; ppH_firstcLog = icArr==ppH_firstcIdx;
ppL_lastb = below.*b1st + ~ext.*s + above.*bEnd;
ppH_firstb = below.*b1st + ~ext.*s + above.*bEnd;
ppL_lastv = below.*v1st + ~ext.*vofsl + above.*vEnd;
ppH_firstv = below.*v1st + ~ext.*vofsl + above.*vEnd;
ppLb(ppL_lastbLog) = ppL_lastb; ppHb(ppH_firstbLog) = ppH_firstb;
ppLv(ppL_lastbLog) = ppL_lastv; ppHv(ppH_firstbLog) = ppH_firstv;
ppL_bNanLog = ibArr >= ppL_bNanIdx; ppH_bNanLog = ibArr <= ppH_bNanIdx;
ppLb(ppL_bNanLog) = nan; ppLv(ppL_bNanLog) = nan;
ppHb(ppH_bNanLog) = nan; ppHv(ppH_bNanLog) = nan; 
ppLb=permshift(reshape(ppLb,szb),md,md-(d-1)); 
ppLv=permshift(reshape(ppLv,szb),md,md-(d-1));
ppHb=permshift(reshape(ppHb,szb),md,md-(d-1)); 
ppHv=permshift(reshape(ppHv,szb),md,md-(d-1));
true;
% DEBUG
%if 1

%sqz_ib1st = squeeze(ib1st); sqz_b1st = squeeze(b1st);
%sqz_bEnd = squeeze(bEnd); sqz_ibEnd = squeeze(ibEnd);
%sqz_seqbEnd = squeeze(seqbEnd); sqz_seqb = squeeze(seqb);
%sqz_intIdx = squeeze(intIdx); sqz_below = squeeze(below);
%sqz_above = squeeze(above); sqz_ppL_lastbIdx = squeeze(ppL_lastbIdx);
%sqz_ppH_firstbIdx = squeeze(ppH_firstbIdx); sqz_ppL_bNanIdx = squeeze(ppL_bNanIdx);
%sqz_ppH_bNanIdx = squeeze(ppH_bNanIdx); sqz_ppL_lastcIdx = squeeze(ppL_lastcIdx);
%sqz_ppH_firstcIdx = squeeze(ppH_firstcIdx); sqz_ppL_cNanIdx = squeeze(ppL_cNanIdx);
%sqz_ppH_cNanIdx = squeeze(ppH_cNanIdx); sqz_ppL_lastb = squeeze(ppL_lastb);
%sqz_ppH_firstb = squeeze(ppH_firstb); sqz_ppL_lastc = squeeze(ppL_lastc);
%sqz_ppH_firstc = squeeze(ppH_firstc);
%vars = whos; vars = {vars.name}; varsIdx = contains(vars,'sqz_'); 
%vars = vars(varsIdx); cellfun(@(var) openvar(var),vars);
%sqz_b = [b(:,1,1);b(:,1,8)]; sqz_sl = squeeze(sl); 
%openvar('sqz_b'); openvar('sqz_sl'); idxArr2 = repmat(ibArr,[2,1]); openvar('idxArr2')

% close all;
% hax(1) = subplot(1,2,1); hold on;
% plot(b(:,1,1),c(:,1,1),'k-','LineWidth',3)
% for ii = 1:7
%   plot(ppL_b(:,1,ii),ppL_c(:,1,ii)+3*ii/10,'b-','LineWidth',3)
%   plot(ppH_b(:,1,ii),ppH_c(:,1,ii)+3*ii/10,'r-','LineWidth',3)
%   scatter(sl(ii),ppL_c(ppL_lastbIdx(ii),1,ii)+3*ii/10,'ko','MarkerFaceColor','g')
% end
% hax(2) = subplot(1,2,2); hold on;
% plot(b(:,1,8),c(:,1,8),'k-','LineWidth',3)
% scatter([b(2,1,8)-.25,b(4,1,8)+.25],c([2,4],1,8),'r*','LineWidth',2)
% for ii = 8:14
%   plot(ppL_b(:,1,ii),ppL_c(:,1,ii)+3*(ii-7)/10,'b-','LineWidth',3)
%   plot(ppH_b(:,1,ii),ppH_c(:,1,ii)+3*(ii-7)/10,'r-','LineWidth',3)
%   scatter(sl(ii),ppL_c(ppL_lastbIdx(ii),1,ii)+3*(ii-7)/10,'ko','MarkerFaceColor','g')
% end
% hax(1).XLim = [0.5,6.5]; hax(2).XLim = [0.5,6.5];
% xlabel(hax(1),'x'); xlabel(hax(2),'x'); ylabel(hax(1),'y(x)'); 
% title(hax(1),'Testing slicing, pp without NaN');
% title(hax(2),'Testing slicing, NaN on pp ends');

% close all; figure; 
% for ii = 1:prod(szb(2:end))
%   hax = subplot(szb(2),szb(4),ii); hold on;
%   plot([slc(ii),slc(ii)],[0,10],'k:','LineWidth',2);
%   plot(b(:,ii),v(:,ii),'k-','LineWidth',5);
%   plot(ppLb(:,ii),ppLv(:,ii),'g--','LineWidth',2);
%   plot(ppHb(:,ii),ppHv(:,ii),'r--','LineWidth',2);
%   hax.XLim = [-.5,10.5]; hax.YLim = [0,10];
% end

%end