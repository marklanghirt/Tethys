function [val] = linppval(ppb,ppv,ppdim,qry)
%[val]=linppval(ppb,ppv,ppdim,qry)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: Evaluates multiple profiles at multiple points. Truncates a profile to 
%  exclude NaNs (to correctly evaluate at the upper limit point) and sets 
%  the evaluation of exterior points to NaN (no interp)

b=ppb; v=ppv; q=qry; szb=size(b); szv=size(v); szq=size(q);
d=ppdim; ndb=numel(szb); ndq=numel(szq); md=max(ndb,ndq);
if sum(szb~=szv); error('ppb and ppv must have same size'); end
szb=[szb,ones([1,md-ndb])]; szq=[szq,ones([1,md-ndq])];
if szq(d) ~= 1; error('qry must be singleton in ppdim'); end
brepsize = ceil(szq./szb); qrepsize=ceil(szb./szq); qrepsize(d)=1;
try b=repmat(b,brepsize); v=repmat(v,brepsize); q=repmat(q,qrepsize);
catch; error('incompatible size(ppb) and size(qry)'); end
b=permshift(b,md,d-1); v=permshift(v,md,d-1); q=permshift(q,md,d-1);
nb=size(b,1); ni=nb-1; szq=size(q); q=q(:,:); b=b(:,:); v=v(:,:); 
nq=size(q,2); b1=b(1:ni,:); v1=v(1:ni,:); dv=(v(2:nb,:)-v1)./(b(2:nb,:)-b1);
allNanLog = sum(~isnan(b),1)==0; intArr = (1:ni).';
[bmin,~]=min(b,[],1,'omitnan'); [bmax,ibmax]=max(b,[],1,'omitnan');
extLog = (q<bmin)|(q>bmax); idx=ones(szq); 
for ii=1:ni; lastIntLog=ii==ibmax-1; intLog=q(1,:)>=b(ii,:) ...
& (q(1,:)<b(ii+1,:)|(q(1,:)==b(ii+1,:)&lastIntLog(1,:)));
idx(intLog)=ii; end; idxLog=intArr==idx; ib=reshape(b1(idxLog),[1,nq]);
iv=reshape(v1(idxLog),[1,nq]); idv=reshape(dv(idxLog),[1,nq]); 
val=(q-ib).*idv+iv; val(allNanLog|extLog)=nan;
val=permshift(reshape(val,szq),md,md-(d-1));