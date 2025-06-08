function [integral] = dimCumTrapz(integrand,dim,nodes,steps)
%[integral]=dimCumTrapz(integrand,dim,nodes,steps)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: Cumulative trapezoidal integration (like cumtrapz) 
%  but with ability to specify integration dimension for nd-arrays

szi = size(integrand); ndi = numel(szi); md = ndi;
if nargin<2; dim = find(szi>1,1,'first'); end
if nargin<3; nodes = []; end
if nargin<4 && isempty(nodes); stps = 1; end; szn = size(nodes);
if nargin<4 && ~isempty(nodes) && szn(dim)==1
if szi(dim)==1; integral = integrand; return
else; error('nodes cannot be singleton in dim'); end; end
if nargin<4 && ~isempty(nodes) && szn(dim)~=1
node = permshift(nodes,md,dim-1); stps = diff(node,1,1); end 
if nargin==4; stps = permshift(steps,md,dim-1); end
itgd = movmean(permshift(integrand,md,dim-1),2,1,'Endpoints','discard');
szi = size(itgd); ndi = numel(szi); szs = size(stps); nds = numel(szs);
szs = [szs,ones([1,md-nds])]; szi = [szi,ones([1,md-ndi])];
srepsize = ceil(szi./szs); srepsize(1) = 1; 
try stps = repmat(stps,srepsize);
catch; error('incompatible size(integrand) and size(nodesOrStep)'); end
szs = size(stps); if ~(szs(1)==1||szs(1)==szi(1))
error('incompatible size(integrand) and size(steps)'); end
integral = permshift([zeros([1,szs(2:end)]);cumsum(itgd.*stps,1,'omitnan')],md,md-dim+1);
end