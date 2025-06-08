function [B] = dimIndex(A,dim,idx)
%[B]=dimIndex(A,dim,idx)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: indexes an nd-array along the specified dimension,
%  a different index value can be specified for each tuple of 
%  coordinates excluding the dimension to index along


sza = size(A); nda = numel(sza);
szi = size(idx); ndi = numel(szi); 
md = max(ndi,nda); szi = [szi,ones([1,md-ndi])]; 
irep = ceil(sza./szi); irep(dim) = 1;
if szi(dim)~=1; error('idx must be singleton in dim'); end
try idx = repmat(idx,irep); catch; error('incompatible size(A) and size(idx)'); end
szb = sza; szb(dim) = 1;
% idxref = repmat(permute((1:sza(dim)).',dim:-1:1),szb);
idxref = repmat(permshift((1:sza(dim)).',md,1-dim),szb);
idxlog = idxref==idx; B=reshape(A(idxlog),szb);
end

