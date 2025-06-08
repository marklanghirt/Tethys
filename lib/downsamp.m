function [B] = downsamp(A,dim,nidx,fidx,lidx)
%[B]=downsamp(A,dim,nidx,fidx,lidx)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: downsample/resample an nd-array along the specified dimension

sza = size(A); nda = numel(sza); 
if sum(size(nidx)~=[1,1])>1; error('nIdx must be scalar'); end
if nargin<4; fidx = 1; lidx = sza(dim); end
szf = size(fidx); ndf = numel(szf);
szl = size(lidx); ndl = numel(szl);
md = max([nda,ndf,ndl]); sza = [sza,ones([1,md-nda])]; 
szf = [szf,ones([1,md-ndf])]; szl = [szl,ones([1,md-ndl])];
if szf(dim)~=1|szl(dim)~=1; error('fIdx,lIdx must be singleton in dim'); end
frepl = ceil(szl./szf); lrepf = ceil(szf./szl);
try fidx = repmat(fidx,frepl); lidx = repmat(lidx,lrepf); catch; error('incompatible size(fIdx) and size(lIdx)'); end
szf = size(fidx); frepa = ceil(sza./szf); frepa(dim) = 1; arepf = ceil(szf./sza);
try A = repmat(A,arepf); fidx = repmat(fidx,frepa); lidx = repmat(lidx,frepa);
catch; error('incompatible size(A) and size(fIDx,lIdx)'); end
szf = size(fidx); sza = size(A);
didx = lidx-fidx; didx = permshift(didx,md,dim-1); didx = didx(:,:);
idx = round((lidx-fidx)./(nidx-1).*permute((0:(nidx-1)).',dim:-1:1)+fidx);
A = permshift(A,md,dim-1); idx = permshift(idx,md,dim-1);
sza = size(A); A = A(:,:); idx = idx(:,:);
B = nan(nidx,size(A,2));
for ii = 1:size(A,2)
  if (didx(ii)+1)>nidx
    B(:,ii) = A(idx(:,ii),ii);
  else
    B(1:(didx(ii)+1),ii) = A(fidx(ii):lidx(ii),ii);
  end
end
B = permshift(reshape(B,[nidx,sza(2:end)]),md,md-(dim-1));


