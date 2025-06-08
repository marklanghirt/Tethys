function [idxOut,valOut] = dimFind(A,dim,n,order,type) %#codegen
% [idxArray]=dimFind(A,dim,n,order)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: vectorized find along specified dimension
%
% INPUTS:
%  A: array to find within
%  dim: dimension to find along
%  n: number of elements to find
%  order: 'first' or 'last'
%   'first' = 'ascending' sort of indices
%   'last' = 'descending' sort of indices

% validate input parameters
if n > size(A,dim); error(['<n> exceeds size(A,',num2str(dim),')']); end
if nargin < 2; dim = 1; end
if nargin < 3; n = 1; end
if nargin < 4; order = 'first'; end
if nargin < 5; switch islogical(A); case 1; type = 'log'; case 0; type = 'num'; end; end
if strcmp(order,'first'); direction = 'ascend'; 
elseif strcmp(order,'last'); direction = 'descend'; 
else error('<order> must be either "first" or "last"'); end
if ~(strcmp(type,'num')||strcmp(type,'log')); error('<type> must be either "log" or "num"'); end
szA = size(A); md = numel(szA); B = permshift(A,md,dim-1); szB = size(B);
% create index stencil array in direction of specified dim find
idxStencil = (1:szA(dim)).';
% copy indices from stencil for elements of A not equal to 0 or not equal to NaN
switch type; case 'num'; idxArray = idxStencil.*~isnan(B);
  case 'log'; idxArray = idxStencil.*~(B==0|isnan(B)); end
% for elements of A that are equal to 0 or NaN, set the indices to NaN or 0
idxArray(idxArray == 0) = nan; 
switch type; case 'num'; B(isnan(idxArray)) = nan;
  case 'log'; B(isnan(idxArray)) = 0; end
% sort the indexing array in the specified order for 'first' or 'last'
[idxArray,sortIdx] = sort(idxArray,1,direction,'MissingPlacement','last');
% take only the first n indices along the specified dim
idxOut = reshape(idxArray(1:n,:),[n,szB(2:end)]);
for ii = 1:prod(szB(2:end)); C(1:szA(dim),ii) = B(sortIdx(:,ii),ii); end
valOut = reshape(C(1:n,:),[n,szB(2:end)]);
idxOut = permshift(idxOut,md,md-(dim-1));
valOut = permshift(valOut,md,md-(dim-1));
end