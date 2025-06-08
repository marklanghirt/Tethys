function [DdA] = dimDiff(A,dim)
%[DdA]=dimDiff(A,dim)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: Computes finite differences of array (like diff)
%  but with ability to specify differencing dimension for nd-arrays

szA = size(A); nd = numel(szA); n = szA(dim);
B = permshift(A,nd,dim-1); szB = size(B);
zz = zeros(1,prod(szB(2:end)));
hw = ones(n-1,1);
mw = [-1;zeros(n-2,1);1];
lw = -ones(n-1,1);
hf = [1;.5.*ones(n-2,1);1];
C = hf.*([hw.*B(2:n,:);zz] ...
        + mw.*B(:,:) ...
       + [zz;lw.*B(1:n-1,:)]);
DdA = permshift(reshape(C,szB),nd,nd-dim+1);
end

