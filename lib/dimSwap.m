function [B] = dimSwap(A,dim1,dim2)
%[B]=dimSwap(A,dim1,dim2)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: swap two dimensions of an nd-array

nd = max([numel(size(A)),dim1,dim2]);
permArr1 = 1:nd; permArr2 = permArr1;
permArr2(permArr1==dim1) = dim2;
permArr2(permArr1==dim2) = dim1;
B = permute(A,permArr2);
end

