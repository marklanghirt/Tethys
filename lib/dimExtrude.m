function [B] = dimExtrude(A,dim)
%[B]=dimExtrude(A,dim)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: extrudes any nd-array along column-major order, 
%  and reshapes to orient the output vector along the specified dimension

B = shiftdim(A(:),1-dim);
end