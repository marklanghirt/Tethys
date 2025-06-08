function [B] = permshift(A,ndim,shift)
%[B]=permshift(A,ndim,shift)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: a version of dimshift using the permute function instead

B = permute(A,circshift(1:ndim,-shift,2));
end

