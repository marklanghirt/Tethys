function [subs] = allsubs(arrayOrSize)
%[subs]=allsubs(arrayOrSize)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: returns a full list of indexing tuples using column-major order 
%  IF arrayOrSize is a vector, then sz = arrayOrSize
%  IF arrayOrSize is not a vector, then sz = size(arrayOrSize)

if isvector(arrayOrSize); sz=arrayOrSize; else; sz=size(arrayOrSize); end
nd = numel(sz); idx = (1:prod(sz)).'; s=cell(1,nd);
[s{:}] = ind2sub(sz,idx); subs = cell2mat(s);