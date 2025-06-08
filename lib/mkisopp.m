function [x,y] = mkisopp(indVar,depVar,funMat,lvls)
%[x,y]=mkisopp(indVar,depVar,funMat,lvls)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: uses matlab's contourc routine to find iso-funMat-evaluated
%  contours, and constructs and returns piecewise polynomial profiles
%
% INPUTS:
%  indVar: vector of independent variable
%  depVar: vector of dependent variable
%  funMat: matrix of function eval of indVar and depVar
%  lvls: vector of levels of desired isocontours


indVar=indVar(:); depVar=depVar(:).'; lvls=lvls(:);
cm = contourc(indVar,depVar,funMat.',lvls);
szcm = size(cm); ncol = szcm(2);
iC = 1; head(iC) = 1; C(iC).n = cm(2,head);
while head(iC) < ncol
  nx(iC) = cm(2,head(iC));
  head(iC+1) = head(iC) + nx(iC) + 1;
  iC = iC + 1;
end; head = head(1:end-1)+1; nC = iC - 1;
mnx = max(nx); tail = head + nx - 1;
x = nan(nC,mnx);
y = nan(nC,mnx);
flipSwitch = cm(1,3)-cm(1,2)<0; dir = 1;
if flipSwitch; temp = head; head = tail; tail = temp; dir = -1; end
for iC = 1:nC
  x(iC,1:nx(iC)) = cm(1,head(iC):dir:tail(iC));
  y(iC,1:nx(iC)) = cm(2,head(iC):dir:tail(iC));
end
end


