function [nVal,arrNum,arrStr] = arrUnfold(nVal,absMin,absMax,opt)
%[nVal,arrNum,arrStr]=arrUnfold(nVal,absMin,absMax[,axis=0,scale=1,fmtspc='%+g'])
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: Creates symmetric array of evenly spaced values symmetric about <axis>:
%    i.e. [ axis-absMax ... axis-absMin , axis+absMin ... axis+absMax ] 

arguments; nVal = 180; absMin = 1; absMax = 90; opt.axis = 0; opt.scale = 1; ...
		opt.fmtspc = '%+g'; end; ov = {'nVal','arrNum','arrStr'}; nao = nargout;
if mod(nVal,2)==1; nVal = nVal + 1; end; lims = abs([absMin,absMax] - opt.axis) .* opt.scale;
arrPstvVals = linspace(lims(1),lims(2),nVal/2); arrNgtvVals = -arrPstvVals(end:-1:1);
arrNum = [arrNgtvVals,arrPstvVals] + opt.axis; arrStr = sprintf(opt.fmtspc,arrNum);
clearvars('-except',ov{1:nao}); end % arrUnfold()