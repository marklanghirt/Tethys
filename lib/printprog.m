function printprog(n,N,t,opt)
%printprog(n,N,[t,opt])
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: routine called by checkprog to print progress of loops 

arguments
  n {mustBeNumeric}
  N {mustBeNumeric}
  t {mustBeNumeric} = []
  opt.label char  = ''
  opt.first {mustBeNumericOrLogical} = false;
end; lbl = opt.label;
% intErr = MException('printprog:notint','n and N must be integers');
% if (n-fix(n))~=0 || (N-fix(N))~=0; throw(intErr); end
if n==0 || opt.first; frst = 1; else; frst = 0; end
[~,~,~,~,dhms]=S2DHMS(t); tlen=16;
llen = numel(lbl); if llen~=0; lbl = [lbl,'  ']; llen = llen+2; end
nd = num2str(fix(log10(N))+1); nfmt = ['%',nd,'.0f'];
prog = sprintf([nfmt,'/',nfmt,'=%2.0f%%  '],n,N,n./N.*100); 
plen = numel(prog); prog = strrep(prog,'%','%%');

if frst; bs = '\n'; else; bs = repmat('\b',1,llen+plen+tlen+1); end
if n>=N; nl = '\n'; else; nl = ''; end
fprintf([bs,lbl,prog,dhms,nl])
end

%% sub S2DHMS
function [D,H,M,S,dhmsStr] = S2DHMS(SS)
%[D,H,M,S,dhmsStr] = S2DHMS(SS)
%   Converts total input seconds to days hours minutes and seconds
D = floor(SS./86400); SS = SS-D.*86400;
H = floor(SS./3600); SS = SS-H.*3600;
M = floor(SS./60); S = SS-M.*60;
dhmsStr = sprintf('(%02.0fd:%02.0fh:%02.0fm:%02.0fs)',D,H,M,S); end