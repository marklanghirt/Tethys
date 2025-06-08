function [D,H,M,S,dhmsStr] = S2DHMS(SS)
%[D,H,M,S,dhmsStr] = S2DHMS(SS)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: Converts total input seconds to days hours minutes and seconds 


D = floor(SS(1)./86400);
SS(2) = SS(1)-D.*86400;
H = floor(SS(2)./3600);
SS(3) = SS(2)-H.*3600;
M = floor(SS(3)./60);
S = SS(3)-M.*60;
dhmsStr = sprintf('%0i:%0i:%0i:%0.1f',D,H,M,S);
end

