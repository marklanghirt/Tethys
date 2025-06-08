function [integral] = compSimpQuad(integrand,dim) %#codegen
%[integral]=compSimpQuad(integrand,dim)
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: This function uses composite simpson's rule with simple adaptations to
%   handle odd numbers of subintervals. The step size is assumed unity (1)
%   and uniform, SO YOU HAVE TO MULTIPLY YOUR OWN STEP SIZE AGAINST THE 
%   OUTPUT OUTSIDE OF THE FUNCTION CALL. 
% Here are the use cases:
%   N = numel(integrand)  <--- number of nodes
%   n = N - 1             <--- number of subintervals
%   N[0] n[ ]      <--- emptyBehavior
%   N[1] n[0]      <--- midpoint rule
%   N[2] n[1]      <--- trapezoidal rule
%   N[3] n[2]      <--- simpson 1/3 rule
%   N[4] n[3]      <--- simpson 3/8 rule
%   N[odd] n[even] <--- [1,N] comp. simpson 1/3 rule
%   N[even] n[odd] <--- [1,N-3] comp. simpson 1/3 rule + [N-3,N] simpson 3/8 rule

N = size(integrand,dim);
permArr = dim:-1:1;

if N == 0 % empty
  integral = 0;
elseif N == 1 % use midpoint rule
  integral = integrand;
elseif N == 2 % use trapezoidal rule
  integral = sum(integrand,dim,'omitnan')./2;
elseif N == 3 % use simpson 1/3 rule
  weights = [1;4;1]./3;
  if dim>1; weights = permute(weights,permArr); end
  integral = sum(weights.*integrand,dim,'omitnan');
elseif N == 4 % use simpson 3/8 rule
  weights = [1;3;3;1].*3./8;
  if dim>1; weights = permute(weights,permArr); end
  integral = sum(weights.*integrand,dim,'omitnan');
else
  if mod(N,2) % N is odd: use comp. simpson 1/3 rule
    weights = ones(N,1);
    weights(2:2:N-1) = 4;
    weights(3:2:N-2) = 2;
    weights = weights./3;
    if dim>1; weights = permute(weights,permArr); end
    integral = sum(weights.*integrand,dim,'omitnan');
  else % N is even: use comp. simpson 1/3 rule with 3/8 on end
    M = N - 3;
    weights = zeros(N,1);
    weights1 = ones(M,1);
    weights1(2:2:M-1) = 4;
    weights1(3:2:M-2) = 2;
    weights1 = weights1./3;
    weights2 = [1;3;3;1].*3./8;
    weights(1:M) = weights(1:M) + weights1;
    weights(M:N) = weights(M:N) + weights2;
    if dim>1; weights = permute(weights,permArr); end
    integral = sum(weights.*integrand,dim,'omitnan');
  end
end
end

