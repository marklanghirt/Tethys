function checkprog(opt)
%[]=checkprog([N=[],label='',deltan=[],deltapercent=[],deltatime='',reset=false])
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: routine called repeatedly to update progress of loops

arguments
  opt.N {mustBeNumeric} = []
  opt.label {mustBeTextScalar} = ''
  opt.deltan {mustBeNumeric} = []
  opt.deltapercent {mustBeNumeric} = []
  opt.deltatime {mustBeTextScalar} = ''
  opt.reset {mustBeNumericOrLogical} = false;
end; 
persistent prog_n prog_N prog_dn prog_dp prog_dt prog_prvtoc prog_tic prog_lbl

% intErr = MException('printprog:notint','n and N must be integers');
% if (n-fix(n))~=0 || (N-fix(N))~=0; throw(intErr); end
if ~isempty(opt.N); prog_N=opt.N; end; prt = 0;
if opt.reset; if isempty(opt.N); error("must specify 'N' on reset"); end
  prog_lbl = opt.label; prog_n=0; prog_dn=[]; prog_dp=[]; prog_dt=[]; prog_prvtoc=[]; prog_tic=tic; 
  prt = 1; d = {opt.deltan,opt.deltapercent,opt.deltatime};
  if ~isempty(d{1}); prog_dn = d{1}; end; if ~isempty(d{2}); prog_dp = d{2}; end
  if ~isempty(d{3}); du = sscanf(d{3},'%f',1); u = sscanf(d{3},'%c',1); 
  mset=[86400,3600,60,1]; m=mset(strcmpi(u,{'d','h','m','s'})); prog_dt=du.*m; end
  if isempty(prog_dn) && isempty(prog_dp) && isempty(prog_dt); prog_dp = 1; end
else; prog_n = prog_n + 1; end

dpp = prog_N.*prog_dp./100; t=toc(prog_tic); lt = prog_prvtoc;
if ~isempty(prog_dn); if fix(prog_n./prog_dn)-fix((prog_n-1)./prog_dn)>0; prt=1; end; end
if ~isempty(prog_dp); if fix(prog_n./dpp)-fix((prog_n-1)./dpp)>0; prt=1; end; end
if ~isempty(prog_dt); if fix(t./prog_dt)-fix(lt./prog_dt)>0; prt=1; end; end
if prt; printprog(prog_n,prog_N,t,'label',prog_lbl,'first',~prog_n); end; prog_prvtoc = t;