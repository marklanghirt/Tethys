function [pool,type,nW] = chspool(opt)
%[pool,type,nW]=chspool(['type','omitcore'])
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: routine to check and choose a parallel pool option based on input options


arguments
  opt.type {mustBeMember(opt.type,{'proc','Proc','PROC','process','Process','processes','Processes',...
		'thrd','Thrd','THRD','thread','Thread','threads','Threads','gpu','Gpu','GPU','bkgd',...
		'Bkgd','BKGD','background','Background','sngl','Sngl','SNGL','single','Single','client','Client',...
		'clnt','Clnt','CLNT','none','None','NONE','off','Off','OFF','false','False','FALSE',...
		'vec','Vec','VEC','debug','Debug','DEBUG','dbg','Dbg','DBG','for','For','FOR'})} = 'threads'
  opt.omitcore {mustBeNumericOrLogical} = 0;
end; type = opt.type; omitcore = opt.omitcore; 
try psp = get(parallel.Settings,'Pool'); 
	initVal = psp.AutoCreate; set(psp,'AutoCreate',true); catch; end
ncore = feature('numcores'); cup=canUseParallelPool; cug=canUseGPU;
try set(psp,'AutoCreate',initVal); catch; end
%t=type; if sum(t==["vec","Vec","VEC","debug","Debug","DEBUG","dbg","Dbg","DBG"]); 
t=type; if sum(t==["vec","Vec","VEC"]); 
	isVec = 1; t="off"; else; isVec = 0; end
	t = find([sum(t==["sngl","Sngl","SNGL","single","Single","client","Client", ...
	"clnt","Clnt","CLNT","none","None","NONE","off","Off","OFF","false","False","FALSE",...
	"vec","Vec","VEC","debug","Debug","DEBUG","dbg","Dbg","DBG","for","For","FOR"]);
  sum(t==["bkgd","Bkgd","BKGD","background","Background"]);
  sum(t==["proc","Proc","PROC","process","Process","processes","Processes"]);
  sum(t==["thrd","Thrd","THRD","thread","Thread","threads","Threads"]);
  sum(t==["gpu","Gpu","GPU"])]);
pool=gcp('nocreate'); if ~isempty(pool); delete(pool); end 
nW = 1 + (cup&(prod(t~=[1,5],2)|((t==5)&~cug))).*(ncore-omitcore-1);
if t==5 % GPU
  if cug; gtbl = gpuDeviceTable(); gAvail = gtbl.('DeviceAvailable');
  gpu = gtbl.('Index')(gAvail); gnm = gtbl.('Name')(gAvail); ng = sum(gAvail); 
  nW = ng; pool = parpool('processes',ng); spmd; for ig=1:ng; gpuDevice(gpu(ig)); 
    fprintf('Worker(%.0f) using GPU(%.0f) %s\n',ig,gpu(ig),gnm(ig)); end; end
  else; warning('GPU pool not available, reverting to threaded pool...'); t=4; end
end 
if t==4 % Threads
  if cup; pool = parpool('threads',nW);
    fprintf('Thread pool started with %.0f workers\n',nW);
  else; t=1; end
end
if t==3 % Processes
  if cup; pool = parpool('processes',nW);
    fprintf('Process pool started with %.0f workers\n',nW);
  else; t=1; end
end
if t==2 % Background;
  if ~cup; warning('Parallel pool not available, reverting to single background process'); end
  pool = parpool('background',nW);
  fprintf('Thread pool started in the background with %.0f workers\n',nW);
end % ⇩ Single Client Process ⇩
if t==1; fprintf(['Parallel features are disabled...\n',...
    'remaining in client session on a single process\n']); nW = 0; end
types = ["off","bkg","prc","thr","gpu"]; type = types(t);
if isVec; type = "vec"; end