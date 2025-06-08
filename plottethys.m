function [] = plottethys(name,opt)
% []=plottethys(name,[suffix="",clim=[40,100]])
%
% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT
%
% DESCRIPTION: a routine to load from .mat-file and plot a computed two-dimensional 
%  slice of the transmission loss field computed by the "Tethys" energy flux model


arguments (Input); name; opt.suffix=""; opt.clim=[40,100]; end; close all;
%% directories and paths
mflnm = mfilename("fullpath"); 
tetdir = strrep(mflnm,"plottethys","");
datdir = tetdir + "dat" + filesep;
figdir = tetdir + "fig" + filesep;
picdir = tetdir + "pic" + filesep;
datpth = datdir + name + ".mat"; fldstr = "_fld2d";
if opt.suffix == "dt"; sfxstr = "_"+dtnow; else; sfxstr=opt.suffix; end
%% extract data from matfile
if ~exist(datpth,"file"); fprintf("missing: %s \n",datpth); return; 
  else load(datpth,"src","rcv","fld"); end
fr = src.f; xs = src.x; ys = src.y; zs = src.z; tl = fld.tl;
xr = rcv.x(:); yr = rcv.y(:)'; zr = permute(rcv.z(:),[3,2,1]);
rr = sqrt((xr-xs).^2 + (yr-ys).^2); tr = atan2d(yr-ys,xr-xs);

%% determine extents of receiver array
lens = [rcv.xN,rcv.yN,rcv.zN];
dims = lens>1; dmcd = sum(dims.*[4,2,1]);
xlims = [min(xr,[],"all"),max(xr,[],"all")];
ylims = [min(yr,[],"all"),max(yr,[],"all")];
zlims = [min(zr,[],"all"),max(zr,[],"all")];
rlims = [min(rr,[],"all"),max(rr,[],"all")]; 
tlims = [min(tr,[],"all"),max(tr,[],"all")];

%% choose plotting based on receiver array
switch dmcd
  case 3 % YZ
    ttl="$$\textrm{frequency}="+fr+"\left(\textrm{Hz}\right);\quad\textrm{x-range}="+xr+"\left(\textrm{m}\right)$$";
    xlbl="$$\textrm{y-range:}\quad y\quad\left(\textrm{m}\right)$$";
    ylbl="$$\textrm{depth:}\quad z\quad\left(\textrm{m}\right)$$";
    tl2 = permute(tl,[3,2,1]); hsf = surf(yr(:),zr(:),tl2,tl2,"EdgeAlpha",0);
  case 5 % XZ
    ttl="$$\textrm{frequency}="+fr+"\left(\textrm{Hz}\right);\quad\textrm{y-range}="+yr+"\left(\textrm{m}\right)$$";
    xlbl="$$\textrm{x-range:}\quad x\quad\left(\textrm{m}\right)$$";
    ylbl="$$\textrm{depth:}\quad z\quad\left(\textrm{m}\right)$$";
    tl2 = permute(tl,[3,1,2]); hsf = surf(xr(:),zr(:),tl2,tl2,"EdgeAlpha",0);
  case 6 % XY
    ttl="$$\textrm{frequency}="+fr+"\left(\textrm{Hz}\right);\quad\textrm{depth}="+zr+"\left(\textrm{m}\right)$$";
    xlbl="$$\textrm{x-range:}\quad x\quad\left(\textrm{m}\right)$$";
    ylbl="$$\textrm{y-range:}\quad y\quad\left(\textrm{m}\right)$$";
    tl2 = permute(tl,[2,1,3]); hsf = surf(xr,yr,tl2,tl2,"EdgeAlpha",0);
  case 7 % XYZ
    slcvar = input("choose a dimension to slice (r,t,x,y,z):  ","s");
    switch slcvar
      case "r"; slclim = rlims;
      case "t"; slclim = tlims;
      case "x"; slclim = xlims;
      case "y"; slclim = ylims;
      case "z"; slclim = zlims; end
    slcval = input("enter a value to slice at ("+slclim(1)+"..."+slclim(2)+"):  ");
    fldstr = fldstr+"_"+slcvar+round(slcval);
    ttl="$$\textrm{frequency}="+fr+"\left(\textrm{Hz}\right);\quad\textrm{"+slcvar+"-range}="+slcval+"\left(\textrm{m}\right)$$";
    % create gridded interpolant and slice <slcval> at <slcvar>
    gi = griddedInterpolant({xr(:),yr(:),zr(:)},tl,"linear","none");
    ntrplen = numel(xr)+numel(yr)+numel(zr); 
    r2 = linspace(rlims(1),rlims(2),ntrplen).';
    t2 = linspace(tlims(1),tlims(2),ntrplen);
    switch slcvar
      case "r"; 
        xlbl="$$\textrm{bearing:}\quad\theta\quad\left(\textrm{deg}\right)$$";
        ylbl="$$\textrm{depth:}\quad z\quad\left(\textrm{m}\right)$$";
        [R,T,Z]=ndgrid(slcval,t2(:),zr); X=R.*cosd(T)+xs; Y=R.*sind(T)+ys;
        tl2=permute(gi(X,Y,Z),[3,2,1]); hsf=surf(t2(:),zr(:),tl2,tl2,EdgeAlpha=0);
      case "t";
        xlbl="$$\textrm{range:}\quad r\quad\left(\textrm{m}\right)$$";
        ylbl="$$\textrm{depth:}\quad z\quad\left(\textrm{m}\right)$$";
        [R,T,Z]=ndgrid(r2(:),slcval,zr); X=R.*cosd(T)+xs; Y=R.*sind(T)+ys;
        tl2=permute(gi(X,Y,Z),[3,1,2]); hsf=surf(r2(:),zr(:),tl2,tl2,EdgeAlpha=0);
      case "z";
        xlbl="$$\textrm{x-range:}\quad x\quad\left(\textrm{m}\right)$$";
        ylbl="$$\textrm{y-range:}\quad y\quad\left(\textrm{m}\right)$$";
        tl2=permute(gi({xr,yr,slcval}),[2,1,3]); hsf=surf(xr(:),yr(:),tl2,tl2,EdgeAlpha=0);
      case "x";
        xlbl="$$\textrm{y-range:}\quad y\quad\left(\textrm{m}\right)$$";
        ylbl="$$\textrm{depth:}\quad z\quad\left(\textrm{m}\right)$$";
        tl2=permute(gi({slcval,yr,zr}),[3,2,1]); hsf=surf(yr(:),zr(:),tl2,tl2,EdgeAlpha=0);
      case "y";
        xlbl="$$\textrm{x-range:}\quad x\quad\left(\textrm{m}\right)$$";
        ylbl="$$\textrm{depth:}\quad z\quad\left(\textrm{m}\right)$$";
        tl2=permute(gi({xr,slcval,zr}),[3,1,2]); hsf=surf(xr(:),zr(:),tl2,tl2,EdgeAlpha=0);
    end%switch
end%switch
%% Common plotting configuration
hfg = gcf; axis tight;
colormap(flipud(turbo(256))); hcb=colorbar(Direction='reverse'); hcb.Label.String="TL"; 
if ~isstring(opt.clim); clim(opt.clim); end; view(0,90); hax = gca; hax.YDir='reverse'; hax.TickDir="both";
xlabel(xlbl,interpreter='latex'); ylabel(ylbl,interpreter='latex'); title(ttl,interpreter='latex')
drawnow();

%% Save plots
filnam=name+fldstr+sfxstr; figpth=figdir+filnam+".fig"; picpth=picdir+filnam+".png";
savefig(hfg,figpth); fprintf("Figure saved to:  %s\n",figpth);
exportgraphics(hfg,picpth); fprintf("Image saved to:  %s\n",picpth);
end%function