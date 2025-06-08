% COPYRIGHT: ©2025, Mark Adam Langhirt
% LICENSE: This research code is distributed under the MIT license,
%  a copy of which is included in the code's repository 
%  or alternatively can be viewed at https://opensource.org/license/MIT

clear variables; close all; tau = 2.*pi; hpi = pi./2; veps=1e-3;
[matDir,~,~] = fileparts(which('isoWedge_quadRay_cycleTest2.m'));
ssty = '\scriptstyle '; dsty = '\displaystyle '; degstr=repmat('^{\circ}',4,1); 
str1 = [ssty,'1^{\mathrm{st}}\,',dsty]; fstxofy = [str1,'x(\widetilde{y})'];
supopnstr = repmat('{}^{',4,1); supclsstr = repmat('}',4,1);
% procedural flags + params
runflag = [1,1,1,1,1]; subrunflag = [1,1,1,1]; svfig = 0; svpng = 0;
viewList = [-245,38];
scrn2pos = [1.0000, 0.3667, 0.6667, 0.5829];

%% Env|Src|Rcv params
x0 = 0; y0 = 4000; z0 = 100; th0 = 12; ph0 = 55; 
thb = 2.8624; bet = tand(thb); H = @(y) bet.*y; 
xmax = 15000; ymax = 10000; zmax = bet.*ymax; 
N = 10000; n = N+1; dx = xmax./N; x = [0:dx:xmax,nan]; 
c = 1500; f = 25; k = tau.*f./c; 
%% duped strings
mthprestr = repmat('|\theta_0|=',4,1); sthprestr=repmat('\theta_0=',4,1);
mphprestr = repmat('|\phi_0|=',4,1); sphprestr=repmat('\phi_0=',4,1); 
mthstr=repmat(sprintf('%.0f',th0),4,1); mphstr=repmat(sprintf('%.0f',ph0),4,1); 
mpheqstr = mat2cell([mphprestr,mphstr,degstr],[1,1,1,1]);
mtheqstr = mat2cell([mthprestr,mthstr,degstr],[1,1,1,1]);
%% Sgn and strings
sg = [-1,+1,-1,+1];
sth0 = [-1;+1;-1;+1]; th0 = sth0.*th0; 
sthstr=[ssty,'-',dsty;ssty,'+',dsty;ssty,'-',dsty;ssty,'+',dsty]; 
sph0 = [-1;-1;+1;+1]; ph0 = sph0.*ph0;
sphstr=[ssty,'-',dsty;ssty,'-',dsty;ssty,'+',dsty;ssty,'+',dsty];
sthsphstr=compose('\\left(%s\\,|\\,%s\\right)',sthstr,sphstr);
spheqstr=mat2cell([sphprestr,supopnstr,sphstr,supclsstr,mphstr,degstr],[1,1,1,1]);
stheqstr=mat2cell([sthprestr,supopnstr,sthstr,supclsstr,mthstr,degstr],[1,1,1,1]);
%% Src k-vec-comps
H0 = H(y0);
kxy0 = k.*cosd(th0);
kz0 = k.*sind(th0); mkz0 = abs(kz0);
kx = kxy0.*cosd(ph0);
ky0 = kxy0.*sind(ph0); mky0 = abs(ky0);
kv0 = [kx,ky0,kz0]; alf0 = kv0./k;
kyz = sqrt(k.^2-kx.^2);
%% Gen k-vec-comps
mkz = @(y) y0./y.*mkz0;
ylwr = y0.*mkz0./kyz;
mkzmax = mkz(ylwr);
kxy = @(y) sqrt(k.^2-(y0./y.*mkz0).^2);
mky = @(y) sqrt(kyz.^2-(y0./y.*mkz0).^2);
mkymax = mky(ymax);
%% Z-cyc-ints
mLz = @(y,z) z.*y./y0./mkz0;
inv_mLz = @(Lz,y) Lz./y.*y0.*mkz0;
mtLz = @(y) bet.*y.^2./y0./mkz0;
m0Lz = @(z) z./mkz0;
mtLz0 = m0Lz(H(y0));
Lz0 = sth0.*m0Lz(z0);
%% Y-cyc-ints
mLy = @(y) y.*mky(y)./kyz.^2;
inv_mLy = @(Ly) sqrt(kyz.^2.*Ly.^2 + y0.^2.*mkz0.^2./kyz.^2);
mtLy = mLy(ymax);
mLy0 = mLy(y0);
Ly0 = sph0.*mLy0;
%% YZ-cyc-ints
mLyz = @(y) atan(mky(y)./mkz(y))./bet;
Lyz = @(y,sy) sy.*mLyz(y);
mLyz0 = mLyz(y0);
Lyz0 = Lyz(y0,sph0);
mtLyz = mLyz(ymax);
t0Lyz = Lyz(ymax,sph0);
%% 1st turning points (refraction/reflection)
y_tp1 = (sph0<=0).*(y0.*mkz0./kyz) + (sph0>0).*ymax;
x_tp1 = (sph0<=0).*kx.*mLy0 + (sph0>0).*kx.*(mtLy-mLy0);
%% Y-cycles
dCPy = pi.*x./sg(2)./kx./mtLy;
CPy0 = - pi.*sg(1)./sg(2).*(Ly0./mtLy);
CPy = dCPy + CPy0;
mpCPy = (pi-abs(mod(CPy,tau)-pi));
spCPy = sign(pi-mod(CPy,tau)); spCPy(spCPy==0)=1;
mLy_x = mtLy./pi.*mpCPy;
y_x = inv_mLy(mLy_x);
mky_x = mky(y_x);
ky_x = spCPy.*mky_x;
%% Z-cycles
mtLz_x = mtLz(y_x);

% zCycMethod1 = intgr8 along y(x)
for ii = 1:size(x,2); if ii == 1; cumLz(:,ii) = [0;0;0;0];
else; cumLz(:,ii) = compSimpQuad(1./mtLz_x(:,1:ii),2);
end; end
dCPz1 = pi./kx./sg(4).*dx.*cumLz;

% zCycMethod2 = pre-intgr8d y-cyc offsets z-cyc
yzCorr = @(y,sgy) sgy.*mLyz(y)./mtLyz - sgy.*mLy(y)./mtLy;
dCPz21 = pi.*x.*mtLyz./sg(4)./kx./mtLy;
dCPz22 = pi.*sg(1)./sg(4).*mtLyz.*yzCorr(y0,sph0);
dCPz23 = pi.*sg(2)./sg(4).*mtLyz.*yzCorr(y_x,spCPy);
dCPz2 = dCPz21 + dCPz22 + dCPz23;
CPz0 = - pi .* sg(3)./sg(4) .* Lz0 ./ mtLz0;
dCPz = cat(3,dCPz1,dCPz2);
CPz = dCPz + repmat(CPz0,[1,1,2]);
mpCPz = (pi-abs(mod(CPz,tau)-pi));
spCPz = sign(pi-mod(CPz,tau)); spCPz(spCPz==0)=1;
mLz_x = mtLz_x./pi.*mpCPz;
z_x = inv_mLz(mLz_x,y_x);
mkz_x = mkz(y_x);
kz_x = spCPz.*mkz_x;
%% metrics
% dz = z_x(:,2:end) - z_x(:,1:end-1);
% dx = repmat(x(1,2:end) - x(1,1:end-1),4,1);
% kz = (mkz_x(:,1:end-1) + mkz_x(:,2:end)) / 2;
% dzdx = dz./dx;
% kzkx = kz./kx;
% hfg = figure; 
% subplot(5,2,1); plot(x,cumtrapz(x,1./mtLz_x(1,:),2)); ylabel('$$\int_0^x\frac{\mathrm{d}x}{|\tilde{L}_{z}(x)|}$$')
% subplot(5,2,2); plot(x,CPz(1,:)); ylabel('$$\Xi_z(x)$$')
% subplot(5,2,3); plot(x,mpCPz(1,:)); ylabel('$$\Xi^{|\angle|}_z(x)$$')
% subplot(5,2,4); plot(x,mtLz_x(1,:)); ylabel('$$|\tilde{L}_z(x)|$$')
% subplot(5,2,5); plot(x,mLz_x(1,:)); ylabel('$$|L_z(x)|$$')
% subplot(5,2,6); plot(x,z_x(1,:)); ylabel('$$z(x)$$')
% subplot(5,2,7); plot(x(1:end-1),dz(1,:)); ylabel('$$\mathrm{d}z(x)$$')
% subplot(5,2,8); plot(x(1:end-1),kz(1,:)); ylabel('$$k_z(x)$$')
% subplot(5,2,9:10); plot(x(1:end-1),abs(dzdx(1,:)));
% hold on; plot(x(1:end-1),abs(kzkx(1,:)));
% ylabel('$$\left|\frac{\mathrm{d}z}{\mathrm{d}x}\right| \quad\mathrm{or}\quad \left|\frac{k_z}{k_x}\right|$$')
% legend({'$$\left|\frac{\mathrm{d}z}{\mathrm{d}x}\right|$$','$$\left|\frac{k_z}{k_x}\right|$$'})



%% cell packaging arrays
X = repmat({x},[4,1,2]); 
Y = repmat(mat2cell(y_x,[1,1,1,1]),1,1,2);
Z = mat2cell(z_x,[1,1,1,1],size(z_x,2),[1,1]);
DZDX = mat2cell((z_x(:,2:end,:)-z_x(:,1:end-1,:))./(x(:,2:end,:)-x(:,1:end-1,:)),[1,1,1,1],size(z_x,2)-1,[1,1]);
X2 = repmat({x(1:end-1)},[4,1,2]);
C = repmat({[1:n,nan]},[4,1,2]);
KX = repmat({kx.*[ones([1,n]),nan]},[4,1,2]);
KY = repmat(mat2cell(ky_x,[1,1,1,1]),[1,1,2]);
KZ = mat2cell(kz_x,[1,1,1,1],size(kz_x,2),[1,1]);
KZKX = mat2cell(kz_x./kx,[1,1,1,1],size(kz_x,2),[1,1]);
MKY = repmat(mat2cell(mky_x,[1,1,1,1]),[1,1,2]);
MKZ = repmat(mat2cell(mkz_x,[1,1,1,1]),[1,1,2]);
PHY = repmat(mat2cell(CPy,[1,1,1,1]),[1,1,2]);
PHZ = mat2cell(CPz,[1,1,1,1],size(CPz,2),[1,1]);

%% 2D plots (fig per ray)
% METHOD 1
if runflag(1); for ii=1:4; if subrunflag(ii); hfg(ii) = figure;
% yCycPhs(x)
hax(ii,1) = subplot(4,2,1); hold on;
title(['$\mathbf{y}\mathrm{-cycle}\qquad',spheqstr{ii},'$']);
hln(ii,1,1) = plot(X{ii,1,1},cos(PHY{ii,1,1}),'.-','Color',matcolors(1),'LineWidth',.5);
hln(ii,1,2) = plot([x_tp1(ii),x_tp1(ii)],[min(cos(PHY{ii,1,1})),max(cos(PHY{ii,1,1}))],':','Color',[.8,.8,.8],'LineWidth',.5);
xlim([0,xmax]); ylim([-1.05,1.05]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$\cos(\Xi_y(x))$');
legend({'$\cos(\Xi_y(x))$',['$',fstxofy,'$']},'Location','northeast');
% y(x)
hax(ii,3) = subplot(4,2,3); hold on;
hln(ii,3,1) = plot(X{ii,1,1},Y{ii,1,1},'Color',matcolors(3),'LineWidth',1);
hln(ii,3,2) = plot([x_tp1(ii),x_tp1(ii)],[0,max(Y{ii,1,1})],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-.05.*ymax,1.05.*ymax]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$y(x)\;[\mathrm{m}]$'); 
legend({'$y(x)$',['$',fstxofy,'$']},'Location','northeast');
% ky(x)
hax(ii,5) = subplot(4,2,5); hold on;
hln(ii,5,1) = plot(X{ii,1,1},KY{ii,1,1},'Color',.6.*matcolors(5),'LineWidth',1);
hln(ii,5,2) = plot(X{ii,1,1},MKY{ii,1,1},'--','Color',matcolors(5),'LineWidth',1);
hln(ii,5,3) = plot([x_tp1(ii),x_tp1(ii)],[min(KY{ii,1,1}),max(KY{ii,1,1})],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-1.05.*mkymax(1),1.05.*mkymax(1)]); 
xlabel('$x\;[\mathrm{m}]$'); ylabel('$k_y(x)$');
legend({'$k_y(x)$','$|k_y(x)|$',['$',fstxofy,'$']},'Location','northeast'); 
% zCycPhs(x)
hax(ii,2) = subplot(4,2,2); hold on; 
title(['zCycMthd:1$\quad\mathbf{z}\mathrm{-cycle}\qquad',stheqstr{ii},'$']);
hln(ii,2,1) = plot(X{ii,1,1},cos(PHZ{ii,1,1}),'Color',matcolors(2),'LineWidth',1);
hln(ii,2,2) = plot([x_tp1(ii),x_tp1(ii)],[min(cos(PHZ{ii,1,1})),max(cos(PHZ{ii,1,1}))],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-1.05,1.05]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$\cos(\Xi_z(x))$');
legend({'$\cos(\Xi_z(x))$',['$',fstxofy,'$']},'Location','northeast');
% z(x)
hax(ii,4) = subplot(4,2,4); hold on;
hln(ii,4,1) = plot(X{ii,1,1},Z{ii,1,1},'Color',matcolors(4),'LineWidth',1);
hln(ii,4,2) = plot([x_tp1(ii),x_tp1(ii)],[0,zmax],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-.05.*zmax,1.05.*zmax]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$z(x)\;[\mathrm{m}]$');
legend({'$z(x)$',['$',fstxofy,'$']},'Location','northeast');
% kz(x)
hax(ii,6) = subplot(4,2,6); hold on;
hln(ii,6,1) = plot(X{ii,1,1},KZ{ii,1,1},'Color',.6.*matcolors(6),'LineWidth',1);
hln(ii,6,2) = plot(X{ii,1,1},MKZ{ii,1,1},'--','Color',matcolors(6),'LineWidth',1); 
hln(ii,6,3) = plot([x_tp1(ii),x_tp1(ii)],[min(KZ{ii,1,1}),max(KZ{ii,1,1})],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-1.05.*mkzmax(1),1.05.*mkzmax(1)]); 
xlabel('$x\;[\mathrm{m}]$'); ylabel('$k_z(x)$');
legend({'$k_z(x)$','$|k_z(x)|$',['$',fstxofy,'$']},'Location','northeast');
% dz/dx and kz/kx
hax(ii,7) = subplot(4,2,7:8); hold on;
hln(ii,7,1) = plot(X2{ii,1,1},abs(DZDX{ii,1,1}),'Color',matcolors(4),'LineWidth',1);
hln(ii,7,2) = plot(X{ii,1,1},abs(KZKX{ii,1,1}),'Color',matcolors(6),'LineWidth',1);
hln(ii,7,3) = plot([x_tp1(ii),x_tp1(ii)],[0,max(KZ{ii,1,1})],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-.05.*mkzmax(1)./kx(1),1.05.*mkzmax(1)./kx(1)]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$\left|\frac{\mathrm{d}z}{\mathrm{d}x}\right| \quad\mathrm{or}\quad \left|\frac{k_z}{k_x}\right|$');
legend({'$\left|\mathrm{d}z/\mathrm{d}x\right|$','$\left|k_z/k_x\right|$',['$',fstxofy,'$']});
% draw + save
drawnow; hfg(ii).Units = 'normalized'; hfg(ii).OuterPosition = scrn2pos; 
figname = ['isoWedgeCycleTest_zCycMthd1_fig',num2str(ii)];
if svfig; savefig(hfg(ii),fullfile(matDir,'fig',[figname,'.fig']),'compact'); end
if svpng; saveas(hfg(ii),fullfile(matDir,'img',[figname,'.png'])); end
end; end; end


% METHOD 2
if runflag(2); for ii=1:4; if subrunflag(ii); hfg(ii) = figure;
      htl(ii) = tiledlayout(hfg(ii),4,2,TileSpacing='tight',Padding='tight');
% yCycPhs(x)
hax(ii,1) = nexttile(1); hold on;
title(['$\mathbf{y}\mathrm{-cycle}\qquad',spheqstr{ii},'$']);
hln(ii,1,1) = plot(X{ii,1,2}./1e3,cos(PHY{ii,1,2}),'.-','Color',matcolors(1),'LineWidth',.5);
hln(ii,1,2) = plot([x_tp1(ii),x_tp1(ii)]./1e3,[min(cos(PHY{ii,1,2})),max(cos(PHY{ii,1,2}))],':','Color',[.1,.1,.1],'LineWidth',.5);
xlim([0,xmax]./1e3); ylim([-1.05,1.05]);
%xlabel('range: $\;x\;$ (km)'); 
ylabel('$\cos(\Xi_y(x))$');
legend({'$\cos(\Xi_y(x))$',['$',fstxofy,'$']},'Location','northeast');
% y(x)
hax(ii,3) = nexttile(3); hax(ii,3).YDir='reverse'; hold on;
hln(ii,3,1) = plot(X{ii,1,2}./1e3,Y{ii,1,2}./1e3,'Color',matcolors(3),'LineWidth',1);
hln(ii,3,2) = plot([x_tp1(ii),x_tp1(ii)]./1e3,[0,max(Y{ii,1,2})]./1e3,':','Color',[.1,.1,.1],'LineWidth',1);
xlim([0,xmax]./1e3); ylim([-.05.*ymax,1.05.*ymax]./1e3);
%xlabel('range: $\;x\;$ (km)');
ylabel('$y(x)\;[\mathrm{m}]$'); 
legend({'$y(x)$',['$',fstxofy,'$']},'Location','northeast');
% ky(x)
hax(ii,5) = nexttile(5); hax(ii,5).YDir = 'reverse'; hold on;
hln(ii,5,1) = plot(X{ii,1,2}./1e3,KY{ii,1,2},'Color',.6.*matcolors(5),'LineWidth',1);
hln(ii,5,2) = plot(X{ii,1,2}./1e3,MKY{ii,1,2},'--','Color',matcolors(5),'LineWidth',1);
hln(ii,5,3) = plot([x_tp1(ii),x_tp1(ii)]./1e3,[min(KY{ii,1,2}),max(KY{ii,1,2})],':','Color',[.1,.1,.1],'LineWidth',1);
xlim([0,xmax]./1e3); ylim([-1.05.*mkymax(1),1.05.*mkymax(1)]); 
xlabel('range: $\;x\;$ (km)'); ylabel('$k_y(x)$');
legend({'$k_y(x)$','$|k_y(x)|$',['$',fstxofy,'$']},'Location','northeast'); 
% zCycPhs(x)
hax(ii,2) = nexttile(2); hold on; 
title(['zCycMthd:2$\quad\mathbf{z}\mathrm{-cycle}\qquad',stheqstr{ii},'$']);
hln(ii,2,1) = plot(X{ii,1,2}./1e3,cos(PHZ{ii,1,2}),'Color',matcolors(2),'LineWidth',1);
hln(ii,2,2) = plot([x_tp1(ii),x_tp1(ii)]./1e3,[min(cos(PHZ{ii,1,2})),max(cos(PHZ{ii,1,2}))],':','Color',[.1,.1,.1],'LineWidth',1);
xlim([0,xmax]./1e3); ylim([-1.05,1.05]);
%xlabel('range: $\;x\;$ (km)');
ylabel('$\cos(\Xi_z(x))$');
legend({'$\cos(\Xi_z(x))$',['$',fstxofy,'$']},'Location','northeast');
% z(x)
hax(ii,4) = nexttile(4); hax(ii,4).YDir='reverse'; hold on;
hln(ii,4,1) = plot(X{ii,1,2}./1e3,Z{ii,1,2},'Color',matcolors(4),'LineWidth',1);
hln(ii,4,2) = plot([x_tp1(ii),x_tp1(ii)]./1e3,[0,zmax],':','Color',[.1,.1,.1],'LineWidth',1);
xlim([0,xmax]./1e3); ylim([-.05.*zmax,1.05.*zmax]);
%xlabel('range: $\;x\;$ (km)');
ylabel('$z(x)\;[\mathrm{m}]$');
legend({'$z(x)$',['$',fstxofy,'$']},'Location','northeast');
% kz(x)
hax(ii,6) = nexttile(6); hax(ii,6).YDir = 'reverse'; hold on;
hln(ii,6,1) = plot(X{ii,1,2}./1e3,KZ{ii,1,2},'Color',.6.*matcolors(6),'LineWidth',1);
hln(ii,6,2) = plot(X{ii,1,2}./1e3,MKZ{ii,1,2},'--','Color',matcolors(6),'LineWidth',1); 
hln(ii,6,3) = plot([x_tp1(ii),x_tp1(ii)]./1e3,[min(KZ{ii,1,2}),max(KZ{ii,1,2})],':','Color',[.1,.1,.1],'LineWidth',1);
xlim([0,xmax]./1e3); ylim([-1.05.*mkzmax(1),1.05.*mkzmax(1)]); 
xlabel('range: $\;x\;$ (km)'); ylabel('$k_z(x)$');
legend({'$k_z(x)$','$|k_z(x)|$',['$',fstxofy,'$']},'Location','northeast');
% dz/dx and kz/kx
hax(ii,7) = nexttile(7,[1 2]); hold on;
hln(ii,7,1) = plot(X2{ii,1,2}./1e3,abs(DZDX{ii,1,2}),'Color',matcolors(4),'LineWidth',1);
hln(ii,7,2) = plot(X{ii,1,2}./1e3,abs(KZKX{ii,1,2}),'Color',matcolors(6),'LineWidth',1);
hln(ii,7,3) = plot([x_tp1(ii),x_tp1(ii)]./1e3,[0,max(KZ{ii,1,2})],':','Color',[.1,.1,.1],'LineWidth',1);
xlim([0,xmax]./1e3); ylim([-.05.*mkzmax(1)./kx(1),1.05.*mkzmax(1)./kx(1)]);
xlabel('range: $\;x\;$ (km)'); ylabel('$\left|\frac{\mathrm{d}z}{\mathrm{d}x}\right| \quad\mathrm{or}\quad \left|\frac{k_z}{k_x}\right|$');
legend({'$\left|\mathrm{d}z/\mathrm{d}x\right|$','$\left|k_z/k_x\right|$',['$',fstxofy,'$']});
% draw + save
drawnow; hfg(ii).Units = 'normalized'; hfg(ii).OuterPosition = scrn2pos; 
figname = ['isoWedgeCycleTest_zCycMthd2_fig',num2str(ii)];
if svfig; savefig(hfg(ii),fullfile(matDir,'fig',[figname,'.fig']),'compact'); end
if svpng; saveas(hfg(ii),fullfile(matDir,'img',[figname,'.png'])); end
end; end; end

% METHOD 1 V 2
if runflag(3); for ii=1:4; if subrunflag(ii); hfg(ii) = figure;
% zCycPhs(x) Mthd1
hax(ii,1) = subplot(4,2,1); hold on; 
title(['zCycMthd:1$\quad\mathbf{z}\mathrm{-cycle}\qquad',stheqstr{ii},'$']);
hln(ii,1,1) = plot(X{ii,1,1},cos(PHZ{ii,1,1}),'Color',matcolors(1),'LineWidth',1);
hln(ii,1,2) = plot([x_tp1(ii),x_tp1(ii)],[min(cos(PHZ{ii,1,1})),max(cos(PHZ{ii,1,1}))],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-1.05,1.05]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$\cos(\Xi_z(x))$');
legend({'$\cos(\Xi_z(x))$',['$',fstxofy,'$']},'Location','northeast');
% z(x) Mthd1
hax(ii,3) = subplot(4,2,3); hold on;
hln(ii,3,1) = plot(X{ii,1,1},Z{ii,1,1},'Color',matcolors(3),'LineWidth',1);
hln(ii,3,2) = plot([x_tp1(ii),x_tp1(ii)],[0,zmax],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-.05.*zmax,1.05.*zmax]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$z(x)\;[\mathrm{m}]$');
legend({'$z(x)$',['$',fstxofy,'$']},'Location','northeast');
% kz(x) Mthd1
hax(ii,5) = subplot(4,2,5); hold on;
hln(ii,5,1) = plot(X{ii,1,1},KZ{ii,1,1},'Color',.6.*matcolors(5),'LineWidth',1);
hln(ii,5,2) = plot(X{ii,1,1},MKZ{ii,1,1},'--','Color',matcolors(5),'LineWidth',1); 
hln(ii,5,3) = plot([x_tp1(ii),x_tp1(ii)],[min(KZ{ii,1,1}),max(KZ{ii,1,1})],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-1.05.*mkzmax(1),1.05.*mkzmax(1)]); 
xlabel('$x\;[\mathrm{m}]$'); ylabel('$k_z(x)$');
legend({'$k_z(x)$','$|k_z(x)|$',['$',fstxofy,'$']},'Location','northeast');
% dz/dx and kz/kx Mthd1
hax(ii,7) = subplot(4,2,7); hold on;
hln(ii,7,1) = plot(X2{ii,1,1},abs(DZDX{ii,1,1}),':','Color',matcolors(3),'LineWidth',1);
hln(ii,7,2) = plot(X{ii,1,1},abs(KZKX{ii,1,1}),'--','Color',matcolors(5),'LineWidth',1);
hln(ii,7,3) = plot([x_tp1(ii),x_tp1(ii)],[0,max(KZ{ii,1,1})],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-.05.*mkzmax(1)./kx(1),1.05.*mkzmax(1)./kx(1)]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$\left|\frac{\mathrm{d}z}{\mathrm{d}x}\right| \quad\mathrm{or}\quad \left|\frac{k_z}{k_x}\right|$');
legend({'$\left|\mathrm{d}z/\mathrm{d}x\right|$','$\left|k_z/k_x\right|$',['$',fstxofy,'$']});
% zCycPhs(x) Mthd2
hax(ii,2) = subplot(4,2,2); hold on; 
title(['zCycMthd2$\quad\mathbf{z}\mathrm{-cycle}\qquad',stheqstr{ii},'$']);
hln(ii,2,1) = plot(X{ii,1,2},cos(PHZ{ii,1,2}),'Color',matcolors(2),'LineWidth',1);
hln(ii,2,2) = plot([x_tp1(ii),x_tp1(ii)],[min(cos(PHZ{ii,1,2})),max(cos(PHZ{ii,1,2}))],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-1.05,1.05]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$\cos(\Xi_z(x))$');
legend({'$\cos(\Xi_z(x))$',['$',fstxofy,'$']},'Location','northeast');
% z(x) Mthd2
hax(ii,4) = subplot(4,2,4); hold on;
hln(ii,4,1) = plot(X{ii,1,2},Z{ii,1,2},'Color',matcolors(4),'LineWidth',1);
hln(ii,4,2) = plot([x_tp1(ii),x_tp1(ii)],[0,zmax],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-.05.*zmax,1.05.*zmax]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$z(x)\;[\mathrm{m}]$');
legend({'$z(x)$',['$',fstxofy,'$']},'Location','northeast');
% kz(x) Mthd2
hax(ii,6) = subplot(4,2,6); hold on;
hln(ii,6,1) = plot(X{ii,1,2},KZ{ii,1,2},'Color',.6.*matcolors(6),'LineWidth',1);
hln(ii,6,2) = plot(X{ii,1,2},MKZ{ii,1,2},'--','Color',matcolors(6),'LineWidth',1); 
hln(ii,6,3) = plot([x_tp1(ii),x_tp1(ii)],[min(KZ{ii,1,2}),max(KZ{ii,1,2})],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-1.05.*mkzmax(1),1.05.*mkzmax(1)]); 
xlabel('$x\;[\mathrm{m}]$'); ylabel('$k_z(x)$');
legend({'$k_z(x)$','$|k_z(x)|$',['$',fstxofy,'$']},'Location','northeast');
% dz/dx and kz/kx Mthd2
hax(ii,8) = subplot(4,2,8); hold on;
hln(ii,8,1) = plot(X2{ii,1,2},abs(DZDX{ii,1,2}),':','Color',matcolors(4),'LineWidth',1);
hln(ii,8,2) = plot(X{ii,1,2},abs(KZKX{ii,1,2}),'--','Color',matcolors(6),'LineWidth',1);
hln(ii,8,3) = plot([x_tp1(ii),x_tp1(ii)],[0,max(KZ{ii,1,2})],':','Color',[.8,.8,.8],'LineWidth',1);
xlim([0,xmax]); ylim([-.05.*mkzmax(1)./kx(1),1.05.*mkzmax(1)./kx(1)]);
xlabel('$x\;[\mathrm{m}]$'); ylabel('$\left|\frac{\mathrm{d}z}{\mathrm{d}x}\right| \quad\mathrm{or}\quad \left|\frac{k_z}{k_x}\right|$');
legend({'$\left|\mathrm{d}z/\mathrm{d}x\right|$','$\left|k_z/k_x\right|$',['$',fstxofy,'$']});
% draw + save
drawnow; hfg(ii).Units = 'normalized'; hfg(ii).OuterPosition = scrn2pos; 
figname = ['isoWedgeCycleTest_','zCycMthd1v2','_ray',num2str(ii)];
if svfig; savefig(hfg(ii),fullfile(matDir,'fig',[figname,'.fig']),'compact'); end
if svpng; saveas(hfg(ii),fullfile(matDir,'img',[figname,'.png'])); end
end; end; end



%% 3D plot (all rays)
for ii = 1:4
xSq{ii} = repmat(x_tp1(ii),[2,2]);
ySq{ii} = [repmat(y_tp1(ii),[2,1]),[0;0]];
zSq{ii} = [0,0;H(y_tp1(ii)),0];
end


%% METHOD 1
if runflag(4)
hfg(5) = figure; hax(1,8,1) = subplot(1,1,1); colormap(permute(g2r(n),[2,3,1])); lgdstrs = {};
ls={'-','--','-.',':'}; msh={'o','square','v','pentagram'}; hold on;
% plot rays
for ii = 1:4; if subrunflag(ii); lgdstrs = [lgdstrs; {['$\vec{x}$-ray$\,',sthsphstr{ii},'$']}];
hln(ii,8,1) = patch('XData',X{ii,1,1}.','YData',Y{ii,1,1}.','ZData',Z{ii,1,1}.','CData',C{ii,1,1}.', ...
  'FaceColor','none','EdgeColor','interp','LineStyle',ls{ii}); end; end
% plot 1st turn loctn
for ii = 1:4; if subrunflag(ii); lgdstrs = [lgdstrs; {['$',fstxofy,'\,',sthsphstr{ii},'$']}];
    hln(ii,9,1) = surf(xSq{ii},ySq{ii},zSq{ii},repmat(.8,[2,2,3]),'FaceColor','none',...
  'EdgeAlpha',.5,'FaceAlpha',0,'LineWidth',.5,'LineStyle',ls{ii},'Marker','.','MarkerSize',6, ...
  'MarkerEdgeColor',[.8,.8,.8]); end; end
% plot src ray dir markers
for ii = 1:4; a0=100.*alf0(ii,:); if subrunflag(ii); lgdstrs = [lgdstrs; {['$\vec{k}$-dir$\,',sthsphstr{ii},'$']}];
hln(ii,10,1) = scatter3(x0+a0(1),y0+a0(2),z0+a0(3),'Marker',msh{ii},'LineWidth',.5,...
  'MarkerEdgeColor','c','MarkerFaceColor','c','MarkerEdgeAlpha',.8,'MarkerFaceAlpha',.5,'SizeData',10); end; end
% plot src|srf|bot
hscat = scatter3(x0-veps,y0,z0,'Marker','o','LineWidth',.5,'MarkerEdgeColor','r','MarkerEdgeAlpha',.8,...
  'MarkerFaceColor','r','MarkerFaceAlpha',.5,'SizeData',10);
hsf(1) = surf([0,xmax;0,xmax],[0,0;ymax,ymax],[0,0;0,0],repmat(cat(3,.5,.5,.9),2,2,1),'FaceAlpha',.05);
hsf(2) = surf([0,xmax;0,xmax],[0,0;ymax,ymax],[0,0;zmax,zmax],'CData',[cat(3,.5,.5,.9),cat(3,.5,.5,.9);cat(3,0,.1,.2),cat(3,0,.1,.2)],...
  'AlphaData',[.1,.1;.6,.6],'FaceAlpha','interp');
% addtnl format
xlabel('$x\;[\mathrm{m}]$'); ylabel('$y\;[\mathrm{m}]$'); zlabel('$z\;[\mathrm{m}]$');
hlg = legend([lgdstrs;{'source';'surface';'seafloor'}],'Location','northeastoutside');
title(['Ray Trajectories (zCycMthd:1)$\qquad',mtheqstr{1},'\qquad',mpheqstr{1},'$']); hax(1,8,1).ZDir = 'reverse'; 
drawnow; hfg(5).Units = 'normalized'; hfg(5).OuterPosition = scrn2pos;
hsf(1).EdgeColor=[1,1,1]; hsf(1).EdgeAlpha=.5; hsf(2).EdgeColor=[1,1,1]; hsf(2).EdgeAlpha=.5;
% loop views + draw + save
for ii = 1:size(viewList,1)
  figname = ['isoWedgeCycleTest_zCycMthd1_fig5_view',num2str(ii)]; 
  view(viewList(ii,1),viewList(ii,2)); hlg.Location='northeastoutside'; drawnow;
  if svpng; saveas(hfg(5),fullfile(matDir,'img',[figname,'.png'])); end
end
figname = ['isoWedgeCycleTest_zCycMthd1_fig5'];
if svfig; savefig(hfg(5),fullfile(matDir,'fig',[figname,'.fig']),'compact'); end
end

%% METHOD 2
if runflag(5)
hfg(5) = figure; hax(1,8,1) = subplot(1,1,1); lgdstrs = {}; 
%colormap(permute(g2r(n),[2,3,1])); 
colormap([0,0,0;.2,.2,.2])
msh={'square','o','v','diamond'}; 
msz=[30,20,20,30];
lw = [1,1,1,2]; ls = {'-','--','-.',':'}; 
mfc = {[.1,.1,1],[1,.1,.1],[.1,1,.1],[.2,.2,0]};
mec = {[0,0,.5],[.5,0,0],[0,.5,0],[.1,.1,0]};
rlc = {[0,0,1],[1,0,0],[0,1,0],[0,0,0]};
hold on;
% plot rays
for ii = 1:4; if subrunflag(ii); lgdstrs = [lgdstrs; {['$\vec{x}$-ray$\,',sthsphstr{ii},'$']}];
%hln(ii,8,1) = patch('XData',X{ii,1,2}.'./1e3,'YData',Y{ii,1,2}.'./1e3,'ZData',Z{ii,1,2}.', ...
%  'CData',repmat(shiftdim(rlc{ii},-1),[1,n+1,1]),'LineWidth',4, ...
%  'FaceColor','none','EdgeColor','interp','LineStyle',ls{ii}); 
hln(ii,8,1) = plot3(X{ii,1,2}.'./1e3,Y{ii,1,2}.'./1e3,Z{ii,1,2}.','Color',rlc{ii},'LineWidth',lw(ii),'LineStyle',ls{ii}); 
end; end
% plot 1st turn loctn (vert surfaces)
for ii = 1:4; if subrunflag(ii); lgdstrs = [lgdstrs; {['$',fstxofy,'\,',sthsphstr{ii},'$']}];
    hln(ii,9,1) = surf(xSq{ii}./1e3,ySq{ii}./1e3,zSq{ii},repmat(.8,[2,2,3]),'FaceColor','none',...
  'EdgeAlpha',.5,'FaceAlpha',0,'LineWidth',.5,'LineStyle',ls{ii},'Marker','.','MarkerSize',6, ...
  'MarkerEdgeColor',[.8,.8,.8]); end; end
% plot src ray dir markers
for ii = 1:4; a0=100.*alf0(ii,:); if subrunflag(ii); lgdstrs = [lgdstrs; {['$\vec{k}$-dir$\,',sthsphstr{ii},'$']}];
hln(ii,10,1) = scatter3((x0+a0(1))./1e3,(y0+a0(2))./1e3,z0+a0(3),'Marker',msh{ii},'LineWidth',.5,...
  'MarkerEdgeColor',mec{ii},'MarkerFaceColor',mfc{ii},'MarkerEdgeAlpha',.9,'MarkerFaceAlpha',.8,'SizeData',msz(ii)); end; end
% plot source position
hscat = scatter3((x0-veps)./1e3,y0./1e3,z0,'Marker','pentagram','LineWidth',1,'MarkerEdgeColor','k','MarkerEdgeAlpha',.8,...
  'MarkerFaceColor','w','MarkerFaceAlpha',.5,'SizeData',60);
% plot sea surface
hsf(1) = surf([0,xmax;0,xmax]./1e3,[0,0;ymax,ymax]./1e3,[0,0;0,0],repmat(cat(3,.5,.5,.9),2,2,1),'FaceAlpha',.1);
% plot sea bottom
hsf(2) = surf([0,xmax;0,xmax]./1e3,[0,0;ymax,ymax]./1e3,[0,0;zmax,zmax],repmat(cat(3,.1,.05,.0),2,2,1),'FaceAlpha',.2);
% plot ymax wall
%hsf(3) = surf([0,xmax;0,xmax],[ymax,ymax;ymax,ymax],[0,0;zmax,zmax],repmat(cat(3,.1,.1,.1),2,2,1),FaceAlpha=.05,EdgeAlpha=.2);
%hsf(2) = surf([0,xmax;0,xmax],[0,0;ymax,ymax],[0,0;zmax,zmax],'CData',[cat(3,.5,.5,.9),cat(3,.5,.5,.9);cat(3,0,.1,.2),cat(3,0,.1,.2)],...
%  'AlphaData',[.1,.1;.6,.6],'FaceAlpha','interp');
% addtnl format
xlabel('range: $\;x\;$ (km)'); ylabel('range: $\;y\;$ (km)'); zlabel('depth: $\;z\;$ (m)');
hlg = legend([lgdstrs;{'source';'surface';'seafloor'}],'Location','northeastoutside');
title(['Ray Trajectories$\qquad',mtheqstr{1},'\qquad',mpheqstr{1},'$']); 
hax(1,8,1).ZDir = 'reverse'; hax(1,8,1).XDir = 'reverse';
drawnow; hfg(5).Units = 'normalized'; hfg(5).OuterPosition = scrn2pos;
hsf(1).EdgeColor=[0,0,0]; hsf(1).EdgeAlpha=.2; hsf(2).EdgeColor=[0,0,0]; hsf(2).EdgeAlpha=.2;
% loop views + draw + save
for ii = 1:size(viewList,1)
  figname = ['isoWedgeCycleTest_zCycMthd2_fig5_view',num2str(ii)]; 
  view(viewList(ii,1),viewList(ii,2)); hlg.Location='northeastoutside'; drawnow;
  if svpng; saveas(hfg(5),fullfile(matDir,'img',[figname,'.png'])); end
end
figname = ['isoWedgeCycleTest_zCycMthd2_fig5'];
if svfig; savefig(hfg(5),fullfile(matDir,'fig',[figname,'.fig']),'compact'); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [clr] = matcolors(n,fmt)

% n = color index
% fmt = 'rgb' or 'hex'
% clr = matlab color in 'rgb' or 'hex' format

if nargin < 2; fmt = 'rgb'; end
if nargin < 1; n = 0; end

rgbList = {[.001,.447,.741];
            [.85,.325,.098];
            [.929,.694,.125];
            [.494,.184,.556];
            [.466,.674,.188];
            [.301,.745,.933];
            [.635,.078,.184]};
hexList = {'#0072BD';
           '#D95319';
           '#EDB120';
           '#7E2F8E';
           '#77AC30';
           '#4DBEEE';
           '#A2142F'};

switch fmt
case 'rgb'
  outList = rgbList;
case 'hex'
  outList = hexList;
end % switch

switch n
case 0
  clr = cell2mat(outList);
otherwise
  clr = outList{n};
end % switch
end % function

function [rgbMat] = g2r(n)
%[rgbMat]=g2r(n)
arguments; n (1,1) {mustBeInteger(n)} = 8; end
res = 1001; ibeg = 440; iend = 880;
cm = turbo(res); cm=cm(ibeg:iend,:);
idx=((ibeg:iend).'-ibeg)./(iend-ibeg);
rgbMat = interp1(idx,cm,(0:n-1)./(n-1));
rgbMat = permute(rgbMat,[3 1 2]); % for dimensions of patch CData
end

