% % Code developed by Jerome Pinti
% % For questions and remarks, please email pintijerome@gmail.com
% % This code is based on the transport matrix method developed by Tim
% % DeVries from UCSB. The transport matrix can be downloaded there: https://tdevries.eri.ucsb.edu/models-and-data-products/
% % References to the transport matrix are DeVries and Primeau 2011 and DeVries 2014.

%% Load the transport matrix
% load('91x180x48_omega3_Kvvariable_geoLuc_AH2e5.mat') %high resolution - ~10 minutes to inverse

load('CTL.mat') % lower resolution - less than a minute

%% Preparation of the transport matrix
grid = output.grid;
TR = output.TR; % yr^-1
msk = output.msk;
M3d = output.M3d; % land = 0. ocean = 1


VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d;
V = VOL(msk.pkeep);

m = size(TR,1);
sink = zeros(m,1);
sink(1:length(msk.hkeep)) = 1e10; % instantaneous Surface SINK
SSINK = spdiags(sink,0,m,m);
A = TR-SSINK; % transport + sink in surface

% for load run-compiled-236
[lonq,latq,zq] = meshgrid(x,y,z);
Q = permute(Q, [2 1 3]);
q_ocim = interp3(lonq,latq,zq,Q,grid.XT3d,grid.YT3d,grid.ZT3d);
q_ocim = q_ocim(msk.pkeep);
q_ocim(isnan(q_ocim)) = 0;
export = V'*q_ocim / 1e15; % [PgC / yr]
tic
cseq = -A \q_ocim;
toc
TotCseq = V'*cseq / 1e15
seqime = TotCseq / export




% Re interpolate to Q
qocim = 0*M3d+NaN; % make a 3-d array of NaN
qocim(msk.pkeep) = q_ocim;
qocim(isnan(qocim)) = 0;

C_eq = 0*M3d+NaN;
C_eq(msk.pkeep) = cseq;
C_eq(isnan(C_eq)) = 0;

%% A simple plot
LAT = repmat(y,nx+1,1);
LON = repmat([x,360]',1,ny);
cc = sum(C_eq.*grid.DZT3d,3); cc = [cc, cc(:,end)]; % gC/m2/day
%North Pacific
xred = [140 180 180 140 140];
yred = [37 37 50 50 37];

% %North Atlantic
% xred = [320 350 350 320 320];
% yred = [53 53 66 66 53];

figure
cc(cc<0) = 0;
ax = axesm('mollweid','Frame','on','MapLatLimit',[-90 90],'Origin', [0 -160 0],'FLineWidth',0.5);
ax.XTick = [-120 -60 0 60 120 180];
ax.YTick = [-40 -20 0 20 40];
% objects = [handlem('grid'); handlem('mlabel'); handlem('plabel')];
% set(objects,'Clipping','on');
box off
axis off
load coast
geoshow(lat, long,'Color','k')

surfm(LAT', LON', cc,'AlphaData',~isnan(cc));%,'EdgeColor','none')
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
plotm(yred,xred,'r')
colorbar

%% Figure from the poles
figure
subplot(121)
ax = axesm('eqdazim','Frame','on','MapLatLimit',[10 90],'Origin', [90 0 0],'FLineWidth',0.5);
% eqdazim
box off
axis off
load coast
% geoshow(lat, long,'Color','k')
surfm(LAT', LON', cc,'AlphaData',~isnan(cc));
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
colorbar
plotm(yred,xred,'r')


subplot(122)
ax = axesm('eqdazim','Frame','on','MapLatLimit',[-90 -10],'Origin', [-90 -180 0],'FLineWidth',0.5);
% eqdazim
box off
axis off
load coast
% geoshow(lat, long,'Color','k')
surfm(LAT', LON', cc,'AlphaData',~isnan(cc));
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
colorbar