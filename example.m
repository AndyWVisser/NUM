% use the OCIM circulation to find the fraction of water ventilated from
% south of 65oS, and the water age

% load the model output
load CTL.mat output % the "control" version of the model -- see DeVries (2014) for
		    % different versions
M3d = output.M3d; % 1 = ocean, 0 = land
grid = output.grid; % grid metrics
MSKS = output.MSKS; % basin masks
iocn = find(M3d(:)==1);
m = length(iocn);
TR = output.TR; % transport operator [yr^-1]

% set surface boundary condition to track waters ventilated in S. Ocean
% south of 65oS
% REG = 1 everywhere south of 65oS
% REG = 0 elsewhere
REG = 0*M3d;
REG(grid.YT3d<=-65) = 1; % latitudes south of 65oS
F = eq_wmfrac(TR,REG,M3d,grid); % ventilation fractions

% zonally average for Atlantic and Pacific basins
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d.*M3d;
Fatl = F.*MSKS.ATL.*VOL;
Fatl = squeeze(nansum(Fatl,2)./nansum(MSKS.ATL.*VOL,2));
Fpac = F.*MSKS.PAC.*VOL;
Fpac = squeeze(nansum(Fpac,2)./nansum(MSKS.PAC.*VOL,2));

% plot
figure(1)
contourf(grid.yt,-grid.zt,Fpac',[0:.05:1])
colormap(flipud(hot(20)))
set(gca,'CLim',[0 1])
title('Fraction of AABW in Pacific Ocean')
colorbar

% calculate water age
age = eqage(TR,grid,M3d); % age in years

% zonally average for Atlantic and Pacific basins
ageatl = age.*MSKS.ATL.*VOL;
ageatl = squeeze(nansum(ageatl,2)./nansum(MSKS.ATL.*VOL,2));
agepac = age.*MSKS.PAC.*VOL;
agepac = squeeze(nansum(agepac,2)./nansum(MSKS.PAC.*VOL,2));

% plot
figure(2)
contourf(grid.yt,-grid.zt,agepac',[0:100:1500])
colormap(flipud(spring(30).^.5).*flipud(bone(30).^.5))
set(gca,'CLim',[0 1500])
title('Age of waters in Pacific Ocean')
colorbar

