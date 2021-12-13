% % Code developed by Jerome Pinti
% % For questions and remarks, please email pintijerome@gmail.com


%% Make base grid
x = 0:1:359; nx = size(x,2); % [degrees] longitude 
y = -90:1:90; ny = size(y,2); % [degrees] latitude 
z = 0:100:8000; nz = size(z,2); % [m] depth
dz = 100; % [m] vertical resolution


%% Make input for rough finmarchicus calculation ATLANTIC
a = 2.3; % [gC / m2 / day] Export 

Q = zeros(nx,ny,nz);
s = gaussmf(z,[50 1000]); %first is sd and then mean
s = s/sum(s); %integral of the gaussian is 1
s = a * s / dz; % [gC / m3 / yr] Respiration rate at each depth
s = reshape(s,[1 1 size(s,2)]);

longres = (x>(360-40)) & (x<(360-10));
latres = (y>53) & (y<66);
Q(longres,latres,:) = repmat(s, sum(longres), sum(latres), 1); % [gC / m3/ yr]

%% Make input for rough Neocalanus calculation PACIFIC
a = 4.3; % [gC / m2 / day] Export 

Q = zeros(nx,ny,nz);
s = gaussmf(z,[50 1000]); %first is sd and then mean
s = s/sum(s); %integral of the gaussian is 1
s = a * s / dz; % [gC / m3 / yr] Respiration rate at each depth
s = reshape(s,[1 1 size(s,2)]);

longres = (x>(140)) & (x<(180));
latres = (y>37) & (y<50);
Q(longres,latres,:) = repmat(s, sum(longres), sum(latres), 1); % [gC / m3/ yr]