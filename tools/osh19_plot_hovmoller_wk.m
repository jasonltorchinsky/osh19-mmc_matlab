function osh19_plot_hovmoller_wk(params, q_mode, window_ndays, ...
    ovlp_ndays, percent_remove, alt, labels)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Make Hovmoller plot in first tile
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Get path to output data.
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, params.component_name);

addpath(out_path);
addpath(exp_path);
addpath(component_path);

% Create directory to hold plots
plot_path = fullfile(exp_path, 'plots');
if ~(isfolder(plot_path))
    mkdir(plot_path);
end

addpath(plot_path);

% Read in run parameters
params_file = fullfile(component_path, 'params.nc');

nx          = ncread(params_file, 'nx');
ny          = ncread(params_file, 'ny');
nz          = ncread(params_file, 'nz');
H           = ncread(params_file, 'H');
P_E         = ncread(params_file, 'P_E');
sim_days    = ncread(params_file, 'sim_days');
out_freq    = ncread(params_file, 'out_freq');

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

xx  = ncread(grid_file, 'xx');
yy  = ncread(grid_file, 'yy');
zzW = ncread(grid_file, 'zzW');
dy  = ncread(grid_file, 'dy');
dz  = ncread(grid_file, 'dz');

lons = 360 * (xx - xx(1)) / (xx(end) - xx(1));
lats = yy / 110.567;

% Scale meridional, vertical coordinate for basis functions
L = 1490; % Equatorial meridional length scale (km)
yy_norm = yy / L;
dy_norm = dy / L;
zzW_norm = pi * zzW / H;

% Get first baroclinic mode, zeroth parabolic cylinder function, standard
% deviation scale of u for projections
parab_cyl_0   = parab_cyl(yy_norm, 0);
q_clin_mode_1 = q_clin_mode(zzW_norm, 1);
q_clin_mode_2 = q_clin_mode(zzW_norm, 2);

% Get the normalization constants due to discretized norms
merid_norm  = dy_norm * (parab_cyl_0.' * parab_cyl_0);
vert_norm_1 = (1/(nz+1)) * (q_clin_mode_1.' * q_clin_mode_1);
vert_norm_2 = (1/(nz+1)) * (q_clin_mode_2.' * q_clin_mode_2);

% Get moisture, potential temperature anomaly at desired altitude and latitude
q     = zeros(floor(sim_days/out_freq) + 1, nx);
t     = zeros(1, floor(sim_days/out_freq) + 1);

n_outfiles = floor(sim_days/out_freq)+1;
out_idxs  = 0:n_outfiles-1;

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    t(1, out_idx + 1) = ncread(state_file, 't');
    q_temp     = ncread(state_file, 'q');
    
    q1 = zeros([ny, nx]);
    q2 = zeros([ny, nx]);
    for jj = 1:ny
       for ii = 1:nx
           q1(jj, ii) = (1/(nz+1)) ...
               * squeeze(q_clin_mode_1.'*squeeze(q_temp(jj, ii, :))) ...
               * (1/vert_norm_1);
           q2(jj, ii) = (1/(nz+1)) ...
               * squeeze(q_clin_mode_2.'*squeeze(q_temp(jj, ii, :))) ...
               * (1/vert_norm_2);
       end
    end
    
    % Pick Q_mid or Q_up REFURB, ACTUALLY LOW
    if strcmpi(q_mode, 'mid') % is actually low here
        Q = 1 / sqrt(3) * (q1 * q_clin_mode(pi/3, 1) - q2 * q_clin_mode(pi/3, 2));
    else
        Q = 1 / sqrt(3) * (q1 * q_clin_mode(2*pi/3, 1) + q2 * q_clin_mode(2*pi/3, 2));
    end
    
    % Project Q onto zeroth parabolic cylinder function
    for ii = 1:nx
        q(out_idx+1, ii) = dy_norm * squeeze(parab_cyl_0.'*Q(:, ii)) ...
            * (1/merid_norm);
    end
end

days_to_secs = 3600*24;
t = t / days_to_secs;

% Create plot

tlo = tiledlayout(1,3);
tlo.TileSpacing = 'compact';
tlo.Padding = 'compact';

tlo_hovmoller = tiledlayout(tlo, 1, 1);
tlo_hovmoller.Layout.Tile = 1;
tlo_hovmoller.Layout.TileSpan = [1 1];
tlo_hovmoller.TileSpacing = 'compact';
tlo_hovmoller.Padding = 'compact';


h(1) = nexttile(tlo_hovmoller, 1);

hold on;
[~, q_plt] = contourf(lons, t, q, ...
    'edgecolor', 'none');

% Colorbar
max_q = max(abs(q), [], 'all');
min_q = -max_q;
cmap = load('rb.mat').rb;
colormap(cmap);
caxis([min_q, max_q]);

cb = colorbar();
cb.Label.String = 'Moisture Anomaly (kg kg^{-1})';

% Title
title_str = 'Hovmoller Diagram';
title(tlo_hovmoller, title_str);

% Axis limits
xlim([0, 360]);
ylim([0, sim_days]);

% Axis labels
xlabel('Longitude');
ylabel('Time (d)');

% Ticks
xticks(0:45:360);
xticklabels({'180', '135W', '90W', '45W', '0', '45E', '90E', '135E', '180'});

% Example zonal windspeeds
speed_3 = 3.1 * (1 / 1000) * (360 / (2 * pi * P_E)) * days_to_secs;
plot(lons, (1./speed_3).*lons, 'k-.');
if labels
    text(10, 0.02*sim_days, '3.1 m s^{-1}',...
    'Color', 'black',...
    'FontWeight', 'bold',...
    'LineStyle', '-',...
    'LineWidth', 0.25,...
    'Margin', 0.25);
end

speed_4 = 4.6 * (1 / 1000) * (360 / (2 * pi * P_E)) * days_to_secs;
plot(lons, (1./speed_4).*lons + sim_days/2, 'k--');
if labels
    text(10, sim_days/2+(0.02*sim_days), '4.6 m s^{-1}',...
    'Color', 'black',...
    'FontWeight', 'bold',...
    'LineStyle', '-',...
    'LineWidth', 0.25,...
    'Margin', 0.25);
end

hold off;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Wheeler-Kiladis diagram in next two tiles. Just copy-paste the WK script
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Read the input parameters to find out how long the run is and how
% frequently output was created.
params_file = fullfile(component_path, 'params.nc');

nx       = ncread(params_file, 'nx');
P_Y      = ncread(params_file, 'P_Y');
P_E      = ncread(params_file, 'P_E');
beta     = ncread(params_file, 'beta')*(3600*24);
g        = ncread(params_file, 'g')*((3600*24).^2); % convert to km day^(-2)

sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');


% Other parameters
L           = 1490; % Equatorial meridional length scale (km)
Fs          = cast(out_freq, 'double');  % Sampling frequency (day^{-1})
circ_earth  = 2*pi*P_E; % Equatorial circumeference of earth (km)
max_wavenum = 10;
max_freq    = 0.5;

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

yy  = ncread(grid_file, 'yy');
zzW = ncread(grid_file, 'zzW');

dy   = ncread(grid_file, 'dy');

% Get correct index for desired altitude and latitude.
[~, altW_idx] = min(abs(zzW-alt));

altW_true = zzW(altW_idx);

% Find the indices for the channel quarter-widths
[~, min_yy_idx] = min(abs(yy - (-P_Y/2)));
[~, max_yy_idx] = min(abs(yy - P_Y/2));

% Trim data to only middle half of channel
yy      = yy(min_yy_idx:max_yy_idx);
yy_norm = yy/L;
ny      = size(yy, 1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Create power spectra for symmetric and anti-symmetric parts
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

qanom = zeros(ny, nx, floor(sim_days/out_freq) + 1);
days = 0:out_freq:sim_days;
half = ceil(nx/2);

for day = days
    state_file_name = strcat(['state_', num2str(day,'%04u'),'.nc']);
    state_file = fullfile(component_path, state_file_name);
    q_day = ncread(state_file, 'q');
    qanom(:, :, day + 1) = q_day(min_yy_idx:max_yy_idx, :, altW_idx)/1000;
    % Units? Want g kg^(-1).
end

% Get symmetric and anti-symmetric parts of qanom
qanom_sym   = detrend(qanom, 0);
qanom_asym  = qanom_sym;

parab_cyl_0 = parab_cyl(yy_norm, 0);
parab_cyl_1 = parab_cyl(yy_norm, 1);

for ii = 1:nx
    for day = days
        qanom_sym(:,ii,day+1)  = qanom(:,ii,day+1).*parab_cyl_0;
        qanom_asym(:,ii,day+1) = qanom(:,ii,day+1).*parab_cyl_1;
    end
end

qanom_sym  = squeeze(sum(qanom_sym,  1)) * dy;
qanom_asym = squeeze(sum(qanom_asym, 1)) * dy;

% Remove running means
for day = window_ndays:size(days,2)
    qanom_sym(:,day) = qanom_sym(:,day) ...
        - mean(qanom_sym(:,day-window_ndays+1:day),2);
    qanom_asym(:,day) = qanom_asym(:,day) ...
        - mean(qanom_asym(:,day-window_ndays+1:day),2);
end

% Remove data from spin-up time
qanom_sym  = qanom_sym(:,window_ndays+1:end);
qanom_asym = qanom_asym(:,window_ndays+1:end);

% Get windows of data to take
ndays_trim = size(qanom_sym, 2);
nwindows = floor((ndays_trim - window_ndays) / ovlp_ndays) + 1;

qanom_sym_windows = zeros(nwindows, nx, window_ndays);
qanom_asym_windows = zeros(nwindows, nx, window_ndays);
for window = 1:nwindows
    qanom_sym_windows(window,:,:) = qanom_sym(:, ...
        ovlp_ndays*(window-1)+1:window_ndays + ovlp_ndays*(window-1));
    qanom_asym_windows(window,:,:) = qanom_asym(:, ...
        ovlp_ndays*(window-1)+1:window_ndays + ovlp_ndays*(window-1));
end

qanom_sym_windows_power = zeros(nwindows, window_ndays, 2*max_wavenum);
qanom_asym_windows_power = zeros(nwindows, window_ndays, 2*max_wavenum);
    
% Calculate the power spectra of each window
cos_taper = ones(nx, window_ndays);
for ii = 1:nx
   low_percent_idx = round(window_ndays*percent_remove);
   high_percent_idx = window_ndays - low_percent_idx;
   cos_taper(ii, 1:low_percent_idx) = 0.5 ...
       * ( 1.0 - cos(pi*[0:low_percent_idx-1]/(low_percent_idx-1)) );
   cos_taper(ii,high_percent_idx+1:end) = 0.5 ...
       * ( 1.0 + cos(pi*[0:low_percent_idx-1]/(low_percent_idx-1)) );
end
for window = 1:nwindows
    qanom_sym = cos_taper.*squeeze(qanom_sym_windows(window,:,:));
    qanom_sym_fft1 = zeros(size(qanom_sym));
    qanom_sym_fft2 = zeros(size(qanom_sym));
    for day = 1:window_ndays
        qanom_sym_fft1(:,day) = fftshift(fft(squeeze(qanom_sym(:,day))));
    end
    for ii = 1:nx
        qanom_sym_fft2(ii,:)  = fftshift(fft(squeeze(qanom_sym_fft1(ii,:))));
    end
    qanom_sym_power = transpose(qanom_sym_fft2.*conj(qanom_sym_fft2));
    qanom_sym_power = log(fliplr(squeeze( ...
        qanom_sym_power(:,half-max_wavenum+1:half+max_wavenum) )));
    qanom_sym_windows_power(window,:,:) = qanom_sym_power;
    
    qanom_asym = cos_taper.*squeeze(qanom_asym_windows(window,:,:));
    qanom_asym_fft1 = zeros(size(qanom_asym));
    qanom_asym_fft2 = zeros(size(qanom_asym));
    for day = 1:window_ndays
        qanom_asym_fft1(:,day) = fftshift(fft(squeeze(qanom_asym(:,day))));
    end
    for ii = 1:nx
        qanom_asym_fft2(ii,:)  = fftshift(fft(squeeze(qanom_asym_fft1(ii,:))));
    end
    qanom_asym_power = transpose(qanom_asym_fft2.*conj(qanom_asym_fft2));
    qanom_asym_power = log(fliplr(squeeze( ...
        qanom_asym_power(:,half-max_wavenum+1:half+max_wavenum) )));
    
    qanom_asym_windows_power(window,:,:) = qanom_asym_power;
end

qanom_sym_power = squeeze(mean(qanom_sym_windows_power, 1));
qanom_asym_power = squeeze(mean(qanom_asym_windows_power, 1));

power_min = min([qanom_sym_power, qanom_asym_power], [], 'all');
power_max = max([qanom_sym_power, qanom_asym_power], [], 'all');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Calculate the dispersion curves
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wavenums = -max_wavenum+1:max_wavenum;
freqs    = ((-window_ndays/2+1):(window_ndays/2)) * (Fs / window_ndays);

nwvs = nx;
disp_wavenums = linspace(-max_wavenum+1, max_wavenum, nwvs);
half = nwvs/2;
syms disp_omega;

% Different Equivalent Depths

% h_e = 12 m
h_e = 12*10^(-3); % Equivalent depth (km)

n    = -1; % Kelvin waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
kelv12 = roots(1, :); % Index 1 is Kelvin Waves

n    = 0; % Mixed Rossby-Gravity Waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
rossbygrav12  = [roots(2, 1:half-11) roots(3, half-10:end)]; 
     % Index 2, 3 is Mixed Rossby-Gravity Waves

n    = 1; % Inertio-gravity and Rossby waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
rossby12  = roots(2, :); % Index 2 is Equatorial Rossby Waves
inertio12 = roots(3, :); % Index 3 is Inertio-Gracity Waves

% h_e = 25 m
h_e = 25*10^(-3); % Equivalent depth (km)

n    = -1; % Kelvin waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
kelv25 = roots(1, :); % Index 3 is Kelvin Waves

n    = 0; % Mixed Rossby-Gravity Waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
rossbygrav25  = [roots(2, 1:half-10) roots(3, half-9:end)]; 
     % Index 2, 3 is Mixed Rossby-Gravity Waves

n    = 1; % Inertio-gravity and Rossby waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
rossby25  = roots(2, :); % Index 2 is Equatorial Rossby Waves
inertio25 = roots(3, :); % Index 3 is Inertio-Gravity Waves


% h_e = 50 m
h_e = 50*10^(-3); % Equivalent depth (km)

n    = -1; % Kelvin waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
kelv50 = roots(1, :); % Index 1 is Kelvin Waves

n    = 0; % Mixed Rossby-Gravity Waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
rossbygrav50  = [roots(2, 1:half-8) roots(3, half-7:end)]; 
     % Index 2, 3 is Mixed Rossby-Gravity Waves

n    = 1; % Inertio-gravity and Rossby waves
eqns = sqrt(g * h_e) / beta ...
    * ( (4 * pi^2 * disp_omega^2 / (g * h_e)) ...
        - (4 * pi^2 * disp_wavenums.^2 / circ_earth^2) ...
        - (disp_wavenums/(disp_omega * circ_earth) * beta) ) == 2*n+1;
roots = zeros(3, nwvs);
for ii = 1:nwvs
    roots(:,ii) = vpasolve(eqns(ii), disp_omega);
end
rossby50  = roots(2, :); % Index 2 is Equatorial Rossby Waves
inertio50 = roots(3, :); % Index 3 is Inertio-Gracity Waves

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example periods
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

per3  = 1/3 * ones(size(wavenums));
per6  = 1/6 * ones(size(wavenums));
per30 = 1/30 * ones(size(wavenums));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Left tile - Wheeler-Kiladis Diagram of symmetric part.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tlo_wk = tiledlayout(tlo, 1, 2);
tlo_wk.Layout.Tile = 2;
tlo_wk.Layout.TileSpan = [1 2];
tlo_wk.TileSpacing = 'compact';
tlo_wk.Padding = 'compact';

h(2) = nexttile(tlo_wk, 1);

hold on;
[~, symplt] = contourf(wavenums, freqs, qanom_sym_power);

% Contour style
set(symplt, 'LineStyle', 'none');

% Axes range
xlim([-max_wavenum+1, max_wavenum]);
ylim([0, max_freq]);

% Example periods
per3_plt = plot(wavenums, per3);
set(per3_plt, 'LineStyle', '--',...
    'Color', 'k',...
    'Marker', 'none');

per6_plt = plot(wavenums, per6);
set(per6_plt, 'LineStyle', '--',...
    'Color', 'k');

per30_plt = plot(wavenums, per30);
set(per30_plt, 'LineStyle', '--',...
    'Color', 'k');

% Vertical lines
wavenum_0_plt = plot(0.*freqs, freqs);
set(wavenum_0_plt, 'LineStyle', '--',...
    'Color', 'k');

% Dispersion curves
he12_plt = plot(disp_wavenums, kelv12, ...
    disp_wavenums, rossby12, ...
    disp_wavenums, inertio12);
set(he12_plt, 'LineStyle', '-',...
    'Color', 'k',...
    'Marker', 'none');

he25_plt = plot(disp_wavenums, kelv25, ...
    disp_wavenums, rossby25, ...
    disp_wavenums, inertio25);
set(he25_plt, 'LineStyle', '-',...
    'Color', 'k');

he50_plt = plot(disp_wavenums, kelv50, ...
    disp_wavenums, rossby50, ...
    disp_wavenums, inertio50);
set(he50_plt, 'LineStyle', '-',...
    'Color', 'k');

% Text labels
if labels
    text(max_wavenum-0.2, 1/3, '3 Days',...
        'FontSize', 9,...
        'HorizontalAlignment', 'right',...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(max_wavenum-0.2, 1/6, '6 Days',...
        'FontSize', 9,...
        'HorizontalAlignment', 'right',...
        'Color', 'black',...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(max_wavenum-0.2, 1/30, '30 Days',...
        'FontSize', 9,...
        'HorizontalAlignment', 'right',...
        'Color', 'black',...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    labels_12_x = [disp_wavenums(30) - 0.4, disp_wavenums(4),   disp_wavenums(12)];
    labels_12_y = [kelv12(30),              rossby12(4) - 0.01, inertio12(12)    ];
    
    text(labels_12_x, labels_12_y, '12 m',...
        'FontSize', 9,...
        'HorizontalAlignment', 'center',...
        'BackgroundColor', 'white',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    labels_25_x = [disp_wavenums(30) - 0.4, disp_wavenums(4), disp_wavenums(12)];
    labels_25_y = [kelv25(30),              rossby25(4),      inertio25(12)    ];
    
    text(labels_25_x, labels_25_y, '25 m',...
        'FontSize', 9,...
        'HorizontalAlignment', 'center',...
        'BackgroundColor', 'white',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    labels_50_x = [disp_wavenums(30) - 0.4, disp_wavenums(4)  ];
    labels_50_y = [kelv50(30),              rossby50(4) + 0.01];
    
    text(labels_50_x, labels_50_y, '50 m',...
        'FontSize', 9,...
        'HorizontalAlignment', 'center',...
        'BackgroundColor', 'white',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(-5.25, 0.09, 'ER',...
        'FontSize', 9,...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(5.5, 0.1, 'Kelvin',...
        'FontSize', 9,...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(5, 0.45, 'EIG',...
        'FontSize', 9,...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(-8, 0.425, 'WIG',...
        'FontSize', 9,...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
end

title('Symmetric Part');

ax = gca;
set(ax, 'Layer', 'top');
xticks(-max_wavenum+2:2:max_wavenum);
xlabels = string(ax.XAxis.TickLabels);
ax.XAxis.TickLabels = xlabels; % set
xtickangle(0)

hold off;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Right tile - Wheeler-Kiladis Diagram of asymmetric part.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

h(3) = nexttile(tlo_wk, 2);

hold on;

[~, asymplt] = contourf(wavenums, freqs, qanom_asym_power);

% Contour style
set(asymplt, 'LineStyle', 'none');

% Axes range
xlim([-max_wavenum+1, max_wavenum]);
ylim([0, max_freq]);

% Example periods
per3_plt = plot(wavenums, per3);
set(per3_plt, 'LineStyle', '--',...
    'Color', 'k');

per6_plt = plot(wavenums, per6);
set(per6_plt, 'LineStyle', '--',...
    'Color', 'k');

per30_plt = plot(wavenums, per30);
set(per30_plt, 'LineStyle', '--',...
    'Color', 'k');

% Vertical lines
wavenum_0_plt = plot(0.*freqs, freqs);
set(wavenum_0_plt, 'LineStyle', '--',...
    'Color', 'k');


% Dispersion curves
he12_plt = plot(disp_wavenums, rossbygrav12, ...
    disp_wavenums, rossby12, ...
    disp_wavenums, inertio12);
set(he12_plt, 'LineStyle', '-',...
    'Color', 'k');

he25_plt = plot(disp_wavenums, rossbygrav25, ...
    disp_wavenums, rossby25, ...
    disp_wavenums, inertio25);
set(he25_plt, 'LineStyle', '-',...
    'Color', 'k');

he50_plt = plot(disp_wavenums, rossbygrav50, ...
    disp_wavenums, rossby50, ...
    disp_wavenums, inertio50);
set(he50_plt, 'LineStyle', '-',...
    'Color', 'k');

% Text labels
if labels
    text(max_wavenum-0.2, 1/3, '3 Days',...
        'FontSize', 9,...
        'HorizontalAlignment', 'right',...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(max_wavenum-0.2, 1/6, '6 Days',...
        'FontSize', 9,...
        'HorizontalAlignment', 'right',...
        'Color', 'black',...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(max_wavenum-0.2, 1/30, '30 Days',...
        'FontSize', 9,...
        'HorizontalAlignment', 'right',...
        'Color', 'black',...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    labels_12_x = [disp_wavenums(18), disp_wavenums(4),   disp_wavenums(12)];
    labels_12_y = [rossbygrav12(18),  rossby12(4) - 0.01, inertio12(12)    ];
    
    text(labels_12_x, labels_12_y, '12 m',...
        'FontSize', 9,...
        'HorizontalAlignment', 'center',...
        'BackgroundColor', 'white',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    labels_25_x = [disp_wavenums(18), disp_wavenums(4), disp_wavenums(12)];
    labels_25_y = [rossbygrav25(18),  rossby25(4),      inertio25(12)    ];
    
    text(labels_25_x, labels_25_y, '25 m',...
        'FontSize', 9,...
        'HorizontalAlignment', 'center',...
        'BackgroundColor', 'white',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    labels_50_x = [disp_wavenums(18), disp_wavenums(4)  ];
    labels_50_y = [rossbygrav50(18),  rossby50(4) + 0.01];
    
    text(labels_50_x, labels_50_y, '50 m',...
        'FontSize', 9,...
        'HorizontalAlignment', 'center',...
        'BackgroundColor', 'white',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(-5.25, 0.09, 'ER',...
        'FontSize', 9,...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(5.5, 0.25, 'MRG',...
        'FontSize', 9,...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(5, 0.45, 'EIG',...
        'FontSize', 9,...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
    
    text(-8, 0.425, 'WIG',...
        'FontSize', 9,...
        'BackgroundColor', 'white',...
        'EdgeColor', 'black',...
        'Color', 'black',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 1);
end

title('Asymmetric Part');

ax = gca;
set(ax, 'Layer', 'top');
xticks(-max_wavenum+2:2:max_wavenum);
xlabels = string(ax.XAxis.TickLabels);
ax.XAxis.TickLabels = xlabels; % set
xtickangle(0)

hold off;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% General figure styling.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Colorbar
set(h(2:3), 'Colormap', spring, ...
    'CLim', [power_min, power_max])
cbh = colorbar(h(3));
cbh.Layout.Tile = 'east';

set(cbh, 'LineWidth', 0.25,...
    'Color', 'black');
set(cbh.Label, ...
    'String', 'Log(Moisture (g kg^{-1}) Power)',...
    'Color', 'black');


% Axes labels and title

titleStr = 'Wheeler-Kiladis Diagram';
title(tlo_wk, titleStr, ...
    'Color', 'black');

xlabel(tlo_wk, 'Zonal Wavenumber', ...
    'Color', 'black');
ylabel(tlo_wk, 'Frequency (d^{-1})', ...
    'Color', 'black');


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Save plot.
alt_str = sprintf(['%3.1f'], abs(altW_true));

%~ Figure size
set(gcf, 'Units', 'inches');
figWidth  = 9; % Figure width in inches.
figHeight = 6.5; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight],...
    'PaperOrientation', 'portrait');


file_name = strcat([params.component_name, '_hovmoller_' , q_mode, '_wk_', ...
    alt_str, '.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters', '-fillpage');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end
