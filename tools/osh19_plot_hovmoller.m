function osh19_plot_hovmoller(params, q_mode)

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

n_outfiles = floor(sim_days/out_freq);
out_idxs  = 0:n_outfiles;

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    t(1, out_idx + 1) = ncread(state_file, 't');
    q_temp     = ncread(state_file, 'q');
    theta_temp = ncread(state_file, 'theta');
    
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
title(title_str);

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
text(10, 0.02*sim_days, '3.1 m s^{-1}',...
    'Color', 'black',...
    'FontWeight', 'bold',...
    'LineStyle', '-',...
    'LineWidth', 0.25,...
    'Margin', 0.25);

speed_4 = 4.6 * (1 / 1000) * (360 / (2 * pi * P_E)) * days_to_secs;
plot(lons, (1./speed_4).*lons + sim_days/2, 'k--');
text(10, sim_days/2+(0.02*sim_days), '4.6 m s^{-1}',...
    'Color', 'black',...
    'FontWeight', 'bold',...
    'LineStyle', '-',...
    'LineWidth', 0.25,...
    'Margin', 0.25);

hold off;


%~ Figure size
set(gcf, 'Units', 'inches');
figWidth  = 4.5; % Figure width in inches.
figHeight = 8; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight],...
    'PaperOrientation', 'portrait');


file_name = strcat([params.component_name, '_hovmoller_' , q_mode, '.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters', '-fillpage');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end

