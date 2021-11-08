function osh19_plot_hovmoller(params, alt, lat)

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
P_E         = ncread(params_file, 'P_E'); % Convert to meters
sim_days    = ncread(params_file, 'sim_days');
out_freq    = ncread(params_file, 'out_freq');

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

xx  = ncread(grid_file, 'xx');
yy  = ncread(grid_file, 'yy');
zzW = ncread(grid_file, 'zzW');

lons = 360 * (xx - xx(1)) / (xx(end) - xx(1));
lats = yy / 110.567;

% Calculate indices for desired altitude, latitude
[~, lat_idx]  = min(abs(lats-lat));
[~, altW_idx] = min(abs(zzW-alt));

lat_true  = lats(lat_idx);
altW_true = zzW(altW_idx);

% Get moisture, potential temperature anomaly at desired altitude and latitude
q     = zeros(nx, floor(sim_days/out_freq) + 1);
theta = zeros(nx, floor(sim_days/out_freq) + 1);
t     = zeros(1, floor(sim_days/out_freq) + 1);

n_outfiles = floor(sim_days/outfreq);
out_idxs  = 0:n_outfiles;

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    t(1, out_idx + 1) = ncread(state_file, 't');
    q_temp     = ncread(state_file, 'q');
    theta_temp = ncread(state_file, 'theta');
    
    q(:, out_idx + 1) = q_temp(lat_idx, :, altW_idx);
    theta(:, out_idx + 1) = theta_temp(lat_idx, :, altW_idx);
end

days_to_secs = 3600*24;
t = t / days_to_secs;

% Create plot

hold on;
[~, q_plt] = contourf(lons, t, q.', ...
    'edgecolor', 'none');

% Colorbar
max_q = max(abs(q), [], 'all');
min_q = -max_q;
cmap = load('rb.mat').rb;
colormap(cmap);
caxis([min_q, max_q]);

cb = colorbar();
cb.Label.String = 'Moisture Anomaly (kg kg^{-1})';


% Potential temperature anomaly
pos_theta = max(theta, 0);
neg_theta = min(theta, 0);
[~, pos_th_plt] = contour(lons, out_idxs, pos_theta.', ...
    'k-');
[~, neg_th_plt] = contour(lons, out_idxs, neg_theta.', ...
    'k--');

% Contour levels
max_theta = max(abs(theta));
th_levs = 0.1:0.2:0.9 * max_theta;

pos_th_plt.LevelListMode = 'manual';
pos_th_plt.LevelList     = th_levs;

neg_th_plt.LevelListMode = 'manual';
neg_th_plt.LevelList     = th_levs;

% Title
title_base_str = sprintf(['Hovmoller Diagram\nAltitude %3.1f km\n'],...
    altW_true);
if (lat_true < 0)
    title_str = strcat([title_base_str, sprintf(['Latitude %3.1fS'],...
        abs(lat_true))]);
elseif (lat_true > 0)
    title_str = strcat([title_base_str, sprintf(['Latitude %3.1fN'],...
        abs(lat_true))]);
elseif (trueLate == 0)
    title_str = strcat([title_base_str, 'Equator']);
end
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

% Save plot.
if (lat_true < 0)
    lat_str = sprintf(['%3.1fS'], abs(lat_true));
elseif (lat_true > 0)
    lat_str = sprintf(['%3.1fN'], abs(lat_true));
else
    lat_str = 'EQ';
end
alt_str = sprintf(['%3.1f'], abs(altW_true));

file_name = strcat([params.component_name, '_hovmoller_', ...
    alt_str, '_', lat_str, '.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters', '-fillpage');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end

