function osh19_plot_evo(params, alt, lat)

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

H        = ncread(params_file, 'H');
sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

xx  = ncread(grid_file, 'xx');
yy  = ncread(grid_file, 'yy');
zzU = ncread(grid_file, 'zzU');
zzW = ncread(grid_file, 'zzW');

lons = 360 * (xx - xx(1)) / (xx(end) - xx(1));
lats = yy / 110.567;

% Calculate indices for desired altitude, latitude
[~, lat_idx]  = min(abs(lats-lat));
[~, altU_idx] = min(abs(zzU-alt));
[~, altW_idx] = min(abs(zzW-alt));

lat_true  = lats(lat_idx);
altU_true = zzU(altU_idx);
altW_true = zzW(altW_idx);

% Store maximal theta
max_theta = 0;
max_q = 0;

% Calculate output indices for times
% 0*sim_days, 1/3*sim_days, 2/3*sim_days, sim_days
out_idxs = zeros(1,4);
out_idxs(1) = 0;
out_idxs(2) = floor((1.0 / 3.0) * sim_days / out_freq);
out_idxs(3) = floor((2.0 / 3.0) * sim_days / out_freq);
out_idxs(4) = floor(sim_days / out_freq);

% Create tiled layout of plots of evolution
tlo = tiledlayout(4,2); % Outer layout

tlo_horz = tiledlayout(tlo, 4, 1); % Horizontal profiles
tlo_horz.Layout.Tile = 1;
tlo_horz.Layout.TileSpan = [4, 1];

tlo_vert = tiledlayout(tlo, 4, 1); % Vertical profiles
tlo_vert.Layout.Tile = 2;
tlo_vert.Layout.TileSpan = [4, 1];

tile_idx = 1;

for out_idx = out_idxs
    h(1, tile_idx) = nexttile(tlo_horz);
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    % Get winds, moisture
    t = ncread(state_file, 't');
    u = ncread(state_file, 'u');
    v = ncread(state_file, 'v');
    w = 150*ncread(state_file, 'w');
    theta = ncread(state_file, 'theta');
    q = ncread(state_file, 'q');
    
    
    days_to_secs = 3600*24;
    t_str = sprintf(['Day %d'], round(t/days_to_secs, 0));
    
    % Horizontal profile at desired altitude
    u_horz = squeeze(u(:, :, altU_idx));
    v_horz = squeeze(v(:, :, altU_idx));
    q_horz = squeeze(q(:, :, altW_idx));
    theta_horz = squeeze(theta(:, :, altW_idx));
    pos_theta_horz = max(theta_horz, 0);
    neg_theta_horz = min(theta_horz, 0);
    
    hold on;
    contourf(lons, lats, q_horz, ...
        'edgecolor', 'none');
    [~, horz_th_plt(1,tile_idx)] = contour(lons, lats, pos_theta_horz, ...
        'k-');
    [~, horz_th_plt(2,tile_idx)] = contour(lons, lats, neg_theta_horz, ...
        'k--');
    quiver(lons, lats, u_horz', v_horz, ...
        'k', ...
        'Autoscale', 'off');
    text(4.5, 0.9*90 - 45, t_str, ...
        'Color', 'black', ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'black', ...
        'LineStyle', '-', ...
        'LineWidth', 0.25,...
        'Margin', 0.25);
    hold off;
    
    h(2, tile_idx) = nexttile(tlo_vert);
    
    % Vertical profile at desired latitude
    u_vert = squeeze(u(lat_idx, :, :));
    w_vert = squeeze(w(lat_idx, :, :));
    q_vert = squeeze(q(lat_idx, :, :));
    theta_vert = squeeze(theta(lat_idx, :, :));
    pos_theta_vert = max(theta_vert, 0);
    neg_theta_vert = min(theta_vert, 0);
    
    hold on;
    contourf(lons, zzW, q_vert', ...
        'edgecolor', 'none');
    [~, vert_th_plt(1,tile_idx)] = contour(lons, zzW, pos_theta_vert', ...
        'k-');
    [~, vert_th_plt(2,tile_idx)] = contour(lons, zzW, neg_theta_vert', ...
        'k--');
    quiver(lons, zzU, u_vert', w_vert', ...
        'k', ...
        'Autoscale', 'off');
    text(4.5, 0.9*H, t_str, ...
        'Color', 'black', ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'black', ...
        'LineStyle', '-', ...
        'LineWidth', 0.25,...
        'Margin', 0.25);
    hold off;
    
    tile_idx = tile_idx + 1;
    
    % Update maximum theta
    max_theta_horz = max(abs(theta_horz), [], 'all');
    max_theta_vert = max(abs(theta_vert), [], 'all');
    max_theta = max([max_theta, max_theta_horz, max_theta_vert]);
    
    % Update maximum q
    max_q_horz = max(abs(q_horz), [], 'all');
    max_q_vert = max(abs(q_vert), [], 'all');
    max_q = max([max_q, max_q_horz, max_q_vert]);
end

% Colorbar
min_q = -max_q;
cmap = load('rb.mat').rb;
set(h(:,:), ...
    'Colormap', cmap, ...
    'CLim', [min_q, max_q])

cbh = colorbar(h(end));
cbh.Layout.Tile = 'east';
cbh.Label.String = 'Moisture Anomaly (kg kg^{-1})';

% Contour levels for theta plots
th_levs = 0.1:0.2:0.9 * max_theta;
for th_plt = [horz_th_plt(1,:) vert_th_plt(1,:)]
    th_plt.LevelListMode = 'manual';
    th_plt.LevelList = th_levs;
end

for th_plt = [horz_th_plt(2,:) vert_th_plt(2,:)]
    th_plt.LevelListMode = 'manual';
    th_plt.LevelList = -th_levs;
end

% Titles
if lat_true < 0
    lat_str = sprintf(['%.2fS'], ...
        round(abs(lat_true), 2));
elseif lat_true > 0
    lat_str = sprintf(['%.2fN'], ...
        round(abs(lat_true), 2));
else
    lat_str = sprintf(['EQ']);
end

title_str = sprintf(['State Evolution at ' lat_str ' and %.2f km'], ...
    round(altW_true, 2));
title(tlo, title_str)

% Axis limits
xlim(h, [0, 360]);

ylim(h(1,:), [-45, 45]);
ylim(h(2,:), [0, H]);

% Axis labels
xlabel(tlo_horz, 'Longitude');
xlabel(tlo_vert, 'Longitude');

ylabel(tlo_horz, 'Latitude');
ylabel(tlo_vert, 'Altitude (km)');

% Ticks
xticks(h, 0:45:360);
xticklabels(h, {'180', '135W', '90W', '45W', '0', '45E', '90E', '135E', '180'});

yticks(h(1,:), -30:30:30);
yticklabels(h(1,:), {'30S', 'EQ', '30N'});

yticks(h(2,:), 0:4:16);


%~ Figure size
set(gcf, 'Units', 'inches');
figWidth  = 11; % Figure width in inches.
figHeight = 8.5; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight],...
    'PaperOrientation', 'portrait');

% Save plot.
file_name = strcat([params.component_name, '_osh19_evo.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters', '-fillpage');

end

