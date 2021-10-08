function osh19_plot_evo(params)

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
sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');
xx  = ncread(grid_file, 'xx');
yy  = ncread(grid_file, 'yy');
zzU = ncread(grid_file, 'zzU');
zzW = ncread(grid_file, 'zzW');

lons = 360 * (xx - xx(1) / (xx(end) - xx(1));
lats = yy / 110.567;

% Calculate indices for desired altitude, latitude
lat = 0.0;
alt = 9.0;

[lat_true, lat_idx]   = min(lats);
[altU_true, altU_idx] = min(zzU);
[altW_true, altW_idx] = min(zzW);

% Calculate output indices for times
% 0*sim_days, 1/3*sim_days, 2/3*sim_days, sim_days
out_idxs = zeros(1,4);
out_idxs(1) = 0;
out_idxs(2) = floor((1.0 / 3.0) * sim_days / out_freq);
out_idxs(3) = floor((2.0 / 3.0) * sim_days / out_freq);
out_idxs(4) = floor(sim_days / out_freq);

% Create tiled layout of plots of evolution
tiledlayout(4,2);

for out_idx = out_idxs
    nexttile;
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    % Get winds, moisture
    u = ncread(state_file, 'u');
    v = ncread(state_file, 'v');
    w = ncread(state_file, 'w');
    q = ncread(state_file, 'q');
    
    % Horizontal profile at desired altitude
    u_horz = u(:, :, altU_idx);
    v_horz = v(:, :, altU_idx);
    q_horz = u(:, :, altW_idx);
    
    hold on;
    quiver(lons, lats, u_horz, v_horz);
    contourf(lons, lats, q_horz);
    hold off;
    
    nexttile;
    
    % Vertical profile at desired altitude
    u_vert = u(lat_idx, :, :);
    w_vert = w(lat_idx, :, :);
    q_vert = q(lat_idx, :, :);
    
    hold on;
    quiver(lons, zzU, u_vert, w_vert);
    contourf(lons, zzW, q_vert);
    hold off;
    
end

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

