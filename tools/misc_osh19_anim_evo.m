function misc_osh19_anim_evo(params_1, params_2, alt, lat)

% Get path to output data.
out_path_1       = params_1.out_path;
exp_path_1       = fullfile(out_path_1, params_1.exp_name);
component_path_1 = fullfile(exp_path_1, params_1.component_name);

addpath(out_path_1);
addpath(exp_path_1);
addpath(component_path_1);

out_path_2       = params_2.out_path;
exp_path_2       = fullfile(out_path_2, params_2.exp_name);
component_path_2 = fullfile(exp_path_2, params_2.component_name);

addpath(out_path_2);
addpath(exp_path_2);
addpath(component_path_2);

% Create directory to hold plots
plot_path = fullfile(exp_path_1, 'plots');
if ~(isfolder(plot_path))
    mkdir(plot_path);
end

addpath(plot_path);

% Read in run parameters (ASSUME SAME FOR BOTH RUNS)
params_file = fullfile(component_path_1, 'params.nc');

nx       = ncread(params_file, 'nx');
ny       = ncread(params_file, 'ny');
nz       = ncread(params_file, 'nz');

H        = ncread(params_file, 'H');
sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');

% Read in grids (ASSUME SAME FOR BOTH RUNS)
grid_file = fullfile(component_path_1, 'grid.nc');

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


out_idxs = 0:out_freq:floor(sim_days/out_freq);
%out_idxs = 0:2*out_freq:100;

% Make animation of system output
file_name = 'misc_osh19_comp_evo.avi';
plot_file = fullfile(plot_path, file_name);
v = VideoWriter(plot_file);
v.FrameRate = 9;
open(v);

tlo = tiledlayout(2, 1, 'TileSpacing', 'Compact');
h(1) = nexttile(tlo,1);
h(2) = nexttile(tlo, 2);
set(h(:), 'nextplot', 'replacechildren');
pbaspect(h(1), [4 1 1]);
pbaspect(h(2), [4 1 1]);

% Set things common to all plots
days_to_secs = 3600*24;

% Titles
if lat_true < 0
    lat_str = sprintf('%.2fS', ...
        round(abs(lat_true), 2));
elseif lat_true > 0
    lat_str = sprintf('%.2fN', ...
        round(abs(lat_true), 2));
else
    lat_str = sprintf('EQ');
end

title_str = sprintf(['State Evolution at ' lat_str ' and %.2f km'], ...
    round(altW_true, 2));
title(tlo, title_str)

% Axis limits
xlim(h, [0, 360]);
ylim(h(1), [-45, 45]);
ylim(h(2), [0, H]);


% Axis labels
xlabel(h(:), 'Longitude');
ylabel(h(1), 'Latitude');
ylabel(h(2), 'Altitude (km)');

% Ticks
xticks(h, 0:45:360);
xticklabels(h(2), {'180', '135W', '90W', '45W', '0', '45E', '90E', '135E', '180'});

yticks(h(1), -30:30:30);
yticks(h(2),   0:4:16);

yticklabels(h(1), {'30S', 'EQ', '30N'});

% Get moisture limits by chekcing moisture at final day
state_file_name = strcat(['state_', num2str(out_idxs(end) ,'%04u'), '.nc']);
state_file_1 = fullfile(component_path_1, state_file_name);
state_file_2 = fullfile(component_path_2, state_file_name);
q_max = max(abs(ncread(state_file_2, 'q') - ncread(state_file_1, 'q')), [], 'all');
    
% Colors
cmap = load('rb.mat').rb;
set(h, 'Colormap', cmap, ...
    'CLim', [-1, 1] * q_max)
cbh = colorbar(h(1));
cbh.Layout.Tile = 'east';
cbh.Label.String = 'Moisture Anomaly (kg kg^{-1})';

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file_1 = fullfile(component_path_1, state_file_name);
    state_file_2 = fullfile(component_path_2, state_file_name);
    
    % Get winds, moisture
    t = ncread(state_file_1, 't');
    u_horz_1 = ncread(state_file_1, 'u', [1 1 altU_idx], [ny nx 1], [1 1 1]);
    v_horz_1 = ncread(state_file_1, 'v', [1 1 altU_idx], [ny nx 1], [1 1 1]);
    q_horz_1 = ncread(state_file_1, 'q', [1 1 altW_idx], [ny nx 1], [1 1 1]);
    
    u_horz_2 = ncread(state_file_2, 'u', [1 1 altU_idx], [ny nx 1], [1 1 1]);
    v_horz_2 = ncread(state_file_2, 'v', [1 1 altU_idx], [ny nx 1], [1 1 1]);
    q_horz_2 = ncread(state_file_2, 'q', [1 1 altW_idx], [ny nx 1], [1 1 1]);
    
    t_str = sprintf('Day %d', round(t/days_to_secs, 0));
    
    
    cla(h(1));
    
    hold(h(1), 'on');
    contourf(h(1), lons, lats, (q_horz_2 - q_horz_1), 32, ...
        'edgecolor', 'none');
    quiver(h(1), lons, lats, (u_horz_2 - u_horz_1)', (v_horz_2 - v_horz_1)', ...
        'k');
    text(h(1), 4.5, 0.9*90 - 45, t_str, ...
        'Color', 'black', ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'black', ...
        'LineStyle', '-', ...
        'LineWidth', 0.25,...
        'Margin', 0.25);
    hold(h(1), 'off');
    
    
    % Vertical profile at desired latitude
    u_vert_1 = squeeze(ncread(state_file_1, 'u', [lat_idx 1 1], [1 nx nz+1], [1 1 1]));
    w_vert_1 = squeeze(ncread(state_file_1, 'w', [lat_idx 1 1], [1 nx nz+1], [1 1 1]));
    q_vert_1 = squeeze(ncread(state_file_1, 'q', [lat_idx 1 1], [1 nx nz+1], [1 1 1]));
    
    u_vert_2 = squeeze(ncread(state_file_2, 'u', [lat_idx 1 1], [1 nx nz+1], [1 1 1]));
    w_vert_2 = squeeze(ncread(state_file_2, 'w', [lat_idx 1 1], [1 nx nz+1], [1 1 1]));
    q_vert_2 = squeeze(ncread(state_file_2, 'q', [lat_idx 1 1], [1 nx nz+1], [1 1 1]));
    
    cla(h(2));
    
    hold(h(2), 'on');
    contourf(h(2), lons, zzW, (q_vert_2 - q_vert_1)', 32, ...
       'edgecolor', 'none');
    quiver(h(2), lons, zzU, (u_vert_2 - u_vert_1)', 150*(w_vert_2 - w_vert_1)', ...
       'k');
    text(h(2), 4.5, 0.9*H, t_str, ...
       'Color', 'black', ...
       'BackgroundColor', 'white', ...
       'EdgeColor', 'black', ...
       'LineStyle', '-', ...
       'LineWidth', 0.25,...
       'Margin', 0.25);
    hold(h(2), 'off')
    
    
    frame = getframe(gcf);
    writeVideo(v, frame);
end


close(v);

end

