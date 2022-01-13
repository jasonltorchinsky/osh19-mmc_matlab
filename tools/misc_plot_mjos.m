function misc_plot_mjos(eof_params, mjoo_params, alt, lat)

% Get path to output data.
out_path       = eof_params.out_path;
exp_path       = fullfile(out_path, eof_params.exp_name);
truth_path     = fullfile(exp_path, eof_params.component_name);
mjoo_path      = fullfile(exp_path, mjoo_params.component_name);
eof_path       = fullfile(exp_path, 'eofs');

addpath(out_path);
addpath(exp_path);
addpath(truth_path);
addpath(eof_path);

% Create directory to hold plots
plot_path = fullfile(exp_path, 'plots');
if ~(isfolder(plot_path))
    mkdir(plot_path);
end

addpath(plot_path);

% Read in run parameters
truth_params_file = fullfile(truth_path, 'params.nc');

nx          = ncread(truth_params_file, 'nx');
P_E         = ncread(truth_params_file, 'P_E'); % Convert to meters
sim_days    = ncread(truth_params_file, 'sim_days');
out_freq    = ncread(truth_params_file, 'out_freq');

% Read in grids
truth_grid_file = fullfile(truth_path, 'grid.nc');

xx  = ncread(truth_grid_file, 'xx');
yy  = ncread(truth_grid_file, 'yy');
zzW = ncread(truth_grid_file, 'zzW');

lons = 360 * (xx - xx(1)) / (xx(end) - xx(1));
lats = yy / 110.567;

% Calculate indices for desired altitude, latitude
[~, lat_idx]  = min(abs(lats-lat));
[~, altW_idx] = min(abs(zzW-alt));

lat_true  = lats(lat_idx);
altW_true = zzW(altW_idx);

% True moisture time-series
q = zeros(nx, floor(sim_days/out_freq) + 1);
t = zeros(1, floor(sim_days/out_freq) + 1);

n_outfiles = floor(sim_days/out_freq);
out_idxs  = 0:n_outfiles;

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file = fullfile(truth_path, state_file_name);
    
    t(1, out_idx + 1) = ncread(state_file, 't');
    q_temp     = ncread(state_file, 'q');
    
    q(:, out_idx + 1) = q_temp(lat_idx, :, altW_idx);
end

days_to_secs = 3600*24;
t = t / days_to_secs;

% Get moisture parts of EOFs at altitude and latitude
eof_file_name = 'eofs.nc';
eof_file = fullfile(eof_path, eof_file_name);
q_mjo1 = ncread(eof_file, 'q_mjo1');
q_mjo1 = squeeze(q_mjo1(lat_idx, :, altW_idx));
q_mjo2 = ncread(eof_file, 'q_mjo2');
q_mjo2 = squeeze(q_mjo2(lat_idx, :, altW_idx));

% Reconstruct MJO from expansion coefficients
exp1 = ncread(eof_file, 'exp1');
exp2 = ncread(eof_file, 'exp2');

q_mjo_eof = zeros(nx, floor(sim_days/out_freq) + 1);

for out_idx = out_idxs
    q_mjo_eof(:, out_idx + 1) = exp1(out_idx + 1) * q_mjo1 ...
        + exp2(out_idx + 1) * q_mjo2;
end

% Reconstruct MJO from MJOO model
u_1  = zeros(size(out_idxs));
u_2  = zeros(size(out_idxs));

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'),'.nc']);
    state_file = fullfile(mjoo_path, state_file_name);
    u_1(out_idx + 1) = ncread(state_file, 'u_1');
    u_2(out_idx + 1) = ncread(state_file, 'u_2');
end

q_mjo_mjoo = zeros(nx, floor(sim_days/out_freq) + 1);

for out_idx = out_idxs
    q_mjo_mjoo(:, out_idx + 1) = 2*(u_1(out_idx + 1) * q_mjo1 ...
        + u_2(out_idx + 1) * q_mjo2);
end

% Colorbar limits
max_q = max(abs([q_mjo_eof q_mjo_mjoo]), [], 'all');
min_q = -max_q;
cmap = load('rb.mat').rb;

% Create the three Hovmoller plots
tlo = tiledlayout(1,2);

% % First plot - true moisture
% h(1,1) = nexttile(tlo,1);
% 
% hold on;
% contourf(lons, t, q.', ...
%     'edgecolor', 'none');
% colormap(cmap);
% caxis([min_q, max_q]);
% title('Raw');
% hold off;

% Second plot - MJO moisture, expansion coefficients from EOF decomp.
h(1,1) = nexttile(tlo,1);

hold on;
contourf(lons, t, q_mjo_eof.', ...
    'edgecolor', 'none');
colormap(cmap);
caxis([min_q, max_q]);
title('MJO (EOF)');
hold off;

% Third plot - MJO moisture, expansion coefficients from MJOO model
h(1,2) = nexttile(tlo,2);

hold on;
contourf(lons, t, q_mjo_mjoo.', ...
    'edgecolor', 'none');
colormap(cmap);
caxis([min_q, max_q]);
title('MJO (MJOO)');
hold off;


cb = colorbar();
cb.Label.String = 'Moisture Anomaly (kg kg^{-1})';

% Title
title_base_str = sprintf('Hovmoller Diagram\nAltitude %3.1f km\n',...
    altW_true);
if (lat_true < 0)
    title_str = strcat([title_base_str, sprintf('Latitude %3.1fS',...
        abs(lat_true))]);
elseif (lat_true > 0)
    title_str = strcat([title_base_str, sprintf('Latitude %3.1fN',...
        abs(lat_true))]);
elseif (trueLate == 0)
    title_str = strcat([title_base_str, 'Equator']);
end
    
title_str = "Hovmoller Diagram";
title(tlo, title_str);

% Axis limits
xlim(h, [0, 360]);
ylim(h, [0, sim_days]);

% Axis labels
xlabel(tlo, 'Longitude');
ylabel(tlo, 'Time (d)');

% Ticks
xticks(h, 0:90:360);
xticklabels(h, {'180', '90W', '0', '90E', '180'});

% Example zonal windspeeds
for ii = 1:2
    hold(h(1,ii), 'on');
    speed_3 = 3.1 * (1 / 1000) * (360 / (2 * pi * P_E)) * days_to_secs;
    plot(h(1,ii), lons, (1./speed_3).*lons, 'k-.');
    text(h(1,ii), 10, 0.02*sim_days, '3.1 m s^{-1}',...
        'Color', 'black',...
        'FontWeight', 'bold',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 0.25);
    
    speed_4 = 4.6 * (1 / 1000) * (360 / (2 * pi * P_E)) * days_to_secs;
    plot(h(1,ii), lons, (1./speed_4).*lons + sim_days/2, 'k--');
    text(h(1,ii), 10, sim_days/2+(0.02*sim_days), '4.6 m s^{-1}',...
        'Color', 'black',...
        'FontWeight', 'bold',...
        'LineStyle', '-',...
        'LineWidth', 0.25,...
        'Margin', 0.25);
    
    hold(h(1,ii), 'off');
end

%~ Figure size
set(gcf, 'Units', 'inches');
figWidth  = 6.5; % Figure width in inches.
figHeight = 8; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight],...
    'PaperOrientation', 'portrait');

% Save plot.
if (lat_true < 0)
    lat_str = sprintf('%3.1fS', abs(lat_true));
elseif (lat_true > 0)
    lat_str = sprintf('%3.1fN', abs(lat_true));
else
    lat_str = 'EQ';
end
alt_str = sprintf('%3.1f', abs(altW_true));

file_name = strcat(['misc_hovmoller_', ...
    alt_str, '_', lat_str, '.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters', '-fillpage');

end

