function eof_plot_eofs(eof_params, alt)

% Get path to output data.
out_path       = eof_params.out_path;
exp_path       = fullfile(out_path, eof_params.exp_name);
component_path = fullfile(exp_path, eof_params.component_name);
eof_path       = fullfile(exp_path, 'eofs');

addpath(out_path);
addpath(exp_path);
addpath(component_path);
addpath(eof_path);

% Create directory to hold plots
plot_path = fullfile(exp_path, 'plots');
if ~(isfolder(plot_path))
    mkdir(plot_path);
end

addpath(plot_path);

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

xx  = ncread(grid_file, 'xx');
yy  = ncread(grid_file, 'yy');
zzU = ncread(grid_file, 'zzU');
zzW = ncread(grid_file, 'zzW');

lons = 360 * (xx - xx(1)) / (xx(end) - xx(1));
lats = yy / 110.567;

% Calculate indices for desired altitude, latitude
[~, altU_idx] = min(abs(zzU-alt));
[~, altW_idx] = min(abs(zzW-alt));

altU_true = zzU(altU_idx);
altW_true = zzW(altW_idx);

% Read in EOFs, reconstructed MJOs
eof_file_name = 'eofs.nc';
eof_file = fullfile(eof_path, eof_file_name);

u_eof1 = ncread(eof_file, 'u_eof1');
u_eof2 = ncread(eof_file, 'u_eof2');

q_eof1 = ncread(eof_file, 'q_eof1');
q_eof2 = ncread(eof_file, 'q_eof2');

u_mjo1 = ncread(eof_file, 'u_mjo1');
u_mjo2 = ncread(eof_file, 'u_mjo2');

q_mjo1 = ncread(eof_file, 'q_mjo1');
q_mjo2 = ncread(eof_file, 'q_mjo2');

exp1 = ncread(eof_file, 'exp1');
exp2 = ncread(eof_file, 'exp2');

% Normalize EOFs, MJOs to have physical units ~ SCALE BY VARIANCE OR SOMETHING
u_eof1 = u_eof1 * max(abs(exp1), [], 'all');
u_eof2 = u_eof2 * max(abs(exp2), [], 'all');

q_eof1 = q_eof1 * max(abs(exp1), [], 'all');
q_eof2 = q_eof2 * max(abs(exp2), [], 'all');

u_mjo1 = u_mjo1 * max(abs(exp1), [], 'all');
u_mjo2 = u_mjo2 * max(abs(exp2), [], 'all');

q_mjo1 = q_mjo1 * max(abs(exp1), [], 'all');
q_mjo2 = q_mjo2 * max(abs(exp2), [], 'all');


% Create tiled layout 
tlo = tiledlayout(2,3); % Outer layout

tlo_eof = tiledlayout(tlo, 2, 1);
tlo_eof.Layout.Tile = 1;
tlo_eof.Layout.TileSpan = [2, 1];

tlo_mjo = tiledlayout(tlo, 2, 1);
tlo_mjo.Layout.Tile = 2;
tlo_mjo.Layout.TileSpan = [2, 2];

% Create plots of EOFs

% Zonal wind speed
h(1,1) = nexttile(tlo_eof);

hold on;

plot(lons, u_eof1, 'k-', ...
    'DisplayName', 'EOF 1');
plot(lons, u_eof2, 'k--', ...
    'DisplayName', 'EOF 2');

legend();

ylabel('Zonal wind speed (m s^{-1})');

hold off;

% Moisture anomaly
h(2,1) = nexttile(tlo_eof);

hold on;

plot(lons, q_eof1, 'k-', ...
    'DisplayName', 'EOF 1');
plot(lons, q_eof2, 'k--', ...
    'DisplayName', 'EOF 2');

legend();

ylabel('Moisture Anomaly (kg kg^{-1})');

hold off;

% Options for EOFs

% Axis labels
xlabel(tlo_eof, 'Longitude');

% Create plots of MJOs

% MJO 1
h(1,2) = nexttile(tlo_mjo);

hold on;

u_mjo1_horz = squeeze(u_mjo1(:, :, altU_idx));
q_mjo1_horz = squeeze(q_mjo1(:, :, altW_idx));

contourf(lons, lats, q_mjo1_horz, ...
        'edgecolor', 'none');
quiver(lons, lats, u_mjo1_horz, 0*u_mjo1_horz, ...
        'k', ...
        'Autoscale', 'off');

title_str = sprintf(['MJO 1 - Altitude %.2f (km)'], round(altW_true, 2));

title(title_str);


hold off;

h(2,2) = nexttile(tlo_mjo);

hold on;

u_mjo2_horz = squeeze(u_mjo2(:, :, altU_idx));
q_mjo2_horz = squeeze(q_mjo2(:, :, altW_idx));

contourf(lons, lats, q_mjo2_horz, ...
        'edgecolor', 'none');
quiver(lons, lats, u_mjo2_horz, 0*u_mjo2_horz, ...
        'k', ...
        'Autoscale', 'off');

title_str = sprintf(['MJO 2 - Altitude %.2f (km)'], round(altW_true, 2));

title(title_str);
    
ylim([-45, 45]);
yticks(-30:30:30);

hold off;

% Options for MJOs
xlabel(tlo_mjo, 'Longitude');
ylabel(h(:,2), 'Latitude');
ylim(h(:,2), [-45, 45]);
yticks(h(:,2), -30:30:30);
yticklabels(h(:,2), {'30S', 'EQ', '30N'});

% Colorbar
max_q = max(abs([q_mjo1_horz, q_mjo2_horz]), [], 'all');
min_q = -max_q;
cmap = load('rb.mat').rb;
set(h(:,2), ...
    'Colormap', cmap, ...
    'CLim', [min_q, max_q])

cbh = colorbar(h(end));
cbh.Layout.Tile = 'east';
cbh.Label.String = 'Moisture Anomaly (kg kg^{-1})';

% Options common to all plots
xlim(h, [0, 360]);
xticks(h, 0:45:360);
xticklabels(h, {'180', '135W', '90W', '45W', '0', '45E', '90E', '135E', '180'});


%~ Figure size
set(gcf, 'Units', 'inches');
figWidth  = 11; % Figure width in inches.
figHeight = 8.5; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight],...
    'PaperOrientation', 'portrait');

% Save plot.
file_name = strcat([eof_params.component_name, '_eofs.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters', '-fillpage');

end

