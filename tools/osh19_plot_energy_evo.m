function osh19_plot_energy_evo(params, alt)

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
max_wavenum = 10;

wavenums    = -nx/2:nx/2-1;

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

% Get power
qanom_sym_fft   = zeros(size(qanom_sym));
qanom_asym_fft  = zeros(size(qanom_asym));
for day = days
    qanom_sym_fft(:,day+1)  = fftshift(fft(qanom_sym(:,day+1)));
    qanom_asym_fft(:,day+1) = fftshift(fft(qanom_asym(:,day+1)));
end

qanom_sym_power   = abs(qanom_sym_fft).^2;
qanom_asym_power  = abs(qanom_asym_fft).^2;

power_min = min([qanom_sym_power, qanom_asym_power], [], 'all');
power_max = max([qanom_sym_power, qanom_asym_power], [], 'all');


% Create master tiled layout
tlo = tiledlayout(1, 2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Left tile - Wheeler-Kiladis Diagram of symmetric part.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

h(1) = nexttile(tlo, 1);

hold on;
[~, symplt] = contourf(wavenums, days, qanom_sym_power');

% Contour style
set(symplt, 'LineStyle', 'none');

% Vertical lines
wavenum_0_plt = plot(0.*wavenums, sim_days.*wavenums);
set(wavenum_0_plt, 'LineStyle', '--',...
    'Color', 'k');

% Axes range
xlim([-max_wavenum+1, max_wavenum]);
ylim([0, sim_days]);


title('Symmetric Part');

hold off;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Right tile - Wheeler-Kiladis Diagram of asymmetric part.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

h(2) = nexttile(tlo, 2);

hold on;

[~, asymplt] = contourf(wavenums, days, qanom_asym_power');

% Contour style
set(asymplt, 'LineStyle', 'none');

% Vertical lines
wavenum_0_plt = plot(0.*wavenums, sim_days.*wavenums);
set(wavenum_0_plt, 'LineStyle', '--',...
    'Color', 'k');

% Axes range
xlim([-max_wavenum+1, max_wavenum]);
ylim([0, sim_days]);

title('Asymmetric Part');

hold off;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% General figure styling.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Page and figure size
set(gcf, 'PaperUnits', 'inches');
figWidth  = 7.5; % Figure width in inches.
figHeight = 6; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight])

%Colorbar
set(h, 'Colormap', spring, ...
    'CLim', [power_min, power_max])
cbh = colorbar(h(end));
cbh.Layout.Tile = 'east';

set(cbh, 'LineWidth', 0.25,...
    'Color', 'black');
set(cbh.Label, ...
    'String', 'Log(Moisture (g kg^{-1}) Power)',...
    'Color', 'black');


% Axes labels and title

titleStr = strcat(['Energy-Time Diagram' newline, 'Altitude ',...
    num2str(altW_true, '%3.1f'), ' km']);
title(tlo, titleStr, ...
    'Color', 'black');

xlabel(tlo, 'Wavenumber', ...
    'Color', 'black');
ylabel(tlo, 'Time (d)', ...
    'Color', 'black');


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Save plot.
alt_str = sprintf(['%3.1f'], abs(altW_true));

file_name = strcat([params.component_name, '_energy_evo_', ...
    alt_str, '.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end