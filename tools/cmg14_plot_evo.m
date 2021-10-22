function cmg14_plot_evo(params)

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

sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Get arrays for each state variable
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out_idxs = 0:out_freq:floor(sim_days/out_freq);
t    = zeros(size(out_idxs));
u_1  = zeros(size(out_idxs));
u_2  = zeros(size(out_idxs));
v    = zeros(size(out_idxs));
w_u  = zeros(size(out_idxs));

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'),'.nc']);
    state_file = fullfile(component_path, state_file_name);
    t(out_idx + 1)   = ncread(state_file, 't');
    u_1(out_idx + 1) = ncread(state_file, 'u_1');
    u_2(out_idx + 1) = ncread(state_file, 'u_2');
    v(out_idx + 1)   = ncread(state_file, 'v');
    w_u(out_idx + 1) = ncread(state_file, 'w_u');
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Make plots
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create master tiled layout
tlo = tiledlayout(3, 1);

% First plot - MJO indices

h(1) = nexttile(tlo, 1);

hold on;
plot(t, u_1, 'k-');
plot(t, u_2, 'k--');

hold off;

% Second plot - stochastic damping

h(2) = nexttile(tlo, 2);

hold on;
plot(t, v, 'k-');

hold off;

% Third plot - stochastic phase

h(3) = nexttile(tlo, 3);

hold on;
plot(t, w_u, 'k-');

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



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Save plot.

file_name = strcat([params.component_name, '_evo.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end