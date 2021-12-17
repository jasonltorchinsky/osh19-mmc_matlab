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

d_u      = ncread(params_file, 'd_u');
gamma    = ncread(params_file, 'gamma');
sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Get arrays for each state variable
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n_outfiles = floor(sim_days/out_freq);
out_idxs = 0:n_outfiles;
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

% Convert seconds to days
days_to_secs = 3600 * 24;
t = t / days_to_secs;
w_u = w_u * days_to_secs;


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
plot(t, 0*t, 'k:');

legend(h(1), 'u_{1}', 'u_{2}');

hold off;

% Second plot - stochastic damping

h(2) = nexttile(tlo, 2);

hold on;
plot(t, v, 'k-');
plot(t, 0*t+(d_u/gamma), 'k:');

legend(h(2), 'v');

hold off;

% Third plot - stochastic phase

h(3) = nexttile(tlo, 3);

hold on;
plot(t, w_u, 'k-');
plot(t, 0*t, 'k:');

ylabel(h(3), '(d^{-1})');

legend(h(3), '\omega_{u}');

hold off;


% Common across all tiles
xlabel(tlo, 'Time (d)');

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