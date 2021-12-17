function osh19_plot_qtot_vert_grad(params)

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
B_vs        = ncread(params_file, 'B_vs');
sim_days    = ncread(params_file, 'sim_days');
out_freq    = ncread(params_file, 'out_freq');

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

dz  = ncread(grid_file, 'dz');


n_outfiles = floor(sim_days/out_freq)+1;
out_idxs  = 0:n_outfiles-1;

t = zeros([1, n_outfiles]);
mean_ddz_qtot = zeros([1, n_outfiles]);

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    t(1, out_idx+1) = ncread(state_file, 't');
    q_tot = ncread(state_file, 'q_tot');
    
    ddz_q_tot = zeros([ny, nx, nz-1]);
    for kk = 2:nz
        ddz_q_tot(:,:,kk-1) = (q_tot(:,:,kk+1) - q_tot(:,:,kk-1)) ...
            / (2 * dz);
    end
    
    mean_ddz_qtot(1, out_idx+1) = mean(ddz_q_tot, 'all');
    
end

days_to_secs = 3600*24;
t = t / days_to_secs;

tlo = tiledlayout(1, 1);

% Plot non-dimensionalize mean vertical moisture gradient
h = nexttile(tlo, 1);

plot(t, mean_ddz_qtot / B_vs, 'k-');

xlabel(h, 'Time (d)');
ylabel(h, '\partial_{z} q_{tot} / B_{vs}');
title(h, 'Non-Dimensional Mean Vertical Moisture Gradient');

% Page and figure size
set(gcf, 'PaperUnits', 'inches');
figWidth  = 7.5; % Figure width in inches.
figHeight = 6; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight])

% Save plot.

file_name = strcat([params.component_name, '_mean_vert_q_grad.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end

