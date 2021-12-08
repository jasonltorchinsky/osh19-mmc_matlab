function [q_proj, q_std] = eof_proj_q(params)

% Get path to output data.
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, params.component_name);

addpath(out_path);
addpath(exp_path);
addpath(component_path);

% Read in run parameters
params_file = fullfile(component_path, 'params.nc');

nx       = ncread(params_file, 'nx');
ny       = ncread(params_file, 'ny');
nz       = ncread(params_file, 'nz');

H        = ncread(params_file, 'H');
sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');

n_outfiles = floor(sim_days/out_freq) + 1;

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

yy  = ncread(grid_file, 'yy');
zzW = ncread(grid_file, 'zzW');

dy  = ncread(grid_file, 'dy');

L = 1490; % Equatorial meridional length scale (km)
yy_norm  = yy / L;
dy_norm  = dy / L;
zzW_norm = pi * zzW / H;

% Get the first baroclinic mode, zeroth parabolic cylinder function for the
% projections.
parab_cyl_0   = parab_cyl(yy_norm, 0);
q_clin_mode_1 = q_clin_mode(zzW_norm, 1);
q_clin_mode_2 = q_clin_mode(zzW_norm, 2);

Q_mode = params.Q_mode;

% Get the normalization constants due to discretized norms
merid_norm  = dy_norm * (parab_cyl_0.' * parab_cyl_0);
vert_norm_1 = (1/(nz+1)) * (q_clin_mode_1.' * q_clin_mode_1);
vert_norm_2 = (1/(nz+1)) * (q_clin_mode_2.' * q_clin_mode_2);

% Set up u_proj, and begin calculating it
q_proj = zeros([n_outfiles, nx]);

out_idxs = 0:(n_outfiles-1);

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'),'.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    q = ncread(state_file, 'q');
    
    % Project q onto first two baroclinic modes
    q1 = zeros([ny, nx]);
    q2 = zeros([ny, nx]);
    for jj = 1:ny
        for ii = 1:nx
           q1(jj, ii) = (1/(nz+1)) ...
               * squeeze(q_clin_mode_1.'*squeeze(q(jj, ii, :))) ...
               * (1/vert_norm_1);
           q2(jj, ii) = (1/(nz+1)) ...
               * squeeze(q_clin_mode_2.'*squeeze(q(jj, ii, :))) ...
               * (1/vert_norm_2);
        end
    end
    
    % Pick Q_mid or Q_up
    if strcmpi(Q_mode, 'mid')
        Q = 1 / sqrt(3) * (q1 * q_clin_mode(pi/3, 1) + q2 * q_clin_mode(pi/3, 2));
    else
        Q = 1 / sqrt(3) * (q1 * q_clin_mode(2*pi/3, 1) - q2 * q_clin_mode(2*pi/3, 2));
    end
    
    % Project Q onto zeroth parabolic cylinder function
    for ii = 1:nx
       q_proj(out_idx+1, ii) = dy_norm * squeeze(Q(:, ii).'*parab_cyl_0) ...
        * (1/merid_norm);
    end
end

% Get standard deviation from q_proj, after taking out time-series mean, and
% non-dimensionalize it
q_proj = detrend(q_proj, 0);
q_std  = std(q_proj, 0, 'all');
q_proj = q_proj / q_std;

end

