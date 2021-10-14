function u_proj = eof_proj_u(params)

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

H        = ncread(params_file, 'H');
sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');

n_outfiles = floor(sim_days/out_freq);

% Read in grids
grid_file = fullfile(component_path, 'grid.nc');

yy  = ncread(grid_file, 'yy');
zzU = ncread(grid_file, 'zzU');

dy  = ncread(grid_file, 'dy');

L = 1490; % Equatorial meridional length scale (km)
yy_norm  = yy / L;
dy_norm  = dy / L;
zzU_norm = pi * zzU / H;

% Get the first baroclinic mode, zeroth parabolic cylinder function for the
% projections.

parab_cyl_0   = parab_cyl(yy_norm, 0);
u_clin_mode_1 = u_clin_mode(zzU_norm, 1);

% Set up u_proj, and begin calculating it
u_proj = zeros([n_outfiles, nx]);

out_idxs = 0:out_freq:n_outfiles;

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'),'.nc']);
    state_file = fullfile(component_path, state_file_name);
    
    u = ncread(state_file, 'u');
    
    % Project u onto first baroclinic mode
    u1 = zeros([ny, nx]);
    for jj = 1:ny
        for ii = 1:nx
           u1(jj, ii) = mean(squeeze(squeeze(u(jj, ii, :)).*u_clin_mode_1));
        end
    end
    
    % Project u1 onto zeroth parabolic cylinder function
    for ii = 1:nx
       u_proj(out_idx+1, ii) = dy_norm * sum(squeeze(u1(:, ii)).*parab_cyl_0); 
    end
end

end

