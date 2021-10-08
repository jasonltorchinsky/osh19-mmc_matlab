function osh19_output_state(params, bg_profs, state, out_idx, time)

% Ensure the directory structure for output is created
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, params.component_name);

if ~(isfolder(out_path))
    mkdir(out_path);
end
if ~(isfolder(exp_path))
    mkdir(exp_path);
end
if ~(isfolder(component_path))
    mkdir(component_path);
end

addpath(out_path);
addpath(exp_path);
addpath(component_path);

% Unpack some common parameters
nx = params.nx;
ny = params.ny;
nz = params.nz;

% Create output file, delete the current one if present
state_file_name = strcat(['state_', num2str(out_idx,'%04u'), '.nc']);
state_file = fullfile(component_path, state_file_name);
if isfile(state_file)
    delete(state_file)
end

% Current time in simulation
nccreate(state_file, 't',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(state_file, 't', 'description', ...
    'Time of current state');
ncwriteatt(state_file, 't', 'units', 's');
ncwrite(state_file, 't', time);

% Diagnostic variables
nccreate(state_file, 'zeta_tau',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx});
ncwriteatt(state_file, 'zeta_tau', 'description', ...
    'Vertical component of barotropic vorticity');
ncwriteatt(state_file, 'zeta_tau', 'units', 's^(-1)');
ncwrite(state_file, 'zeta_tau', state.zeta_tau);

nccreate(state_file, 'u_psi',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'u_psi', 'description', ...
    'Baroclinic zonal wind speed');
ncwriteatt(state_file, 'u_psi', 'units', 'm s^(-1)');
ncwrite(state_file, 'u_psi', 1000*state.u_psi);

nccreate(state_file, 'v_psi',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'v_psi', 'description', ...
    'Baroclinic meridional wind speed');
ncwriteatt(state_file, 'v_psi', 'units', 'm s^(-1)');
ncwrite(state_file, 'v_psi', 1000*state.v_psi);

nccreate(state_file, 'theta',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'theta', 'description', ...
    'Potential temperature anomaly');
ncwriteatt(state_file, 'theta', 'units', 'K');
ncwrite(state_file, 'theta', state.theta);

nccreate(state_file, 'q',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'q', 'description', ...
    'Moisture anomaly');
ncwriteatt(state_file, 'q', 'units', 'kg kg^(-1)');
ncwrite(state_file, 'q', state.q);

% Prognostic variables
nccreate(state_file, 'u_tau',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx});
ncwriteatt(state_file, 'u_tau', 'description', ...
    'Barotropic zonal wind speed');
ncwriteatt(state_file, 'u_tau', 'units', 'm s^(-1)');
ncwrite(state_file, 'u_tau', 1000*state.u_tau);

nccreate(state_file, 'v_tau',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx});
ncwriteatt(state_file, 'v_tau', 'description', ...
    'Barotropic meridional wind speed');
ncwriteatt(state_file, 'v_tau', 'units', 'm s^(-1)');
ncwrite(state_file, 'v_tau', 1000*state.v_tau);

nccreate(state_file, 'u',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'u', 'description', ...
    'Total zonal wind speed');
ncwriteatt(state_file, 'u', 'units', 'm s^(-1)');
ncwrite(state_file, 'u', 1000*state.u);

nccreate(state_file, 'v',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'v', 'description', ...
    'Total meridional wind speed');
ncwriteatt(state_file, 'v', 'units', 'm s^(-1)');
ncwrite(state_file, 'v', 1000*state.v);

nccreate(state_file, 'w',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'w', 'description', ...
    'Total vertical wind speed');
ncwriteatt(state_file, 'w', 'units', 'm s^(-1)');
ncwrite(state_file, 'w', 1000*state.w);

nccreate(state_file, 'p',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'p', 'description', ...
    'Pressure anomalies');
ncwriteatt(state_file, 'p', 'units', 'Pa');
ncwrite(state_file, 'p', state.p);

% Total moisture and temperature
nccreate(state_file, 'theta_tot',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'theta_tot', 'description', ...
    'Total potential temperature');
ncwriteatt(state_file, 'theta_tot', 'units', 'K');
ncwrite(state_file, 'theta_tot', state.theta + bg_profs.theta_bg_mat);

nccreate(state_file, 'q_tot',...
    'Datatype', 'double',...
    'Format', 'netcdf4',...
    'Dimensions', {'y', ny, 'x', nx, 'z', nz + 1});
ncwriteatt(state_file, 'q_tot', 'description', ...
    'Total mositure');
ncwriteatt(state_file, 'q_tot', 'units', 'kg kg^(-1)');
ncwrite(state_file, 'q_tot', state.q + bg_profs.q_bg_mat);

end

