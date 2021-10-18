function cmg14_output_params(params)

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

% Create output file, delete the current one if present
params_file_name = 'params.nc';
params_file = fullfile(component_path, params_file_name);
if isfile(params_file)
    delete(params_file)
end

nccreate(params_file, 'd_u',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'd_u', 'description', ...
    'Damping for MJO modes');
ncwriteatt(params_file, 'd_u', 'units', 's^(-1)');
ncwrite(params_file, 'd_u', params.d_u);

nccreate(params_file, 'd_v',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'd_v', 'description', ...
    'Damping for stochastic damping');
ncwriteatt(params_file, 'd_v', 'units', 's^(-1)');
ncwrite(params_file, 'd_v', params.d_v);

nccreate(params_file, 'd_w',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'd_w', 'description', ...
    'Damping for stochastic phase');
ncwriteatt(params_file, 'd_w', 'units', 's^(-1)');
ncwrite(params_file, 'd_w', params.d_w);

nccreate(params_file, 'gamma',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'gamma', 'description', ...
    'Strength of non-linear interaction');
ncwriteatt(params_file, 'gamma', 'units', 'N/A');
ncwrite(params_file, 'gamma', params.gamma);

nccreate(params_file, 'a',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'a', 'description', ...
    'Background state phase of MJO modes');
ncwriteatt(params_file, 'a', 'units', 's^(-1)');
ncwrite(params_file, 'a', params.a);

nccreate(params_file, 'w_u_hat',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'w_u_hat', 'description', ...
    'Background mean state of stochastic phase');
ncwriteatt(params_file, 'w_u_hat', 'units', 's^(-1)');
ncwrite(params_file, 'w_u_hat', params.w_u_hat);

nccreate(params_file, 'sigma_u',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'sigma_u', 'description', ...
    'Strength of stochastic forcing for MJO modes');
ncwriteatt(params_file, 'sigma_u', 'units', 's^(-1)');
ncwrite(params_file, 'sigma_u', params.sigma_u);

nccreate(params_file, 'sigma_v',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'sigma_v', 'description', ...
    'Strength of stochastic forcing for stochastic damping');
ncwriteatt(params_file, 'sigma_v', 'units', 's^(-1)');
ncwrite(params_file, 'sigma_v', params.sigma_v);

nccreate(params_file, 'sigma_w',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'sigma_w', 'description', ...
    'Strength of stochastic forcing for stochastic phase');
ncwriteatt(params_file, 'sigma_w', 'units', 's^(-1)');
ncwrite(params_file, 'sigma_w', params.sigma_w);

nccreate(params_file, 'f_0',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'f_0', 'description', ...
    'Mean time-periodic damping');
ncwriteatt(params_file, 'f_0', 'units', 's^(-1)');
ncwrite(params_file, 'f_0', params.f_0);

nccreate(params_file, 'f_t',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'f_t', 'description', ...
    'Amplitude of time-periodic damping');
ncwriteatt(params_file, 'f_t', 'units', 's^(-1)');
ncwrite(params_file, 'f_t', params.f_t);

nccreate(params_file, 'w_f',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'w_f', 'description', ...
    'Frequency of time-periodic damping');
ncwriteatt(params_file, 'w_f', 'units', 's^(-1)');
ncwrite(params_file, 'w_f', params.w_f);

nccreate(params_file, 'phi',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'phi', 'description', ...
    'Phase-shift of time-periodic damping');
ncwriteatt(params_file, 'phi', 'units', 'N/A');
ncwrite(params_file, 'phi', params.phi);

nccreate(params_file, 'dt',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'dt', 'description', ...
    'Time-step size');
ncwriteatt(params_file, 'dt', 'units', 's^(-1)');
ncwrite(params_file, 'dt', params.dt);

days_to_secs = 3600*24;

nccreate(params_file, 'sim_days',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'sim_days', 'description', ...
    'Number of days in simulation');
ncwriteatt(params_file, 'sim_days', 'units', 'd');
ncwrite(params_file, 'sim_days', params.sim_days/days_to_secs);

nccreate(params_file, 'out_freq',...
    'Datatype', 'double',...
    'Format', 'netcdf4');
ncwriteatt(params_file, 'out_freq', 'description', ...
    'Number of days between output files');
ncwriteatt(params_file, 'out_freq', 'units', 'd');
ncwrite(params_file, 'out_freq', params.out_freq/days_to_secs);


end

